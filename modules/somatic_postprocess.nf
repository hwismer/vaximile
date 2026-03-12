process INDEX_FINAL_VCF {
    /*

    tabix index a vcf file

    */

    cpus 2
    memory "8GB"

    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/${meta.somatic_name}/variants", mode: "copy"

    input:
        tuple val(meta), path(vcf)
        

    output:
        tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
        """
        bgzip $vcf
        tabix -p vcf "${vcf}.gz"
        """

}


process VCF_TO_TABLE {
    cpus 1
    memory "8GB"

    container "broadinstitute/gatk:4.6.1.0"

    publishDir "${params.outdir}/${meta.somatic_name}/variants", mode: "copy"

    input:
        tuple val(meta), path(vcf), path(vcf_index)

    output:
        path "${meta.somatic_name}_variants.tsv"

    script:
        """

        gatk VariantsToTable \
            -V $vcf \
            -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC -F AF -F set -F FILTER -F CSQ \
            -GF AD -GF DP -GF GT -GF AF \
            -GF RDP -GF RAF -GF RAD -GF RADF -GF RADR -GF TX -GF GX \
            -O "${meta.somatic_name}_variants.tsv"

        """


}

process MERGE_SOMATIC_VCFS {

    /*

    Use the deprecated CombineVariants from GATK 3.6.0 to combine vcf files. 
    
    Currently just from mutect + strelka, but I should probably change this to
    take in a dynamic number of vcf files.

    */

    cpus 4
    memory "32GB"
    cache "lenient"
    
    container "broadinstitute/gatk3:3.6-0"

    publishDir "${params.outdir}/${merged_meta.somatic_name}/somatic/", mode: "copy"

    input:
        tuple val(somatic_name), val(call_sets) 
        path(reference_fa)
        path(reference_index_dir)
    
    output:
        tuple val(merged_meta), path("${merged_meta.somatic_name}_variants.vcf.gz")


    script:

        def call1 = call_sets[0]
        def call2 = call_sets[1]

        def vcf1_meta = call1[0]
        def vcf1  = call1[1]
        def tbi1  = call1[2]

        def vcf2_meta = call2[0]
        def vcf2  = call2[1]
        def tbi2  = call2[2]

        merged_meta = [
            somatic_name: vcf1_meta.somatic_name,
            somatic_caller: 'mutect_strelka',
            tumor_metamap: vcf1_meta.tumor_metamap,
            normal_metamap: vcf1_meta.normal_metamap
        ]

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            -genotypeMergeOptions PRIORITIZE \
            --rod_priority_list mutect,strelka \
            -V:${vcf1_meta.somatic_caller} $vcf1 \
            -V:${vcf2_meta.somatic_caller} $vcf2 \
            -o "${merged_meta.somatic_name}_variants.vcf.gz"

        """
    
}

process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */
    
    cpus 4
    memory "16GB"
    
    container "biocontainers/bcftools:v1.9-1-deb_cv1"
    
    publishDir "${params.outdir}/${somatic_meta.somatic_name}/somatic/${somatic_meta.somatic_caller}/", mode:"copy"

    
    input:
        tuple val(somatic_meta), path(somatic_vcf)
        
    output:
        tuple val(somatic_meta), 
            path("${somatic_meta.somatic_name}_variants.vcf.gz"), 
            path("${somatic_meta.somatic_name}_variants.vcf.gz.tbi"), emit: filter_vcf
        
        tuple val(somatic_meta), 
            path("${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz"), 
            path("${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
        
    script:
        """
        cp $somatic_vcf "${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz"
        bcftools index -t "${somatic_meta.somatic_name}_variants_unfiltered.vcf.gz"
        bcftools view -f PASS -Oz -o "${somatic_meta.somatic_name}_variants.vcf.gz" $somatic_vcf
        bcftools index -t "${somatic_meta.somatic_name}_variants.vcf.gz"
        """
        


}


process ADD_VCF_GT_FIELD {

    /*

    Adds the GT fields to a vcf lacking it. This is currently used for Strelka/Manta vcfs and the field is populated by 0/1 by default.

    */

    cpus 2
    memory "16GB"
    
    container "griffithlab/vatools:5.2.0"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
         
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${somatic_meta.somatic_caller}_gt.vcf.gz"), emit: gt_vcf

    script:
        """
        vcf-genotype-annotator $somatic_vcf \
            "${somatic_meta.tumor_metamap.sample_name}" \
            0/1 \
            -o "${somatic_meta.somatic_name}_${somatic_meta.somatic_caller}_gt.vcf.gz"

        """

}

process VEP_FILTER {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 2
    memory "16GB"
    cache "lenient"


    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(vep_vcf)
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("${meta.somatic_name}_vep.vcf"), emit: vep_vcf

    script:
        """

        filter_vep -i $vep_vcf \
            -o "${meta.somatic_name}_vep.vcf" \
            --format vcf \
            --filter "gnomADe_AF < 0.001 or not gnomADe_AF"

        """
}

process VT_SOMATIC_POSTPROCESS {

    /*

        Postprocess somatic variants from a vcf. This is used to normalize variants and decompose biallelics
        into different entries. The vcf is then sorted and duplicate entries removed. The final vcf is then indexed.

    */

    cpus 4
    memory "16GB"
    cache "lenient"

    container "staphb/bcftools:1.23"

    publishDir "${params.outdir}/${meta.somatic_name}/somatic/${meta.somatic_caller}/", mode:"copy"

    input:
        tuple val(meta), path(somatic_vcf), path(somatic_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(meta),
            path("${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz"),
            path("${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz.tbi"), emit: vt_vcf

    script:

        """
        bcftools norm -m -any -f $reference_fa $somatic_vcf -Oz -o norm_vcf.vcf.gz
        bcftools sort norm_vcf.vcf.gz -Oz -o norm_sort.vcf.gz
        bcftools norm -d exact norm_sort.vcf.gz -Oz -o "${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz"
        bcftools index -t "${meta.somatic_name}_${meta.somatic_caller}_normalized_variants.vcf.gz"
        """
}
