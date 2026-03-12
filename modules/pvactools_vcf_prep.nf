process BAMREADCOUNT {

    /*

    Run bamreadcount to add coverage information to a vcf using a BAM file. The sample names in the VCF
    must match the read groups and sample names in the BAM.

    */

    cpus 4
    memory "32GB"
    cache "lenient"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    publishDir "${params.outdir}/${vcf_meta.somatic_name}/coverage/", mode: "copy"

    input:
        tuple val(vcf_meta), path(vcf), val(bam_meta), path(bam), path(bai)
        path reference_fa

    output:

        tuple val(vcf_meta), path(vcf), val(bam_meta),
            path("${bam_meta.sample_name}_${bam_meta.molecule}_bamrc_helper/*indel.tsv"), 
            path("${bam_meta.sample_name}_${bam_meta.molecule}_bamrc_helper/*snv.tsv"), emit: brc_files
    script:
        
        """
        mkdir ${bam_meta.sample_name}_${bam_meta.molecule}_bamrc_helper
        bam_readcount_helper.py \
            $vcf \
            ${bam_meta.sample_name} \
            $reference_fa \
            $bam \
            ${bam_meta.molecule} \
            ${bam_meta.sample_name}_${bam_meta.molecule}_bamrc_helper
        """

}



process VEP_ANNOTATE {

    /*

    Use VEP to annotate a vcf files. Requires path to an installed cache
    as well as a vep plugin directory. If using PVAC later on, the vep plugins
    will need to includet those specified by pvac in their docs.

    */

    cpus 8
    memory "16GB"
    cache "lenient"

    
    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(somatic_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("somatic_vep.vcf"), emit: vep_vcf

    script:
        """
        vep \
            --input_file $somatic_vcf  \
            --output_file somatic_vep.vcf \
            --everything \
            --format vcf --vcf --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs \
            --fasta $reference_fa  \
            --offline --cache \
            --plugin Frameshift --plugin Wildtype \
            --pick \
            --dir_plugins $vep_plugins \
            --dir_cache $vep_cache \
            --transcript_version
        
        """
}

process ANNOTATE_VCF_COVERAGE {

    /*

        Use the output of bamreadcount to annotate coverage given a particular sample.
    */

    cpus 4
    memory "32GB"

    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/${somatic_meta.somatic_name}/coverage/", mode: "copy"

    input:
        tuple val(id), val(somatic_meta), path(vcf), path(tumor_bamrc), path(normal_bamrc), path(tumor_rna_bamrc)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_annotated.vcf")

    script:
        def tumor_indels = tumor_bamrc[0]
        def tumor_snvs = tumor_bamrc[1]

        def normal_indels = normal_bamrc[0]
        def normal_snvs = normal_bamrc[1]

        def tumor_rna_indels = tumor_rna_bamrc[0]
        def tumor_rna_snvs = tumor_rna_bamrc[1]

        """

        echo "${somatic_meta.tumor_metamap.sample_name}"
        echo "${somatic_meta.tumor_metamap.molecule}"

        vcf-readcount-annotator \
            $vcf \
            $tumor_indels \
            DNA \
            -s ${somatic_meta.tumor_metamap.sample_name} \
            -t indel \
            -o vcf1.vcf

        vcf-readcount-annotator \
            vcf1.vcf \
            $tumor_snvs \
            DNA \
            -s ${somatic_meta.tumor_metamap.sample_name} \
            -t snv \
            -o vcf2.vcf


        vcf-readcount-annotator \
            vcf2.vcf \
            $normal_snvs \
            DNA \
            -s ${somatic_meta.normal_metamap.sample_name} \
            -t snv \
            -o vcf3.vcf

        vcf-readcount-annotator \
            vcf3.vcf \
            $normal_indels \
            DNA \
            -s ${somatic_meta.normal_metamap.sample_name} \
            -t indel \
            -o vcf4.vcf

        vcf-readcount-annotator \
            vcf4.vcf \
            $tumor_rna_snvs \
            RNA \
            -s ${somatic_meta.tumor_metamap.sample_name} \
            -t snv \
            -o vcf5.vcf

        vcf-readcount-annotator \
            vcf5.vcf \
            $tumor_rna_indels \
            RNA \
            -s ${somatic_meta.tumor_metamap.sample_name} \
            -t indel \
            -o ${somatic_meta.somatic_name}_annotated.vcf

        """

}

process ANNOTATE_VCF_EXPRESSION {

    /*

    Use the abundance estimates from kallist to annotate transcript expression in a vcf file.

    */
    cpus 2
    memory "16GB"

    container "griffithlab/vatools:5.2.0"

    publishDir "${params.outdir}/${vcf_meta.somatic_name}/coverage/", mode: "copy"

    input:
        tuple val(vcf_meta), path(vcf), val(kallisto_meta), path(tx_abundance), path(gene_abundance)


    output:
        tuple val(vcf_meta), path("${vcf_meta.somatic_name}_cov_expr_annotated.vcf")

    script:

        """
        vcf-expression-annotator \
            $vcf \
            -s ${kallisto_meta.sample_name} \
            $tx_abundance \
            kallisto transcript \
            -o tx.vcf

        vcf-expression-annotator \
            tx.vcf \
            $gene_abundance \
            custom gene \
            -i ENSEMBLID \
            -e TPM \
            -s ${kallisto_meta.sample_name} \
            --ignore-ensembl-id-version \
            -o "${vcf_meta.somatic_name}_cov_expr_annotated.vcf"

        """

}

process PHASE_VCF_SELECT_VARIANTS {

    /*

    Part of creating a phased germline vcf
    Takes a somatic vcf and extracts just the tumor sample.

    */

    cpus 2
    memory "32GB"
    cache "lenient"


    container "broadinstitute/gatk:4.6.1.0"

    input:

        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path("tumor_only.vcf.gz"), path("tumor_only.vcf.gz.tbi")

    script:
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_fa}" \
            --sample-name ${somatic_meta.tumor_metamap.sample_name} \
            -O tumor_only.vcf.gz

        gatk IndexFeatureFile \
            -I tumor_only.vcf.gz

        """

}

process PHASE_VCF_COMBINE_VARIANTS {

    /*

    Part of creating a phased germline vcf.
    Combines the variants from the tumor-only vcf and the germline_vcf (which has been renamed).


    */

    cpus 2
    memory "32GB"

    container "broadinstitute/gatk3:3.6-0"

    input:

        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path(tumor_only_vcf), path(tumor_only_vcf_index)

        path(reference_fa)
        path(reference_index_files)



    output:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path("combined_somatic_plus_germline.vcf")

    script:
        """

        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_fa \
                --variant $germline_vcf \
                --variant $tumor_only_vcf \
                -o combined_somatic_plus_germline.vcf \
                --assumeIdenticalSamples

        """

}

process PHASE_VCF_SORT_VCF {

    /*

    Part of creating a phased vcf.
    Sort the somatic + germline combined vcf.

    */

    cpus 2
    memory "32GB"

    container 'broadinstitute/picard:3.4.0'

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path(combined_vcf)
        path(reference_dict)

    output:
        tuple val(somatic_meta), path("combined.sorted.vcf"), emit: vcf

    script:
        """
        java -jar /usr/picard/picard.jar \
            SortVcf \
                -I $combined_vcf \
                -O combined.sorted.vcf \
                -SD $reference_dict

        """

}

process PHASE_VCF_RENAME {

    /*

    Part of creating a phase vcf file.
    Germline sample name will be the name of the NORMAL sample, but to combine variants with the tumor sample, the names must match.
    Here the sample name in the germline vcf is renamed to the name of the tumor sample.

    */

    cpus 2
    memory "16GB"

    container "biocontainers/bcftools:v1.9-1-deb_cv1"

    input:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path(germline_vcf), path(germline_vcf_index),
            path(tumor_only_vcf), path(tumor_only_vcf_index)

    output:
        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index),
            val(germline_meta), path("vt_germline_rename.vcf.gz" ), path("vt_germline_rename.vcf.gz.tbi"),
            path(tumor_only_vcf), path(tumor_only_vcf_index)

    script:
        """
        echo ${somatic_meta.tumor_metamap.sample_name} > new_names.txt

        bcftools reheader \
            --samples new_names.txt \
            --output vt_germline_rename.vcf.gz \
            $germline_vcf

        bcftools index -t vt_germline_rename.vcf.gz
        """


}

process PHASE_VCF_RBPHASING {

    /*

    Use deprecated ReadBackedPhasing from GATK 3.6.0 to phase variants
    in the somatic + germline combined and sorted vcf

    */

    cpus 4
    memory "32GB"


    container "broadinstitute/gatk3:3.6-0"

    input:
        tuple val(somatic_meta), path(combined_sorted_vcf), val(tumor_meta), path(tumor_reads), path(tumor_reads_index)
        path(reference_fa)
        path(reference_index_dir)

    output:
        tuple val(somatic_meta), path("phased.vcf")

    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_fa \
                -I $tumor_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o phased.vcf

        """
}

process PHASE_VCF_VEP {

    /*

        Use VEP to annotated the phased vcf.

    */
    cpus 8
    memory "32GB"


    container "ensemblorg/ensembl-vep:release_115.0"

    input:
        tuple val(meta), path(phased_vcf)
        path reference_fa
        path vep_cache
        path vep_plugins

    output:
        tuple val(meta), path("phased_vcf_vep.vcf"), emit: vcf

    script:
        """
        vep \
            --input_file $phased_vcf  \
            --output_file "phased_vcf_vep.vcf" \
            --format vcf --vcf --symbol --terms SO --tsl --biotype \
            --hgvs --fasta $reference_fa  \
            --offline --cache $vep_cache \
            --plugin Frameshift --plugin Wildtype --plugin Downstream \
            --pick \
            --fork ${task.cpus} \
            --dir_plugins $vep_plugins \
            --transcript_version

        """

}

process PHASE_VCF_INDEX {

    /*

        Index the final phased vcf.

    */

    cpus 1
    memory "16GB"

    conda "bioconda::tabix=0.2.6"

    publishDir "${params.outdir}/${meta.somatic_name}/variants/phased_variants/", mode: "copy"

    input:
        tuple val(meta), path(phased_vcf)

    output:
        tuple val(meta), path("${meta.somatic_name}_phased_annotated.vcf.gz"), path("${meta.somatic_name}_phased_annotated.vcf.gz.tbi"), emit: phased_vcf


    script:
        """
        echo ${meta.somatic_name}
        bgzip -c $phased_vcf > ${meta.somatic_name}_phased_annotated.vcf.gz

        tabix -p vcf ${meta.somatic_name}_phased_annotated.vcf.gz
        """

}
