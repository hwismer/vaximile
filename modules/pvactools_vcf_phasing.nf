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

        tuple val(somatic_meta), path(somatic_vcf), path(somatic_vcf_index)
        tuple path(reference_fa), path(reference_fai), path(reference_dict)

    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_tumor_only.vcf.gz"), path("${somatic_meta.tumor_meta.sample_name}_tumor_only.vcf.gz.tbi")

    script:
        """
        gatk SelectVariants \
            -V $somatic_vcf \
            -R "${reference_fa}" \
            --sample-name ${somatic_meta.tumor_meta.sample_name} \
            --create-output-variant-index \
            -O ${somatic_meta.tumor_meta.sample_name}_tumor_only.vcf.gz
        """

}

process PHASE_VCF_COMBINE_VARIANTS {

    /*

    Part of creating a phased germline vcf.
    Combines the variants from the tumor-only vcf and the germline_vcf (which has been renamed).


    */

    cpus 2
    memory "16GB"

    container "broadinstitute/gatk3:3.6-0"

    input:
        tuple val(somatic_meta), path(tumor_only_vcf), path(tumor_only_index), path(germline_vcf), path(germline_index)  
        tuple path(reference_fa), path(reference_fai), path(reference_dict)


    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_combined_somatic_plus_germline.vcf")

    script:
        """
        java -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
                -R $reference_fa \
                --variant $germline_vcf \
                --variant $tumor_only_vcf \
                -o ${somatic_meta.tumor_meta.sample_name}_combined_somatic_plus_germline.vcf \
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
        tuple val(somatic_meta), path(combined_vcf)
        tuple path(reference_fa), path(reference_fai), path(reference_dict)

    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_combined.sorted.vcf"), emit: sorted_vcf

    script:
        """
        java -jar /usr/picard/picard.jar \
            SortVcf \
                -I $combined_vcf \
                -O ${somatic_meta.tumor_meta.sample_name}_combined.sorted.vcf \
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

    container "staphb/bcftools:1.23.1" 

    input:
        tuple val(normal_meta), val(somatic_meta), path(germline_vcf), path(germline_vcf_index)

    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz"), path("${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz.tbi")

    script:
        """
        cat > sample_map.txt <<EOF
        ${normal_meta.sample_name} ${somatic_meta.tumor_meta.sample_name}
        EOF
        
        bcftools reheader \
            -N sample_map.txt \
            --threads $task.cpus \
            -o ${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz \
            $germline_vcf

        bcftools index -t --threads $task.cpus ${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz 
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
        tuple val(tumor_meta), val(somatic_meta), path(combined_sorted_vcf), path(tumor_reads), path(tumor_reads_index)
        tuple path(reference_fa), path(reference_index), path(reference_dict)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_phased.vcf")

    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_fa \
                -I $tumor_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o ${somatic_meta.somatic_name}_phased.vcf

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
