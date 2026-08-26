process MUTECT2_FILTER_MUTECT_CALLS {
    
    cpus 2
    memory "16GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Filtering Mutect calls for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(unfiltered_vcf), path(orientation_model), path(stats), path(contamination_table)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_mutect_filtered.vcf.gz"), path("${somatic_meta.somatic_name}_mutect_filtered.vcf.gz.tbi"), emit: filtered_vcf


    script:
    """
    gatk FilterMutectCalls \
        -R $reference_fa \
        -V $unfiltered_vcf \
        --orientation-bias-artifact-priors $orientation_model \
        --contamination-table $contamination_table \
        -stats $stats \
        --create-output-variant-index \
        -O "${somatic_meta.somatic_name}_mutect_filtered.vcf.gz"
    """
}
