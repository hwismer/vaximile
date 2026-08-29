process MUTECT2_FILTER_MUTECT_CALLS {
    
    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Filtering Mutect calls for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(unfiltered_vcf), path(orientation_model), path(stats), path(contamination_table)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(somatic_meta), path("*_mutect_filtered.vcf.gz"), path("*_mutect_filtered.vcf.gz.tbi"), emit: filtered_vcf
        path "versions.yml", topic: versions


    script:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    """
    gatk FilterMutectCalls \
        -R $reference_fa \
        -V $unfiltered_vcf \
        --orientation-bias-artifact-priors $orientation_model \
        --contamination-table $contamination_table \
        -stats $stats \
        --create-output-variant-index \
        -O "${prefix}_mutect_filtered.vcf.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    """
    touch ${prefix}_mutect_filtered.vcf.gz
    touch ${prefix}_mutect_filtered.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
