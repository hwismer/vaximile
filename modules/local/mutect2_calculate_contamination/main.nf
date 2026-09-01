process MUTECT2_CALCULATE_CONTAMINATION {

    label 'process_low_memory'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Calculalting contamination for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(tumor_pileups), path(normal_pileups)

    output:
        tuple val(somatic_meta), path("*_contamination.table"), emit: table
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"

    """
    gatk CalculateContamination \
        -I $tumor_pileups \
        -matched $normal_pileups \
        -O "${prefix}_contamination.table"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"

    """
    touch "${prefix}_contamination.table"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """

}
