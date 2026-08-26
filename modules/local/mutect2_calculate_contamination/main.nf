process MUTECT2_CALCULATE_CONTAMINATION {

    cpus 4
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Calculalting contamination for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(tumor_pileups), path(normal_pileups)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_contamination.table")

    script:

    """
    gatk CalculateContamination \
        -I $tumor_pileups \
        -matched $normal_pileups \
        -O "${somatic_meta.somatic_name}_contamination.table"
    """

}
