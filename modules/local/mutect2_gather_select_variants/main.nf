process MUTECT2_GATHER_SELECT_VARIANTS {

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Selecting Variants within interval for ${somatic_meta.somatic_name} in $interval_shard"

    input:
        tuple val(somatic_meta), val(somatic_name), path(vcf), path(vcf_index),  path(interval_shard)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_${interval_shard}.vcf.gz")

    script:
    
    """
    gatk SelectVariants \
        -V $vcf \
        -L $interval_shard \
        -O "${somatic_name}_${interval_shard}.vcf.gz"
    """

    stub:
    """
    touch "${somatic_name}_${interval_shard}.vcf.gz"
    """

}
