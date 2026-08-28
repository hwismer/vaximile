process HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS {

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Select variants from ${sample_meta.sample_name} within ${interval_shard}"

    input:
        tuple val(sample_meta), val(sample_name), path(vcf), path(vcf_index), path(interval_shard)

    output:
        tuple val(sample_meta), path("${sample_meta.sample_name}_${interval_shard}.vcf.gz")

    script:
    
    """
    gatk SelectVariants \
        -V $vcf \
        -L $interval_shard \
        -O "${sample_name}_${interval_shard}.vcf.gz"
    """
}
