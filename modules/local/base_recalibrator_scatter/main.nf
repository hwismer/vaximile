process BASE_RECALIBRATOR_SCATTER {

    // BaseRecalibrator on one interval shard.

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "BaseRecalibrator on ${meta.sample_name} ${interval_shard}"

    input:
        tuple val(meta), path(markdup_bam), path(markdup_bam_bai), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path(reference_dict)
        tuple path(known_sites_dbsnp), path(known_sites_dbsnp_index)
        tuple path(known_sites_1000g_snps), path(known_sites_1000g_snps_index)
        tuple path(known_indels), path(known_indels_index)
        tuple path(mills), path(mills_index)

    output:
        tuple val(meta), path("*_recal_table.table"), emit: table

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}_${interval_shard}"
    """
    gatk BaseRecalibrator \
        -I $markdup_bam \
        -O "${prefix}_recal_table.table" \
        -R $reference_fa \
        --known-sites $known_sites_dbsnp \
        --known-sites $known_sites_1000g_snps \
        --known-sites $known_indels \
        --known-sites $mills \
        -L $interval_shard

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}_${interval_shard}"
    """
    touch "${prefix}_recal_table.table"
    """
}
