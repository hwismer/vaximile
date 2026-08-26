process APPLY_BQSR_SCATTER {

    /*
    Apply base quality score recalibration on a provided interval.
    */
    
    cpus 4
    memory "12GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "ApplyBQSR on ${meta.sample_name} ${interval_shard}"

    input:
        tuple val(meta), path(markdup_bam), path(markdup_bam_bai), path(recal_table), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_${interval_shard}_bqsr.bam")
    script:
    """
    gatk ApplyBQSR \
        -R $reference_fa \
        -I $markdup_bam \
        -L $interval_shard \
        --bqsr-recal-file $recal_table \
        -O "${meta.sample_name}_${meta.molecule}_${interval_shard}_bqsr.bam"
    """

}
