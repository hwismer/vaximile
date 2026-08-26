process APPLY_BQSR_SCATTER {

    /*
    Apply base quality score recalibration on a provided interval.
    */
    
    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "ApplyBQSR on ${meta.sample_name} ${interval_shard}"

    input:
        tuple val(meta), path(markdup_bam), path(markdup_bam_bai), path(recal_table), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(meta), path("*_bqsr.bam")
    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}_${interval_shard}"
    """
    gatk ApplyBQSR \
        -R $reference_fa \
        -I $markdup_bam \
        -L $interval_shard \
        --bqsr-recal-file $recal_table \
        -O "${prefix}_bqsr.bam"
    """

}
