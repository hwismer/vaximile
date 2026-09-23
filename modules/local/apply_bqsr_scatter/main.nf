process APPLY_BQSR_SCATTER {

    /*
    Apply base quality score recalibration on a provided interval.
    */
    
    // process_low, matching nf-core's gatk4/applybqsr. ApplyBQSR just rewrites quality
    // scores and does not need the 32 GB of process_low_memory; it also runs once per
    // sample per interval, so the saving is multiplied by scatter_count.
    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "ApplyBQSR on ${meta.sample_name} ${interval_shard}"

    input:
        tuple val(meta), path(markdup_bam), path(markdup_bam_bai), path(recal_table), path(interval_shard)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(meta), path("*_bqsr.bam"), emit: bam
    script:
    // Explicit heap, as nf-core's GATK4 modules do. Without --java-options the JVM picks
    // its own maximum, which is a fraction of whatever memory it believes it has - not
    // necessarily the amount the scheduler granted. Sizing it from task.memory keeps the
    // heap inside the reservation, which matters more now the reservation is smaller.
    def avail_mem = (task.memory.mega * 0.8).intValue()
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}_${interval_shard}"
    """
    gatk --java-options "-Xmx${avail_mem}M" ApplyBQSR \
        -R $reference_fa \
        -I $markdup_bam \
        -L $interval_shard \
        --bqsr-recal-file $recal_table \
        -O "${prefix}_bqsr.bam"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}_${interval_shard}"
    """
    touch ${prefix}_bqsr.bam
    """

}
