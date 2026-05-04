process COMBINE_FASTQS {

    cpus 2
    memory "8GB"
    
    tag "Combing FASTQs from ${meta.sample_name}"

    input:
        tuple val(meta), path(fastqs_r1), path(fastqs_r2)
    output:
        tuple val(meta), path("${meta.somatic_name}_R1.merged.fastq.gz"), path("${meta.somatic_name}_R2.merged.fastq.gz")


    script:
    """
    cat ${fastqs_r1.join(' ')} > ${meta.somatic_name}_R1.merged.fastq.gz
    cat ${fastqs_r2.join(' ')} > ${meta.somatic_name}_R2.merged.fastq.gz
    """


}
process SPLIT_INTERVALS {

    /*

        Given a file of genomic intervals, split into scatter_count number of shards.

    */

    cpus 2
    memory "8GB"
    cache "lenient"
    
    tag "Splitting ${intervals_file} into ${scatter_count} shards w/ ${interval_padding} bp padding"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple path(reference_fa), path(reference_fa_index), path(reference_dict)
        path intervals_file
        val scatter_count
        val interval_padding

    output:
        tuple val(interval_padding), path("*-scattered.interval_list"), emit: interval_shards

    script:
    """
    gatk SplitIntervals \
        -R $reference_fa \
        -L $intervals_file \
        --scatter-count $scatter_count \
        --interval-padding $interval_padding \
        -O .

    """
}
