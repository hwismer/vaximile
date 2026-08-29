process SPLIT_INTERVALS {

    /*
        Takes a file of genomic intervals such as a picard interval file or BED file splits into scatter_count number of shards for parallel processing.

    */

    label 'process_low'
    cache "lenient"
    
    tag "Splitting ${intervals_file} into ${scatter_count} shards w/ ${interval_padding} bp padding"

    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple path(reference_fa), path(reference_fa_index)
        path reference_dict
        tuple val(capture_kit), path(intervals_file)
        val scatter_count
        val interval_padding

    output:
         tuple val(capture_kit), path("*-scattered.interval_list"), emit: interval_shards
         path "versions.yml", topic: versions

    script:
    """
    gatk SplitIntervals \
        -R $reference_fa \
        -L $intervals_file \
        --scatter-count $scatter_count \
        --interval-padding $interval_padding \
        -O .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    """
    for i in \$(seq 1 ${scatter_count}); do
        touch "\$(printf '%04d' \$((i - 1)))-scattered.interval_list"
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}

//***************************************************************************************************************************
// RESOURCE PULLING FROM WEB SOURCES
