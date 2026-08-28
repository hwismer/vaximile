process APPLY_BQSR_GATHER {

    /*
    Gathers all bqsr games to create a final merged bam with adjusted base quality scores.
    */

    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GatherBams on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), path(bams)
    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_bqsr.bam")

    script:
    // GatherBamFiles concatenates without re-sorting, so the shards have to arrive in
    // coordinate order. groupTuple does not preserve order, so sort explicitly here.
    // SplitIntervals names shards with a zero-padded index (0000-scattered.interval_list,
    // 0001-...), which makes sorting by filename equivalent to coordinate order.
    // MUTECT2_GATHER_VCFS and HAPLOTYPE_CALLER_GATHER_VCFS do the same thing.
    def sorted_bams = bams.toSorted { a, b -> a.name <=> b.name }
    """
    gatk GatherBamFiles \
        ${sorted_bams.collect { "-I ${it}" }.join(' ')} \
        -O "${sample_name}_${molecule}_bqsr.bam"
    """

}
