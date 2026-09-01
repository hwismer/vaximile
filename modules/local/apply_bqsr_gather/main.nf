process APPLY_BQSR_GATHER {

    /*
    Gathers all bqsr shards into a final coordinate-sorted, indexed BAM with adjusted base
    quality scores.

    `samtools merge`, not `gatk GatherBamFiles`. GatherBamFiles concatenates, and
    concatenation cannot produce a sorted BAM here: APPLY_BQSR_SCATTER passes -L <shard>,
    and GATK emits every read *overlapping* the interval rather than only those starting
    inside it, so a read spanning a boundary is written near the start of the later shard
    while beginning before the end of the earlier one. Concatenating left positions out of
    order at every seam - by less than a read length, so the header still claimed
    SO:coordinate and only `samtools index` noticed:

        [E::hts_idx_push] Unsorted positions on sequence #1: 52263835 followed by 52263738

    That used to be papered over by re-sorting the whole gathered BAM afterwards. A k-way
    merge of already-sorted shards gets the same result in a single streaming pass, so the
    separate sort is gone rather than merely hidden.

    -c is not optional. Every shard carries the same @RG, and without it merge treats those
    as colliding IDs and renames all but one - "S1" becomes "S1-7A2E5CD9" - which would split
    one sample's reads across two read groups and break the GATK steps downstream. -p does
    the same for @PG.

    Shard order does not matter to merge, so the filename sort the GatherBamFiles version
    needed is gone too.
    */

    // process_high (8 CPU) rather than process_medium (4): -@ feeds samtools merge, and
    // this took over the compression work the removed SORT_BAM used to do. Not higher - the
    // k-way merge itself is serial and the threads only go to BAM (de)compression, which
    // flattens out around 8.
    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "GatherBams on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), path(bams)
    output:
        tuple val(meta),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam"),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam.bai"), emit: bam
        path "versions.yml", topic: versions

    script:
    def prefix = "${sample_name}_${molecule}_bqsr"
    """
    samtools merge \\
        -c -p \\
        -@ $task.cpus \\
        --write-index \\
        -o "${prefix}.bam##idx##${prefix}.bam.bai" \\
        ${bams.join(' ')}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch "${sample_name}_${molecule}_bqsr.bam"
    touch "${sample_name}_${molecule}_bqsr.bam.bai"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """

}
