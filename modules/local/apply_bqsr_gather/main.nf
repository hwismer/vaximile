process APPLY_BQSR_GATHER {

    // Merge the BQSR shards into one sorted, indexed BAM. samtools merge rather than concatenation, since
    // shards overlap at their boundaries; -c/-p keep the shared @RG/@PG from being renamed.

    // Threads help BAM compression; the merge itself is serial.
    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "GatherBams on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), path(bams)
    output:
        tuple val(meta),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam"),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam.bai"), emit: bam

    script:
    def prefix = "${sample_name}_${molecule}_bqsr"
    """
    samtools merge \\
        -c -p \\
        -@ $task.cpus \\
        --write-index \\
        -o "${prefix}.bam##idx##${prefix}.bam.bai" \\
        ${bams.join(' ')}
    """

    stub:
    """
    touch "${sample_name}_${molecule}_bqsr.bam"
    touch "${sample_name}_${molecule}_bqsr.bam.bai"
    """

}
