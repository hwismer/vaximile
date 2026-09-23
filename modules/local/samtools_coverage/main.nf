process SAMTOOLS_COVERAGE {
    
    // `samtools coverage` has no -@/--threads option at all, so process_high was
    // reserving 8 CPUs it could never use. nf-core's samtools/coverage is process_single.
    label 'process_single'
    conda "bioconda::samtools=1.23.1"
    
    tag "Samtools coverage on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("*_coverage.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    samtools coverage $bam > ${prefix}_coverage.tsv
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}_coverage.tsv
    """

}
