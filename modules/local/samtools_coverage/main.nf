process SAMTOOLS_COVERAGE {
    
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
