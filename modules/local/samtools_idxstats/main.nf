process SAMTOOLS_IDXSTATS {
    
    label 'process_medium'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Samtools idxstats on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("*_idxstats.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    samtools idxstats --threads $task.cpus $bam > ${prefix}_idxstats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}_idxstats.tsv
    """
}
