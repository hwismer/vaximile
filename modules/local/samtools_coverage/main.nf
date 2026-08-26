process SAMTOOLS_COVERAGE {
    
    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Samtools coverage on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("*_coverage.tsv")

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    samtools coverage $bam > ${prefix}_coverage.tsv
    """
}
