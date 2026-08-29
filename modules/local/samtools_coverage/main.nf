process SAMTOOLS_COVERAGE {
    
    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Samtools coverage on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("*_coverage.tsv"), emit: tsv
       path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    samtools coverage $bam > ${prefix}_coverage.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}_coverage.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """

}
