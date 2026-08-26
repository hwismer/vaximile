process SAMTOOLS_COVERAGE {
    
    cpus 8
    memory "24GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Samtools coverage on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}_coverage.tsv")

    script:
    """
    samtools coverage $bam > ${meta.sample_name}_${meta.molecule}_coverage.tsv
    """
}
