process SAMTOOLS_IDXSTATS {
    
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Samtools idxstats on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("${meta.sample_name}_${meta.molecule}_idxstats.tsv")

    script:
    """
    samtools idxstats $bam > ${meta.sample_name}_${meta.molecule}_idxstats.tsv
    """
}
