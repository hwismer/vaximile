process INDEX_BAM {

    // Index a sorted BAM as .bai

    label 'process_low'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "Indexing ${bam.name}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path(bam), path("${bam}.bai"), emit: bam

    script:
        """
        samtools index -@ $task.cpus $bam
        """

    stub:
        """
        touch ${bam}.bai
        """
}
