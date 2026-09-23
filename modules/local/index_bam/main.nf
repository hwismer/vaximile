process INDEX_BAM {

    /*
    Index a coordinate-sorted BAM, producing a .bai.

    Used twice: after duplicate marking, and after BQSR gather in place of a re-sort.

    SAMTOOLS_SORMADUP writes only the BAM. Its --write-index option produces a .csi,
    which is not interchangeable here: every downstream module in this pipeline declares
    a path(bai) input, and Strelka and Manta read .bai specifically. So the index is made
    as a separate step rather than by configuring the nf-core module, which also keeps
    that module unpatched and updatable.
    */

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
