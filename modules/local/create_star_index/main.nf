process CREATE_STAR_INDEX {

    // Build a STAR index from the reference FASTA and GTF.

    label 'process_max'
    cache 'lenient'

    conda "bioconda::star=2.7.11b"

    tag "Creating STAR index with ${reference_fa} and ${gtf}"

    publishDir "./resources/star/", mode: "copy"

    input:
        tuple path(reference_fa), path(reference_index)
        path(gtf)

    output:
        path("STARGenomeDir"), type: "dir", emit: star_index

    script:
        """
        gzip -d -c $gtf > gtf.gtf

        STAR \
            --runThreadN $task.cpus \
            --runMode genomeGenerate \
            --genomeDir STARGenomeDir \
            --genomeFastaFiles $reference_fa \
            --sjdbGTFfile gtf.gtf

        """

    stub:
        """
        mkdir -p STARGenomeDir
        """


}
