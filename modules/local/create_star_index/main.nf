process CREATE_STAR_INDEX {

    /*
    
    Use a reference fasta and gtf to create a star index for star 2.7.10

    */

    cpus 24
    memory "64GB"
    cache 'lenient'

    container "alexdobin/star:2.7.10a_alpha_220506"

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


}
