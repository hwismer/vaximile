process CREATE_KALLISTO_INDEX {

    /*

    Create a kallisto index using a reference genome.

    */
    
    label 'process_high'
    cache 'lenient'
    
    publishDir "./resources/kallisto", mode: "copy"

    conda "bioconda::kallisto=0.51.1"

    tag "Creating kallisto index on ${transcriptome_fa}"

    input:
        path(transcriptome_fa)

    output:
        path("*_kallisto_index.idx")

    script:
        def prefix = task.ext.prefix ?: "${transcriptome_fa}"
        """
        kallisto index -i "${prefix}_kallisto_index.idx" $transcriptome_fa
        """

}
