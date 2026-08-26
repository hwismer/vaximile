process CREATE_KALLISTO_INDEX {

    /*

    Create a kallisto index using a reference genome.

    */
    
    cpus 8
    memory "32GB"
    cache 'lenient'
    
    publishDir "./resources/kallisto", mode: "copy"

    conda "bioconda::kallisto=0.51.1"

    tag "Creating kallisto index on ${transcriptome_fa}"

    input:
        path(transcriptome_fa)

    output:
        path("${transcriptome_fa}_kallisto_index.idx")

    script:
        """
        kallisto index -i "${transcriptome_fa}_kallisto_index.idx" $transcriptome_fa
        """

}
