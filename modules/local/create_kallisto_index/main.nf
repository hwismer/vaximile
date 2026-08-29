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
        path("*_kallisto_index.idx"), emit: index
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${transcriptome_fa}"
        """
        kallisto index -i "${prefix}_kallisto_index.idx" $transcriptome_fa
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: \$(kallisto version 2>&1 | sed 's/kallisto, version //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${transcriptome_fa}"
        """
        touch ${prefix}_kallisto_index.idx
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallisto: 0.51.1
        END_VERSIONS
        """

}
