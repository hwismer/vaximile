process CREATE_SALMON_INDEX {

    cpus 8
    memory "32GB"
    
    conda "bioconda::salmon=1.11.4"
    
    publishDir "./resources/salmon/", mode: "copy"

    tag "Creating salmon index on ${transcripts_fa}"

    input:
        path(transcripts_fa)

    output:
        path("${transcripts_fa}_salmon_index")

    script:
    """
    salmon index \
        -t $transcripts_fa \
        -i "${transcripts_fa}_salmon_index" \
        -k 31
    """


}
