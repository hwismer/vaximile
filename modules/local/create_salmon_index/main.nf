process CREATE_SALMON_INDEX {

    label 'process_high'
    
    conda "bioconda::salmon=1.11.4"
    
    publishDir "./resources/salmon/", mode: "copy"

    tag "Creating salmon index on ${transcripts_fa}"

    input:
        path(transcripts_fa)

    output:
        path("*_salmon_index"), emit: dir
        path "versions.yml", topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${transcripts_fa}"
    """
    salmon index \
        --threads $task.cpus \
        -t $transcripts_fa \
        -i "${prefix}_salmon_index" \
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version 2>&1 | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${transcripts_fa}"
    """
    mkdir -p "${prefix}_salmon_index"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: 1.11.4
    END_VERSIONS
    """


}
