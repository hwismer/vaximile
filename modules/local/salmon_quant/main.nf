process SALMON_QUANT {

    cpus 12
    memory "32GB"

    conda "bioconda::salmon=1.11.4"

    tag "Running salmon quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(salmon_index)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_salmon_quant"), emit: quant

    script:
    """
    salmon quant \
        -i $salmon_index \
        --libType A \
        -1 $fastq1 \
        -2 $fastq2 \
        --validateMappings \
        -p $task.cpus \
        -o "${meta.sample_name}_${meta.molecule}_salmon_quant"
    """

}
