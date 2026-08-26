process SALMON_QUANT {

    label 'process_very_high'

    conda "bioconda::salmon=1.11.4"

    tag "Running salmon quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(salmon_index)

    output:
        tuple val(meta), path("*_salmon_quant"), emit: quant

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    salmon quant \
        -i $salmon_index \
        -1 $fastq1 \
        -2 $fastq2 \
        $args \
        -p $task.cpus \
        -o "${prefix}_salmon_quant"
    """

}
