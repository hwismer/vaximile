process SALMON_QUANT {

    label 'process_very_high'

    conda "bioconda::salmon=1.11.4"

    tag "Running salmon quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(salmon_index)
        path(gtf)

    output:
        tuple val(meta), path("*_salmon_quant"), emit: quant
        tuple val(meta), path("*_salmon_quant/quant.sf"), emit: tsv
        tuple val(meta), path("*.gene_tpm.tsv"), emit: genes_tsv

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    gzip -d -c $gtf > genes.gtf

    salmon quant \
        -i $salmon_index \
        $args \
        -g genes.gtf \
        -1 $fastq1 \
        -2 $fastq2 \
        -p $task.cpus \
        -o "${prefix}_salmon_quant"

    cp "${prefix}_salmon_quant/quant.genes.sf" "${prefix}.gene_tpm.tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    mkdir -p "${prefix}_salmon_quant"
    touch "${prefix}_salmon_quant/quant.sf"
    touch "${prefix}_salmon_quant/quant.genes.sf"
    touch "${prefix}.gene_tpm.tsv"
    """

}
