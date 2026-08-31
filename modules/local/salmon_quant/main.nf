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
        // quant.sf on its own as well as the directory. GET_RNA_STRANDEDNESS and MultiQC
        // want the whole run directory; SALMON_TXIMPORT and the VCF expression annotators
        // want the transcript table, and vcf-expression-annotator takes a file, not a dir.
        tuple val(meta), path("*_salmon_quant/quant.sf"), emit: tsv
        // Gene-level TPM, aggregated by salmon itself from --geneMap rather than by a
        // separate tximport step. Same columns as quant.sf, with gene IDs in Name.
        //
        // Copied out under the sample name because this one is published: salmon calls it
        // quant.genes.sf for every sample, so publishing it in place would either collide
        // or drag the run directory along to disambiguate.
        tuple val(meta), path("*.gene_tpm.tsv"), emit: genes_tsv
        path "versions.yml", topic: versions

    script:
    // $args MUST stay ahead of -1/-2. It carries --libType, and salmon rejects a
    // library type that appears after the read files:
    //   "The (--libType/-l) option must precede the input files."
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    # salmon reads the map itself; it will not take the .gz, so decompress as STAR does.
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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version 2>&1 | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    mkdir -p "${prefix}_salmon_quant"
    touch "${prefix}_salmon_quant/quant.sf"
    touch "${prefix}_salmon_quant/quant.genes.sf"
    touch "${prefix}.gene_tpm.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: 1.11.4
    END_VERSIONS
    """

}
