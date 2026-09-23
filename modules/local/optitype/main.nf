process OPTITYPE {

    // HLA class I typing with OptiType. Memory intensive.

    conda "bioconda::optitype=1.5.0"
    container "quay.io/biocontainers/optitype:1.5.0--pyhdfd78af_1"
    label 'process_high'

    tag "Optitype calls for ${meta.sample_name}"

    input:
        tuple val(meta), val(molecule), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("optitype_out/*_result.tsv"), path("optitype_out/*_coverage_plot.pdf"), emit: hla_calls

    script:

        def molecule_flag = molecule.toLowerCase()
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${molecule}"
        """
        optitype run \\
            -i $fastq1 \\
            -i $fastq2 \\
            --$molecule_flag \\
            --prefix "${prefix}" \\
            --outdir optitype_out \\
            --solver cbc \\
            --threads ${task.cpus} \\
            --ilp-threads ${task.cpus}
        """

    stub:

        def prefix = task.ext.prefix ?: "${meta.sample_name}_${molecule}"
        """
        mkdir -p optitype_out
        touch optitype_out/${prefix}_result.tsv
        touch optitype_out/${prefix}_coverage_plot.pdf
        """
}
