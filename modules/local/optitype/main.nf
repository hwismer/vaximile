process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive.
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */

    // OptiType 1.5.0, pinned. Both a conda spec and the matching biocontainer, so this
    // module runs under either provisioning method.
    //
    // 1.5.0 is a CLI rewrite, not a drop-in bump. The entry point is a click group -
    // `optitype run`, not `OptiTypePipeline.py` - and the OptiType.ini file this module
    // used to write is gone: mapper, solver and thread counts are now flags. The previous
    // image was `fred2/optitype:latest`, a floating tag last pushed in 2018.
    conda "bioconda::optitype=1.5.0"
    container "quay.io/biocontainers/optitype:1.5.0--pyhdfd78af_1"
    label 'process_high'

    tag "Optitype calls for ${meta.sample_name}"

    input:
        tuple val(meta), val(molecule), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("optitype_out/*_result.tsv"), path("optitype_out/*_coverage_plot.pdf"), emit: hla_calls
        path "versions.yml", topic: versions

    script:

        def molecule_flag = molecule.toLowerCase()
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${molecule}"
        """
        # -i is `multiple=True` in 1.5.0, so each read file needs its own flag; passing
        # both after a single -i silently drops the second.
        #
        # --prefix matters beyond naming: without it OptiType invents a timestamp prefix
        # AND nests results in <outdir>/<timestamp>/, which would not match the output
        # globs above. With it, results land directly in optitype_out/ as
        # <prefix>_result.tsv and <prefix>_coverage_plot.pdf.
        optitype run \\
            -i $fastq1 \\
            -i $fastq2 \\
            --$molecule_flag \\
            --prefix "${prefix}" \\
            --outdir optitype_out \\
            --solver cbc \\
            --threads ${task.cpus} \\
            --ilp-threads ${task.cpus}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            optitype: 1.5.0
        END_VERSIONS
        """

    stub:

        def prefix = task.ext.prefix ?: "${meta.sample_name}_${molecule}"
        """
        mkdir -p optitype_out
        touch optitype_out/${prefix}_result.tsv
        touch optitype_out/${prefix}_coverage_plot.pdf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            optitype: 1.5.0
        END_VERSIONS
        """
}
