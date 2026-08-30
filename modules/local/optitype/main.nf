process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    // Pinned to the image's newest published tag rather than the floating `latest`,
    // which that repo last pushed in 2018. Deliberately NOT bumped to OptiType 1.5.0:
    // 1.5.0 is a CLI rewrite - a click group (`optitype run`) whose flags replace the
    // OptiType.ini config this script writes - so moving to it is a migration that needs
    // validating against real data, not a version pin. See docs/usage.md.
    container "fred2/optitype:release-v1.3.1"
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
        cat << EOF > OptiType.ini
        [mapping]
        razers3=/usr/local/bin/razers3
        threads=${task.cpus}
        [ilp]
        solver=cbc
        threads=${task.cpus}
        [behavior]
        deletebam=true
        unpaired_weight=0
        use_discordant=false
        EOF

        # No `which` probes here: the script runs under `set -e`, and
        # `which OptiTypePipeline.py` exits non-zero in this image (the tool lives at an
        # absolute path, not on PATH), which would abort the task before OptiType ran.
        python /usr/local/bin/OptiType/OptiTypePipeline.py \
            -i $fastq1 $fastq2 \
            --$molecule_flag \
            -c OptiType.ini \
            --prefix "${prefix}" \
            --outdir optitype_out
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            optitype: 1.3.1
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
            optitype: 1.3.1
        END_VERSIONS
        """
}
