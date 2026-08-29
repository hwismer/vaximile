process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    // Pinned to OptiType 1.5.0. The previous image, fred2/optitype:latest, was a floating
    // tag last pushed in 2018 (that repo's newest tag is release-v1.3.1), so runs were
    // neither reproducible nor current. The biocontainer puts OptiTypePipeline.py and
    // razers3 on PATH rather than under /usr/local/bin/OptiType, hence the script change.
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
        cat << EOF > OptiType.ini
        [mapping]
        razers3=\$(which razers3)
        threads=${task.cpus}
        [ilp]
        solver=cbc
        threads=${task.cpus}
        [behavior]
        deletebam=true
        unpaired_weight=0
        use_discordant=false
        EOF

        OptiTypePipeline.py \
            -i $fastq1 $fastq2 \
            --$molecule_flag \
            -c OptiType.ini \
            --prefix "${prefix}" \
            --outdir optitype_out
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
