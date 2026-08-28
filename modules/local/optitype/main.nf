process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    container "fred2/optitype:latest"
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
        
        which python
        which OptiTypePipeline.py
        
        pwd
        ls -lh

        python /usr/local/bin/OptiType/OptiTypePipeline.py \
            -i $fastq1 $fastq2 \
            --$molecule_flag \
            -c OptiType.ini \
            --prefix "${prefix}" \
            --outdir optitype_out
        """
}
