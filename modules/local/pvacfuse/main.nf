process PVACFUSE {
    
    label 'process_high_memory'
    // Container only: pVACtools runs the IEDB predictors from /opt/iedb, which only the image has.
    container "griffithlab/pvactools:7.1.5"

    
    tag "pVACfuse on ${somatic_name}"
    input:
        tuple val(somatic_name), val(arriba_meta), path(arriba_fusions), val(star_meta), path(starfusion_calls), val(hla_meta), val(hla_alleles)
        path human_ref_peptides

    output:
        tuple val(star_meta), path("*_pvacfuse"), emit: pvacfuse_dir

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        pvacfuse run \
            $arriba_fusions \
            "${somatic_name}" \
            "${hla_alleles}" \
            all \
            "${prefix}_pvacfuse" \
            --starfusion-file $starfusion_calls \
            $args \
            --iedb-install-directory /opt/iedb \
            --run-reference-proteome-similarity \
            --peptide-fasta $human_ref_peptides \
            -t $task.cpus
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_name}"
        """
        mkdir -p ${prefix}_pvacfuse
        """

}
