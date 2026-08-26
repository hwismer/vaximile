process PVACFUSE {
    
    label 'process_high_memory'
    conda "bioconda::pvactools=7.0.1"
    container "griffithlab/pvactools:7.0.1"

    
    tag "pVACfuse on ${somatic_name}"
    input:
        tuple val(somatic_name), val(arriba_meta), path(arriba_fusions), val(star_meta), path(starfusion_calls), val(hla_meta), path(hla_pvac_input)
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
            \$(head $hla_pvac_input -n 1) \
            all \
            "${prefix}_pvacfuse" \
            --starfusion-file $starfusion_calls \
            $args \
            --iedb-install-directory /opt/iedb \
            --run-reference-proteome-similarity \
            --peptide-fasta $human_ref_peptides \
            -t $task.cpus
        """

}
