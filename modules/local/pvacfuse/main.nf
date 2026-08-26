process PVACFUSE {
    
    cpus 8
    memory "64GB"
    container "griffithlab/pvactools:7.0.1"

    
    tag "pVACfuse on ${somatic_name}"
    input:
        tuple val(somatic_name), val(arriba_meta), path(arriba_fusions), val(star_meta), path(starfusion_calls), val(hla_meta), path(hla_pvac_input)
        path human_ref_peptides

    output:
        tuple val(star_meta), path("${somatic_name}_pvacfuse"), emit: pvacfuse_dir
    
    script:
        """
        pvacfuse run \
            $arriba_fusions \
            "${somatic_name}" \
            \$(head $hla_pvac_input -n 1) \
            all \
            "${somatic_name}_pvacfuse" \
            --starfusion-file $starfusion_calls \
            -e1 8,9,10,11 \
            -e2 12,13,14,15,16,17,18 \
            --iedb-install-directory /opt/iedb \
            --run-reference-proteome-similarity \
            --peptide-fasta $human_ref_peptides \
            -m median \
            --top-score-metric2 'combined_percentile','ic50' \
            -t $task.cpus
        """

}
