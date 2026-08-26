process STAR_FUSION {

    /*
        Run star-fusion to detect RNA fusion events.

    */

    label 'process_very_high'
    conda "bioconda::star-fusion=1.15.0"
    container "trinityctat/starfusion:1.15.0"

    tag "Running STARfusion on ${meta.sample_name}"

    input:
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2)
        path ctat_resource_lib

    output:
        tuple val(meta), path("*_starfusion/*.fusion_predictions.tsv"), emit: fusion_preds
        tuple val(meta), path("*_starfusion/*.fusion_predictions.abridged.tsv"), emit: abridged_preds
        //tuple val(meta), path("${meta.sample_name}_starfusion/*.coding_effect.tsv"), emit: coding_effect
        //tuple val(meta), path("${meta.sample_name}_starfusion/"), emit: all_output


    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --output_dir "./${prefix}_starfusion"

        """

}
