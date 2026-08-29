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
        path "versions.yml", topic: versions


    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --left_fq $fastq1 \
             --right_fq $fastq2 \
             --CPU $task.cpus \
             --output_dir "./${prefix}_starfusion"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star-fusion: \$(STAR-Fusion --version 2>&1 | grep -Eo '[0-9]+\.[0-9.]+' | head -1)
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        mkdir -p ${prefix}_starfusion
        touch ${prefix}_starfusion/star-fusion.fusion_predictions.tsv
        touch ${prefix}_starfusion/star-fusion.fusion_predictions.abridged.tsv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star-fusion: 1.15.0
        END_VERSIONS
        """

}
