process STAR_FUSION {

    // Call RNA fusions with STAR-Fusion.

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


    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --left_fq $fastq1 \
             --right_fq $fastq2 \
             --CPU $task.cpus \
             --output_dir "./${prefix}_starfusion"

        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        mkdir -p ${prefix}_starfusion
        touch ${prefix}_starfusion/star-fusion.fusion_predictions.tsv
        touch ${prefix}_starfusion/star-fusion.fusion_predictions.abridged.tsv
        """

}
