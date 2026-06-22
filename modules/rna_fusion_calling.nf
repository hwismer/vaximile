process STAR_FUSION {

    /*
        Run star-fusion to detect RNA fusion events.

    */

    cpus 12
    memory "90GB"
    container "trinityctat/starfusion:1.15.0"

    tag "Running STARfusion on ${meta.sample_name}"

    input:
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2)
        path ctat_resource_lib

    output:
        tuple val(meta), path("${meta.sample_name}_starfusion/*.fusion_predictions.tsv"), emit: fusion_preds
        tuple val(meta), path("${meta.sample_name}_starfusion/*.fusion_predictions.abridged.tsv"), emit: abridged_preds
        //tuple val(meta), path("${meta.sample_name}_starfusion/*.coding_effect.tsv"), emit: coding_effect
        //tuple val(meta), path("${meta.sample_name}_starfusion/"), emit: all_output


    script:
        """
        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --output_dir "./${meta.sample_name}_starfusion"

        """

}


process ARRIBA_FUSION {

    cpus 6
    memory "350GB"
    conda "bioconda::arriba=2.5.1"

    tag "Running Arriba fusion calling on ${meta.sample_name}"

    input:
        tuple val(meta), path(star_bam), path(star_bam_index)
        tuple path(reference_fa), path(reference_index)
        path gtf
        tuple path(arriba_blacklist), path(arriba_known_fusions), path(arriba_protein_domains)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_arriba_fusions.tsv"), emit: arriba_fusions
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_arriba_fusions.discarded.tsv"), emit: discarded_fusions

    script:
        """
        arriba -x $star_bam \
            -g $gtf \
            -a $reference_fa \
            -b $arriba_blacklist \
            -k $arriba_known_fusions \
            -p $arriba_protein_domains \
            -o "${meta.sample_name}_${meta.molecule}_arriba_fusions.tsv" \
            -O "${meta.sample_name}_${meta.molecule}_arriba_fusions.discarded.tsv"
        """

}
