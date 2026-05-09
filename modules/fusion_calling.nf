process STAR_FUSION {

    /*
        NOT IMPLEMENTED YET
        Run star-fusion to detect RNA fusion events.

    */

    cpus 4
    memory "64GB"

    container "trinityctat/starfusion:1.15.0"

    publishDir "${params.outdir}/${meta.somatic_sample}/fusions/", mode: "copy"

    input:
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2)
        path ctat_resource_lib

    output:
        tuple val(meta), path("${meta.sample_name}_starfusion/*.fusion_predictions.tsv"), emit: fusion_preds
        tuple val(meta), path("${meta.sample_name}_starfusion/*.fusion_predictions.abridged.tsv"), emit: abridged_preds
        tuple val(meta), path("${meta.sample_name}_starfusion/*.coding_effect.tsv"), emit: coding_effect
        tuple val(meta), path("${meta.sample_name}_starfusion/"), emit: all_output


    script:
        """
        STAR-Fusion --genome_lib_dir $ctat_resource_lib \
             -J $chimeric_out \
             --examine_coding_effect \
             --FusionInspector validate \
             --left_fq $fastq1 \
             --right_fq $fastq2 \
             --denovo_reconstruct \
             --output_dir "./${meta.sample_name}_starfusion"

        """

}


process ARRIBA_FUSION {

    cpus 4
    memory "32GB"

    conda "bioconda::arriba=2.5.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/fusions/arriba", mode: "copy"

    input:
        tuple val(meta), path(star_bam), path(star_bam_index)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
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
