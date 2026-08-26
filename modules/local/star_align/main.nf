process STAR_ALIGN {

    /*

    Align RNA reads with STAR. Parameters included from star-fusion to be able to use the output of this process
    in a downstream star-fusion or arriba process without having to re-map.

    */

    label 'process_very_high'

    container "alexdobin/star:2.7.10a_alpha_220506"

    tag "Aligning ${meta.sample_name} with STAR"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(star_index_dir)
        path(gtf)

    output:
        tuple val(meta), path("*_Aligned.out.bam"), emit: star_bam
        tuple val(meta), path("*_ReadsPerGene.out.tab"), emit: gene_quant
        tuple val(meta), path("*_Log.final.out"), emit:final_log
        tuple val(meta), path("*_SJ.out.tab"), emit: sj_out
        tuple val(meta), path("*_Chimeric.out.junction"), path(fastq1), path(fastq2), emit: chimeric_out
        tuple val(meta), path("*"), emit: tutto
    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """

        gzip -d -c $gtf > gencode.gtf

        STAR \
            --runThreadN $task.cpus \
            --genomeDir $star_index_dir \
            --readFilesIn $fastq1 $fastq2 \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped None \
            --twopassMode Basic \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --outFilterMultimapNmax 50 \
            --chimOutJunctionFormat 1 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:"${meta.sample_name}" SM:"${meta.sample_name}" \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --chimNonchimScoreDropMin 10 \
            --chimOutType Junctions WithinBAM HardClip \
            --chimScoreDropMax 30 \
            --peOverlapNbasesMin 10 \
            --peOverlapMMp 0.1 \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSplicedMateMapLmin 30 \
            --outFileNamePrefix ./${prefix}_ \
            --quantMode GeneCounts \
            --sjdbGTFfile gencode.gtf

        """

}
