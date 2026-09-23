process STAR_ALIGN {

    // Align RNA reads with STAR, using STAR-Fusion's parameters so Arriba and STAR-Fusion can reuse the output.

    // process_max: STAR takes all available threads; drop to process_very_high if queueing is slow.
    label 'process_max'
    
    conda "bioconda::star=2.7.11b"

    tag "Aligning ${meta.sample_name} with STAR"

    input:
        tuple val(meta), val(sample_name), path(fastq1), path(fastq2)
        path(star_index_dir)
        path(gtf)

    output:
        tuple val(meta), path("*_Aligned.out.bam"), emit: star_bam
        tuple val(meta), path("*_Log.final.out"), emit:final_log
        tuple val(meta), path("*_SJ.out.tab"), emit: sj_out
        tuple val(meta), path("*_Chimeric.out.junction"), path(fastq1), path(fastq2), emit: chimeric_out
        // Removed: `path("*")` globbed the whole work directory, so staged inputs
        // (FASTQs, index dirs, the decompressed GTF) and versions.yml were emitted
        // as results. Nothing consumed it.
    script:
        def prefix = task.ext.prefix ?: "${sample_name}_${meta.molecule}"
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
            --outSAMattrRGline ID:"${sample_name}" SM:"${sample_name}" \
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
            --sjdbGTFfile gencode.gtf

        """

    stub:
        def prefix = task.ext.prefix ?: "${sample_name}_${meta.molecule}"
        """
        touch ${prefix}_Aligned.out.bam
            touch ${prefix}_Log.final.out
        touch ${prefix}_SJ.out.tab
        touch ${prefix}_Chimeric.out.junction
        """

}
