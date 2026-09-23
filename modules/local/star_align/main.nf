process STAR_ALIGN {

    /*

    Align RNA reads with STAR. Parameters included from star-fusion to be able to use the output of this process
    in a downstream star-fusion or arriba process without having to re-map.

    */

    // process_max, not process_very_high: STAR is the slowest step in the RNA path and
    // --runThreadN takes whatever the tier gives it. Scaling is sub-linear past ~16 threads,
    // so 32 buys well under 2x, and on a busy cluster the wider reservation may cost more in
    // queue time than it saves. Drop back to process_very_high if that trade goes the wrong
    // way. Memory is 96 GB in both tiers, which is what GRCh38 plus two-pass needs.
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
