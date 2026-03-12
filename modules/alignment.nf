process STAR_ALIGN {

    /*

    Align RNA reads with STAR. Parameters included from star-fusion to be able to use the output of this process
    in a downstream star-fusion or arriba process without having to re-map.

    */

    cpus 32

    memory "80GB"

    container "alexdobin/star:2.7.10a_alpha_220506"

    //publishDir "${params.outdir}/alignment/star_raw/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(star_index_dir)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Aligned.out.bam"), emit: star_bam 
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Log.final.out"), emit:final_log
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_SJ.out.tab"), emit: sj_out
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Chimeric.out.junction"),path(fastq1), path(fastq2), emit: chimeric_out
    script:
        """
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
            --outFileNamePrefix ./${meta.sample_name}_${meta.molecule}_

        """

}


process CREATE_STAR_INDEX {

    /*
    
    Use a reference fasta and gtf to create a star index for star 2.7.10

    */

    cpus 32
    memory "64GB"
    cache 'lenient'

    container "alexdobin/star:2.7.10a_alpha_220506"

    input:
        path(reference_fa)
        path(reference_index_files)
        path(gtf)

    output:
        path("./STARGenomeDir"), emit: star_index

    script:
        """
        gzip -d -c $gtf > gencode.gtf

        STAR \
            --runThreadN $task.cpus \
            --runMode genomeGenerate \
            --genomeDir ./STARGenomeDir \
            --genomeFastaFiles $reference_fa \
            --sjdbGTFfile gencode.gtf

        """


}

process STAR_INDEX_BAM {

    /*

        Index the BAM file from a star process.

    */

    cpus 8
    memory "32GB"

    container "biocontainers/samtools:v1.9-4-deb_cv1"

    publishDir "${params.outdir}/${meta.somatic_sample}/alignment/star/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(bam)
        tuple val(meta), path(final_log)
        tuple val(meta), path(sj_out)
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam.bai"), emit: star_bam
        tuple val(meta), path(final_log), emit: final_log
        tuple val(meta), path(sj_out), emit: sj_out
        tuple val(meta), path(chimeric_out), path(fastq1), path(fastq2), emit: chimeric_out

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"

        """

}

process BWA_MAP {

    /*
        Map fastq files using BWA. Outputs a sorted BAM file and its index
        Reads groups are created using metadata information and currently are basically just the same name.

    */

    cpus 16
    memory "32GB"
    cache "lenient"

    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path reference_fa
        path reference_index_dir
        path bwa_index

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.sam")

    script:
    """
    NEW_RG="@RG\\tID:${meta.sample_name}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:ILLUMINA"

    bwa-mem2 mem -t $task.cpus -R \$NEW_RG $reference_fa ${fastq1} ${fastq2} > "${meta.sample_name}_${meta.molecule}.sam"

    """
}

process BWA_POSTPROCESS {

    cpus 32
    memory "64GB"
    cache "lenient"

    conda "bioconda::samtools"

    publishDir "${params.outdir}/${meta.somatic_sample}/alignment/bwa/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), val(bwa_sam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai"), emit: mapped_bam

    script:
    """

    samtools sort --threads $task.cpus $bwa_sam -o "${meta.sample_name}_${meta.molecule}_sorted.bam"
    samtools index "${meta.sample_name}_${meta.molecule}_sorted.bam"

    """
}

process CREATE_BWA_INDEX {

    /*
        Create a bwa-mem2 index for bwa mem mapping.
    */

    cpus 16
    memory "32GB"
    
    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    cache 'lenient'

    input:
        path reference_fa
        path reference_fa_index

    output:
        path "*{.bwt.2bit.64,.sa,.pac,.amb,.ann,.0123}", emit: bwa_index

    script:
        """
        bwa-mem2 index $reference_fa
        """
}

process PREPROCESS_BAM {

    /*

    IN THE FUTURE THESE STEPS COULD BE SPLIT UP AND PARALLELIZED.

    Pre-process mapped BAM files according to GATK best practices.
    Currently this calls:
        MarkDuplicatesSpark
        BaseRecalibrator
        ApplyBQSR
        GetPileupSummaries

    */



    cpus 16
    memory "48GB"
    cache "lenient"

    container "broadinstitute/gatk:4.6.1.0"

    publishDir "${params.outdir}/${meta.somatic_sample}/preprocess_bam/${meta.sample_name}_${meta.molecule}", mode: "copy"

    input:
        tuple val(meta), path(reads), path(reads_index)
        path reference_fa
        path reference_fa_index
        path known_sites_dbsnp
        path known_sites_dbsnp_index
        path known_sites_1000g_snps
        path known_sites_1000g_snps_index
        path known_indels
        path known_indels_index
        path mills
        path mills_index
        path common_germline
        path common_germline_index


    output:
        tuple val(meta),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bam"),
            path("${meta.sample_name}_${meta.molecule}_bqsr.bai"),
        emit: preproc_bams

        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_pileups.table"), emit: preproc_bams_pileups


    script:
        """
        gatk MarkDuplicatesSpark \
            -I $reads \
            -O "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            --tmp-dir "\${PWD}"

        gatk BaseRecalibrator \
            -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            -O "${meta.sample_name}_${meta.molecule}_recal_table.table" \
            -R $reference_fa \
            --known-sites $known_sites_dbsnp \
            --known-sites $known_sites_1000g_snps \
            --known-sites $known_indels \
            --known-sites $mills

        gatk ApplyBQSR \
            -R $reference_fa \
            -I "${meta.sample_name}_${meta.molecule}_dedup.bam" \
            --bqsr-recal-file "${meta.sample_name}_${meta.molecule}_recal_table.table" \
            -O "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
            --create-output-bam-index

        gatk GetPileupSummaries \
            -I "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
            -V "${common_germline}" \
            -L "${common_germline}" \
            -O "${meta.sample_name}_${meta.molecule}_pileups.table"
        """
}


