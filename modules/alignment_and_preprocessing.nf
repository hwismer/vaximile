
process BWA_MAP {

    /*
        Map fastq files using BWA. Outputs a sorted BAM file and its index
        Reads groups are created using metadata information and currently are basically just the same name.

    */

    cpus 16
    memory "64GB"
    cache "lenient"

    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    tag "BWA Alignment on ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
        path bwa_index

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.sam")

    script:
    """
    NEW_RG="@RG\\tID:${meta.sample_name}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:${meta.molecule}_${meta.sequencing_type}"

    bwa-mem2 mem -t $task.cpus -R \$NEW_RG $reference_fa $fastq1 $fastq2 > "${meta.sample_name}_${meta.molecule}.sam"

    """
}

process CREATE_BWA_INDEX {

    /*
        Create a bwa-mem2 index for bwa mem mapping.
    */

    cpus 24
    memory "100GB"
    
    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    tag "Creating BWA index for $reference_fa"

    publishDir "./resources/bwa/bwa_mem2_index_${reference_fa}"

    cache 'lenient'

    input:
        tuple path(reference_fa), path(reference_index), path(reference_dict)

    output:
        path "*{.bwt.2bit.64,.sa,.pac,.amb,.ann,.0123}", emit: bwa_index

    script:
        """
        bwa-mem2 index $reference_fa
        """
}


process MARK_DUPLICATES_SPARK {

    cpus 16
    memory "64GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "MarkDuplicatesSpark on ${meta.sample_name}"

    input:
        tuple val(meta), path(aligned_sam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_markdup.bam"), path("${meta.sample_name}_${meta.molecule}_markdup.bam.bai")
    
    script:

    def tmp = task.workDir

    """
    gatk MarkDuplicatesSpark \
        -I $aligned_sam \
        -O "${meta.sample_name}_${meta.molecule}_markdup.bam" \
        --create-output-bam-index \
        --tmp-dir ${tmp} \
        --spark-master local[${task.cpus}] \
        --conf spark.local.dir=${tmp} \
        --conf spark.sql.shuffle.partitions=${task.cpus * 3} \
        --conf spark.executor.memory=${(task.memory.toGiga() * 0.8) as int}g \
        --conf spark.driver.memory=8g
    """

}

process BASE_RECALIBRATOR_SCATTER {

    cpus 2
    memory "12GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "BaseRecalibrator on ${meta.sample_name} ${interval_shard}"

    input:
        tuple val(meta), path(markdup_bam), path(markdup_bam_bai), val(index), path(interval_shard)
        tuple path(reference_fa), path(reference_index), path(reference_dict)
        tuple path(known_sites_dbsnp), path(known_sites_dbsnp_index)
        tuple path(known_sites_1000g_snps), path(known_sites_1000g_snps_index)
        tuple path(known_indels), path(known_indels_index)
        tuple path(mills), path(mills_index)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_${interval_shard}_recal_table.table")
    
    script:
    """
    gatk BaseRecalibrator \
        -I $markdup_bam \
        -O "${meta.sample_name}_${meta.molecule}_${interval_shard}_recal_table.table" \
        -R $reference_fa \
        --known-sites $known_sites_dbsnp \
        --known-sites $known_sites_1000g_snps \
        --known-sites $known_indels \
        --known-sites $mills \
        -L $interval_shard

    """
}

process BASE_RECALIBRATOR_GATHER {
    
    cpus 2
    memory "812GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GatherBQSRReports on ${meta.sample_name}"

    input:
        tuple val(meta), path(recal_tables)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_recal_table.table")
    
    script:
    """
    gatk GatherBQSRReports \
        ${recal_tables.collect { "-I ${it}" }.join(' ')} \
        -O "${meta.sample_name}_${meta.molecule}_recal_table.table"
    """

}



process APPLY_BQSR_SCATTER {
    
    cpus 4
    memory "16GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "ApplyBQSR on ${meta.sample_name} ${interval_index}"

    input:
        tuple val(meta), path(markdup_bam), path(markdup_bam_bai), path(recal_table), val(interval_index), path(interval_shard)
        tuple path(reference_fa), path(reference_index), path(reference_dict)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_${interval_shard}_bqsr.bam")
    script:
    """
    gatk ApplyBQSR \
        -R $reference_fa \
        -I $markdup_bam \
        -L $interval_shard \
        --bqsr-recal-file $recal_table \
        -O "${meta.sample_name}_${meta.molecule}_${interval_shard}_bqsr.bam"
    """

}

process APPLY_BQSR_GATHER {

    cpus 4
    memory "32GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GatherBams on ${meta.sample_name}"

    input:
        tuple val(meta), path(bams)
    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_bqsr.bam"), path("${meta.sample_name}_${meta.molecule}_bqsr.bam.bai")
    
    script:
    def sorted_bams = bams.sort { it.name }
    """
    gatk GatherBamFiles \
        ${bams.collect { "-I ${it}" }.join(' ')} \
        -O "${meta.sample_name}_${meta.molecule}_bqsr.bam"

    gatk BuildBamIndex \
        -I "${meta.sample_name}_${meta.molecule}_bqsr.bam" \
        -O "${meta.sample_name}_${meta.molecule}_bqsr.bam.bai"

    """

}

process GET_PILEUP_SUMMARIES {
    
    cpus 4
    memory "24GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GetPileupSummaries on ${meta.sample_name}"

    input:
        tuple val(meta), path(bqsr_bam), path(bqsr_bai)
        tuple path(common_germline), path(common_germline_index)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_pileups.table")
    
    script:
    """
    gatk GetPileupSummaries \
        -I $bqsr_bam \
        -V $common_germline \
        -L $common_germline \
        -O "${meta.sample_name}_${meta.molecule}_pileups.table"

    """



}

