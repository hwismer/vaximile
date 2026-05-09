include { BWA_MAP;  CREATE_BWA_INDEX; } from "../modules/alignment_and_preprocessing.nf"
include { MARK_DUPLICATES_SPARK; BASE_RECALIBRATOR_SCATTER; BASE_RECALIBRATOR_GATHER; APPLY_BQSR_SCATTER; APPLY_BQSR_GATHER; GET_PILEUP_SUMMARIES } from "../modules/alignment_and_preprocessing.nf"
include { SORT_BAM } from "../modules/utilities.nf"

workflow DNA_ALIGN_AND_PREPROC {

    take:
        fastqs // (metamap, fastq1, fastq2)
        reference_genome // (reference fasta, reference fasta index, reference dict)
        bwa_index // Either null if not supplied in main workflow or path to bwa index files
        known_sites_dbsnp // (vcf, vcf index)
        known_sites_1000g_snps // (vcf, vcf_index)
        known_indels // (vcf, vcf_index)
        mills // (vcf, vcf_index)
        common_germline // (vcf, vcf_index)
        intervals // interval files created from the SPLIT_INTERVALS process
        num_intervals // Integer for how many intervals are present
        
    main:

        // If bwq_index is null (not explicitly defined as a param), then generate a bew index using the supplied reference genome
        if ( bwa_index ) {
           bwa_index_ch = Channel.fromPath(bwa_index).collect()
        } else {
            bwa_index_ch = CREATE_BWA_INDEX(reference_genome)
        }
        //bwa_index.view()
        
        // Map with BWA-mem2
        bwa_sam = BWA_MAP(fastqs, reference_genome, bwa_index_ch)


        // **********************************************************
        // GATK PRE-PROCESSING BEST PRACTICES

        // MarkDuplicatesSpark
        mark_dup = MARK_DUPLICATES_SPARK(bwa_sam)
       
        // BASE RECALIBRATION
        // Call BaseRecalibrator on each interval
        base_recal_input = mark_dup.combine(intervals)
        base_recal = BASE_RECALIBRATOR_SCATTER(
            base_recal_input,
            reference_genome,
            known_sites_dbsnp,
            known_sites_1000g_snps,
            known_indels,
            mills,
        )
        base_recal_gather_input = base_recal.groupTuple(size: num_intervals)
        base_recal_gathered = BASE_RECALIBRATOR_GATHER(base_recal_gather_input)
        
        // Scatter BQSR calls
        bqsr_input = mark_dup.join(base_recal_gathered).combine(intervals)
        bqsr = APPLY_BQSR_SCATTER(bqsr_input, reference_genome)
        bqsr_scattered = bqsr.groupTuple(size: num_intervals)
        bqsr_gather = APPLY_BQSR_GATHER(bqsr_scattered) // Get final BQSR bams
        bqsr_sort = SORT_BAM(bqsr_gather)

        pileup_summaries = GET_PILEUP_SUMMARIES(bqsr_gather, common_germline)

    emit:
        preproc_bams = bqsr_sort // For somatic calling
        markdup_bams = mark_dup // Non-recalibrated BAMS for callers like Strelka that don't expect recalibrated scores
        base_recal = base_recal_gathered // Recalibration metrics for MultiQC report
        pileup_summaries = pileup_summaries

}

