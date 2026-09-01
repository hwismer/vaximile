include { BWA_MAP } from "../../../modules/local/bwa_map/main"
include { BAM_MARKDUPLICATES } from "../bam_markduplicates/main"
include { BASE_RECALIBRATOR_SCATTER } from "../../../modules/local/base_recalibrator_scatter/main"
include { BASE_RECALIBRATOR_GATHER } from "../../../modules/local/base_recalibrator_gather/main"
include { APPLY_BQSR_SCATTER } from "../../../modules/local/apply_bqsr_scatter/main"
include { APPLY_BQSR_GATHER } from "../../../modules/local/apply_bqsr_gather/main"
include { GET_PILEUP_SUMMARIES } from "../../../modules/local/get_pileup_summaries/main"
include { SAMTOOLS_FLAGSTAT } from "../../../modules/local/samtools_flagstat/main"
include { SAMTOOLS_COVERAGE } from "../../../modules/local/samtools_coverage/main"
include { SAMTOOLS_IDXSTATS } from "../../../modules/local/samtools_idxstats/main"

workflow DNA_ALIGN_AND_PREPROC {

    take:
        fastqs // (metamap, fastq1, fastq2)
        reference_genome // (reference fasta, reference fasta index)
        reference_dict // reference_dict file
        bwa_index // Either null if not supplied in main workflow or path to bwa index files
        known_sites_dbsnp // (vcf, vcf index)
        known_sites_1000g_snps // (vcf, vcf_index)
        known_indels // (vcf, vcf_index)
        mills // (vcf, vcf_index)
        common_germline // (vcf, vcf_index)
        intervals // (kit name, tuple(interval_files)) interval files created from the SPLIT_INTERVALS process
        num_intervals // Integer for how many intervals are present
        
    main:
        
        // Map with minibwa
        bwa_map_input = fastqs
            .map { meta, fastq1, fastq2 ->
                tuple(meta, meta.sample_name, meta.molecule, meta.sequencing_type, fastq1, fastq2)
            }

        bwa_bam = BWA_MAP(bwa_map_input, reference_genome, bwa_index).bam

        // **********************************************************
        // GATK PRE-PROCESSING BEST PRACTICES

        // Coordinate-sort and mark duplicates (samtools collate/fixmate/sort/markdup,
        // via the nf-core SAMTOOLS_SORMADUP module) in place of MarkDuplicatesSpark.
        markduplicates = BAM_MARKDUPLICATES(bwa_bam, reference_genome)
        mark_dup = markduplicates.bam
        markdup_metrics = markduplicates.metrics

        // BASE RECALIBRATION
        // Call BaseRecalibrator on each interval
        markdup_kit = mark_dup.map{meta, bam, bai ->
            tuple(meta.capture_kit, meta, bam, bai)
        }
        
        intervals_map = intervals.flatMap{ kit, interval_list ->
            interval_list.collect { interval ->
                tuple(kit, interval)
            }
        }
        

        base_recal_input = markdup_kit.combine(intervals_map, by: 0).map{kit, meta, bam, bai, interval -> tuple(meta, bam, bai, interval)}
    
        base_recal = BASE_RECALIBRATOR_SCATTER(
            base_recal_input,
            reference_genome,
            reference_dict,
            known_sites_dbsnp,
            known_sites_1000g_snps,
            known_indels,
            mills
        ).table
        base_recal_gather_input = base_recal.groupTuple(size: num_intervals)
            .map { meta, recal_tables ->
                tuple(meta, meta.sample_name, meta.molecule, recal_tables)
            }
        base_recal_gathered = BASE_RECALIBRATOR_GATHER(base_recal_gather_input).table
        
        // Scatter BQSR calls
        bqsr_input = mark_dup.join(base_recal_gathered)
        .map{meta, bam, bai, bqsr ->
            tuple(meta.capture_kit, meta, bam, bai, bqsr)
        }
        .combine(intervals_map, by:0 )
        .map{kit, meta, bam, bai, bqsr, interval -> tuple(meta, bam, bai, bqsr, interval) }

        bqsr = APPLY_BQSR_SCATTER(bqsr_input, reference_genome, reference_dict).bam
        bqsr_scattered = bqsr.groupTuple(size: num_intervals)
            .map { meta, bams ->
                tuple(meta, meta.sample_name, meta.molecule, bams)
            }
        bqsr_gather = APPLY_BQSR_GATHER(bqsr_scattered).bam // Get final BQSR bams
        // APPLY_BQSR_GATHER now emits (meta, bam, bai) directly: it merges the shards
        // rather than concatenating them, so the result is coordinate-sorted and indexed in
        // one streaming pass and there is no separate sort or index step to do here.
        bqsr_sort = bqsr_gather

        pileup_summaries = GET_PILEUP_SUMMARIES(bqsr_sort, common_germline).table
        
        flagstat = SAMTOOLS_FLAGSTAT(mark_dup).flagstat
        coverage = SAMTOOLS_COVERAGE(mark_dup).tsv
        idxstats = SAMTOOLS_IDXSTATS(mark_dup).tsv

        
    emit:
        preproc_bams = bqsr_sort // For somatic calling
        markdup_bams = mark_dup // Non-recalibrated BAMS for callers like Strelka that don't expect recalibrated scores
        base_recal = base_recal_gathered // Recalibration metrics for MultiQC report
        markdup_metrics = markdup_metrics // Duplicate stats from samtools markdup for MultiQC
        pileup_summaries = pileup_summaries
        flagstat = flagstat
        coverage = coverage
        idxstats = idxstats

}
