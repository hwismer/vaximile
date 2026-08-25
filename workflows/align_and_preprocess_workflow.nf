include { BWA_MAP;  CREATE_BWA_INDEX;  } from "../modules/local/dna_processing.nf"
include { PREPARE_FASTA; INDEX_FASTA; MAKE_FASTA_DICT } from "../modules/local/utilities.nf"
include { MARK_DUPLICATES_SPARK; BASE_RECALIBRATOR_SCATTER; BASE_RECALIBRATOR_GATHER; APPLY_BQSR_SCATTER; APPLY_BQSR_GATHER; GET_PILEUP_SUMMARIES } from "../modules/local/dna_processing.nf"
include { SORT_BAM } from "../modules/local/utilities.nf"
include { SAMTOOLS_FLAGSTAT; SAMTOOLS_COVERAGE; SAMTOOLS_IDXSTATS } from "../modules/local/quality_control.nf"


workflow BWA_INDEX {
    
    take:
        reference_genome
        bwa_index

    main:
        if ( bwa_index ) {
           bwa_index_ch = Channel.fromPath("${bwa_index}/*").collect()
        } else {
            bwa_index_ch = CREATE_BWA_INDEX(reference_genome)
        }


    emit:
        bwa_index = bwa_index_ch

}


workflow PREPARE_REFERENCE_FASTA {
    
    take:
        fasta
    main:
        
        fasta_proc = PREPARE_FASTA(fasta)
        fasta_plus_fai = INDEX_FASTA(fasta_proc).fai
        dict = MAKE_FASTA_DICT(fasta_plus_fai).dict

    emit:
        fa_fai_pair = fasta_plus_fai
        dict  = dict
        
}


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
        
        // Map with BWA-mem2
        bwa_sam = BWA_MAP(fastqs, reference_genome, bwa_index)

        // **********************************************************
        // GATK PRE-PROCESSING BEST PRACTICES

        // MarkDuplicatesSpark
        mark_dup = MARK_DUPLICATES_SPARK(bwa_sam)
       
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
        )
        base_recal_gather_input = base_recal.groupTuple(size: num_intervals)
        base_recal_gathered = BASE_RECALIBRATOR_GATHER(base_recal_gather_input)
        
        // Scatter BQSR calls
        bqsr_input = mark_dup.join(base_recal_gathered)
        .map{meta, bam, bai, bqsr ->
            tuple(meta.capture_kit, meta, bam, bai, bqsr)
        }
        .combine(intervals_map, by:0 )
        .map{kit, meta, bam, bai, bqsr, interval -> tuple(meta, bam, bai, bqsr, interval) }

        bqsr = APPLY_BQSR_SCATTER(bqsr_input, reference_genome, reference_dict)
        bqsr_scattered = bqsr.groupTuple(size: num_intervals)
        bqsr_gather = APPLY_BQSR_GATHER(bqsr_scattered) // Get final BQSR bams
        bqsr_sort = SORT_BAM(bqsr_gather)

        pileup_summaries = GET_PILEUP_SUMMARIES(bqsr_sort, common_germline)
        
        flagstat = SAMTOOLS_FLAGSTAT(mark_dup)
        coverage = SAMTOOLS_COVERAGE(mark_dup)
        idxstats = SAMTOOLS_IDXSTATS(mark_dup)

        
    emit:
        preproc_bams = bqsr_sort // For somatic calling
        markdup_bams = mark_dup // Non-recalibrated BAMS for callers like Strelka that don't expect recalibrated scores
        base_recal = base_recal_gathered // Recalibration metrics for MultiQC report
        pileup_summaries = pileup_summaries
        flagstat = flagstat
        coverage = coverage
        idxstats = idxstats

}

