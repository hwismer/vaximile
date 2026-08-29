include { HAPLOTYPE_CALLER_SCATTER } from "../../../modules/local/haplotype_caller_scatter/main"
include { HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS } from "../../../modules/local/haplotype_caller_cnn_score_variants/main"
include { HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS } from "../../../modules/local/haplotype_caller_gather_select_variants/main"
include { HAPLOTYPE_CALLER_GATHER_VCFS } from "../../../modules/local/haplotype_caller_gather_vcfs/main"
include { HAPLOTYPE_CALLER_FILTER_VARIANTS } from "../../../modules/local/haplotype_caller_filter_variants/main"

workflow HAPLOTYPECALLER_WORKFLOW {

    take:
        sample_bams
        intervals
        num_intervals
        interval_padding
        reference_genome
        reference_dict
        hapmap
        mills

    main:
        
        intervals_map = intervals.flatMap{ kit, interval_list ->
            interval_list.collect { interval ->
                tuple(kit, interval)
            }
        }

        sample_bams_kit = sample_bams.map{meta, bam, bai ->
            tuple(meta.capture_kit, meta, bam, bai)
        }

        bams_intervals = sample_bams_kit.combine(intervals_map, by:0).map{ kit, meta, bam, bai, interval -> tuple(meta, bam, bai, interval) }
        haplotype_scatter = HAPLOTYPE_CALLER_SCATTER(bams_intervals, reference_genome, reference_dict, interval_padding,).vcf
        cnn_score = HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS(haplotype_scatter, reference_genome, reference_dict, interval_padding).vcf
        select_variants_input = cnn_score
            .map { meta, vcf, vcf_index, interval_shard ->
                tuple(meta, meta.sample_name, vcf, vcf_index, interval_shard)
            }

        select_variants = HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS(select_variants_input).vcf
        gathered_vcf = HAPLOTYPE_CALLER_GATHER_VCFS(select_variants.groupTuple(size: num_intervals)).vcf
        filtered_vcf = HAPLOTYPE_CALLER_FILTER_VARIANTS(gathered_vcf, hapmap, mills).germline_vcf
    
    emit:
        vcf = filtered_vcf

}
