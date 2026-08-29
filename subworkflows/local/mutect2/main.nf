include { MUTECT2_SCATTER } from "../../../modules/local/mutect2_scatter/main"
include { MUTECT2_GATHER_SELECT_VARIANTS } from "../../../modules/local/mutect2_gather_select_variants/main"
include { MUTECT2_GATHER_VCFS } from "../../../modules/local/mutect2_gather_vcfs/main"
include { MUTECT2_CALCULATE_CONTAMINATION } from "../../../modules/local/mutect2_calculate_contamination/main"
include { MUTECT2_LEARN_READ_ORIENTATION } from "../../../modules/local/mutect2_learn_read_orientation/main"
include { MUTECT2_MERGE_STATS } from "../../../modules/local/mutect2_merge_stats/main"
include { MUTECT2_FILTER_MUTECT_CALLS } from "../../../modules/local/mutect2_filter_mutect_calls/main"

workflow MUTECT2 {

    take:
        somatic_pairs // (somatic metamap, tumor_bam, tumor_bai, normal_bam, normal_bai)
        somatic_pileups // (somatic metamap , tumor pileups, normal pileups)
        intervals // (shard name, interval_shard)
        num_intervals // int
        interval_padding
        reference_genome // (fasta, fasta.fai)
        reference_dict
        gnomad // (vcf, tbi)
        pon // (vcf, tbi)

    main:
       
        intervals_map = intervals.flatMap{ kit, interval_list ->
            interval_list.collect { interval ->
                tuple(kit, interval)
            }
        }

        somatic_pairs_kit = somatic_pairs.map{ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }


        somatic_pair_interval = somatic_pairs_kit.combine(intervals_map, by:0).map{kit, meta, tb, tbai, nb, nbai, interval -> tuple(meta, meta.normal_meta.sample_name, tb, tbai, nb, nbai, interval) }
        mutect2_scatter = MUTECT2_SCATTER(somatic_pair_interval, reference_genome, reference_dict, gnomad, pon, interval_padding)

        mutect_vcfs = mutect2_scatter.vcf
        mutect_f1r2s = mutect2_scatter.f1r2
        mutect_stats = mutect2_scatter.stats

        gather_select_variants_input = mutect_vcfs
            .map { somatic_meta, vcf, vcf_index, interval_shard ->
                tuple(somatic_meta, somatic_meta.somatic_name, vcf, vcf_index, interval_shard)
            }

        select_variants = MUTECT2_GATHER_SELECT_VARIANTS(gather_select_variants_input).vcf
        select_variants_grouped = select_variants.groupTuple(size: num_intervals)
        gather_vcfs = MUTECT2_GATHER_VCFS(select_variants_grouped).vcf

        contamination = MUTECT2_CALCULATE_CONTAMINATION(somatic_pileups).table
        read_orientation = MUTECT2_LEARN_READ_ORIENTATION(mutect_f1r2s.groupTuple(size: num_intervals)).tar
        stats = MUTECT2_MERGE_STATS(mutect_stats.groupTuple(size: num_intervals)).stats

        mutect2_filtering_input = gather_vcfs.join(read_orientation).join(stats).join(contamination)
        filtered_calls = MUTECT2_FILTER_MUTECT_CALLS(mutect2_filtering_input, reference_genome, reference_dict)

    emit:
        mutect2_vcf = filtered_calls.filtered_vcf

}
