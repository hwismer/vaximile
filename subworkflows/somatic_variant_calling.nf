// Somatic variant calling: Mutect2, Strelka/Manta and DeepSomatic, and their consensus.
include { ADD_VCF_GT_FIELD } from "../modules/local/add_vcf_gt_field/main"
include { FILTER_VCF } from "../modules/local/filter_vcf/main"
include { POSTPROCESS_VCF } from "../modules/local/postprocess_vcf/main"
include { MERGE_SOMATIC_VCFS } from "../modules/local/merge_somatic_vcfs/main"
include { MUTECT2_SCATTER } from "../modules/local/mutect2_scatter/main"
include { MUTECT2_GATHER_SELECT_VARIANTS } from "../modules/local/mutect2_gather_select_variants/main"
include { MUTECT2_GATHER_VCFS } from "../modules/local/mutect2_gather_vcfs/main"
include { MUTECT2_CALCULATE_CONTAMINATION } from "../modules/local/mutect2_calculate_contamination/main"
include { MUTECT2_LEARN_READ_ORIENTATION } from "../modules/local/mutect2_learn_read_orientation/main"
include { MUTECT2_MERGE_STATS } from "../modules/local/mutect2_merge_stats/main"
include { MUTECT2_FILTER_MUTECT_CALLS } from "../modules/local/mutect2_filter_mutect_calls/main"
include { MANTA } from "../modules/local/manta/main"
include { STRELKA } from "../modules/local/strelka/main"
include { POSTPROCESS_STRELKA } from "../modules/local/postprocess_strelka/main"
include { DEEPSOMATIC } from "../modules/local/deepsomatic/main"

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

workflow STRELKA_WORKFLOW {

    take:
        somatic_pairs // (somatic metamap, tumor_bam, tumor_bai, normal_bam, normal_bai)
        reference_genome // (fasta, fasta.fai, dict)
        strelka_bed

    main:
        
        somatic_pairs_kit = somatic_pairs.map{ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }
        .combine(strelka_bed, by:0)
        .map{ capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bed, bed_tbi ->
            tuple(meta, meta.tumor_meta.sequencing_type, tumor_bam, tumor_bai, normal_bam, normal_bai, bed, bed_tbi)
        }
        
        manta = MANTA(somatic_pairs_kit, reference_genome).dir

        strelka_input = somatic_pairs.join(manta).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir)
        }
        .combine(strelka_bed, by:0)
        .map { capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir, bed, bed_tbi ->
            tuple(meta, meta.tumor_meta.sequencing_type, tumor_bam, tumor_bai, normal_bam, normal_bai, manta_dir, bed, bed_tbi)
        }

        strelka = STRELKA(strelka_input, reference_genome)

        postprocess_strelka_input = strelka.strelka_vcfs
            .map { meta, strelka_snvs, strelka_snvs_index, strelka_indels, strelka_indels_index ->
                tuple(meta, meta.tumor_meta.sample_name, meta.normal_meta.sample_name, strelka_snvs, strelka_snvs_index, strelka_indels, strelka_indels_index)
            }

        strelka_postprocess = POSTPROCESS_STRELKA(postprocess_strelka_input)

    emit:
        strelka_vcf = strelka_postprocess.vcf
}

workflow DEEPSOMATIC_WORKFLOW {

    take:
        somatic_pairs // (somatic metamap, tumor_bam, tumor_bai, normal_bam, normal_bai)
        reference_genome // (fasta, fasta.fai, dict)
        capture_kits // (kit_name, bed file)

    main:

        somatic_pairs_kit = somatic_pairs.map{ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(meta.capture_kit, meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }
        .combine(capture_kits, by:0)
        .map { kit, meta, tb, tbai, nbam, nbai, bed ->
            tuple(meta, meta.tumor_meta.sample_name, meta.normal_meta.sample_name, meta.tumor_meta.sequencing_type, tb, tbai, nbam, nbai, bed)
        }

        deepsomatic = DEEPSOMATIC(somatic_pairs_kit, reference_genome).vcf

    emit:
        deepsomatic_vcf = deepsomatic
}


// 2-of-3 consensus: add GT to Strelka calls, keep PASS, normalise, then CombineVariants --minimumN 2.
workflow SOMATIC_CONSENSUS {

    take:
        mutect_vcf
        strelka_vcf
        deepsomatic_vcf
        reference_genome
        reference_dict

    main:

    // Add GT to strelka calls
    add_vcf_gt_field_input = strelka_vcf
        .map { meta, vcf, tbi ->
            tuple(meta, meta.tumor_meta.sample_name, vcf, tbi)
        }

    strelka_gt = ADD_VCF_GT_FIELD(add_vcf_gt_field_input).vcf

    strelka = strelka_gt.map{meta, vcf -> tuple(meta, "strelka", vcf)}
    mutect = mutect_vcf.map{meta, vcf, _tbi -> tuple(meta, "mutect", vcf)}
    deepsomatic = deepsomatic_vcf.map{meta, vcf, _tbi -> tuple(meta, "deepsomatic", vcf) }

    vcfs = mutect.mix(deepsomatic).mix(strelka)
    vcfs_filtered = FILTER_VCF(vcfs).filtered_vcf
    postprocess_vcf_input = vcfs_filtered
        .map { meta, caller, vcf, tbi ->
            tuple(meta, meta.somatic_name, caller, vcf, tbi)
        }

    vcfs_normalized = POSTPROCESS_VCF(postprocess_vcf_input, reference_genome).vt_vcf
    
    callers = vcfs_normalized.branch{ meta, caller, vcf, tbi ->
        mutect: caller == "mutect"
        strelka: caller == "strelka"
        deepsomatic: caller == "deepsomatic"
    }
    
    merged_callers = callers.mutect.join(callers.deepsomatic).join(callers.strelka)

    merge_somatic_vcfs_input = merged_callers
        .map { meta, vcf1_caller, vcf1, vcf1_index, vcf2_caller, vcf2, vcf2_index, vcf3_caller, vcf3, vcf3_index ->
            tuple(meta, meta.somatic_name, vcf1_caller, vcf1, vcf1_index, vcf2_caller, vcf2, vcf2_index, vcf3_caller, vcf3, vcf3_index)
        }

    merged_vcf = MERGE_SOMATIC_VCFS(merge_somatic_vcfs_input, reference_genome, reference_dict).vcf

    emit:
        vcf = merged_vcf
}
