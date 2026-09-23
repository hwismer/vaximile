/*
    germline variant calling: haplotypecaller_workflow, germline_workflow

    One file per pipeline step. Each workflow keeps the take/emit signature it had
    as its own subworkflow directory, so callers are unchanged.
*/
include { HAPLOTYPE_CALLER_SCATTER } from "../modules/local/haplotype_caller_scatter/main"
include { HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS } from "../modules/local/haplotype_caller_cnn_score_variants/main"
include { HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS } from "../modules/local/haplotype_caller_gather_select_variants/main"
include { HAPLOTYPE_CALLER_GATHER_VCFS } from "../modules/local/haplotype_caller_gather_vcfs/main"
include { HAPLOTYPE_CALLER_FILTER_VARIANTS } from "../modules/local/haplotype_caller_filter_variants/main"
include { DEEPVARIANT } from "../modules/local/deepvariant/main"
include { STRELKA_GERMLINE } from "../modules/local/strelka_germline/main"
include { POSTPROCESS_VCF } from "../modules/local/postprocess_vcf_germline/main"
include { MERGE_GERMLINE_VCFS } from "../modules/local/merge_germline_vcfs/main"
include { FILTER_VCF } from "../modules/local/filter_vcf_germline/main"
include { VEP_ANNOTATE } from "../modules/local/vep_annotate_germline/main"
include { INDEX_VCF } from "../modules/local/index_vcf/main"
include { VCF_TO_TABLE } from "../modules/local/vcf_to_table/main"

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

workflow GERMLINE_WORKFLOW {

    take:
        sample_bams // (sample metamap, bam, bai)
        capture_kits // (kit_name, bed)
        processed_regions // (kit_name, bed, tbi)
        intervals // (capture_kit, [interval_shards])
        num_intervals // int
        interval_padding
        reference_genome // (fasta, fasta.fai)
        reference_dict
        hapmap
        mills
        vep_cache
        vep_plugins

    main:

        haplotype_caller = HAPLOTYPECALLER_WORKFLOW(sample_bams, intervals, num_intervals, interval_padding, reference_genome, reference_dict, hapmap, mills)
        haplotype_caller_vcf = haplotype_caller.vcf
        
        
        strelka_input = sample_bams.map { meta, bam, bai ->
            tuple(meta.capture_kit, meta, bam, bai)
        }
        .combine(processed_regions, by:0)
        .map { capture_kit, meta, bam, bai, bed, bed_tbi ->
            tuple(meta, meta.sequencing_type, bam, bai, bed, bed_tbi)
        }
        strelka_germline = STRELKA_GERMLINE(strelka_input, reference_genome).vcf
        
        deepsomatic_input = sample_bams.map { meta, bam, bai ->
            tuple(meta.capture_kit, meta, bam, bai)
        }
        .combine(capture_kits, by:0)
        .map { capture_kit, meta, bam, bai, bed ->
            tuple(meta, meta.sequencing_type, bam, bai, bed)
        }
        
        deepvariant = DEEPVARIANT(deepsomatic_input, reference_genome)
        deepvariant_vcf = deepvariant.vcf

        vcfs = haplotype_caller_vcf.mix(strelka_germline).mix(deepvariant_vcf)
        filter_vcfs = FILTER_VCF(vcfs).filtered_vcf
        postprocess_vcf_input = filter_vcfs
            .map { meta, caller, vcf, tbi ->
                tuple(meta, meta.sample_name, caller, vcf, tbi)
            }

        vcfs_norm = POSTPROCESS_VCF(postprocess_vcf_input, reference_genome).postproc_vcf

        callers = vcfs_norm.branch{ meta, caller, vcf, tbi ->
            strelka: caller == "strelka"
            deepvariant: caller == "deepvariant"
            haplotypecaller: caller == "haplotypecaller"
        }
        
        merged_callers = callers.deepvariant.join(callers.haplotypecaller).join(callers.strelka)
        
        merged_vcf = MERGE_GERMLINE_VCFS(merged_callers, reference_genome, reference_dict).vcf

        merged_vep = VEP_ANNOTATE(merged_vcf, reference_genome, vep_cache, vep_plugins)
        vep_vcf = merged_vep.vcf
        vep_report = merged_vep.report

        
        merged_vcf_input = vep_vcf.map{meta, vcf ->
            tuple(meta.sample_name, meta, vcf)
        }
        merged_vcf_indexed = INDEX_VCF(merged_vcf_input, "germline").vcf

        vcf_table_name = merged_vcf_indexed.map{ meta, vcf, tbi -> tuple(meta, meta.sample_name + "_germline", vcf, tbi) }
        vcf_table = VCF_TO_TABLE(vcf_table_name).tsv

    
    emit:
        germline_vcf = merged_vcf_indexed
        germline_vcf_table = vcf_table
        germline_vep = vep_report

}
