include { HAPLOTYPE_CALLER_SCATTER; HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS; 
    HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS; HAPLOTYPE_CALLER_GATHER_VCFS; 
    HAPLOTYPE_CALLER_FILTER_VARIANTS;  } from "../modules/germline_variant_calling.nf"

include { DEEPVARIANT; STRELKA_GERMLINE } from "../modules/germline_variant_calling.nf"

include { POSTPROCESS_VCF; MERGE_GERMLINE_VCFS; FILTER_VCF; VEP_ANNOTATE } from "../modules/germline_variant_calling.nf"

include { INDEX_VCF; } from "../modules/utilities.nf"

include { VCF_TO_TABLE } from "../modules/vcf_postprocessing.nf"



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
        haplotype_scatter = HAPLOTYPE_CALLER_SCATTER(bams_intervals, reference_genome, reference_dict, interval_padding,)
        cnn_score = HAPLOTYPE_CALLER_CNN_SCORE_VARIANTS(haplotype_scatter, reference_genome, reference_dict, interval_padding)
        select_variants = HAPLOTYPE_CALLER_GATHER_SELECT_VARIANTS(cnn_score)
        gathered_vcf = HAPLOTYPE_CALLER_GATHER_VCFS(select_variants.groupTuple(size: num_intervals))
        filtered_vcf = HAPLOTYPE_CALLER_FILTER_VARIANTS(gathered_vcf, hapmap, mills)
    
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
            tuple(meta, bam, bai, bed, bed_tbi)
        } 
        strelka_germline = STRELKA_GERMLINE(strelka_input, reference_genome)
        
        deepsomatic_input = sample_bams.map { meta, bam, bai ->
            tuple(meta.capture_kit, meta, bam, bai)
        }
        .combine(capture_kits, by:0)
        .map { capture_kit, meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }
        
        deepvariant = DEEPVARIANT(deepsomatic_input, reference_genome)
        deepvariant_vcf = deepvariant.vcf

        vcfs = haplotype_caller_vcf.mix(strelka_germline).mix(deepvariant_vcf)
        filter_vcfs = FILTER_VCF(vcfs)
        vcfs_norm = POSTPROCESS_VCF(filter_vcfs, reference_genome)

        callers = vcfs_norm.branch{ meta, caller, vcf, tbi ->
            strelka: caller == "strelka"
            deepvariant: caller == "deepvariant"
            haplotypecaller: caller == "haplotypecaller"
        }
        
        merged_callers = callers.deepvariant.join(callers.haplotypecaller).join(callers.strelka)
        
        merged_vcf = MERGE_GERMLINE_VCFS(merged_callers, reference_genome, reference_dict)

        merged_vep = VEP_ANNOTATE(merged_vcf, reference_genome, vep_cache, vep_plugins)
        vep_vcf = merged_vep.vcf
        vep_report = merged_vep.report

        
        merged_vcf_input = vep_vcf.map{meta, vcf ->
            tuple(meta.sample_name, meta, vcf)
        }
        merged_vcf_indexed = INDEX_VCF(merged_vcf_input, "germline")

        vcf_table_name = merged_vcf_indexed.map{ meta, vcf, tbi -> tuple(meta, meta.sample_name + "_germline", vcf, tbi) }
        vcf_table = VCF_TO_TABLE(vcf_table_name)

    
    emit:
        germline_vcf = merged_vcf_indexed
        germline_vcf_table = vcf_table
        germline_vep = vep_report

}

