include { DEEPVARIANT } from "../../../modules/local/deepvariant/main"
include { STRELKA_GERMLINE } from "../../../modules/local/strelka_germline/main"
include { POSTPROCESS_VCF } from "../../../modules/local/postprocess_vcf_germline/main"
include { MERGE_GERMLINE_VCFS } from "../../../modules/local/merge_germline_vcfs/main"
include { FILTER_VCF } from "../../../modules/local/filter_vcf_germline/main"
include { VEP_ANNOTATE } from "../../../modules/local/vep_annotate_germline/main"
include { INDEX_VCF } from "../../../modules/local/index_vcf/main"
include { VCF_TO_TABLE } from "../../../modules/local/vcf_to_table/main"

include { HAPLOTYPECALLER_WORKFLOW } from "../haplotypecaller_workflow/main"

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
