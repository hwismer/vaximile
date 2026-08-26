include { VEP_ANNOTATE } from "../../../modules/local/vep_annotate/main"
include { INDEX_VCF } from "../../../modules/local/index_vcf/main"
include { PHASE_VCF_SELECT_VARIANTS } from "../../../modules/local/phase_vcf_select_variants/main"
include { PHASE_VCF_RENAME } from "../../../modules/local/phase_vcf_rename/main"
include { PHASE_VCF_COMBINE_VARIANTS } from "../../../modules/local/phase_vcf_combine_variants/main"
include { PHASE_VCF_SORT_VCF } from "../../../modules/local/phase_vcf_sort_vcf/main"
include { PHASE_VCF_RBPHASING } from "../../../modules/local/phase_vcf_rbphasing/main"

workflow PVAC_VCF_PHASING {

    take:
        tumor_bam
        somatic_vcf
        germline_vcf
        reference_genome
        reference_dict
        vep_cache
        vep_plugins
       
    main:

    select_variants = PHASE_VCF_SELECT_VARIANTS(somatic_vcf, reference_genome, reference_dict)
    

    somatic_meta_germline = somatic_vcf.map{meta, vcf, tbi ->
        tuple(meta.normal_meta, meta)
    }.join(germline_vcf, by:0)

    germline_renamed = PHASE_VCF_RENAME(somatic_meta_germline)


    combine_variants_input = select_variants.join(germline_renamed, by:0)
    combined_variants = PHASE_VCF_COMBINE_VARIANTS(combine_variants_input, reference_genome, reference_dict)

    combined_sorted = PHASE_VCF_SORT_VCF(combined_variants, reference_genome, reference_dict)

    combined_sorted_by_tumor_sample = combined_sorted.map{ meta, vcf ->
        tuple(meta.tumor_meta, meta, vcf)
    }.join(tumor_bam)

    rbphased = PHASE_VCF_RBPHASING(combined_sorted_by_tumor_sample, reference_genome, reference_dict)

    phased_vep = VEP_ANNOTATE(rbphased, reference_genome, vep_cache, vep_plugins)

    
    index_input = phased_vep.vcf.map{meta, vcf ->
        tuple(meta.tumor_meta.sample_name, meta, vcf)
    }
    final_phased = INDEX_VCF(index_input, "phased")

    
    emit:
        final_phased = final_phased

}
