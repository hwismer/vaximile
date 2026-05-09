include { ADD_VCF_GT_FIELD; MERGE_SOMATIC_VCFS; FILTER_VCF; POSTPROCESS_VCF } from "../modules/somatic_postprocess.nf"

include { VEP_ANNOTATE; VEP_POPULATION_FILTER; BAMREADCOUNT; ANNOTATE_VCF_TRANSCRIPT_EXPRESSION; ANNOTATE_VCF_GENE_EXPRESSION } from "../modules/pvactools_vcf_prep.nf"

include {
    ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_TUMOR_DNA; 
    ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_NORMAL_DNA;
    ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_TUMOR_RNA
} from "../modules/pvactools_vcf_prep.nf"
include { INDEX_VCF } from "../modules/utilities.nf"

include { 
    PHASE_VCF_SELECT_VARIANTS;
    PHASE_VCF_RENAME;
    PHASE_VCF_COMBINE_VARIANTS;
    PHASE_VCF_SORT_VCF;
    PHASE_VCF_RBPHASING;
    PHASE_VCF_VEP;
    PHASE_VCF_INDEX;
} from "../modules/pvactools_vcf_phasing.nf"

        



workflow PVAC_INPUT_PREP_WORKFLOW {

    take:
        mutect_vcf
        strelka_vcf
        preproc_bams
        star_bam
        kallisto_tx_abundance
        kallisto_gene_abundance
        reference_genome
        vep_cache
        vep_plugins
        germline_vcf

    main:
    
    // Add GT to strelka calls
    strelka_gt = ADD_VCF_GT_FIELD(strelka_vcf)

    strelka = strelka_gt.map{meta, vcf, tbi -> tuple(meta, "strelka", vcf, tbi)}
    mutect = mutect_vcf.map{meta, vcf, tbi -> tuple(meta, "mutect", vcf, tbi)}

    vcfs = mutect.mix(strelka)
    vcfs_filtered = FILTER_VCF(vcfs).filtered_vcf
    vcfs_normalized = POSTPROCESS_VCF(vcfs_filtered, reference_genome)

    callers = vcfs_normalized.branch{ meta, caller, vcf, tbi ->
        mutect: caller == "mutect"
        strelka: caller == "strelka"
    }
    
    merged_callers = callers.mutect.join(callers.strelka)

    merged_vcf = MERGE_SOMATIC_VCFS(merged_callers, reference_genome)
    vep = VEP_ANNOTATE(merged_vcf, reference_genome, vep_cache, vep_plugins)
    vep_filtered = VEP_POPULATION_FILTER(vep.vcf, vep_cache, vep_plugins)

    vep_filtered_somatic_name = vep_filtered.map {meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    all_samples_somatic_name = preproc_bams.mix(star_bam).map {meta, bam, bai -> tuple(meta.somatic_name, meta, bam,bai)}


    bamreadcount_helper_input = vep_filtered_somatic_name.combine(all_samples_somatic_name, by:0)

    brc_helper = BAMREADCOUNT(bamreadcount_helper_input, reference_genome)
    
    brc_helper_branched = brc_helper.branch { somatic_meta, sample_meta, indels, snvs ->
        tumor_dna: sample_meta.sample_type == "Tumor" && sample_meta.molecule == "DNA"
        normal_dna: sample_meta.sample_type == "Normal" && sample_meta.molecule == "DNA"
        tumor_rna: sample_meta.sample_type == "Tumor" && sample_meta.molecule == "RNA"
    }
    
    tumor_dna = ANNOTATE_VCF_COVERAGE_TUMOR_DNA(brc_helper_branched.tumor_dna.join(vep_filtered))
    normal_dna_tdna = ANNOTATE_VCF_COVERAGE_NORMAL_DNA(brc_helper_branched.normal_dna.join(tumor_dna)) 
    tumor_rna_ndna_tdna = ANNOTATE_VCF_COVERAGE_TUMOR_RNA(brc_helper_branched.tumor_rna.join(normal_dna_tdna))

    somatic_name_vcf_coverage = tumor_rna_ndna_tdna.map{meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    somatic_name_tx = kallisto_tx_abundance.map{meta, tx -> tuple(meta.somatic_name, meta, tx) }
    tx_vcf = ANNOTATE_VCF_TRANSCRIPT_EXPRESSION(somatic_name_vcf_coverage.join(somatic_name_tx))

    tx_vcf_somatic_name = tx_vcf.map{ meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    somatic_name_gene = kallisto_gene_abundance.map{ meta, gene -> tuple(meta.somatic_name, meta, gene) }
    gene_vcf = ANNOTATE_VCF_GENE_EXPRESSION(tx_vcf_somatic_name.join(somatic_name_gene))

    final_vcf = INDEX_VCF(gene_vcf) 


    // CREATE PHASED VCF 
    tumor_normal_samples = preproc_bams.branch {meta, bam, bai ->
        tumor:meta.sample_type == "Tumor"
        normal:meta.sample_type == "Normal"
    }
    
    // Phase with germline calls
    phased_vcf = PVAC_VCF_PHASING(
        tumor_normal_samples.tumor,
        final_vcf,
        germline_vcf,
        reference_genome,
        vep_cache,
        vep_plugins
    )

    phased_vcf_final = phased_vcf.final_phased

    emit:
        somatic_vcf = final_vcf
        phased_vcf = phased_vcf_final
}


workflow PVAC_VCF_PHASING {

    take:
        tumor_bam
        somatic_vcf
        germline_vcf
        reference_genome
        vep_cache
        vep_plugins
       
    main:

    select_variants = PHASE_VCF_SELECT_VARIANTS(somatic_vcf, reference_genome)

    meta_germline = somatic_vcf.map{somatic_meta, vcf, tbi ->
        tuple(somatic_meta.normal_meta, somatic_meta)
    }.join(germline_vcf, by: 0)

    germline_renamed = PHASE_VCF_RENAME(meta_germline)


    combine_variants_input = select_variants.join(germline_renamed, by:0)
    combined_variants = PHASE_VCF_COMBINE_VARIANTS(combine_variants_input, reference_genome)

    combined_sorted = PHASE_VCF_SORT_VCF(combined_variants, reference_genome)

    combined_sorted_by_tumor_sample = combined_sorted.map{ meta, vcf ->
        tuple(meta.tumor_meta, meta, vcf)
    }.join(tumor_bam)

    rbphased = PHASE_VCF_RBPHASING(combined_sorted_by_tumor_sample, reference_genome)

    phased_vep = VEP_ANNOTATE(rbphased, reference_genome, vep_cache, vep_plugins)


    final_phased = INDEX_VCF(phased_vep.vcf)

    
    emit:
        final_phased = final_phased

}
