include { ADD_VCF_GT_FIELD } from "../../../modules/local/add_vcf_gt_field/main"
include { MERGE_SOMATIC_VCFS } from "../../../modules/local/merge_somatic_vcfs/main"
include { FILTER_VCF } from "../../../modules/local/filter_vcf/main"
include { POSTPROCESS_VCF } from "../../../modules/local/postprocess_vcf/main"
include { VCF_TO_TABLE } from "../../../modules/local/vcf_to_table/main"
include { VEP_ANNOTATE } from "../../../modules/local/vep_annotate/main"
include { VEP_POPULATION_FILTER } from "../../../modules/local/vep_population_filter/main"
include { BAMREADCOUNT } from "../../../modules/local/bamreadcount/main"
include { ANNOTATE_VCF_TRANSCRIPT_EXPRESSION } from "../../../modules/local/annotate_vcf_transcript_expression/main"
include { ANNOTATE_VCF_GENE_EXPRESSION } from "../../../modules/local/annotate_vcf_gene_expression/main"
include { ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_TUMOR_DNA } from "../../../modules/local/annotate_vcf_coverage/main"
include { ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_NORMAL_DNA } from "../../../modules/local/annotate_vcf_coverage/main"
include { ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_TUMOR_RNA } from "../../../modules/local/annotate_vcf_coverage/main"
include { INDEX_VCF } from "../../../modules/local/index_vcf/main"

include { PVAC_VCF_PHASING } from "../pvac_vcf_phasing/main"

workflow PVAC_INPUT_PREP_WORKFLOW {

    take:
        mutect_vcf
        strelka_vcf
        deepsomatic_vcf
        preproc_bams
        star_bam
        kallisto_tx_abundance
        kallisto_gene_abundance
        reference_genome
        reference_dict
        vep_cache
        vep_plugins
        germline_vcf

    main:
    
    // Add GT to strelka calls
    strelka_gt = ADD_VCF_GT_FIELD(strelka_vcf)

    strelka = strelka_gt.map{meta, vcf, tbi -> tuple(meta, "strelka", vcf, tbi)}
    mutect = mutect_vcf.map{meta, vcf, tbi -> tuple(meta, "mutect", vcf, tbi)}
    deepsomatic = deepsomatic_vcf.map{meta, vcf, tbi -> tuple(meta, "deepsomatic", vcf, tbi) }

    vcfs = mutect.mix(deepsomatic).mix(strelka)
    vcfs_filtered = FILTER_VCF(vcfs).filtered_vcf
    vcfs_normalized = POSTPROCESS_VCF(vcfs_filtered, reference_genome)
    
    callers = vcfs_normalized.branch{ meta, caller, vcf, tbi ->
        mutect: caller == "mutect"
        strelka: caller == "strelka"
        deepsomatic: caller == "deepsomatic"
    }
    
    merged_callers = callers.mutect.join(callers.deepsomatic).join(callers.strelka)


    merged_vcf = MERGE_SOMATIC_VCFS(merged_callers, reference_genome, reference_dict)
    vep = VEP_ANNOTATE(merged_vcf, reference_genome, vep_cache, vep_plugins)
    vep_filtered = VEP_POPULATION_FILTER(vep.vcf, vep_cache, vep_plugins)

    vep_filtered_somatic_name = vep_filtered.map {meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    all_samples_somatic_name = preproc_bams.mix(star_bam).map {meta, bam, bai -> tuple(meta.somatic_name, meta, bam,bai)}
    

    bamreadcount_helper_input = vep_filtered_somatic_name.combine(all_samples_somatic_name, by:0)

    brc_helper = BAMREADCOUNT(bamreadcount_helper_input, reference_genome)
    
    brc_helper_branched = brc_helper.branch { somatic_meta, sample_meta, indels, snvs ->
        tumor_dna: sample_meta.sample_type == "TUMOR" && sample_meta.molecule == "DNA"
        normal_dna: sample_meta.sample_type == "NORMAL" && sample_meta.molecule == "DNA"
        tumor_rna: sample_meta.sample_type == "TUMOR" && sample_meta.molecule == "RNA"
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

    index_input = gene_vcf.map{meta, vcf ->
        tuple(meta.somatic_name, meta, vcf)
    }
    final_vcf = INDEX_VCF(index_input, "variants")
    final_vcf_table_input = final_vcf.map { meta, vcf, tbi -> tuple(meta, meta.somatic_name + "_somatic_variants", vcf, tbi) }
    final_vcf_table = VCF_TO_TABLE(final_vcf_table_input)


    // CREATE PHASED VCF 
    tumor_normal_samples = preproc_bams.branch {meta, bam, bai ->
        tumor:meta.sample_type == "TUMOR"
        normal:meta.sample_type == "NORMAL"
    }
    
    // Phase with germline calls
    phased_vcf = PVAC_VCF_PHASING(
        tumor_normal_samples.tumor,
        final_vcf,
        germline_vcf,
        reference_genome,
        reference_dict,
        vep_cache,
        vep_plugins
    )

    phased_vcf_final = phased_vcf.final_phased

    emit:
        somatic_vcf = final_vcf
        somatic_vcf_table = final_vcf_table
        phased_vcf = phased_vcf_final
        vep_report = vep.report
}
