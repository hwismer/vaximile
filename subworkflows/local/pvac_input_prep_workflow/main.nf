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
include { fan_out_pairs; expand_pairs } from "../utils_nfcore_vaximile_pipeline"

workflow PVAC_INPUT_PREP_WORKFLOW {

    take:
        mutect_vcf
        strelka_vcf
        deepsomatic_vcf
        preproc_bams
        star_bam
        salmon_tx_abundance
        salmon_gene_abundance
        reference_genome
        reference_dict
        vep_cache
        vep_plugins
        germline_vcf

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
    vep = VEP_ANNOTATE(merged_vcf, reference_genome, vep_cache, vep_plugins)
    vep_filtered = VEP_POPULATION_FILTER(vep.vcf, vep_cache, vep_plugins).vcf

    vep_filtered_somatic_name = vep_filtered.map {meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    // preproc_bams and star_bam are sample-level, so their metas carry a somatic_names
    // list rather than a scalar. Fan out to one entry per pair before combining below.
    all_samples_somatic_name = fan_out_pairs(preproc_bams.mix(star_bam))
    

    bamreadcount_helper_input = vep_filtered_somatic_name.combine(all_samples_somatic_name, by:0)
        .map { somatic_name, somatic_meta, vcf, sample_meta, bam, bai ->
            tuple(somatic_name, somatic_meta, sample_meta.sample_name, sample_meta.molecule, vcf, sample_meta, bam, bai)
        }

    brc_helper = BAMREADCOUNT(bamreadcount_helper_input, reference_genome).brc_files
    
    brc_helper_branched = brc_helper.branch { somatic_meta, sample_meta, indels, snvs ->
        tumor_dna: sample_meta.sample_type == "TUMOR" && sample_meta.molecule == "DNA"
        normal_dna: sample_meta.sample_type == "NORMAL" && sample_meta.molecule == "DNA"
        tumor_rna: sample_meta.sample_type == "TUMOR" && sample_meta.molecule == "RNA"
    }
    

    annotate_vcf_coverage_tumor_dna_input = brc_helper_branched.tumor_dna.join(vep_filtered)
        .map { somatic_meta, sample_meta, indels, snvs, vcf ->
            tuple(somatic_meta, sample_meta, sample_meta.sample_name, sample_meta.molecule, somatic_meta.somatic_name, indels, snvs, vcf)
        }

    tumor_dna = ANNOTATE_VCF_COVERAGE_TUMOR_DNA(annotate_vcf_coverage_tumor_dna_input).vcf

    annotate_vcf_coverage_normal_dna_input = brc_helper_branched.normal_dna.join(tumor_dna)
        .map { somatic_meta, sample_meta, indels, snvs, vcf ->
            tuple(somatic_meta, sample_meta, sample_meta.sample_name, sample_meta.molecule, somatic_meta.somatic_name, indels, snvs, vcf)
        }

    normal_dna_tdna = ANNOTATE_VCF_COVERAGE_NORMAL_DNA(annotate_vcf_coverage_normal_dna_input).vcf

    annotate_vcf_coverage_tumor_rna_input = brc_helper_branched.tumor_rna.join(normal_dna_tdna)
        .map { somatic_meta, sample_meta, indels, snvs, vcf ->
            tuple(somatic_meta, sample_meta, sample_meta.sample_name, sample_meta.molecule, somatic_meta.somatic_name, indels, snvs, vcf)
        }

    tumor_rna_ndna_tdna = ANNOTATE_VCF_COVERAGE_TUMOR_RNA(annotate_vcf_coverage_tumor_rna_input).vcf


    somatic_name_vcf_coverage = tumor_rna_ndna_tdna.map{meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    // salmon output is RNA sample-level, so fan out to key it by pair.
    somatic_name_tx = fan_out_pairs(salmon_tx_abundance)
    annotate_vcf_transcript_expression_input = somatic_name_vcf_coverage.join(somatic_name_tx)
        .map { somatic_name, somatic_meta, vcf, sample_meta, tx ->
            tuple(somatic_name, somatic_meta, sample_meta.sample_name, vcf, sample_meta, tx)
        }

    tx_vcf = ANNOTATE_VCF_TRANSCRIPT_EXPRESSION(annotate_vcf_transcript_expression_input).vcf

    tx_vcf_somatic_name = tx_vcf.map{ meta, vcf -> tuple(meta.somatic_name, meta, vcf) }
    somatic_name_gene = fan_out_pairs(salmon_gene_abundance)
    annotate_vcf_gene_expression_input = tx_vcf_somatic_name.join(somatic_name_gene)
        .map { somatic_name, somatic_meta, vcf, sample_meta, gene ->
            tuple(somatic_name, somatic_meta, sample_meta.sample_name, vcf, sample_meta, gene)
        }

    gene_vcf = ANNOTATE_VCF_GENE_EXPRESSION(annotate_vcf_gene_expression_input).vcf

    index_input = gene_vcf.map{meta, vcf ->
        tuple(meta.somatic_name, meta, vcf)
    }
    final_vcf = INDEX_VCF(index_input, "variants").vcf
    final_vcf_table_input = final_vcf.map { meta, vcf, tbi -> tuple(meta, meta.somatic_name + "_somatic_variants", vcf, tbi) }
    final_vcf_table = VCF_TO_TABLE(final_vcf_table_input).tsv


    // CREATE PHASED VCF 
    // expand_pairs, not fan_out_pairs: PVAC_VCF_PHASING joins these BAMs on the meta map
    // held in somatic_meta.tumor_meta, so the meta must match that one exactly.
    tumor_normal_samples = expand_pairs(preproc_bams).branch {meta, bam, bai ->
        tumor:meta.sample_type == "TUMOR"
        normal:meta.sample_type == "NORMAL"
    }
    
    // Phase with germline calls
    phased_vcf = PVAC_VCF_PHASING(
        tumor_normal_samples.tumor,
        final_vcf,
        expand_pairs(germline_vcf),
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
