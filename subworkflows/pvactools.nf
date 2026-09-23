/*
    pvactools: pvac_vcf_phasing, pvac_input_prep_workflow, pvactools_workflow

    One file per pipeline step. Each workflow keeps the take/emit signature it had
    as its own subworkflow directory, so callers are unchanged.
*/
include { VEP_ANNOTATE } from "../modules/local/vep_annotate/main"
include { INDEX_VCF } from "../modules/local/index_vcf/main"
include { PHASE_VCF_SELECT_VARIANTS } from "../modules/local/phase_vcf_select_variants/main"
include { PHASE_VCF_RENAME } from "../modules/local/phase_vcf_rename/main"
include { PHASE_VCF_COMBINE_VARIANTS } from "../modules/local/phase_vcf_combine_variants/main"
include { PHASE_VCF_SORT_VCF } from "../modules/local/phase_vcf_sort_vcf/main"
include { PHASE_VCF_RBPHASING } from "../modules/local/phase_vcf_rbphasing/main"
include { VCF_TO_TABLE } from "../modules/local/vcf_to_table/main"
include { VEP_POPULATION_FILTER } from "../modules/local/vep_population_filter/main"
include { BAMREADCOUNT } from "../modules/local/bamreadcount/main"
include { ANNOTATE_VCF_TRANSCRIPT_EXPRESSION } from "../modules/local/annotate_vcf_transcript_expression/main"
include { ANNOTATE_VCF_GENE_EXPRESSION } from "../modules/local/annotate_vcf_gene_expression/main"
include { ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_TUMOR_DNA } from "../modules/local/annotate_vcf_coverage/main"
include { ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_NORMAL_DNA } from "../modules/local/annotate_vcf_coverage/main"
include { ANNOTATE_VCF_COVERAGE as ANNOTATE_VCF_COVERAGE_TUMOR_RNA } from "../modules/local/annotate_vcf_coverage/main"
include { fan_out_pairs; expand_pairs } from "./pipeline_init.nf"
include { PVACSEQ } from "../modules/local/pvacseq/main"
include { PVACFUSE } from "../modules/local/pvacfuse/main"
include { COMBINE_PVACSEQ_AGGREGATED_REPORT } from "../modules/local/combine_pvacseq_aggregated_report/main"

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

    phase_vcf_select_variants_input = somatic_vcf
        .map { somatic_meta, vcf, vcf_index ->
            tuple(somatic_meta, somatic_meta.tumor_meta.sample_name, vcf, vcf_index)
        }

    select_variants = PHASE_VCF_SELECT_VARIANTS(phase_vcf_select_variants_input, reference_genome, reference_dict).vcf
    

    somatic_meta_germline = somatic_vcf.map{meta, vcf, tbi ->
        tuple(meta.normal_meta, meta)
    }.join(germline_vcf, by:0)

    phase_vcf_rename_input = somatic_meta_germline
        .map { normal_meta, somatic_meta, vcf, vcf_index ->
            tuple(normal_meta, somatic_meta, normal_meta.sample_name, somatic_meta.tumor_meta.sample_name, vcf, vcf_index)
        }

    germline_renamed = PHASE_VCF_RENAME(phase_vcf_rename_input).vcf


    combine_variants_input = select_variants.join(germline_renamed, by:0)
    combined_variants = PHASE_VCF_COMBINE_VARIANTS(combine_variants_input, reference_genome, reference_dict).vcf

    combined_sorted = PHASE_VCF_SORT_VCF(combined_variants, reference_genome, reference_dict).sorted_vcf

    combined_sorted_by_tumor_sample = combined_sorted.map{ meta, vcf ->
        tuple(meta.tumor_meta, meta, vcf)
    }.join(tumor_bam)

    rbphased = PHASE_VCF_RBPHASING(combined_sorted_by_tumor_sample, reference_genome, reference_dict).vcf

    phased_vep = VEP_ANNOTATE(rbphased, reference_genome, vep_cache, vep_plugins)

    
    index_input = phased_vep.vcf.map{meta, vcf ->
        tuple(meta.tumor_meta.sample_name, meta, vcf)
    }
    final_phased = INDEX_VCF(index_input, "phased").vcf

    
    emit:
        final_phased = final_phased

}

workflow PVAC_INPUT_PREP_WORKFLOW {

    take:
        // The consensus callset, not the three caller outputs: assembling it moved to
        // SOMATIC_CONSENSUS, where the callers themselves live.
        somatic_vcf
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
    
    vep = VEP_ANNOTATE(somatic_vcf, reference_genome, vep_cache, vep_plugins)
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

workflow PVACTOOLS_WORKFLOW {

    take:
        pvacseq_input // (somatic_name, somatic_meta, somatic_vcf, somatic_vcf_index, phased_vcf, phase_vcf_index, hla_calls)
        pvacfuse_input // (somatic_name, arriba_meta, arriba_fusions, star_meta, star_fusions, hla_meta, hla_calls)
        proteome_reference

    main:

    pvacseq_ch = pvacseq_input
        .map { somatic_name, somatic_meta, somatic_vcf, somatic_vcf_index, phased_vcf, phased_vcf_index, hla_meta, hla_pvac_input ->
            tuple(somatic_name, somatic_meta, somatic_meta.tumor_meta.sample_name, somatic_meta.normal_meta.sample_name, somatic_vcf, somatic_vcf_index, phased_vcf, phased_vcf_index, hla_meta, hla_pvac_input)
        }

    pvacseq = PVACSEQ(pvacseq_ch, proteome_reference)
    pvacfuse = PVACFUSE(pvacfuse_input, proteome_reference).pvacfuse_dir

    pvacseq_patient = pvacseq.pvaseq_mhc_i_aggr.map{meta, report -> tuple(meta.patient, report)}.groupTuple()
    combined_report = COMBINE_PVACSEQ_AGGREGATED_REPORT(pvacseq_patient).tsv


    emit:
        pvacseq = pvacseq.pvacseq_dir
        pvacseq_mhc_i_combined = combined_report
        pvacfuse = pvacfuse

}
