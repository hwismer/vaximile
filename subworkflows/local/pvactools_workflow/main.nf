include { PVACSEQ } from "../../../modules/local/pvacseq/main"
include { PVACFUSE } from "../../../modules/local/pvacfuse/main"
include { COMBINE_PVACSEQ_AGGREGATED_REPORT } from "../../../modules/local/combine_pvacseq_aggregated_report/main"

workflow PVACTOOLS_WORKFLOW {

    take:
        pvacseq_input // (somatic_name, somatic_meta, somatic_vcf, somatic_vcf_index, phased_vcf, phase_vcf_index, hla_calls)
        pvacfuse_input // (somatic_name, arriba_meta, arriba_fusions, star_meta, star_fusions, hla_meta, hla_calls)
        proteome_reference

    main:

    pvacseq = PVACSEQ(pvacseq_input, proteome_reference)
    pvacfuse = PVACFUSE(pvacfuse_input, proteome_reference)

    pvacseq_patient = pvacseq.pvaseq_mhc_i_aggr.map{meta, report -> tuple(meta.patient, report)}.groupTuple()
    combined_report = COMBINE_PVACSEQ_AGGREGATED_REPORT(pvacseq_patient)


    emit:
        pvacseq = pvacseq.pvacseq_dir
        pvacseq_mhc_i_combined = combined_report
        pvacfuse = pvacfuse

}
