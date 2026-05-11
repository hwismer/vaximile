include { PVACSEQ } from "../modules/pvactools.nf"

workflow PVACTOOLS_WORKFLOW {

    take:
        pvacseq_input // (somatic_name, somatic_meta, somatic_vcf, somatic_vcf_index, phased_vcf, phase_vcf_index, hla_calls)
        proteome_reference

    main:

    pvacseq = PVACSEQ(pvacseq_input, proteome_reference)


    emit:
        pvacseq = pvacseq

}
