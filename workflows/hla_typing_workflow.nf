include { OPTITYPE; HLAHD; HLA_CALLS_PVAC }  from "../modules/hla_typing.nf"

workflow HLA_TYPING_WORKFLOW {

    take:
        fastqs // (metamap, fastq1, fastq2)

    main:
        
        // Run Optitype and HLA-HD on fastqs
        optitype = OPTITYPE(fastqs)
        hlahd = HLAHD(fastqs)
        
        // Combine HLA calls as input for consensus calls
        combined_typing_results = optitype.hla_calls.join(hlahd.final_hla_calls)

        // Final HLA calls
        hla_calls_pvac_input = HLA_CALLS_PVAC(combined_typing_results)

    emit:
        optitype = optitype.hla_calls
        hlahd = hlahd.final_hla_calls
        pvac_input = hla_calls_pvac_input.pvac_calls

}
