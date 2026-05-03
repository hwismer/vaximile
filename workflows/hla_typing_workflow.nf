include { OPTITYPE; HLAHD; HLA_CALLS_PVAC }  from "../modules/hla_typing.nf"

workflow HLA_TYPING_WORKFLOW {

    take:
        fastqs

    main:
        
        optitype = OPTITYPE(fastqs)
        hlahd = HLAHD(fastqs)

        combined_typing_results = optitype.join(hlahd.final_hla_calls)

        hla_calls_pvac_input = HLA_CALLS_PVAC(combined_typing_results)

        
        /*
        hla_optitype_postprocess = POSTPROCESS_OPTITYPE(hla_optitype)

        hla_calls_hlahd = HLAHD_HLA_CALLS(hla_hlahd.hla_calls)
    
        hla_calls = hla_calls_hlahd.map { meta, hla_call ->
                                      def hla_value = meta.hla != "CALL" ? meta.hla : hla_call
                                      return [meta,hla_value]
                                    }.filter { meta, hla_call ->
                                               meta.sample_type == "Normal"
                                        }
        */
    

    emit:
        optitype = optitype
        hlahd = hlahd.final_hla_calls

}
