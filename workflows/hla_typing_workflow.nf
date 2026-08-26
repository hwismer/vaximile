include { OPTITYPE; HLAHD; HLA_CALLS_PVAC; EXTRACT_MHC_REGION; BAM_TO_FASTQ; MHC_REGION_FASTQS; HLAHD_TO_TSV; HLA_BED; MHCFLOW }  from "../modules/local/hla_typing.nf"

workflow HLA_TYPING_WORKFLOW {

    take:
        bams // (metamap, sorted bam, bai)
        chr_prefix
        hla_fasta
        hla_kmers
        hla_freqs

    main:

        hla_bed = HLA_BED(chr_prefix)
        
        mhcflow = MHCFLOW(bams, hla_fasta, hla_bed, hla_kmers, hla_freqs)


        fastqs = MHC_REGION_FASTQS(bams)
        //mhc_region = EXTRACT_MHC_REGION(bams)
        //fastqs = BAM_TO_FASTQ(mhc_region)
        // Run Optitype and HLA-HD on fastqs
        optitype = OPTITYPE(fastqs)
        hlahd = HLAHD(fastqs)
        
        // Combine HLA calls as input for consensus calls
        combined_typing_results = optitype.hla_calls.join(hlahd.final_hla_calls)
        
        hlahd_tsv = HLAHD_TO_TSV(hlahd.final_hla_calls)
        // Final HLA calls
        hla_calls_pvac_input = HLA_CALLS_PVAC(combined_typing_results)

    emit:
        optitype = optitype.hla_calls
        hlahd = hlahd.final_hla_calls
        hlahd_tsv = hlahd_tsv
        pvac_input = hla_calls_pvac_input.pvac_calls

}
