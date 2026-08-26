include { OPTITYPE } from "../../../modules/local/optitype/main"
include { HLAHD } from "../../../modules/local/hlahd/main"
include { HLA_CALLS_PVAC } from "../../../modules/local/hla_calls_pvac/main"
include { EXTRACT_MHC_REGION } from "../../../modules/local/extract_mhc_region/main"
include { BAM_TO_FASTQ } from "../../../modules/local/bam_to_fastq/main"
include { MHC_REGION_FASTQS } from "../../../modules/local/mhc_region_fastqs/main"
include { HLAHD_TO_TSV } from "../../../modules/local/hlahd_to_tsv/main"
include { HLA_BED } from "../../../modules/local/hla_bed/main"
include { MHCFLOW } from "../../../modules/local/mhcflow/main"

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
