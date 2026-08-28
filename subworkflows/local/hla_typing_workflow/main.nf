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
        
        mhcflow_input = bams
            .map { meta, bam, bai ->
                tuple(meta, meta.sample_name, bam, bai)
            }

        mhcflow = MHCFLOW(mhcflow_input, hla_fasta, hla_bed, hla_kmers, hla_freqs)


        fastqs = MHC_REGION_FASTQS(bams)
        //mhc_region = EXTRACT_MHC_REGION(bams)
        //fastqs = BAM_TO_FASTQ(mhc_region)
        // Run Optitype and HLA-HD on fastqs
        optitype_input = fastqs
            .map { meta, fastq1, fastq2 ->
                tuple(meta, meta.molecule, fastq1, fastq2)
            }

        optitype = OPTITYPE(optitype_input)
        hlahd = HLAHD(fastqs)
        
        // Combine HLA calls as input for consensus calls
        combined_typing_results = optitype.hla_calls.join(hlahd.final_hla_calls)
        
        hlahd_to_tsv_input = hlahd.final_hla_calls
            .map { meta, hlahd_result ->
                tuple(meta, meta.sample_name, hlahd_result)
            }

        hlahd_tsv = HLAHD_TO_TSV(hlahd_to_tsv_input)
        // Final HLA calls
        hla_calls_pvac_input = HLA_CALLS_PVAC(combined_typing_results)

    emit:
        optitype = optitype.hla_calls
        hlahd = hlahd.final_hla_calls
        hlahd_tsv = hlahd_tsv
        pvac_input = hla_calls_pvac_input.pvac_calls

}
