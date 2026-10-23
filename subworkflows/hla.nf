// HLA: typing (OptiType, HLA-HD) and loss of heterozygosity (mhcflow, lohhlamod).
include { OPTITYPE } from "../modules/local/optitype/main"
include { HLAHD } from "../modules/local/hlahd/main"
include { HLA_CALLS_PVAC } from "../modules/local/hla_calls_pvac/main"
include { EXTRACT_MHC_REGION } from "../modules/local/extract_mhc_region/main"
include { BAM_TO_FASTQ } from "../modules/local/bam_to_fastq/main"
include { MHC_REGION_FASTQS } from "../modules/local/mhc_region_fastqs/main"
include { HLAHD_TO_TSV } from "../modules/local/hlahd_to_tsv/main"
include { HLA_BED } from "../modules/local/hla_bed/main"
include { MHCFLOW } from "../modules/local/mhcflow/main"
include { NOVOALIGN_HLA_FASTA } from "../modules/local/novoalign_hla_fasta/main"
include { MHCFLOW_REALIGN } from "../modules/local/mhcflow_realign/main"
include { LOHHLAMOD } from "../modules/local/lohhlamod/main"
include { LOHHLAPLOT } from "../modules/local/lohhlaplot/main"
include { LOHHLA_MQC } from "../modules/local/lohhla_mqc/main"
include { fan_out_pairs } from "./pipeline_init.nf"

workflow HLA_TYPING_WORKFLOW {

    take:
        bams // (metamap, sorted bam, bai)
        chr_prefix

    main:

        hla_bed = HLA_BED(chr_prefix).bed

        fastqs = MHC_REGION_FASTQS(bams).reads
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

        hlahd_tsv = HLAHD_TO_TSV(hlahd_to_tsv_input).hlahd_tsv
        // Consensus calls for pVACtools
        hla_consensus = HLA_CALLS_PVAC(combined_typing_results)

    emit:
        hla_bed = hla_bed  // MHC region, for HLA_LOH_WORKFLOW's mhcflow runs
        optitype = optitype.hla_calls
        hlahd = hlahd.final_hla_calls
        hlahd_tsv = hlahd_tsv
        pvac_input = hla_consensus.pvac_calls
        consensus_mqc = hla_consensus.mqc

}

// HLA LOH: mhcflow builds each normal's HLA reference, each tumour is realigned against its own
// normal's reference, and lohhlamod calls loss. Runs only with --novoalign_dir, as mhcflow needs novoalign.
workflow HLA_LOH_WORKFLOW {

    take:
        markdup_bams      // (meta, bam, bai) - every DNA library
        hla_reference     // (hla_fasta, hla_fai)
        purityploidy      // (somatic_meta, purityploidy) - ASCAT, per pair
        hla_bed
        hla_kmers
        hla_freqs
        novoalign_dir     // --novoalign_dir: novoalign and novoindex

    main:

        // Step 1: mhcflow on each normal, giving its HLA reference and its reads realigned to it.
        // Value channel, so every MHCFLOW task reuses the one index.
        hla_reference_indexed = NOVOALIGN_HLA_FASTA(hla_reference, novoalign_dir).out.first()
        mhcflow_input = markdup_bams
            .filter { meta, _bam, _bai -> meta.sample_type == "NORMAL" }
            .map { meta, bam, bai -> tuple(meta, meta.sample_name, bam, bai) }
        mhcflow = MHCFLOW(mhcflow_input, hla_reference_indexed, hla_bed, hla_kmers, hla_freqs, novoalign_dir)
        sample_hla_ref = mhcflow.sample_hla_ref
        realn_bam = mhcflow.realn_bam

        // fan_out_pairs re-keys a library by somatic_name, which is what puts a normal
        // shared by two tumours into both pairs - and what stops a tumour from being
        // realigned against the wrong subject's reference.
        normal_ref = fan_out_pairs(
            sample_hla_ref.filter { meta, _fasta, _nix -> meta.sample_type == "NORMAL" }
        ).map { somatic_name, _meta, fasta, nix -> tuple(somatic_name, fasta, nix) }

        normal_realn = fan_out_pairs(
            realn_bam.filter { meta, _bam, _bai -> meta.sample_type == "NORMAL" }
        ).map { somatic_name, _meta, bam, bai -> tuple(somatic_name, bam, bai) }

        tumor_bams = fan_out_pairs(
            markdup_bams.filter { meta, _bam, _bai -> meta.sample_type == "TUMOR" }
        )

        // Step 2: the tumour against the normal's reference, realignment only.
        realign_input = tumor_bams.join(normal_ref)
        tumor_realn = MHCFLOW_REALIGN(realign_input, hla_bed, hla_kmers, hla_freqs, novoalign_dir).realn_bam

        // remainder: true so pairs without ASCAT results still get an LOH call, using lohhlamod's defaults.
        loh_input = tumor_realn
            .join(normal_realn)
            .join(normal_ref)
            .map { somatic_name, meta, tbam, tbai, nbam, nbai, fasta, _nix ->
                tuple(somatic_name, meta, tbam, tbai, nbam, nbai, fasta)
            }
            .join(
                purityploidy.map { meta, pp -> tuple(meta.somatic_name, pp) },
                remainder: true
            )
            .filter { items -> items[1] != null }
            .map { items ->
                def pp = items[7]
                items[0..6] + [pp ?: []]
            }

        loh = LOHHLAMOD(loh_input)

        // Per-gene LOH plots, joined to their results by somatic_name.
        plot_input = loh.loh_res
            .map { somatic_name, meta, res -> tuple(somatic_name, meta, res) }
            .join(loh.loh_rds.map { somatic_name, _meta, rds -> tuple(somatic_name, rds) })

        plots = LOHHLAPLOT(plot_input).plots

        // The plots are PDFs, which MultiQC cannot embed, so a rasterised copy goes to the
        // report and the PDFs stay the published artefact. The result table goes the same
        // way, re-keyed by pair.
        mqc_input = plots.join(loh.loh_res.map { somatic_name, _meta, res -> tuple(somatic_name, res) })
        loh_mqc = LOHHLA_MQC(mqc_input)

    emit:
        loh_dir = loh.loh_dir
        loh_res = loh.loh_res
        loh_plots = plots
        loh_plots_png = loh_mqc.png
        loh_res_mqc = loh_mqc.tsv
}
