/*
    hla: hla_typing_workflow, hla_loh_workflow

    One file per pipeline step. Each workflow keeps the take/emit signature it had
    as its own subworkflow directory, so callers are unchanged.
*/
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
        hla_fasta
        hla_kmers
        hla_freqs

    main:

        hla_bed = HLA_BED(chr_prefix).bed
        
        mhcflow_input = bams
            .map { meta, bam, bai ->
                tuple(meta, meta.sample_name, bam, bai)
            }

        // .first() to make the indexed reference an explicit value channel, matching how
        // the un-indexed reference reaches this subworkflow (main.nf .first()s each HLA
        // reference param). Nextflow would convert this single-item queue implicitly, but
        // relying on that means the reference silently stops being reusable the day it
        // emits more than one item.
        hla_reference_indexed = NOVOALIGN_HLA_FASTA(hla_fasta).out.first()

        mhcflow = MHCFLOW(mhcflow_input, hla_reference_indexed, hla_bed, hla_kmers, hla_freqs)


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
        // Final HLA calls
        hla_calls_pvac_input = HLA_CALLS_PVAC(combined_typing_results)

    emit:
        // The subject-specific reference and realignment the LOH workflow builds on, and
        // hla_bed, which HLA_LOH_WORKFLOW needs for the tumour's second mhcflow run and
        // which is produced in here rather than by the caller.
        mhcflow_hla_ref = mhcflow.sample_hla_ref
        mhcflow_realn_bam = mhcflow.realn_bam
        hla_bed = hla_bed
        optitype = optitype.hla_calls
        hlahd = hlahd.final_hla_calls
        hlahd_tsv = hlahd_tsv
        pvac_input = hla_calls_pvac_input.pvac_calls

}

/*
    HLA loss of heterozygosity, following the three steps mhcflow documents.

    Step 1 is already done by the time this runs: MHCFLOW types every DNA library, and its
    finalizer stage leaves each sample realigned against the alleles inferred for that
    sample. What is left is to pair a tumour with its own normal, realign the tumour against
    THAT normal's HLA reference, and compare the two.

    Pairing is the whole point. A tumour realigned against its own alleles cannot be compared
    with a normal realigned against the normal's: lohhlamod reads coverage at mismatch sites
    across both BAMs, so both have to sit on one reference, and that reference is the
    normal's - the germline genotype the tumour is being tested for losing.
*/
workflow HLA_LOH_WORKFLOW {

    take:
        markdup_bams      // (meta, bam, bai) - every DNA library, pre-mhcflow
        sample_hla_ref    // (meta, hla_fasta, hla_nix) - MHCFLOW finalizer, per library
        realn_bam         // (meta, realn_bam, realn_bai) - MHCFLOW finalizer, per library
        purityploidy      // (somatic_meta, purityploidy) - ASCAT, per pair
        hla_bed
        hla_kmers
        hla_freqs
        montage_script

    main:

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
        tumor_realn = MHCFLOW_REALIGN(realign_input, hla_bed, hla_kmers, hla_freqs).realn_bam

        // Step 3: both BAMs, the normal's reference, and ASCAT's purity and ploidy.
        //
        // join(remainder: true) for ASCAT, not a plain join: ASCAT is allowed to fail, and
        // an inner join would drop the whole pair from the LOH analysis when it does,
        // silently. The remainder arrives as null, which the module reads as "no estimates"
        // and lohhlamod fills with its own defaults.
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

        // Coverage, logR and BAF profiles per HLA gene. join on somatic_name rather than
        // zipping the two emits: they come from the same process, but order across pairs is
        // not guaranteed, and a mismatch here would plot one pair's result over another's
        // intermediates without failing.
        plot_input = loh.loh_res
            .map { somatic_name, meta, res -> tuple(somatic_name, meta, res) }
            .join(loh.loh_rds.map { somatic_name, _meta, rds -> tuple(somatic_name, rds) })

        plots = LOHHLAPLOT(plot_input).plots

        // The plots are PDFs, which MultiQC cannot embed, so a rasterised copy goes to the
        // report and the PDFs stay the published artefact. The result table goes the same
        // way, re-keyed by pair.
        mqc_input = plots.join(loh.loh_res.map { somatic_name, _meta, res -> tuple(somatic_name, res) })
        loh_mqc = LOHHLA_MQC(mqc_input, montage_script)

    emit:
        loh_dir = loh.loh_dir
        loh_res = loh.loh_res
        loh_plots = plots
        loh_plots_png = loh_mqc.png
        loh_res_mqc = loh_mqc.tsv
}
