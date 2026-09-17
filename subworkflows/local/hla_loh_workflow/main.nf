include { MHCFLOW_REALIGN } from "../../../modules/local/mhcflow_realign/main"
include { LOHHLAMOD } from "../../../modules/local/lohhlamod/main"
include { LOHHLAPLOT } from "../../../modules/local/lohhlaplot/main"
include { LOHHLA_PLOTS_MQC } from "../../../modules/local/lohhla_plots_mqc/main"
include { fan_out_pairs } from "../utils_nfcore_vaximile_pipeline"

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
        // report and the PDFs stay the published artefact.
        plots_png = LOHHLA_PLOTS_MQC(plots).png

    emit:
        loh_dir = loh.loh_dir
        loh_res = loh.loh_res
        loh_plots = plots
        loh_plots_png = plots_png
}
