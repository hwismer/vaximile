// vaximile: tumour neoantigen discovery from paired tumour/normal DNA and tumour RNA sequencing.

include { VAXIMILE } from './workflows/vaximile'
include { PIPELINE_INITIALISATION } from './subworkflows/pipeline_init.nf'
include { PIPELINE_COMPLETION } from './subworkflows/pipeline_init.nf'

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

workflow {

    main:

    // --help is printed by nf-schema before the workflow runs; exit here rather than continue.
    if (params.help || params.containsKey('helpFull')) {
        System.exit(0)
    }

    PIPELINE_INITIALISATION(
        params.samplesheet,
        params.capture_kits
    )

    VAXIMILE(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.capture_kits
    )

    PIPELINE_COMPLETION(
        VAXIMILE.out.multiqc_reports
    )

    publish:
    multiqc_reports = VAXIMILE.out.multiqc_reports
    markdup_bams = VAXIMILE.out.markdup_bams
    preproc_bams = VAXIMILE.out.preproc_bams
    star_bam = VAXIMILE.out.star_bam
    germline_vcf_table = VAXIMILE.out.germline_vcf_table
    ascat_results = VAXIMILE.out.ascat_results
    hla_loh = VAXIMILE.out.hla_loh
    hla_loh_plots = VAXIMILE.out.hla_loh_plots
    somatic_vcf = VAXIMILE.out.somatic_vcf
    somatic_vcf_table = VAXIMILE.out.somatic_vcf_table
    germline_vcf = VAXIMILE.out.germline_vcf
    optitype_calls = VAXIMILE.out.optitype_calls
    hlahd_calls = VAXIMILE.out.hlahd_calls
    hla_pvac_input = VAXIMILE.out.hla_pvac_input
    pvacseq = VAXIMILE.out.pvacseq
    pvacseq_mhc_i_combined = VAXIMILE.out.pvacseq_mhc_i_combined
    pvacfuse = VAXIMILE.out.pvacfuse
    salmon_gene = VAXIMILE.out.salmon_gene
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PUBLISH TARGETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

// Publish pair-level outputs under <somatic_name>/ and per-library outputs under samples/<sample_name>/.
def publish_scope(meta) {
    return meta.somatic_name ? "${meta.somatic_name}" : "samples/${meta.sample_name}"
}

output {
    multiqc_reports {
        path { patient, report -> "${params.outdir}/${patient}/multiqc/" }
    }
    somatic_vcf {
        path { meta, vcf, vcf_index -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/variants" }
    }
    somatic_vcf_table {
        path { meta, table -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/variants" }
    }
    optitype_calls {
        path { meta, tsv, pdf -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/hla/optitype/" }
    }
    hlahd_calls {
        path { meta, calls -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/hla/hlahd/" }
    }
    hla_pvac_input {
        path { meta, calls -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/hla/" }
    }
    pvacseq {
        path { meta, pvacseq_dir -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/pvactools/" }
    }
    pvacseq_mhc_i_combined {
        path { patient, report -> "${params.outdir}/${patient}/pvactools_report" }
    }
    pvacfuse {
        path { meta, pvacfuse_dir -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/pvactools" }
    }
    germline_vcf {
        path { meta, vcf, tbi -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/germline/" }
    }
    // Alongside the VCF, matching how somatic_vcf and somatic_vcf_table share variants/.
    germline_vcf_table {
        path { meta, tsv -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/germline/" }
    }
    // All three BAM sets share alignment/; their filenames already distinguish them
    // (_markdup.bam, _bqsr.bam, _STAR_sorted.bam).
    markdup_bams {
        path { meta, bam, bai -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/alignment/" }
    }
    preproc_bams {
        path { meta, bam, bai -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/alignment/" }
    }
    star_bam {
        path { meta, bam, bai -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/alignment/" }
    }
    // Pair-level: ASCAT is called on a tumour/normal pair. Published beside the LOH
    // results rather than in its own directory, because its purity and ploidy estimates are
    // what the LOH copy number inference is built on, and the two are read together.
    ascat_results {
        path { meta, result -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/hla/loh/" }
    }
    salmon_gene {
        path { meta, gene_abundance -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/salmon" }
    }
    // Pair-level, and keyed by the pair rather than by a sample: LOH is a property of a
    // tumour measured against its own normal, so the somatic_meta carries the patient and
    // the two channels arrive as (somatic_name, meta, files).
    hla_loh {
        path { _somatic_name, meta, res -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/hla/loh/" }
    }
    hla_loh_plots {
        path { _somatic_name, meta, plots -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/hla/loh/" }
    }

}
