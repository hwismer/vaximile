/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VAXIMILE
    Tumor neoantigen discovery from paired tumour-normal bulk DNA and RNA sequencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

include { VAXIMILE } from './workflows/vaximile'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_vaximile_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_vaximile_pipeline'

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

workflow {

    main:

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

    // Every module reports to the `versions` topic. A topic channel collects them with no
    // per-process wiring, which is why versions are declared with `topic:` rather than
    // `emit:` - 94 modules would otherwise each need threading through their subworkflow.
    //
    // Two shapes arrive here. Local modules emit a versions.yml path. Unmodified nf-core
    // modules (SAMTOOLS_SORMADUP, ASCAT) instead emit a tuple of
    // (process, tool, version) built with eval(). Handing those tuples straight to
    // collectFile made it read the process name as a *filename*: the version was dropped
    // from the report and an empty file with a process-shaped name appeared in
    // pipeline_info/. So the tuples are rendered to YAML text, and the versions.yml files
    // are read to text alongside them: `sort: true` compares the collected items against
    // each other, and it cannot compare a String to a Path, so both must be the same type.
    versions = channel.topic('versions').branch { item ->
        nf_core: item instanceof Collection
        versions_yml: true
    }

    software_versions = versions.versions_yml
        .map { yml -> yml.text }
        .mix( versions.nf_core.map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" } )
        .collectFile(name: 'software_versions.yml', sort: true, newLine: false)

    publish:
    software_versions = software_versions
    multiqc_reports = VAXIMILE.out.multiqc_reports
    somatic_vcf = VAXIMILE.out.somatic_vcf
    somatic_vcf_table = VAXIMILE.out.somatic_vcf_table
    germline_vcf = VAXIMILE.out.germline_vcf
    optitype_calls = VAXIMILE.out.optitype_calls
    hlahd_calls = VAXIMILE.out.hlahd_calls
    hla_pvac_input = VAXIMILE.out.hla_pvac_input
    pvacseq = VAXIMILE.out.pvacseq
    pvacseq_mhc_i_combined = VAXIMILE.out.pvacseq_mhc_i_combined
    pvacfuse = VAXIMILE.out.pvacfuse
    kallisto_gene = VAXIMILE.out.kallisto_gene
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PUBLISH TARGETS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*
    Where an output belongs in the published tree.

    Pair-level outputs carry a somatic_meta with a scalar somatic_name and go under the
    pair. Sample-level outputs - germline calls, per-sample HLA typing, RNA quantification -
    carry a somatic_names list instead, because a library shared between pairs is processed
    once and no longer belongs to exactly one of them. Those go under the sample.

    Reading meta.somatic_name unconditionally is what published them to a literal "null"
    directory once libraries began being deduplicated.
*/
def publish_scope(meta) {
    return meta.somatic_name ?: meta.sample_name
}

output {
    software_versions {
        path { _v -> "${params.outdir}/pipeline_info/" }
    }
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
    kallisto_gene {
        path { meta, gene_abundance -> "${params.outdir}/${meta.patient}/${publish_scope(meta)}/kallisto" }
    }

}
