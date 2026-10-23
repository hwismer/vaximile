// Quality control: alignment stats, coverage, library damage, relatedness and variant QC, gathered
// with every other step's reports into one MultiQC report per patient.
include { SAMTOOLS_FLAGSTAT } from "../modules/local/samtools_flagstat/main"
include { SAMTOOLS_COVERAGE } from "../modules/local/samtools_coverage/main"
include { SAMTOOLS_IDXSTATS } from "../modules/local/samtools_idxstats/main"
include { MOSDEPTH } from "../modules/local/mosdepth/main"
include { PICARD_ARTIFACT_METRICS } from "../modules/local/picard_artifact_metrics/main"
include { SOMALIER_EXTRACT } from "../modules/local/somalier_extract/main"
include { SOMALIER_RELATE } from "../modules/local/somalier_relate/main"
include { SOMATIC_QC } from "../modules/local/somatic_qc/main"
include { GERMLINE_QC } from "../modules/local/germline_qc/main"
include { VAF_GRIDS } from "../modules/local/vaf_grids/main"
include { MULTIQC } from "../modules/local/multiqc/main"

// The pipeline's git commit: Nextflow knows it when run from GitHub, otherwise ask the local clone.
def pipeline_commit() {
    if ( workflow.commitId ) {
        return workflow.commitId
    }
    try {
        def commit = "git -C ${projectDir} rev-parse --short HEAD".execute().text.trim()
        def dirty = "git -C ${projectDir} status --porcelain --untracked-files=no".execute().text.trim()
        return commit ? (dirty ? "${commit} (with uncommitted changes)" : commit) : 'unknown'
    } catch ( Exception e ) {
        return 'unknown'
    }
}

// How this run was launched and every parameter it used, for the report's Run Information section.
def run_info() {
    def details = [
        "Command": workflow.commandLine,
        "Run name": workflow.runName,
        "Session ID": workflow.sessionId,
        "Started": workflow.start,
        "Resumed": workflow.resume,
        "Nextflow version": nextflow.version,
        "Pipeline commit": pipeline_commit(),
        "Profile": workflow.profile,
        "Launch directory": workflow.launchDir,
        "Work directory": workflow.workDir,
        "Output directory": params.outdir,
    ].collectEntries { k, v -> [k, v.toString()] }
    // Resources left unset are downloaded or built into ./vaximile_resources/; show those paths instead.
    def resources = "${workflow.launchDir}/vaximile_resources"
    def resource_defaults = [
        "bwa_index": "${resources}/bwa_index/ (built)",
        "star_index": "${resources}/star_index/STARGenomeDir (built)",
        "salmon_index": "${resources}/salmon_index/${file(params.transcriptome_reference).name}_salmon_index (built)",
        "ctat_resource_lib": "${resources}/ctat_resource_dir/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir (downloaded)",
        "vep_cache": "${resources}/vep_cache/vep_cache (downloaded)",
        "vep_plugins": "${resources}/vep_pvac_plugins/VEP_plugins (downloaded)",
        "arriba_blacklist": "${resources}/arriba_hg38/arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz (downloaded)",
        "arriba_known_fusions": "${resources}/arriba_hg38/arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz (downloaded)",
        "arriba_protein_domains": "${resources}/arriba_hg38/arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3 (downloaded)",
        "ascat_resources": "${resources}/ascat_hg38_${params.reference_includes_chr_prefix ? 'chr' : 'nochr'} (downloaded)",
    ]
    def run_params = params.sort().collectEntries { k, v -> [k, v != null ? v.toString() : (resource_defaults[k] ?: '')] }
    return [details: details, params: run_params]
}

workflow QC {

    take:
        markdup_bams          // (meta, bam, bai), every DNA library
        capture_kits          // (capture kit, bed)
        reference_genome      // (fasta, fai)
        reference_dict
        somalier_sites
        somatic_caller_vcfs   // (pair meta, caller, vcf, tbi), each caller's PASS calls
        somatic_consensus     // (pair meta, vcf)
        somatic_vcf           // (pair meta, vcf, tbi), after every filter
        contamination         // (pair meta, table)
        pvacseq_aggregated    // (pair meta, report)
        germline_caller_vcfs  // (normal meta, caller, vcf, tbi)
        germline_vcf          // (normal meta, vcf, tbi)
        reports               // (meta, file or files) from other steps, which MultiQC reads as they are

    main:

    // Alignment stats, capture coverage and library damage, per DNA library.
    flagstat = SAMTOOLS_FLAGSTAT(markdup_bams).flagstat
    coverage = SAMTOOLS_COVERAGE(markdup_bams).tsv
    idxstats = SAMTOOLS_IDXSTATS(markdup_bams).tsv

    mosdepth_input = markdup_bams
        .map { meta, bam, bai -> tuple(meta.capture_kit, meta, bam, bai) }
        .combine(capture_kits, by: 0)
        .map { _kit, meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }
    mosdepth = MOSDEPTH(mosdepth_input).dist

    artifact_metrics = PICARD_ARTIFACT_METRICS(markdup_bams, reference_genome, reference_dict).mqc

    // Relatedness of every library from the same patient, to catch sample swaps.
    somalier_extract = SOMALIER_EXTRACT(markdup_bams, reference_genome, somalier_sites).somalier
        .unique { meta, _files -> meta.sample_name }
        .map { meta, extracted -> tuple(meta.patient, meta, extracted) }
        .groupTuple()
    somalier = SOMALIER_RELATE(somalier_extract)

    // Somatic, per pair: each caller's calls, the consensus, the final calls, contamination, the capture BED and pVACseq's report.
    somatic_qc_input = somatic_caller_vcfs
        .map { meta, _caller, vcf, _tbi -> tuple(meta.somatic_name, vcf) }
        .groupTuple(size: 3)
        .join(somatic_consensus.map { meta, vcf -> tuple(meta.somatic_name, meta, vcf) })
        .join(somatic_vcf.map { meta, vcf, _tbi -> tuple(meta.somatic_name, vcf) })
        .join(contamination.map { meta, table -> tuple(meta.somatic_name, table) })
        .map { _name, caller_vcfs, meta, consensus, final_vcf, contamination_table ->
            tuple(meta.capture_kit, meta, caller_vcfs, consensus, final_vcf, contamination_table)
        }
        .combine(capture_kits, by: 0)
        .map { _kit, meta, caller_vcfs, consensus, final_vcf, contamination_table, bed ->
            tuple(meta.somatic_name, meta, caller_vcfs, consensus, final_vcf, contamination_table, bed)
        }
        // remainder: true so a pair whose pVACseq run failed still gets its variant QC.
        .join(pvacseq_aggregated.map { meta, report -> tuple(meta.somatic_name, report) }, remainder: true)
        .filter { items -> items[1] instanceof Map }
        .map { _name, meta, caller_vcfs, consensus, final_vcf, contamination_table, bed, report ->
            tuple(meta, caller_vcfs, consensus, final_vcf, contamination_table, bed, report ?: [])
        }
    somatic_qc = SOMATIC_QC(somatic_qc_input)

    // Germline, per normal: each caller's calls and the consensus.
    germline_qc_input = germline_caller_vcfs
        .map { meta, _caller, vcf, _tbi -> tuple(meta.sample_name, vcf) }
        .groupTuple(size: 3)
        .join(germline_vcf.map { meta, vcf, _tbi -> tuple(meta.sample_name, meta, vcf) })
        .map { _name, caller_vcfs, meta, consensus -> tuple(meta, caller_vcfs, consensus) }
    germline_qc = GERMLINE_QC(germline_qc_input)

    // The VAF histograms and scatter plots, drawn per patient with one panel per pair or normal.
    vaf_tables = somatic_qc.vaf_tables
        .mix(germline_qc.vaf_tables)
        .map { meta, tables -> tuple(meta.patient, tables) }
        .groupTuple()
        .map { patient, tables -> tuple(patient, tables.flatten()) }
    vaf_grids = VAF_GRIDS(vaf_tables).png

    // MultiQC takes (patient, library name, file). The library name lets it rename samples to
    // <sample>_<molecule>; files for a pair or the whole patient have none.
    library_and_pair_files = reports
        .mix(
            flagstat,
            coverage,
            idxstats,
            mosdepth.map { meta, global, region, summary -> tuple(meta, [global, region, summary]) },
            artifact_metrics.map { meta, oxog, table -> tuple(meta, [oxog, table]) },
            somatic_qc.mqc,
            germline_qc.mqc
        )
        .flatMap { meta, files ->
            def library = meta.sample_name ? meta.sample_name + "_" + meta.molecule : null
            [files].flatten().collect { f -> tuple(meta.patient, library, f) }
        }
    patient_files = somalier.pairs
        .mix(somalier.samples, vaf_grids)
        .flatMap { patient, files -> [files].flatten().collect { f -> tuple(patient, null, f) } }

    multiqc = MULTIQC(library_and_pair_files.mix(patient_files).groupTuple(), run_info()).html

    emit:
        multiqc = multiqc
}
