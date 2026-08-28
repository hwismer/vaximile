/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow with functionality specific to the vaximile pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow PIPELINE_INITIALISATION {

    take:
    samplesheet    // string: path to the patient samplesheet
    capture_kits   // string: path to the capture kit samplesheet

    main:

    validateParameters()
    log.info(paramsSummaryLog(workflow))
    validateReferenceInputs()

    ch_capture_kits = Channel.fromPath(capture_kits, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.kit, row.bed) }

    ch_samplesheet = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->

            // Be semi-flexible in parsing the sequencing type
            def molecule = (row.sequencing_type.toLowerCase() == 'rna') ? 'RNA' :
                (row.sequencing_type.toLowerCase() in ['exome', 'genome', "exome_FFPE", "genome_FFPE"]) ? 'DNA' :
                null

            def meta = [
                somatic_name: row.somatic_name,
                patient: row.patient,
                sample_name: row.sample_name,
                sex: row.sex,
                sample_type: row.sample_type.toUpperCase(),
                sequencing_type: row.sequencing_type,
                capture_kit: row.capture_kit,
                molecule: molecule
            ]

            def reads = [
                file(row.fastqr1, checkIfExists: true),
                file(row.fastqr2, checkIfExists: true)
            ]

            return [meta, reads]
        }

    emit:
    samplesheet  = ch_samplesheet
    capture_kits = ch_capture_kits
}

workflow PIPELINE_COMPLETION {

    take:
    multiqc_report // channel: [ patient, report ]

    main:

    // Capture the metadata object here, while the workflow binding is still resolvable.
    // Referencing `workflow` from inside the closure instead resolves to null when the
    // handler actually fires from within a named workflow body, which ended every run
    // with "Failed to invoke `workflow.onComplete` event handler" and an NPE on
    // `workflow.success`.
    def wf = workflow
    def outdir = params.outdir

    wf.onComplete {
        log.info(
            wf.success
                ? "Pipeline completed successfully. Results in ${outdir}"
                : "Pipeline completed with errors. Exit status: ${wf.exitStatus}"
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Checks nextflow_schema.json cannot express: the VEP cache is a user-supplied
// directory whose absence only surfaces deep into annotation otherwise.
//
def validateReferenceInputs() {
    if (params.vep_cache && !file(params.vep_cache).exists()) {
        error("VEP cache directory not found: ${params.vep_cache}")
    }
}
