// Pipeline initialisation: parameter validation, samplesheet parsing and shared helpers.
include { validateParameters; paramsSummaryLog } from "plugin/nf-schema"



workflow PIPELINE_INITIALISATION {

    take:
    samplesheet    // string: path to the patient samplesheet
    capture_kits   // string: path to the capture kit samplesheet

    main:

    validateParameters()
    log.info(paramsSummaryAll())
    validateReferenceInputs()

    // file() so relative BED paths resolve and missing files fail at launch.
    ch_capture_kits = Channel.fromPath(capture_kits, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.kit, file(row.bed, checkIfExists: true)) }

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

    // Capture `workflow` here; it is null inside the onComplete closure.
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

// Functions

// Log every parameter by schema section, marking non-defaults (nf-schema's own summary shows only non-defaults).
def paramsSummaryAll() {
    def schema = new groovy.json.JsonSlurper().parse(file("${projectDir}/nextflow_schema.json").toFile())
    def described = [] as Set
    def lines = ["", "-" * 78, "Parameters for this run", "-" * 78]

    schema['$defs'].each { section_key, section ->
        def props = section.properties ?: [:]
        if (!props) {
            return
        }
        lines << ""
        lines << "  ${section.title ?: section_key}"
        props.keySet().sort().each { name ->
            described << name
            def spec = props[name]
            def value = params.containsKey(name) ? params[name] : null
            def has_default = spec.containsKey('default')
            def overridden = value != null && (!has_default || value != spec.default)
            def shown = value == null ? '(not set)' : value
            lines << "    ${name.padRight(30)} : ${shown}${overridden ? '   *' : ''}"
        }
    }

    // Anything in params that the schema does not describe - profile-only settings,
    // -params-file extras, typos in a --flag that validateParameters() let through.
    def extra = params.keySet().findAll { k -> !(k in described) }.sort()
    if (extra) {
        lines << ""
        lines << "  Not described in nextflow_schema.json"
        extra.each { name -> lines << "    ${name.padRight(30)} : ${params[name]}" }
    }

    lines << ""
    lines << "  * = differs from the schema default"
    lines << "-" * 78
    return lines.join('\n')
}

// Checks the schema can't express, e.g. that --vep_cache exists.
def validateReferenceInputs() {
    if (params.vep_cache && !file(params.vep_cache).exists()) {
        error("VEP cache directory not found: ${params.vep_cache}")
    }
}


// Library deduplication

// A library is identified by patient + sample name + molecule.
def library_key(meta) {
    return [meta.patient, meta.sample_name, meta.molecule]
}

// Collapse rows describing the same library (e.g. a normal listed once per pair) so it is processed once.
// Metas carry a somatic_names list; fan_out_pairs() restores per-pair items where needed.
def dedupe_libraries(reads_ch) {
    return reads_ch
        .map { meta, reads -> tuple(library_key(meta), meta, reads) }
        .groupTuple()
        .map { key, metas, read_sets ->
            def (patient, sample_name, molecule) = key

            // Distinct data must be given distinct sample_names; that is the samplesheet's
            // contract. A contradiction is a typo rather than a decision, so it fails at
            // launch instead of silently processing whichever row happened to sort first.
            def distinct_reads = read_sets.collect { rs -> rs.collect { f -> f.toString() } }.unique()
            if ( distinct_reads.size() > 1 ) {
                error(
                    "Samplesheet rows for sample_name '${sample_name}' (patient '${patient}', " +
                    "${molecule}) list different FastQs:\n\n" +
                    distinct_reads.collect { rs -> "  " + rs.join(', ') }.join('\n') + "\n\n" +
                    "Rows sharing patient + sample_name + molecule are treated as one library and\n" +
                    "processed once. Give genuinely different data a different sample_name."
                )
            }

            def inconsistent = ['sample_type', 'sequencing_type', 'capture_kit', 'sex'].findAll { f ->
                metas.collect { m -> m[f] }.unique().size() > 1
            }
            if ( inconsistent ) {
                error(
                    "Samplesheet rows for sample_name '${sample_name}' (patient '${patient}', " +
                    "${molecule}) disagree on: ${inconsistent.join(', ')}.\n\n" +
                    inconsistent.collect { f ->
                        "  ${f}: " + metas.collect { m -> m[f] }.unique().join(' vs ')
                    }.join('\n') + "\n\n" +
                    "These rows are treated as one library, so these values must match."
                )
            }

            def base = metas[0].findAll { k, _v -> k != 'somatic_name' }
            tuple(base + [somatic_names: metas.collect { m -> m.somatic_name }.unique().sort()], read_sets[0])
        }
}

// Inverse of dedupe_libraries: one item per pair, with a scalar somatic_name on the meta.
def expand_pairs(ch) {
    return ch.flatMap { item ->
        def meta = item[0]
        def payload = item.size() > 1 ? item[1..-1] : []
        meta.somatic_names.collect { sn ->
            [meta.findAll { k, _v -> k != 'somatic_names' } + [somatic_name: sn]] + payload
        }
    }
}

// expand_pairs, with somatic_name prepended as a join key.
def fan_out_pairs(ch) {
    return expand_pairs(ch).map { item -> [item[0].somatic_name] + item }
}


// Pair a sample-level channel into (somatic_meta, tumour payload..., normal payload...).
def pair_tumor_normal(ch) {
    def branched = fan_out_pairs(ch).branch { item ->
        tumor:  item[1].sample_type == "TUMOR"
        normal: item[1].sample_type == "NORMAL"
    }

    return branched.tumor.join(branched.normal).map { item ->
        // join emits [somatic_name, tumour_meta, tumour payload..., normal_meta, normal
        // payload...]; both payloads are the same length, so size = 3 + 2n.
        def n = (item.size() - 3).intdiv(2)
        def somatic_name = item[0]
        def tumor_meta = item[1]
        def tumor_payload = n > 0 ? item[2..(1 + n)] : []
        def normal_meta = item[2 + n]
        def normal_payload = n > 0 ? item[(3 + n)..-1] : []

        def somatic_meta = [
            somatic_name: somatic_name,
            patient: tumor_meta.patient,
            capture_kit: tumor_meta.capture_kit,
            tumor_meta: tumor_meta,
            normal_meta: normal_meta
        ]

        [somatic_meta] + tumor_payload + normal_payload
    }
}
