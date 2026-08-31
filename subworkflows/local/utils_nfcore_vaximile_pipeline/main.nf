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
    log.info(paramsSummaryAll())
    validateReferenceInputs()

    // file(), not the raw string: downstream processes declare `path(bed)`, and a bare
    // string only resolves if it happens to be absolute - a relative path in the CSV
    // fails with "Not a valid path value" once a task is submitted. checkIfExists also
    // turns a wrong BED path into a startup error instead of a mid-run one.
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
// Log EVERY parameter the run will use, grouped by the sections in
// nextflow_schema.json, with a marker on the ones that differ from their schema default.
//
// nf-schema's paramsSummaryLog() deliberately prints only values that differ from the
// defaults ("Only displaying parameters that differ from the pipeline defaults"), which
// hides most of what a run actually depends on - every reference URI, every GATK resource
// VCF, the scatter settings. Those defaults are exactly what you want recorded next to a
// set of results, so print all of them.
//
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

//
// Checks nextflow_schema.json cannot express: the VEP cache is a user-supplied
// directory whose absence only surfaces deep into annotation otherwise.
//
def validateReferenceInputs() {
    if (params.vep_cache && !file(params.vep_cache).exists()) {
        error("VEP cache directory not found: ${params.vep_cache}")
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LIBRARY DEDUPLICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// A library is identified by patient + sample name + molecule.
def library_key(meta) {
    return [meta.patient, meta.sample_name, meta.molecule]
}

/*
    Collapse samplesheet rows that describe the same library.

    The samplesheet carries one row per (library, tumour/normal pair), so a normal shared
    by two tumours is listed twice with a different somatic_name each time. somatic_name is
    part of the meta map, and the meta map is a channel item's identity, so those rows are
    two distinct items - and every sample-level process runs twice on byte-identical data:
    fastp, alignment, markdup, BQSR, germline calling, HLA typing, somalier and the samtools
    QC. With N tumours on one normal, the normal is processed N times, and its germline VCF
    is produced N times into the same published name.

    Collapsing here yields one item per library carrying the list of pairs it belongs to.
    fan_out_pairs() restores the per-pair view at the tumour/normal join sites, so nothing
    downstream of those joins changes shape.
*/
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

/*
    Inverse of dedupe_libraries: emit one copy of each item per tumour/normal pair the
    library belongs to, with the scalar somatic_name put back on the meta in place of the
    somatic_names list.

    Arity-agnostic: element 0 is the meta, the rest is payload carried through untouched.

    Use this form where the meta map itself is the join key - PVAC_VCF_PHASING joins
    `somatic_meta.normal_meta` against the germline VCF channel and `somatic_meta.tumor_meta`
    against the tumour BAM channel, so those channels must carry a meta that is equal to the
    one stored inside somatic_meta, scalar somatic_name and all.
*/
def expand_pairs(ch) {
    return ch.flatMap { item ->
        def meta = item[0]
        def payload = item.size() > 1 ? item[1..-1] : []
        meta.somatic_names.collect { sn ->
            [meta.findAll { k, _v -> k != 'somatic_names' } + [somatic_name: sn]] + payload
        }
    }
}

/*
    expand_pairs with the somatic_name prepended as element 0, for the sites that branch and
    join on the name as a plain string key.
*/
def fan_out_pairs(ch) {
    return expand_pairs(ch).map { item -> [item[0].somatic_name] + item }
}


/*
    Pair a sample-level channel into tumour/normal pairs.

    Every somatic caller needs the same shape - one item per pair, tumour payload then
    normal payload, behind a somatic_meta - and this was hand-rolled three times in
    workflows/vaximile/main.nf with the same five-field somatic_meta copied verbatim. Three
    copies is how they stop agreeing.

    Arity-agnostic. `ch` is (meta, payload...) with any number of payload elements, and the
    result is (somatic_meta, tumour payload..., normal payload...). With (meta, bam, bai) in
    that gives (somatic_meta, tumour_bam, tumour_bai, normal_bam, normal_bai); with
    (meta, pileup) it gives (somatic_meta, tumour_pileup, normal_pileup).
*/
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
