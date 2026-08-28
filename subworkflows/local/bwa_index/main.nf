include { CREATE_BWA_INDEX } from "../../../modules/local/create_bwa_index/main"

/*
    BWA_MAP passes the reference FASTA to bwa-mem2 as the index prefix:

        bwa-mem2 mem ... $reference_fa $fastq1 $fastq2

    so bwa-mem2 looks for "<reference_fa>.0123", "<reference_fa>.amb" and so on. A
    prebuilt index is therefore only usable if its files are named after the PREPARED
    reference - the "<stem>_prc.fa" that PREPARE_FASTA writes - rather than after whatever
    FASTA it happened to be built from elsewhere.

    That name is derived from params.reference_fa, so it is checked here at launch instead
    of letting bwa-mem2 fail per-sample with a missing-index error hours into a run.
*/
def expected_index_prefix() {
    def stem = file(params.reference_fa).name.replaceFirst(/\.(fasta|fa)(\.gz)?$/, '')
    return "${stem}_prc.fa"
}

def prebuilt_bwa_index(index_dir) {
    // The five files `bwa-mem2 index` writes next to the FASTA it indexes.
    // CREATE_BWA_INDEX emits exactly this set, and a prebuilt index must match it.
    def extensions = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']

    def dir = file(index_dir, checkIfExists: true)
    if ( !dir.isDirectory() ) {
        error("--bwa_index must be a directory containing a bwa-mem2 index, but '${index_dir}' is a file.")
    }

    def prefix = expected_index_prefix()
    def wanted = extensions.collect { ext -> file("${dir}/${prefix}${ext}") }
    def missing = wanted.findAll { f -> !f.exists() }

    if ( missing ) {
        def present = dir.list().sort()
        error(
            "--bwa_index '${index_dir}' is not a usable bwa-mem2 index.\n\n" +
            "Missing:\n" + missing.collect { f -> "  ${f.name}" }.join('\n') + "\n\n" +
            "BWA_MAP uses the prepared reference FASTA as the bwa-mem2 index prefix, so the\n" +
            "index files must be named '${prefix}.<ext>'. An index built from a FASTA with a\n" +
            "different name will not be found, even if it is otherwise valid.\n\n" +
            "Directory contains: " + (present ? present.join(', ') : '(empty)') + "\n\n" +
            "To get a directory in the right shape, run once without --bwa_index and reuse\n" +
            "what CREATE_BWA_INDEX publishes to ./resources/bwa/."
        )
    }

    return wanted
}

workflow BWA_INDEX {

    take:
        reference_genome   // tuple(prepared_fasta, fai)
        bwa_index          // params.bwa_index: directory holding a prebuilt index, or null

    main:
        if ( bwa_index ) {
            // MUST be channel.value, not fromList/of.
            //
            // BWA_MAP is invoked once per sample against a queue channel of FASTQs, and
            // this index has to be readable by every one of those tasks. A value channel
            // is read without being consumed; a one-item queue channel is consumed by the
            // first task, so BWA_MAP would align only a single sample and silently skip
            // the rest.
            //
            // The auto branch works because a process output that emits exactly once is
            // treated as a value channel too - so both branches broadcast, which is what
            // makes them interchangeable here.
            bwa_index_ch = channel.value( prebuilt_bwa_index(bwa_index) )
        } else {
            bwa_index_ch = CREATE_BWA_INDEX(reference_genome).bwa_index
        }

    emit:
        bwa_index = bwa_index_ch
}
