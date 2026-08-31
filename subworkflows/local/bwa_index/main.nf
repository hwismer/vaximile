include { CREATE_BWA_INDEX } from "../../../modules/local/create_bwa_index/main"

/*
    BWA_MAP passes the reference FASTA to minibwa as the index prefix:

        minibwa map ... $reference_fa $fastq1 $fastq2

    so minibwa looks for "<reference_fa>.l2b" and "<reference_fa>.mbw". A
    prebuilt index is therefore only usable if its files are named after the PREPARED
    reference - the "<stem>_prc.fa" that PREPARE_FASTA writes - rather than after whatever
    FASTA it happened to be built from elsewhere.

    That name is derived from params.reference_fa, so it is checked here at launch instead
    of letting minibwa fail per-sample with a missing-index error hours into a run.
*/
def expected_index_prefix() {
    def stem = file(params.reference_fa).name.replaceFirst(/\.(fasta|fa)(\.gz)?$/, '')
    return "${stem}_prc.fa"
}

def prebuilt_bwa_index(index_dir) {
    // The two files `minibwa index` writes next to the FASTA it indexes: .l2b holds the
    // 2-bit encoded reference, .mbw the BWT and sampled suffix array.
    // CREATE_BWA_INDEX emits exactly this set, and a prebuilt index must match it.
    def extensions = ['.l2b', '.mbw']

    def dir = file(index_dir, checkIfExists: true)
    if ( !dir.isDirectory() ) {
        error("--bwa_index must be a directory containing a minibwa index, but '${index_dir}' is a file.")
    }

    def prefix = expected_index_prefix()
    def wanted = extensions.collect { ext -> file("${dir}/${prefix}${ext}") }

    // Distinguish "not there" from "there but unusable". A dangling symlink lists by name
    // but fails exists(), so reporting it as merely missing produces the contradictory
    // message of naming a file the very next line shows in the directory.
    def dangling = wanted.findAll { f ->
        !f.exists() && java.nio.file.Files.exists(f, java.nio.file.LinkOption.NOFOLLOW_LINKS)
    }
    def absent = wanted.findAll { f ->
        !f.exists() && !java.nio.file.Files.exists(f, java.nio.file.LinkOption.NOFOLLOW_LINKS)
    }

    if ( dangling || absent ) {
        def present = dir.list().sort()
        def report = "--bwa_index '${index_dir}' is not a usable minibwa index.\n\n"

        // The upgrade case: a directory left over from when this pipeline used bwa-mem2.
        // It is a perfectly valid index, just for the wrong aligner, so the generic
        // "named after the wrong FASTA" advice below would send the reader hunting for a
        // naming bug that is not there. Detect it by its signature extensions and say so.
        def looks_like_bwa_mem2 = present.any { n -> n.endsWith('.bwt.2bit.64') || n.endsWith('.0123') }
        if ( looks_like_bwa_mem2 && !dangling ) {
            error(
                "--bwa_index '${index_dir}' holds a bwa-mem2 index, not a minibwa index.\n\n" +
                "This pipeline switched from bwa-mem2 to minibwa. minibwa uses two files\n" +
                "('.l2b', '.mbw') in place of bwa-mem2's '.0123'/'.amb'/'.ann'/'.bwt.2bit.64'/'.pac',\n" +
                "so an index built by the older pipeline cannot be reused and must be rebuilt.\n\n" +
                "Rebuild it by running once without --bwa_index; the new index is published to\n" +
                "./resources/bwa/. Note that minibwa alignments are not identical to bwa-mem2's,\n" +
                "so BAMs produced before and after the switch should not be mixed in one cohort.\n\n" +
                "Directory contains: " + present.join(', ')
            )
        }
        if ( dangling ) {
            report += "Broken symlinks (present, but their target is gone):\n" +
                dangling.collect { f -> "  ${f.name}" }.join('\n') + "\n\n" +
                "This is what an index published by an older run looks like after work/ has\n" +
                "been cleaned: publishDir used to default to symlinks. Rebuild the index by\n" +
                "running once without --bwa_index; it is now copied rather than linked.\n\n"
        }
        if ( absent ) {
            report += "Missing:\n" + absent.collect { f -> "  ${f.name}" }.join('\n') + "\n\n" +
                "BWA_MAP uses the prepared reference FASTA as the minibwa index prefix, so\n" +
                "the index files must be named '${prefix}.<ext>'. An index built from a FASTA\n" +
                "with a different name will not be found, even if it is otherwise valid.\n\n"
        }
        report += "Directory contains: " + (present ? present.join(', ') : '(empty)')
        error(report)
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
