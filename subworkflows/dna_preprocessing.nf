/*
    dna preprocessing: dna_qc_workflow, bwa_index, dna_align_and_preproc, bam_markduplicates, somalier

    One file per pipeline step. Each workflow keeps the take/emit signature it had
    as its own subworkflow directory, so callers are unchanged.
*/
include { FASTP } from "../modules/local/fastp/main"
include { CREATE_BWA_INDEX } from "../modules/local/create_bwa_index/main"
include { BWA_MAP } from "../modules/local/bwa_map/main"
include { BASE_RECALIBRATOR_SCATTER } from "../modules/local/base_recalibrator_scatter/main"
include { BASE_RECALIBRATOR_GATHER } from "../modules/local/base_recalibrator_gather/main"
include { APPLY_BQSR_SCATTER } from "../modules/local/apply_bqsr_scatter/main"
include { APPLY_BQSR_GATHER } from "../modules/local/apply_bqsr_gather/main"
include { GET_PILEUP_SUMMARIES } from "../modules/local/get_pileup_summaries/main"
include { SAMTOOLS_FLAGSTAT } from "../modules/local/samtools_flagstat/main"
include { SAMTOOLS_COVERAGE } from "../modules/local/samtools_coverage/main"
include { SAMTOOLS_IDXSTATS } from "../modules/local/samtools_idxstats/main"
include { SAMTOOLS_SORMADUP } from "../modules/nf-core/samtools/sormadup/main"
include { INDEX_BAM } from "../modules/local/index_bam/main"
include { SOMALIER_EXTRACT } from "../modules/local/somalier_extract/main"
include { SOMALIER_RELATE } from "../modules/local/somalier_relate/main"

workflow DNA_QC_WORKFLOW {

    take:
        reads_channel
        
    main:
        fastp = FASTP(reads_channel)

    emit:
        fastp_fastqs = fastp.fastqs
        fastp_reports = fastp.reports
}

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

workflow DNA_ALIGN_AND_PREPROC {

    take:
        fastqs // (metamap, fastq1, fastq2)
        reference_genome // (reference fasta, reference fasta index)
        reference_dict // reference_dict file
        bwa_index // Either null if not supplied in main workflow or path to bwa index files
        known_sites_dbsnp // (vcf, vcf index)
        known_sites_1000g_snps // (vcf, vcf_index)
        known_indels // (vcf, vcf_index)
        mills // (vcf, vcf_index)
        common_germline // (vcf, vcf_index)
        intervals // (kit name, tuple(interval_files)) interval files created from the SPLIT_INTERVALS process
        num_intervals // Integer for how many intervals are present
        
    main:
        
        // Map with minibwa
        bwa_map_input = fastqs
            .map { meta, fastq1, fastq2 ->
                tuple(meta, meta.sample_name, meta.molecule, meta.sequencing_type, fastq1, fastq2)
            }

        bwa_bam = BWA_MAP(bwa_map_input, reference_genome, bwa_index).bam

        // **********************************************************
        // GATK PRE-PROCESSING BEST PRACTICES

        // Coordinate-sort and mark duplicates (samtools collate/fixmate/sort/markdup,
        // via the nf-core SAMTOOLS_SORMADUP module) in place of MarkDuplicatesSpark.
        markduplicates = BAM_MARKDUPLICATES(bwa_bam, reference_genome)
        mark_dup = markduplicates.bam
        markdup_metrics = markduplicates.metrics

        // BASE RECALIBRATION
        // Call BaseRecalibrator on each interval
        markdup_kit = mark_dup.map{meta, bam, bai ->
            tuple(meta.capture_kit, meta, bam, bai)
        }
        
        intervals_map = intervals.flatMap{ kit, interval_list ->
            interval_list.collect { interval ->
                tuple(kit, interval)
            }
        }
        

        base_recal_input = markdup_kit.combine(intervals_map, by: 0).map{kit, meta, bam, bai, interval -> tuple(meta, bam, bai, interval)}
    
        base_recal = BASE_RECALIBRATOR_SCATTER(
            base_recal_input,
            reference_genome,
            reference_dict,
            known_sites_dbsnp,
            known_sites_1000g_snps,
            known_indels,
            mills
        ).table
        base_recal_gather_input = base_recal.groupTuple(size: num_intervals)
            .map { meta, recal_tables ->
                tuple(meta, meta.sample_name, meta.molecule, recal_tables)
            }
        base_recal_gathered = BASE_RECALIBRATOR_GATHER(base_recal_gather_input).table
        
        // Scatter BQSR calls
        bqsr_input = mark_dup.join(base_recal_gathered)
        .map{meta, bam, bai, bqsr ->
            tuple(meta.capture_kit, meta, bam, bai, bqsr)
        }
        .combine(intervals_map, by:0 )
        .map{kit, meta, bam, bai, bqsr, interval -> tuple(meta, bam, bai, bqsr, interval) }

        bqsr = APPLY_BQSR_SCATTER(bqsr_input, reference_genome, reference_dict).bam
        bqsr_scattered = bqsr.groupTuple(size: num_intervals)
            .map { meta, bams ->
                tuple(meta, meta.sample_name, meta.molecule, bams)
            }
        bqsr_gather = APPLY_BQSR_GATHER(bqsr_scattered).bam // Get final BQSR bams
        // APPLY_BQSR_GATHER now emits (meta, bam, bai) directly: it merges the shards
        // rather than concatenating them, so the result is coordinate-sorted and indexed in
        // one streaming pass and there is no separate sort or index step to do here.
        bqsr_sort = bqsr_gather

        pileup_summaries = GET_PILEUP_SUMMARIES(bqsr_sort, common_germline).table
        
        flagstat = SAMTOOLS_FLAGSTAT(mark_dup).flagstat
        coverage = SAMTOOLS_COVERAGE(mark_dup).tsv
        idxstats = SAMTOOLS_IDXSTATS(mark_dup).tsv

        
    emit:
        preproc_bams = bqsr_sort // For somatic calling
        markdup_bams = mark_dup // Non-recalibrated BAMS for callers like Strelka that don't expect recalibrated scores
        base_recal = base_recal_gathered // Recalibration metrics for MultiQC report
        markdup_metrics = markdup_metrics // Duplicate stats from samtools markdup for MultiQC
        pileup_summaries = pileup_summaries
        flagstat = flagstat
        coverage = coverage
        idxstats = idxstats

}

/*
    Coordinate-sort the aligned BAM and mark duplicates, replacing MarkDuplicatesSpark.

    SAMTOOLS_SORMADUP is an unmodified nf-core module, so it expects nf-core conventions
    this pipeline does not otherwise follow, and this wrapper is what bridges them:

      - It tags tasks with `meta.id` and derives its default prefix from it. This
        pipeline's meta map has no `id`, so an `id` is added for the duration of the call
        and stripped again from the output.

        Stripping matters. The meta map is the join key in DNA_ALIGN_AND_PREPROC
        (`mark_dup.join(base_recal_gathered)`) and the grouping key for HLA typing and the
        somatic pairs. Leaving an extra key on it would silently fail to match the metas
        carried by every other channel.

      - Its second input is a `tuple(meta2, fasta, fai)`, where this pipeline passes a
        bare `tuple(fasta, fai)`, so a meta is prepended here.

      - It emits the BAM and index on separate channels and, with --write-index, produces
        a .csi. Downstream wants one `tuple(meta, bam, bai)`, so INDEX_BAM makes the .bai
        and re-joins the shape. See that module for why .csi will not do.
*/
workflow BAM_MARKDUPLICATES {

    take:
        aligned_bam        // tuple(meta, bam) from BWA_MAP, unsorted
        reference_genome   // tuple(prepared_fasta, fai)

    main:
        // `id` is what SAMTOOLS_SORMADUP tags with; ext.prefix in conf/modules.config
        // appends the _markdup suffix, keeping the output names the Spark module used.
        sormadup_input = aligned_bam.map { meta, bam ->
            tuple(meta + [id: "${meta.sample_name}_${meta.molecule}"], bam)
        }

        sormadup_reference = reference_genome.map { fasta, fai ->
            tuple([id: fasta.simpleName], fasta, fai)
        }

        sormadup = SAMTOOLS_SORMADUP(sormadup_input, sormadup_reference)

        // Restore the original meta before anything downstream keys on it.
        markdup_bam = sormadup.bam.map { meta, bam ->
            tuple(meta.findAll { k, _v -> k != 'id' }, bam)
        }

        indexed = INDEX_BAM(markdup_bam).bam

        markdup_metrics = sormadup.metrics.map { meta, metrics ->
            tuple(meta.findAll { k, _v -> k != 'id' }, metrics)
        }

    emit:
        bam     = indexed          // tuple(meta, bam, bai)
        metrics = markdup_metrics  // tuple(meta, metrics) - duplicate stats for MultiQC
}

workflow SOMALIER {

    take:
        sample_bams
        reference_genome
        somalier_sites_vcf

    main:
        extract = SOMALIER_EXTRACT(sample_bams, reference_genome, somalier_sites_vcf).somalier
        extract_unique = extract.unique { meta, files -> 
            meta.sample_name
        }
        extract_by_patient = extract_unique.map{meta, extracted ->
            tuple(meta.patient, meta, extracted)
        }
        .groupTuple()

        relate = SOMALIER_RELATE(extract_by_patient)
    
    emit:
        extract = extract
        pairs = relate.pairs
        samples = relate.samples
        groups = relate.groups
        html = relate.html
        

}
