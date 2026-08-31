include { SAMTOOLS_SORMADUP } from "../../../modules/nf-core/samtools/sormadup/main"
include { INDEX_BAM         } from "../../../modules/local/index_bam/main"

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
