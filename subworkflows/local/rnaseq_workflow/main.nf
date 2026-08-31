include { CREATE_STAR_INDEX } from "../../../modules/local/create_star_index/main"
include { STAR_ALIGN } from "../../../modules/local/star_align/main"
include { STAR_SORT_INDEX_BAM } from "../../../modules/local/star_sort_index_bam/main"
include { CREATE_SALMON_INDEX } from "../../../modules/local/create_salmon_index/main"
include { SALMON_QUANT } from "../../../modules/local/salmon_quant/main"
include { GET_RNA_STRANDEDNESS } from "../../../modules/local/get_rna_strandedness/main"

workflow RNASEQ_WORKFLOW {

    take:
        fastqs // (meta, fastq1, fastq2)
        reference_genome // (reference fasta, reference fastq index, reference dict)
        gtf // (gtf file .gz) GTF file corresponding to reference_genome and transcriptome_fa
        star_index // null or (/StarGenomeDir/)
        salmon_index // null or salmon index file
        transcriptome_fa // (transcriptome fasta file) 

    main:
        
        // If null parameter, then generate star index from reference_genome, otherwise use existing
        if ( star_index ) {
            star_index_ch = Channel.fromPath(star_index).collect()
        } else {
            star_index_ch = CREATE_STAR_INDEX(reference_genome, gtf).star_index
        }

        // Align fastqs to genome
        star_align_input = fastqs
            .map { meta, fastq1, fastq2 ->
                tuple(meta, meta.sample_name, fastq1, fastq2)
            }

        star = STAR_ALIGN(star_align_input, star_index_ch, gtf)
        star_sorted = STAR_SORT_INDEX_BAM(star.star_bam)

        // If null parameter, generate salmon index using transcriptome_fa
        if ( salmon_index ) {
            salmon_index_ch = Channel.fromPath(salmon_index).collect()
        } else {
            salmon_index_ch = CREATE_SALMON_INDEX(transcriptome_fa).dir
        }
        // Quantification, strandedness prediction, and - via --geneMap, which takes the
        // GTF - gene-level aggregation, all from the one salmon run. This replaced kallisto
        // plus a tximport step that produced the same gene TPMs.
        salmon = SALMON_QUANT(fastqs, salmon_index_ch, gtf)

        // Get the predicted RNA strandedness from salmon quant
        rna_strand = GET_RNA_STRANDEDNESS(salmon.quant)

        rna_strand_predictions = rna_strand.strand_txt
            .map { meta, pred ->
                tuple(meta, pred.text.trim())
            }

    emit:
        star_bam = star_sorted.bam
        star_chimeric_out = star.chimeric_out
        star_final_log = star.final_log
        salmon_tx = salmon.tsv            // quant.sf, transcript-level
        salmon_gene = salmon.genes_tsv    // quant.genes.sf, gene-level
        salmon_dir = salmon.quant         // whole run directory, for MultiQC
        rna_strand = rna_strand_predictions
}
