include { CREATE_STAR_INDEX } from "../../../modules/local/create_star_index/main"
include { STAR_ALIGN } from "../../../modules/local/star_align/main"
include { STAR_SORT_INDEX_BAM } from "../../../modules/local/star_sort_index_bam/main"
include { CREATE_KALLISTO_INDEX } from "../../../modules/local/create_kallisto_index/main"
include { KALLISTO_QUANT } from "../../../modules/local/kallisto_quant/main"
include { KALLISTO_TXIMPORT } from "../../../modules/local/kallisto_tximport/main"
include { CREATE_SALMON_INDEX } from "../../../modules/local/create_salmon_index/main"
include { SALMON_QUANT } from "../../../modules/local/salmon_quant/main"
include { GET_RNA_STRANDEDNESS } from "../../../modules/local/get_rna_strandedness/main"

workflow RNASEQ_WORKFLOW {

    take:
        fastqs // (meta, fastq1, fastq2)
        reference_genome // (reference fasta, reference fastq index, reference dict)
        gtf // (gtf file .gz) GTF file corresponding to reference_genome and transcriptome_fa
        star_index // null or (/StarGenomeDir/)
        kallisto_index // null or kallisto index file
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

        // If null parameter, then generate kallisto index using transcriptome_fa
        if ( kallisto_index ) {
            kallisto_index_ch = Channel.fromPath(kallisto_index).collect()
        } else {
            kallisto_index_ch = CREATE_KALLISTO_INDEX(transcriptome_fa).index
        }

        // Perform quantification with kallisto
        kallisto = KALLISTO_QUANT(fastqs, kallisto_index_ch)
        // Get gene-level abundance using tximport
        kallisto_gene_quant = KALLISTO_TXIMPORT(kallisto.abundance, gtf)

        // If null parameter, generate salmon index using transcriptome_fa
        if ( salmon_index ) {
            salmon_index_ch = Channel.fromPath(salmon_index).collect()
        } else {
            salmon_index_ch = CREATE_SALMON_INDEX(transcriptome_fa).dir
        }
        // Perform quantification and strandedness prediction with salmon
        salmon = SALMON_QUANT(fastqs, salmon_index_ch)

        // Get the predicted RNA strandedness from salmon quant
        rna_strand = GET_RNA_STRANDEDNESS(salmon.quant)

        rna_strand_predictions = rna_strand.strand_txt
            .map { meta, pred ->
                tuple(meta, pred.text.trim())
            }

    emit:
        star_bam = star_sorted.bam
        star_chimeric_out = star.chimeric_out
        star_gene_quant = star.gene_quant
        star_final_log = star.final_log
        kallisto_tx = kallisto.abundance
        kallisto_gene = kallisto_gene_quant.gene_abundance
        salmon_tx = salmon.quant
        rna_strand = rna_strand_predictions
}
