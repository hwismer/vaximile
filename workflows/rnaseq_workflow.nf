include { CREATE_STAR_INDEX; STAR_ALIGN; STAR_SORT_INDEX_BAM; CREATE_KALLISTO_INDEX; KALLISTO_QUANT; KALLISTO_TXIMPORT; CREATE_SALMON_INDEX; SALMON_QUANT; GET_RNA_STRANDEDNESS } from "../modules/rnaseq.nf"

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
            star_index_ch = CREATE_STAR_INDEX(reference_genome, gtf)
        }

        // Align fastqs to genome
        star = STAR_ALIGN(fastqs, star_index_ch, gtf)
        star_sorted = STAR_SORT_INDEX_BAM(star.star_bam)

        // If null parameter, then generate kallisto index using transcriptome_fa
        if ( kallisto_index ) {
            kallisto_index_ch = Channel.fromPath(kallisto_index).collect()
        } else {
            kallisto_index_ch = CREATE_KALLISTO_INDEX(transcriptome_fa)
        }

        // Perform quantification with kallisto
        kallisto = KALLISTO_QUANT(fastqs, kallisto_index_ch)
        // Get gene-level abundance using tximport
        kallisto_gene_quant = KALLISTO_TXIMPORT(kallisto.abundance, gtf)

        // If null parameter, generate salmon index using transcriptome_fa
        if ( salmon_index ) {
            salmon_index_ch = Channel.fromPath(salmon_index).collect()
        } else {
            salmon_index_ch = CREATE_SALMON_INDEX(transcriptome_fa)
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

