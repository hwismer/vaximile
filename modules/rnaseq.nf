process KALLISTO_TXIMPORT {

    cpus 1
    memory "16GB"
    
    conda "conda-forge::r-base=4.4.3 bioconda::bioconductor-tximport=1.34.0 conda-forge::r-readr=2.1.6 conda-forge::r-dplyr=1.1.4 bioconda::bioconductor-rtracklayer=1.66.0"

    publishDir "${params.outdir}/${meta.somatic_sample}/rnaseq/", mode: "copy"

    input:
        tuple val(meta), path(kallisto_abundance)
        path gtf

    output:
        tuple val(meta), path("${meta.sample_name}.gene_tpm.tsv")

    script:
    """
    #!/usr/bin/env Rscript
    
    library(rtracklayer)
    library(dplyr)
    library(readr)
    library(tximport)
    
    
    gtf <- rtracklayer::import("${gtf}", format = "gtf")
    gtf_df=as.data.frame(gtf)
    tx2gene <- gtf_df[,c("transcript_id","gene_id", "gene_name")]

    
    
    kallisto_tsv <- tximport("${kallisto_abundance}", 
                         type = "kallisto", abundanceCol = "tpm",
                         countsFromAbundance = "no",
                         tx2gene = tx2gene, ignoreAfterBar = T)


    gene_tpm <- as.data.frame(kallisto_tsv\$abundance)
    gene_tpm\$gene_id <- rownames(gene_tpm)

    gene_tpm <- gene_tpm %>%
      select(gene_id, everything())

    gene_annot <- tx2gene %>%
        select(gene_id, gene_name) %>%
            distinct()

    gene_tpm <- left_join(gene_tpm, gene_annot, by = "gene_id")
    colnames(gene_tpm) <- c("ENSEMBLID", "TPM", "Gene")

    write_tsv(gene_tpm, "${meta.sample_name}.gene_tpm.tsv")
    
    """


}
process KALLISTO_QUANT {


    /*

        Get transcript abundance estimations using kallisto quant.

    */

    cpus 16
    memory "32GB"

    conda "bioconda::kallisto=0.51.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/rnaseq/kallisto/", mode: "copy"

    input:
        tuple val(meta), path(read1), path(read2)
        path kallisto_index

    output:
        tuple val(meta), path("${meta.sample_name}_kallisto/abundance.tsv"), emit: abundance
        tuple val(meta), path("${meta.sample_name}_kallisto"), emit: kallisto_dir

    script:
        """
        kallisto quant -i $kallisto_index -o ${meta.sample_name}_kallisto -t ${task.cpus} $read1 $read2 
        """
}


process KALLISTO_INDEX {

    /*

    Create a kallisto index using a reference genome.

    */
    
    cpus 32
    memory "32GB"
    cache 'lenient'
    

    conda "bioconda::kallisto=0.51.1"

    input:
        path(kallisto_reference)

    output:
        path("kallisto_index.idx")

    script:
        """
        kallisto index -i kallisto_index.idx $kallisto_reference

        """

}


