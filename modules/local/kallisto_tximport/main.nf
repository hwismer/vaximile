process KALLISTO_TXIMPORT {

    label 'process_low'

    conda "conda-forge::r-base=4.4.3 bioconda::bioconductor-tximport=1.34.0 conda-forge::r-readr=2.1.6 conda-forge::r-dplyr=1.1.4 bioconda::bioconductor-rtracklayer=1.66.0"

    tag "Getting gene abundance for ${meta.sample_name}"

    input:
        tuple val(meta), path(kallisto_abundance)
        path gtf

    output:
        tuple val(meta), path("*.gene_tpm.tsv"), emit: gene_abundance
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
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

    write_tsv(gene_tpm, "${prefix}.gene_tpm.tsv")
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.3
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch ${prefix}.gene_tpm.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.3
    END_VERSIONS
    """


}
