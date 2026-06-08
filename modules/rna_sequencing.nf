process STAR_ALIGN {

    /*

    Align RNA reads with STAR. Parameters included from star-fusion to be able to use the output of this process
    in a downstream star-fusion or arriba process without having to re-map.

    */

    cpus 16

    memory "48GB"

    container "alexdobin/star:2.7.10a_alpha_220506"

    tag "Aligning ${meta.sample_name} with STAR"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(star_index_dir)
        path(gtf)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Aligned.out.bam"), emit: star_bam
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_ReadsPerGene.out.tab"), emit: gene_quant
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Log.final.out"), emit:final_log
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_SJ.out.tab"), emit: sj_out
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_Chimeric.out.junction"), path(fastq1), path(fastq2), emit: chimeric_out
        tuple val(meta), path("*"), emit: tutto
    script:
        """

        gzip -d -c $gtf > gencode.gtf

        STAR \
            --runThreadN $task.cpus \
            --genomeDir $star_index_dir \
            --readFilesIn $fastq1 $fastq2 \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped None \
            --twopassMode Basic \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --outFilterMultimapNmax 50 \
            --chimOutJunctionFormat 1 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:"${meta.sample_name}" SM:"${meta.sample_name}" \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --chimNonchimScoreDropMin 10 \
            --chimOutType Junctions WithinBAM HardClip \
            --chimScoreDropMax 30 \
            --peOverlapNbasesMin 10 \
            --peOverlapMMp 0.1 \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSplicedMateMapLmin 30 \
            --outFileNamePrefix ./${meta.sample_name}_${meta.molecule}_ \
            --quantMode GeneCounts \
            --sjdbGTFfile gencode.gtf

        """

}


process CREATE_STAR_INDEX {

    /*
    
    Use a reference fasta and gtf to create a star index for star 2.7.10

    */

    cpus 24
    memory "64GB"
    cache 'lenient'

    container "alexdobin/star:2.7.10a_alpha_220506"

    tag "Creating STAR index with ${reference_fa} and ${gtf}"

    publishDir "./resources/star/", mode: "copy"

    input:
        tuple path(reference_fa), path(reference_index)
        path(gtf)

    output:
        path("STARGenomeDir"), type: "dir", emit: star_index

    script:
        """
        gzip -d -c $gtf > gtf.gtf

        STAR \
            --runThreadN $task.cpus \
            --runMode genomeGenerate \
            --genomeDir STARGenomeDir \
            --genomeFastaFiles $reference_fa \
            --sjdbGTFfile gtf.gtf

        """


}

process STAR_SORT_INDEX_BAM {

    /*

        Index the BAM file from a star process.

    */

    cpus 8
    memory "32GB"
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_STAR_sorted.bam.bai"), emit: bam

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_STAR_sorted.bam"
        """

}

process KALLISTO_TXIMPORT {

    cpus 1
    memory "16GB"
    
    conda "conda-forge::r-base=4.4.3 bioconda::bioconductor-tximport=1.34.0 conda-forge::r-readr=2.1.6 conda-forge::r-dplyr=1.1.4 bioconda::bioconductor-rtracklayer=1.66.0"

    tag "Getting gene abundance for ${meta.sample_name}"

    input:
        tuple val(meta), path(kallisto_abundance)
        path gtf

    output:
        tuple val(meta), path("${meta.sample_name}.gene_tpm.tsv"), emit: gene_abundance

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

    cpus 8
    memory "32GB"

    conda "bioconda::kallisto=0.51.1"
    
    tag "Kallisto quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(read1), path(read2)
        path kallisto_index

    output:
        tuple val(meta), path("${meta.sample_name}_kallisto/abundance.tsv"), emit: abundance
        tuple val(meta), path("${meta.sample_name}_kallisto"), emit: kallisto_dir
        tuple val(meta), path("*"), emit: tutto

    script:
        """
        kallisto quant -i $kallisto_index -o ${meta.sample_name}_kallisto -t ${task.cpus} $read1 $read2 > "${meta.sample_name}_kallist_stdout.out" 
        """
}


process CREATE_KALLISTO_INDEX {

    /*

    Create a kallisto index using a reference genome.

    */
    
    cpus 8
    memory "32GB"
    cache 'lenient'
    
    publishDir "./resources/kallisto", mode: "copy"

    conda "bioconda::kallisto=0.51.1"

    tag "Creating kallisto index on ${transcriptome_fa}"

    input:
        path(transcriptome_fa)

    output:
        path("${transcriptome_fa}_kallisto_index.idx")

    script:
        """
        kallisto index -i "${transcriptome_fa}_kallisto_index.idx" $transcriptome_fa
        """

}

process CREATE_SALMON_INDEX {

    cpus 8
    memory "32GB"
    
    conda "bioconda::salmon=1.11.4"
    
    publishDir "./resources/salmon/", mode: "copy"

    tag "Creating salmon index on ${transcripts_fa}"

    input:
        path(transcripts_fa)

    output:
        path("${transcripts_fa}_salmon_index")

    script:
    """
    salmon index \
        -t $transcripts_fa \
        -i "${transcripts_fa}_salmon_index" \
        -k 31
    """


}
process SALMON_QUANT {

    cpus 12
    memory "32GB"

    conda "bioconda::salmon=1.11.4"

    tag "Running salmon quant on ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        path(salmon_index)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_salmon_quant"), emit: quant

    script:
    """
    salmon quant \
        -i $salmon_index \
        --libType A \
        -1 $fastq1 \
        -2 $fastq2 \
        --validateMappings \
        -p $task.cpus \
        -o "${meta.sample_name}_${meta.molecule}_salmon_quant"
    """

}

process GET_RNA_STRANDEDNESS {

    cpus 1
    memory "4GB"

    conda "python=3.10 pandas=2.1"

    tag "Predicting RNA strandedness on ${meta.sample_name}"
    
    input:
        tuple val(meta), path(salmon_quant)

    output:
        tuple val(meta), path("${meta.sample_name}_strandedness.txt"), emit: strand_txt
    
    script:
    """
    #!/usr/bin/env python3

    import json
    import sys
    
    with open("${salmon_quant}/lib_format_counts.json", "r") as f:
        data = json.load(f)

    expected_format = data.get("expected_format")
    print(expected_format)
    if expected_format[1] == "U":
        strandedness = "XS"
    elif expected_format[1] == "S":
        if expected_format[2] == "R":
            strandedness = "RF"
        elif expected_format[2] == "F":
            strandedness = "FR"
    else:
        print("Unable to parse strandedness. Check salmon output")
        sys.exit(1)
        
    with open("${meta.sample_name}_strandedness.txt", "w") as f:
        f.write(strandedness + "\\n")

    """


}
