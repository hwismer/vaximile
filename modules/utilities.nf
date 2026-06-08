//***************************************************************************************************************************
// FILE COMBINING (BAMS, FASTQS)

process MERGE_BAMS {

    /*
    Merges an arbitrary number of bam files with the same metadata.
    */
    
    cpus 6
    memory "18GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Merging Bams from ${meta.sample_name}"
    input:
        tuple val(meta), path(bams), path(bais)

    output:
       tuple val(meta), path("${meta.sample_name}.bam"), path("${meta.sample_name}.bam.bai")

    script:
	"""
    set -euo pipefail

    samtools merge \
        -@ ${task.cpus} \
        -f \
        ${meta.sample_name}.bam \
        ${bams.join(' ')}

    samtools index \
        -@ ${task.cpus} \
        ${meta.sample_name}.bam
    """
}



process COMBINE_FASTQS {
        
    // Combines FASTQS

    cpus 2
    memory "8GB"
    
    tag "Combing FASTQs from ${meta.sample_name}"

    input:
        tuple val(meta), path(fastqs_r1), path(fastqs_r2)
    output:
        tuple val(meta), path("${meta.somatic_name}_R1.merged.fastq.gz"), path("${meta.somatic_name}_R2.merged.fastq.gz")


    script:
    """
    cat ${fastqs_r1.join(' ')} > ${meta.somatic_name}_R1.merged.fastq.gz
    cat ${fastqs_r2.join(' ')} > ${meta.somatic_name}_R2.merged.fastq.gz
    """


}

//***************************************************************************************************************************
// INTERVAL PROCESSING

process SPLIT_INTERVALS {

    /*
        Takes a file of genomic intervals such as a picard interval file or BED file splits into scatter_count number of shards for parallel processing.

    */

    cpus 2
    memory "8GB"
    cache "lenient"
    
    tag "Splitting ${intervals_file} into ${scatter_count} shards w/ ${interval_padding} bp padding"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple path(reference_fa), path(reference_fa_index)
        path reference_dict
        tuple val(capture_kit), path(intervals_file)
        val scatter_count
        val interval_padding

    output:
         tuple val(capture_kit), path("*-scattered.interval_list"), emit: interval_shards

    script:
    """
    gatk SplitIntervals \
        -R $reference_fa \
        -L $intervals_file \
        --scatter-count $scatter_count \
        --interval-padding $interval_padding \
        -O .

    """
}

//***************************************************************************************************************************
// RESOURCE PULLING FROM WEB SOURCES

process PULL_ARRIBA_RESOURCES {

    // Pulls arribra resources from release 2.5.1. Runs locally incase job nodes don't have internet.

    cpus 1
    memory "4GB"
    executor "local"

    tag "Pulling Arriba resources v2.5.1"

    output:
        tuple path("./arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz"), 
            path("./arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz"),
            path("./arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3"), emit: resources
    
    script:

        
    """
        wget https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz
        ls
        echo "done"
        tar -xzf arriba_v2.5.1.tar.gz
        ls

    """

}

process PULL_VEP_PVAC_PLUGINS {

    // Pulls the VEP plugins necessary to run pvactools. Runs locally to ensure internet connection.

    cpus 1
    memory "4GB"
    executor "local"

    container "griffithlab/pvactools:6.0.3"

    tag "Pulling Frameshift and Wildtype VEP plugins"

    output:
        path("VEP_plugins"), emit: plugins

    script:
    """
    git clone https://github.com/Ensembl/VEP_plugins.git

    pvacseq install_vep_plugin VEP_plugins

    """

}

process PULL_CTAT_RESOURCE_BUNDLE {

    // Pulls the hg38 CTAT resource bundle needed for STARfusion and other tools.

    cpus 1
    memory "8GB"
    executor "local"

    tag "Pulling CTAT plug-n-play resource bundle"

    output:
        path("./GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir"), emit: ctat_resource_dir

    script:
    """
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
    tar -xzf GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz

    """

}


//***************************************************************************************************************************
// FILE INDEXING OPERATIONS

process INDEX_VCF {
    
    // TBI index a vcf file


    cpus 2
    memory "8GB"
    container "staphb/bcftools:1.23"

    tag "Indexing $vcf"

    input:
        tuple val(vcf_name), val(meta), path(vcf) // vcf name should be a string that corresponds to the file name ie. Sample1
        val(filename_suffix) // suffix corresponds to name after vcf_name ie somatic_variants -> Sample1_somatic_variants.vcf.gz

    output:
        tuple val(meta), path("${vcf_name}_${filename_suffix}.vcf.gz"), path("${vcf_name}_${filename_suffix}.vcf.gz.tbi")

    script:
        """
        bcftools view $vcf -Oz -o "${vcf_name}_${filename_suffix}.vcf.gz"
        bcftools index -t ${vcf_name}_${filename_suffix}.vcf.gz
        """

}

process SORT_BAM {

    /*
    Sorts a BAM file and indexes it.
    */

    cpus 8
    memory "24GB"

    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "Sorting ${meta.sample_name}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai")

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_sorted.bam"
        """

}
process CAPTURE_KIT_BED_PROCESS {

    // Tools like Strelka require a bgzipped and tbi indexes BED file when running with specified regions.
    // Take a BED file and bgzip and TBI index it

    cpus 2
    memory "8GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Preprocessing $kit BED file $bed"

    input:
        tuple val(kit), path(bed)

    output:
        tuple val(kit), path("${bed.baseName}_sorted.bed.gz"), path("${bed.baseName}_sorted.bed.gz.tbi")

    script:
    """
    bedtools sort -i $bed > ${bed.baseName}_sorted.bed
    bgzip -@ $task.cpus -c ${bed.baseName}_sorted.bed > ${bed.baseName}_sorted.bed.gz
    tabix -p bed ${bed.baseName}_sorted.bed.gz
    """

}

//*******************************************************************************************************************
// REFERENCE FASTA PREPARATION

process PREPARE_FASTA {

    // Unzips a fasta.gz or fa.gz or otherwise renames to a standard name

	cpus 2
	memory "8GB"

    tag "Preprocessing $fasta"

    input:
    	path fasta

    output:
    	path "*_prc.fa", emit: fasta

    script:

    	def is_gz = fasta.name.endsWith('.gz')
        def prefix = fasta.name.replaceFirst(/\.(fasta|fa)(\.gz)?$/, '')

    	"""
    	set -euo pipefail

    	if ${is_gz}; then
        	gunzip -c ${fasta} > ${prefix}_prc.fa
    	else
        	cp ${fasta} ${prefix}_prc.fa
    	fi
    	"""
}

process INDEX_FASTA {

    // Indexes a fasta file with samtools faidx

    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Indexing $fasta"

    input:
    	path fasta

    output:
    	tuple path(fasta), path("${fasta}.fai"), emit: fai

    script:
        """
        set -euo pipefail
    	samtools faidx $fasta
        """
}

process MAKE_FASTA_DICT {

    // Generate picard fasta.dict file for use with GATK tools

    cpus 4
    memory "16GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Creating reference dict for $fasta"

    input:
        tuple path(fasta), path(fai)

    output:
        path("${fasta.baseName}.dict"), emit: dict

    script:
    """
    gatk CreateSequenceDictionary \
        R=${fasta} \
        O=${fasta.baseName}.dict
    """
}

//*******************************************************************************************************************

