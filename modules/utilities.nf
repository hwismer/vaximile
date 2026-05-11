process MULTIQC {

    cpus 2
    memory "16GB"
    conda "bioconda::multiqc=1.34-0"

    tag "Running MultiQC on ${somatic_name}"

    publishDir "./multiqc/"

    input:
        tuple val(somatic_name), path(files)

    output:
        tuple val(somatic_name), path("${somatic_name}_multiqc_report.html")

    script:
    """
    multiqc \
        -n ${somatic_name}_multiqc_report.html \
        .

    """





}


process COMBINE_FASTQS {

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

process SORT_BAM {

    /*

        Index the BAM file from a star process.

    */

    cpus 16
    memory "32GB"

    container "biocontainers/samtools:v1.9-4-deb_cv1"

    tag "Sorting ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_sorted.bam"), path("${meta.sample_name}_${meta.molecule}_sorted.bam.bai")

    script:
        """
        samtools sort --threads $task.cpus  $bam -o "${meta.sample_name}_${meta.molecule}_sorted.bam"
        samtools index -@ $task.cpus  "${meta.sample_name}_${meta.molecule}_sorted.bam"
        """

}
process SPLIT_INTERVALS {

    /*

        Given a file of genomic intervals, split into scatter_count number of shards.

    */

    cpus 2
    memory "8GB"
    cache "lenient"
    
    tag "Splitting ${intervals_file} into ${scatter_count} shards w/ ${interval_padding} bp padding"

    container "broadinstitute/gatk:4.6.1.0"

    input:
        tuple path(reference_fa), path(reference_fa_index), path(reference_dict)
        path intervals_file
        val scatter_count
        val interval_padding

    output:
         path("*-scattered.interval_list"), emit: interval_shards

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


process PULL_ARRIBA_RESOURCES {

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

process INDEX_VCF {
    /*

    tabix index a vcf file

    */

    cpus 2
    memory "8GB"
    
    container "staphb/bcftools:1.23"

    publishDir "./testout/"

    input:
        tuple val(vcf_name), val(meta), path(vcf)
        val(filename_suffix)

    output:
        tuple val(meta), path("${vcf_name}_${filename_suffix}.vcf.gz"), path("${vcf_name}_${filename_suffix}.vcf.gz.tbi")

    script:
        """
        bcftools view $vcf -Oz -o "${vcf_name}_${filename_suffix}.vcf.gz"
        bcftools index -t ${vcf_name}_${filename_suffix}.vcf.gz
        """

}

process BED_BGZIP_INDEX {

    cpus 2
    memory "8GB"

    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    input:
        path(bed)

    output:
        tuple path("${bed.baseName}_sorted.bed.gz"), path("${bed.baseName}_sorted.bed.gz.tbi")

    script:
    """
    bedtools sort -i $bed > ${bed.baseName}_sorted.bed
    bgzip -@ $task.cpus -c ${bed.baseName}_sorted.bed > ${bed.baseName}_sorted.bed.gz
    tabix -p bed ${bed.baseName}_sorted.bed.gz
    """

}
