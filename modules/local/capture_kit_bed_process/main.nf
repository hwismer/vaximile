process CAPTURE_KIT_BED_PROCESS {

    // Tools like Strelka require a bgzipped and tbi indexes BED file when running with specified regions.
    // Take a BED file and bgzip and TBI index it

    label 'process_low'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"

    tag "Preprocessing $kit BED file $bed"

    input:
        tuple val(kit), path(bed)

    output:
        tuple val(kit), path("*_sorted.bed.gz"), path("*_sorted.bed.gz.tbi"), emit: bed
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    """
    bedtools sort -i $bed > ${prefix}_sorted.bed
    bgzip -@ $task.cpus -c ${prefix}_sorted.bed > ${prefix}_sorted.bed.gz
    tabix -p bed ${prefix}_sorted.bed.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    """
    touch ${prefix}_sorted.bed.gz
    touch ${prefix}_sorted.bed.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """

}

//*******************************************************************************************************************
// REFERENCE FASTA PREPARATION
