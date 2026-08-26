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
