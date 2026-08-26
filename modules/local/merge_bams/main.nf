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
