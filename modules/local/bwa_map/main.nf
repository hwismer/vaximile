process BWA_MAP {

    // Align reads with minibwa and pipe straight to BAM. Read groups come from the samplesheet metadata.
    // pipefail is set explicitly so a minibwa crash can't produce a truncated BAM.
    // No container: no published image has both minibwa and samtools.

    label 'process_very_high'
    conda "bioconda::minibwa=0.7 bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "BWA Alignment on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), val(sequencing_type), path(fastq1), path(fastq2)
        tuple path(reference_fa), path(reference_index)
        path bwa_index

    output:
        tuple val(meta), path("*.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
    def bam_threads = Math.max(1, task.cpus.intdiv(4))
    """
    set -euo pipefail

    NEW_RG="@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:${molecule}_${sequencing_type}"

    minibwa map -t $task.cpus -R \$NEW_RG $reference_fa $fastq1 $fastq2 \\
        | samtools view -1 -@ ${bam_threads} -o "${prefix}.bam" -

    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
    """
    touch ${prefix}.bam
    """
}
