process BWA_MAP {

    /*
        Map fastq files using minibwa. Outputs an unsorted BAM.
        Reads groups are created using metadata information and currently are basically just the same name.
        Creates read group solely based on provided metadata from samplesheet. Any readgroup information
        present in the FASTQs is ignored.

        minibwa takes -R in the same '@RG\tID:foo\tSM:bar' form bwa-mem2 did, so the read
        group string below is unchanged. Note that minibwa uses a different algorithm from
        bwa-mem2 and its alignments are NOT identical to it.

        minibwa only writes SAM, but SAMTOOLS_SORMADUP begins with `samtools cat`, which
        reads BAM and CRAM only and fails on SAM with "input is not BAM or CRAM".
        MarkDuplicatesSpark took SAM directly, so nothing needed this before. Converting in
        the pipe rather than in a following process means the SAM is never written at all -
        for WGS that removes a multi-hundred-GB intermediate, so this is faster than the
        bwa-mem2 pipeline was.

        -1 is fast BAM compression: this file is transient, read once by SORMADUP, so the
        cost of a higher level would not be repaid.

        The `set -euo pipefail` below duplicates the global process.shell in
        nextflow.config, deliberately. Without pipefail a minibwa crash mid-pipe would be
        masked by samtools exiting 0, shipping a silently truncated BAM. Nextflow invokes
        .command.sh as `env bash -C -e -u -o pipefail .command.sh`, so the shebang inside
        the file is never used - re-running it by hand while debugging loses those options.
        The explicit line keeps the pipe safe in that case too.

        No container. The published minibwa images (biocontainers, staphb) carry minibwa
        alone, and this now needs samtools in the same task. Under -profile docker or
        singularity this process will fall back to the host, as the other samtools-only
        modules in this pipeline already do. See docs/usage.md.
    */

    label 'process_very_high'
    conda "bioconda::minibwa=0.7 bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    tag "BWA Alignment on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), val(sequencing_type), path(fastq1), path(fastq2)
        tuple path(reference_fa), path(reference_index)
        path bwa_index

    output:
        tuple val(meta), path("*.bam"), emit: bam
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
    // minibwa is the bottleneck; give the compressor a slice rather than a second full set.
    def bam_threads = Math.max(1, task.cpus.intdiv(4))
    """
    set -euo pipefail

    NEW_RG="@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:${molecule}_${sequencing_type}"

    minibwa map -t $task.cpus -R \$NEW_RG $reference_fa $fastq1 $fastq2 \\
        | samtools view -1 -@ ${bam_threads} -o "${prefix}.bam" -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: 0.7
        samtools: 1.23.1
    END_VERSIONS
    """
}
