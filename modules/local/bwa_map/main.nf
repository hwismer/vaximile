process BWA_MAP {

    /*
        Map fastq files using minibwa. Outputs a SAM file.
        Reads groups are created using metadata information and currently are basically just the same name.
        Creates read group solely based on provided metadata from samplesheet. Any readgroup information
        present in the FASTQs is ignored.

        minibwa takes -R in the same '@RG\tID:foo\tSM:bar' form bwa-mem2 did, so the read
        group string below is unchanged. Note that minibwa uses a different algorithm from
        bwa-mem2 and its alignments are NOT identical to it.
    */

    label 'process_high'
    conda "bioconda::minibwa=0.7"
    container "quay.io/biocontainers/minibwa:0.7--h118bc1c_0"

    tag "BWA Alignment on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), val(sequencing_type), path(fastq1), path(fastq2)
        tuple path(reference_fa), path(reference_index)
        path bwa_index

    output:
        tuple val(meta), path("*.sam"), emit: sam
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
    """
    NEW_RG="@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:${molecule}_${sequencing_type}"

    minibwa map -t $task.cpus -R \$NEW_RG $reference_fa $fastq1 $fastq2 > "${prefix}.sam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
    """
    touch ${prefix}.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: 0.7
    END_VERSIONS
    """
}
