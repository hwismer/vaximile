process CREATE_BWA_INDEX {

    /*
        Create a minibwa index for minibwa map mapping.
    */

    label 'process_high_memory'
    conda "bioconda::minibwa=0.7"
    container "quay.io/biocontainers/minibwa:0.7--h118bc1c_0"

    tag "Creating BWA index for $reference_fa"

    // mode: "copy" - without it Nextflow defaults to "symlink", so ./resources/bwa/ ends
    // up holding links into work/. Those dangle as soon as work/ is cleaned, and a later
    // run passing --bwa_index at this directory then sees files that list but do not
    // resolve. The other three index builders already copy.
    publishDir "./resources/bwa/", mode: "copy"

    cache 'lenient'

    input:
        tuple path(reference_fa), path(reference_index)

    output:
        path "*{.l2b,.mbw}", emit: bwa_index
        path "versions.yml", topic: versions

    script:
        // minibwa index needs ~18x the genome size in RAM, so ~56 GB for GRCh38 - within
        // this tier's 96 GB, and less than bwa-mem2 index required. Add -l here if a
        // smaller node is all that is available; it trades build speed for memory.
        """
        minibwa index -t $task.cpus $reference_fa
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minibwa: \$(minibwa version)
        END_VERSIONS
        """

    stub:
        """
        touch ${reference_fa}.l2b
        touch ${reference_fa}.mbw
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minibwa: 0.7
        END_VERSIONS
        """
}
