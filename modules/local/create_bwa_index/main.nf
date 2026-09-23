process CREATE_BWA_INDEX {

    /*
        Create a minibwa index for minibwa map mapping.
    */

    label 'process_high_memory'
    conda "bioconda::minibwa=0.7"
    container "quay.io/biocontainers/minibwa:0.7--h118bc1c_0"

    tag "Creating BWA index for $reference_fa"

    // Copy, not symlink, so the published index survives cleaning work/.
    publishDir "./resources/bwa/", mode: "copy"

    cache 'lenient'

    input:
        tuple path(reference_fa), path(reference_index)

    output:
        path "*{.l2b,.mbw}", emit: bwa_index

    script:
        // minibwa index needs ~18x the genome size in RAM, so ~56 GB for GRCh38 - within
        // this tier's 96 GB, and less than bwa-mem2 index required. Add -l here if a
        // smaller node is all that is available; it trades build speed for memory.
        """
        minibwa index -t $task.cpus $reference_fa
        """

    stub:
        """
        touch ${reference_fa}.l2b
        touch ${reference_fa}.mbw
        """
}
