process CREATE_BWA_INDEX {

    /*
        Create a minibwa index for minibwa map mapping.
    */

    label 'process_high_memory'
    conda "bioconda::minibwa=0.7"
    container "quay.io/biocontainers/minibwa:0.7--h118bc1c_0"

    tag "Creating BWA index for $reference_fa"

    publishDir "./resources/bwa/", mode: "copy"

    cache 'lenient'

    input:
        tuple path(reference_fa), path(reference_index)

    output:
        path "*{.l2b,.mbw}", emit: bwa_index

    script:
        """
        minibwa index -t $task.cpus $reference_fa
        """

    stub:
        """
        touch ${reference_fa}.l2b
        touch ${reference_fa}.mbw
        """
}
