process CREATE_BWA_INDEX {

    /*
        Create a bwa-mem2 index for bwa mem mapping.
    */

    label 'process_high_memory'
    conda "bioconda::bwa-mem2=2.2.1"
    container "iarcbioinfo/bwa-mem2-tools:v1.0"

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
        path "*{.bwt.2bit.64,.sa,.pac,.amb,.ann,.0123}", emit: bwa_index
        path "versions.yml", topic: versions

    script:
        """
        bwa-mem2 index $reference_fa
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -1)
        END_VERSIONS
        """

    stub:
        """
        touch ${reference_fa}.bwt.2bit.64
        touch ${reference_fa}.sa
        touch ${reference_fa}.pac
        touch ${reference_fa}.amb
        touch ${reference_fa}.ann
        touch ${reference_fa}.0123
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa-mem2: 2.2.1
        END_VERSIONS
        """
}
