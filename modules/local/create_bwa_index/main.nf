process CREATE_BWA_INDEX {

    /*
        Create a bwa-mem2 index for bwa mem mapping.
    */

    cpus 8
    memory "80GB"
    container "iarcbioinfo/bwa-mem2-tools:v1.0"

    tag "Creating BWA index for $reference_fa"

    publishDir "./resources/bwa/"

    cache 'lenient'

    input:
        tuple path(reference_fa), path(reference_index)

    output:
        path "*{.bwt.2bit.64,.sa,.pac,.amb,.ann,.0123}", emit: bwa_index

    script:
        """
        bwa-mem2 index $reference_fa
        """
}
