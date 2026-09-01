process BAMREADCOUNT {

    /*

    Run bamreadcount to add coverage information to a vcf using a BAM file. The sample names in the VCF
    must match the read groups and sample names in the BAM.

    */

    // bam-readcount is single-threaded and the helper script does not parallelise,
    // so the extra cores of process_high went unused.
    label 'process_low'
    cache "lenient"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    input:
        tuple val(somatic_name), val(somatic_meta), val(sample_name), val(molecule), path(vcf), val(sample_meta), path(bam), path(bai)
        tuple path(reference_fa), path(reference_index)

    output:

        tuple val(somatic_meta), val(sample_meta), path("*_bamrc_helper/*indel.tsv"), path("*_bamrc_helper/*snv.tsv"), emit: brc_files
        path "versions.yml", topic: versions
    script:
        def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
        """
        mkdir ${prefix}_bamrc_helper
        bam_readcount_helper.py \
            $vcf \
            ${sample_name} \
            $reference_fa \
            $bam \
            ${molecule} \
            ${prefix}_bamrc_helper
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bam-readcount: 1.2.1
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
        """
        mkdir -p ${prefix}_bamrc_helper
        touch ${prefix}_bamrc_helper/${prefix}_indel.tsv
        touch ${prefix}_bamrc_helper/${prefix}_snv.tsv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bam-readcount: 1.2.1
        END_VERSIONS
        """

}
