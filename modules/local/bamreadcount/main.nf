process BAMREADCOUNT {

    // Run bam-readcount at the VCF's sites. VCF sample names must match the BAM read groups.

    label 'process_low'
    cache "lenient"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    input:
        tuple val(somatic_name), val(somatic_meta), val(sample_name), val(molecule), path(vcf), val(sample_meta), path(bam), path(bai)
        tuple path(reference_fa), path(reference_index)

    output:

        tuple val(somatic_meta), val(sample_meta), path("*_bamrc_helper/*indel.tsv"), path("*_bamrc_helper/*snv.tsv"), emit: brc_files
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
        """

    stub:
        def prefix = task.ext.prefix ?: "${sample_name}_${molecule}"
        """
        mkdir -p ${prefix}_bamrc_helper
        touch ${prefix}_bamrc_helper/${prefix}_indel.tsv
        touch ${prefix}_bamrc_helper/${prefix}_snv.tsv
        """

}
