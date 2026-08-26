process BAMREADCOUNT {

    /*

    Run bamreadcount to add coverage information to a vcf using a BAM file. The sample names in the VCF
    must match the read groups and sample names in the BAM.

    */

    label 'process_high'
    cache "lenient"

    container "mgibio/bam_readcount_helper-cwl:1.2.1"

    input:
        tuple val(somatic_name), val(somatic_meta), path(vcf), val(sample_meta), path(bam), path(bai)
        tuple path(reference_fa), path(reference_index)

    output:

        tuple val(somatic_meta), val(sample_meta), path("*_bamrc_helper/*indel.tsv"), path("*_bamrc_helper/*snv.tsv"), emit: brc_files
    script:
        def prefix = task.ext.prefix ?: "${sample_meta.sample_name}_${sample_meta.molecule}"
        """
        mkdir ${prefix}_bamrc_helper
        bam_readcount_helper.py \
            $vcf \
            ${sample_meta.sample_name} \
            $reference_fa \
            $bam \
            ${sample_meta.molecule} \
            ${prefix}_bamrc_helper
        """

}
