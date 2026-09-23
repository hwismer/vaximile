process MERGE_GERMLINE_VCFS {

    // Merge the three germline callsets with GATK3 CombineVariants.

    label 'process_low_memory'
    // GATK3 only: CombineVariants/ReadBackedPhasing have no GATK4 equivalent, and bioconda's gatk3 needs a licensed jar.
    container "broadinstitute/gatk3:3.6-0"

    tag "Combining VCFs: $vcf1 $vcf2 $vcf3"

    input:
        tuple val(sample_meta),
            val(vcf1_caller), path(vcf1), path(vcf1_index),
            val(vcf2_caller), path(vcf2), path(vcf2_index),
            val(vcf3_caller), path(vcf3), path(vcf3_index)
        tuple path(reference_fa), path(reference_index)
        path(reference_dict)

    output:
        tuple val(sample_meta), path("*_germline_variants.vcf.gz"), emit: vcf


    script:

        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"
        """
        java -Xmx${task.memory.toGiga() - 1}g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            $args \
            --rod_priority_list $vcf1_caller,$vcf2_caller,$vcf3_caller \
            -V:$vcf1_caller $vcf1 \
            -V:$vcf2_caller $vcf2 \
            -V:$vcf3_caller $vcf3 \
            -o "${prefix}_germline_variants.vcf.gz"
        """

    stub:

        def prefix = task.ext.prefix ?: "${sample_meta.sample_name}"
        """
        touch ${prefix}_germline_variants.vcf.gz
        """

}
