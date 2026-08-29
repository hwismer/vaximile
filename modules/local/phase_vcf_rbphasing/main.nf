process PHASE_VCF_RBPHASING {

    /*

    Use deprecated ReadBackedPhasing from GATK 3.6.0 to phase variants
    in the somatic + germline combined and sorted vcf

    */

    label 'process_medium'
    // GATK3 ONLY - deliberately has no conda spec, and must not be given a gatk4 one.
    // This tool has no GATK4 equivalent: CombineVariants and ReadBackedPhasing were both
    // dropped in GATK4. bioconda's `gatk` 3.x is only a wrapper that needs the licensed
    // jar registered by hand, so this module stays container-only.
    container "broadinstitute/gatk3:3.6-0"

    tag "ReadBacked Phasing for ${somatic_meta.somatic_name}"

    input:
        tuple val(tumor_meta), val(somatic_meta), path(combined_sorted_vcf), path(tumor_reads), path(tumor_reads_index)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(somatic_meta), path("*_phased.vcf")

    script:

        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_fa \
                -I $tumor_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o ${prefix}_phased.vcf

        """
    stub:

        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        touch ${prefix}_phased.vcf
        """

}
