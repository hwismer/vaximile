process PHASE_VCF_RBPHASING {

    // Phase the combined, sorted VCF with GATK3 ReadBackedPhasing.

    label 'process_low_memory'
    // GATK3 only: CombineVariants/ReadBackedPhasing have no GATK4 equivalent, and bioconda's gatk3 needs a licensed jar.
    container "broadinstitute/gatk3:3.6-0"

    tag "ReadBacked Phasing for ${somatic_meta.somatic_name}"

    input:
        tuple val(tumor_meta), val(somatic_meta), path(combined_sorted_vcf), path(tumor_reads), path(tumor_reads_index)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(somatic_meta), path("*_phased.vcf"), emit: vcf

    script:

        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        java -Xmx${task.memory.toGiga() - 1}g -jar /usr/GenomeAnalysisTK.jar \
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
