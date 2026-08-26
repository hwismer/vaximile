process PHASE_VCF_RBPHASING {

    /*

    Use deprecated ReadBackedPhasing from GATK 3.6.0 to phase variants
    in the somatic + germline combined and sorted vcf

    */

    cpus 4
    memory "32GB"
    container "broadinstitute/gatk3:3.6-0"

    tag "ReadBacked Phasing for ${somatic_meta.somatic_name}"

    input:
        tuple val(tumor_meta), val(somatic_meta), path(combined_sorted_vcf), path(tumor_reads), path(tumor_reads_index)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_phased.vcf")

    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T ReadBackedPhasing \
                -R $reference_fa \
                -I $tumor_reads \
                --variant $combined_sorted_vcf \
                -L $combined_sorted_vcf \
                -o ${somatic_meta.somatic_name}_phased.vcf

        """
}
