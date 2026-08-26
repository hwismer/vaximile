process MERGE_SOMATIC_VCFS {

    /*

    Use the deprecated CombineVariants from GATK 3.6.0 to combine vcf files. 

    */

    cpus 4
    memory "16GB"
    container "broadinstitute/gatk3:3.6-0"

    tag "Merge 3 somatic vcfs from $vcf1 $vcf2 $vcf3"

    input:
        tuple val(somatic_meta), 
            val(vcf1_caller), path(vcf1), path(vcf1_index), 
            val(vcf2_caller), path(vcf2), path(vcf2_index),
            val(vcf3_caller), path(vcf3), path(vcf3_index)
        tuple path(reference_fa), path(reference_index)
        path(reference_dict)
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_variants.vcf.gz")


    script:

        """
        java -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            -genotypeMergeOptions PRIORITIZE \
            --rod_priority_list $vcf1_caller,$vcf2_caller,$vcf3_caller \
            -V:$vcf1_caller $vcf1 \
            -V:$vcf2_caller $vcf2 \
            -V:$vcf3_caller $vcf3 \
            --minimumN 2 \
            -o "${somatic_meta.somatic_name}_variants.vcf.gz"

        """
    
}
