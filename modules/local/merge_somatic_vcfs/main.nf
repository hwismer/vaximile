process MERGE_SOMATIC_VCFS {

    /*

    Use the deprecated CombineVariants from GATK 3.6.0 to combine vcf files. 

    */

    label 'process_medium'
    // GATK3 ONLY - deliberately has no conda spec, and must not be given a gatk4 one.
    // This tool has no GATK4 equivalent: CombineVariants and ReadBackedPhasing were both
    // dropped in GATK4. bioconda's `gatk` 3.x is only a wrapper that needs the licensed
    // jar registered by hand, so this module stays container-only.
    container "broadinstitute/gatk3:3.6-0"

    tag "Merge 3 somatic vcfs from $vcf1 $vcf2 $vcf3"

    input:
        tuple val(somatic_meta), val(somatic_name),
            val(vcf1_caller), path(vcf1), path(vcf1_index),
            val(vcf2_caller), path(vcf2), path(vcf2_index),
            val(vcf3_caller), path(vcf3), path(vcf3_index)
        tuple path(reference_fa), path(reference_index)
        path(reference_dict)
    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_variants.vcf.gz"), emit: vcf
        path "versions.yml", topic: versions


    script:

        def args = task.ext.args ?: ''
        """
        java -Xmx${task.memory.toGiga() - 1}g -jar /usr/GenomeAnalysisTK.jar \
            -T CombineVariants \
            -R $reference_fa \
            -genotypeMergeOptions PRIORITIZE \
            --rod_priority_list $vcf1_caller,$vcf2_caller,$vcf3_caller \
            -V:$vcf1_caller $vcf1 \
            -V:$vcf2_caller $vcf2 \
            -V:$vcf3_caller $vcf3 \
            $args \
            -o "${somatic_name}_variants.vcf.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk3: 3.6-0
        END_VERSIONS
        """

    stub:
        """
        touch ${somatic_meta.somatic_name}_variants.vcf.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk3: 3.6-0
        END_VERSIONS
        """

}
