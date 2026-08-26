process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */

    cpus 4
    memory "16GB"
    container "staphb/bcftools:1.23"
    
    tag "Filtering VCF $vcf"

    input:
        tuple val(meta), val(caller), path(vcf), path(tbi)

    output:
        tuple val(meta), val(caller), path("${meta.sample_name}_${caller}_filtered_variants.vcf.gz"), path("${meta.sample_name}_${caller}_filtered_variants.vcf.gz.tbi"), emit: filtered_vcf


    script:
        """
        bcftools view -f PASS -Oz -o "${meta.sample_name}_${caller}_filtered_variants.vcf.gz" $vcf
        bcftools index -t "${meta.sample_name}_${caller}_filtered_variants.vcf.gz"
        """



}
