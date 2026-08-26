process FILTER_VCF {

    /*

    Filter VCF, retaining only PASS variants.

    */
    
    cpus 4
    memory "16GB"
    container "staphb/bcftools:1.23"

    tag "Filtering non-passing variants from $caller $somatic_vcf"
    
    
    input:
        tuple val(somatic_meta), val(caller), path(somatic_vcf), path(tbi)
        
    output:
        tuple val(somatic_meta), val(caller), path("${somatic_meta.somatic_name}_${caller}_filtered_variants.vcf.gz"), path("${somatic_meta.somatic_name}_${caller}_filtered_variants.vcf.gz.tbi"), emit: filtered_vcf
        
        
    script:
        """
        bcftools view -f PASS -Oz -o "${somatic_meta.somatic_name}_${caller}_filtered_variants.vcf.gz" $somatic_vcf
        bcftools index -t "${somatic_meta.somatic_name}_${caller}_filtered_variants.vcf.gz"
        """
        


}
