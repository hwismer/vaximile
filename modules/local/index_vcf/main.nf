process INDEX_VCF {
    
    // TBI index a vcf file


    cpus 2
    memory "8GB"
    container "staphb/bcftools:1.23"

    tag "Indexing $vcf"

    input:
        tuple val(vcf_name), val(meta), path(vcf) // vcf name should be a string that corresponds to the file name ie. Sample1
        val(filename_suffix) // suffix corresponds to name after vcf_name ie somatic_variants -> Sample1_somatic_variants.vcf.gz

    output:
        tuple val(meta), path("${vcf_name}_${filename_suffix}.vcf.gz"), path("${vcf_name}_${filename_suffix}.vcf.gz.tbi")

    script:
        """
        bcftools view $vcf -Oz -o "${vcf_name}_${filename_suffix}.vcf.gz"
        bcftools index -t ${vcf_name}_${filename_suffix}.vcf.gz
        """

}
