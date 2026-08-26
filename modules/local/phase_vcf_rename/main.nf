process PHASE_VCF_RENAME {

    /*

    Part of creating a phase vcf file.
    Germline sample name will be the name of the NORMAL sample, but to combine variants with the tumor sample, the names must match.
    Here the sample name in the germline vcf is renamed to the name of the tumor sample.

    */

    cpus 2
    memory "16GB"
    container "staphb/bcftools:1.23.1" 

    tag "Renaming germline sample ${normal_meta.sample_name} to ${somatic_meta.tumor_meta.sample_name}"

    input:
        tuple val(normal_meta), val(somatic_meta), path(germline_vcf), path(germline_vcf_index)

    output:
        tuple val(somatic_meta), path("${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz"), path("${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz.tbi")

    script:
        """
        cat > sample_map.txt <<EOF
        ${normal_meta.sample_name} ${somatic_meta.tumor_meta.sample_name}
        EOF
        
        bcftools reheader \
            -N sample_map.txt \
            --threads $task.cpus \
            -o ${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz \
            $germline_vcf

        bcftools index -t --threads $task.cpus ${somatic_meta.tumor_meta.sample_name}_germline_rename.vcf.gz 
        """


}
