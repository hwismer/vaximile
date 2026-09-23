process PHASE_VCF_RENAME {

    /*

    Part of creating a phase vcf file.
    Germline sample name will be the name of the NORMAL sample, but to combine variants with the tumor sample, the names must match.
    Here the sample name in the germline vcf is renamed to the name of the tumor sample.

    */

    label 'process_low'
    conda "bioconda::bcftools=1.23.1"
    container "staphb/bcftools:1.23.1"

    tag "Renaming germline sample ${normal_meta.sample_name} to ${somatic_meta.tumor_meta.sample_name}"

    input:
        tuple val(normal_meta), val(somatic_meta), val(normal_sample_name), val(tumor_sample_name), path(germline_vcf), path(germline_vcf_index)

    output:
        tuple val(somatic_meta), path("*_germline_rename.vcf.gz"), path("*_germline_rename.vcf.gz.tbi"), emit: vcf

    script:
        def prefix = task.ext.prefix ?: "${tumor_sample_name}"
        """
        cat > sample_map.txt <<EOF
        ${normal_sample_name} ${tumor_sample_name}
        EOF
        
        bcftools reheader \
            -N sample_map.txt \
            --threads $task.cpus \
            -o ${prefix}_germline_rename.vcf.gz \
            $germline_vcf

        bcftools index -t --threads $task.cpus ${prefix}_germline_rename.vcf.gz
        """

    stub:
        def prefix = task.ext.prefix ?: "${tumor_sample_name}"
        """
        touch ${prefix}_germline_rename.vcf.gz
        touch ${prefix}_germline_rename.vcf.gz.tbi
        """


}
