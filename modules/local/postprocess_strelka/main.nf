process POSTPROCESS_STRELKA {

    /*

    Combine the SNVs and Indels from a Strelka and Manta run and rename the default TUMOR and NORMAL
    samples to the specified sample names.

    */

    label 'process_medium'

    conda "bioconda::bcftools=1.23.1"
    container "staphb/bcftools:1.23.1"

    input:
        tuple val(somatic_meta), val(tumor_sample_name), val(normal_sample_name),
              path(strelka_snvs), path(strelka_snvs_index),
              path(strelka_indels), path(strelka_indels_index)

    output:
        tuple val(somatic_meta), path("*_strelka_snvs_indels.vcf.gz"), path("*_strelka_snvs_indels.vcf.gz.tbi"), emit: vcf


    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        bcftools concat \
            $args \
            --threads $task.cpus \
            -Oz \
            -o strelka_merged.vcf.gz \
            $strelka_snvs $strelka_indels

        cat > sample_map.txt <<EOF
        NORMAL ${normal_sample_name}
        TUMOR ${tumor_sample_name}
        EOF
        
        bcftools reheader \
            -N sample_map.txt \
            --threads $task.cpus \
            -o ${prefix}_strelka_snvs_indels.vcf.gz \
            strelka_merged.vcf.gz

        bcftools index -t --threads $task.cpus ${prefix}_strelka_snvs_indels.vcf.gz
        """

    stub:
        def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
        """
        touch ${prefix}_strelka_snvs_indels.vcf.gz
        touch ${prefix}_strelka_snvs_indels.vcf.gz.tbi
        """
}
