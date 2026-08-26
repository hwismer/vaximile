process POSTPROCESS_STRELKA {

    /*

    Combine the SNVs and Indels from a Strelka and Manta run and rename the default TUMOR and NORMAL
    samples to the specified sample names.

    */

    label 'process_medium'

    container "staphb/bcftools:1.23.1" 

    input:
        tuple val(somatic_meta),
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
            -W=tbi \
            -o strelka_merged.vcf.gz \
            --threads $task.cpus \
            $strelka_snvs $strelka_indels

        cat > sample_map.txt <<EOF
        NORMAL ${somatic_meta.normal_meta.sample_name}
        TUMOR ${somatic_meta.tumor_meta.sample_name}
        EOF
        
        bcftools reheader \
            -N sample_map.txt \
            --threads $task.cpus \
            -o ${prefix}_strelka_snvs_indels.vcf.gz \
            strelka_merged.vcf.gz

        bcftools index -t --threads $task.cpus ${prefix}_strelka_snvs_indels.vcf.gz
        """
}
