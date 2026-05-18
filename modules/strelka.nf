process POSTPROCESS_STRELKA {

    /*

    Combine the SNVs and Indels from a Strelka and Manta run and rename the default TUMOR and NORMAL
    samples to the specified sample names.

    */

    cpus 4
    memory "16GB"
    
    container "staphb/bcftools:1.23.1" 

    input:
        tuple val(somatic_meta),
              path(strelka_snvs), path(strelka_snvs_index),
              path(strelka_indels), path(strelka_indels_index)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_strelka_snvs_indels.vcf.gz"), path("${somatic_meta.somatic_name}_strelka_snvs_indels.vcf.gz.tbi"), emit: vcf


    script:
        """
        bcftools concat \
            --allow-overlaps \
            --remove-duplicates \
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
            -o ${somatic_meta.somatic_name}_strelka_snvs_indels.vcf.gz \
            strelka_merged.vcf.gz

        bcftools index -t --threads $task.cpus ${somatic_meta.somatic_name}_strelka_snvs_indels.vcf.gz
        """
}

process MANTA {

    /*

    Run Manta Indel Caller

    */

   cpus 8
   memory "24GB"
   tag "Running Manta on ${somatic_meta.somatic_name}"
    
   container "mgibio/manta_somatic-cwl:1.6.0"

   input:
        tuple val(somatic_meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bed), path(bed_index)
        tuple path(reference_fa), path(reference_fai), path(reference_dict)
    output:
        tuple val(somatic_meta), path("./${somatic_meta.somatic_name}_manta")

    script:
    def exome_flag = (somatic_meta.tumor_meta.sequencing_type == "exome" || somatic_meta.tumor_meta.sequencing_type == "exome_ffpe") ? "--exome" : ""

    """
    /usr/bin/manta/bin/configManta.py \
        --normalBam $normal_bam \
        --tumorBam $tumor_bam \
        --referenceFasta $reference_fa \
        ${exome_flag} \
        --callRegions $bed \
        --runDir ./${somatic_meta.somatic_name}_manta/
        
    ./${somatic_meta.somatic_name}_manta/runWorkflow.py -j $task.cpus
    """

}

process STRELKA {
    
    /*

    Run Manta Indel Caller

    */
    
    cpus 8
    memory "24GB"
    
    container "mgibio/strelka:2.9.9"

    tag "Running Strelka on ${somatic_meta.somatic_name}"

   input:
        tuple val(somatic_meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(manta_dir), path(bed), path(bed_index)
        tuple path(reference_fa), path(reference_fai), path(reference_dict)

    output:
        tuple val(somatic_meta), 
        path("./strelka/results/variants/somatic.snvs.vcf.gz"), path("./strelka/results/variants/somatic.snvs.vcf.gz.tbi"), 
        path("./strelka/results/variants/somatic.indels.vcf.gz"), path("./strelka/results/variants/somatic.indels.vcf.gz.tbi"), emit: strelka_vcfs

    script:

    def exome_flag = (somatic_meta.tumor_meta.sequencing_type == "exome" || somatic_meta.tumor_meta.sequencing_type == "exome_ffpe") ? "--exome" : ""

    """
    /opt/strelka/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam $normal_bam \
        --tumorBam $tumor_bam \
        ${exome_flag} \
        --callRegions $bed \
        --referenceFasta $reference_fa \
        --indelCandidates "${manta_dir}/results/variants/candidateSmallIndels.vcf.gz" \
        --runDir "./strelka"

    ./strelka/runWorkflow.py -m local -j $task.cpus
    """

}

