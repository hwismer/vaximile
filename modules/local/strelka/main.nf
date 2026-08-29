process STRELKA {
    
    /*

    Run Manta Indel Caller

    */
    
    label 'process_high'

    container "mgibio/strelka:2.9.9"

    tag "Running Strelka on ${somatic_meta.somatic_name}"

   input:
        tuple val(somatic_meta), val(tumor_sequencing_type), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(manta_dir), path(bed), path(bed_index)
        tuple path(reference_fa), path(reference_fai)

    output:
        tuple val(somatic_meta), 
        path("./strelka/results/variants/somatic.snvs.vcf.gz"), path("./strelka/results/variants/somatic.snvs.vcf.gz.tbi"), 
        path("./strelka/results/variants/somatic.indels.vcf.gz"), path("./strelka/results/variants/somatic.indels.vcf.gz.tbi"), emit: strelka_vcfs
        path "versions.yml", topic: versions

    script:

    def exome_flag = (tumor_sequencing_type == "exome" || tumor_sequencing_type == "exome_ffpe") ? "--exome" : ""

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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: 2.9.9
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ./strelka/results/variants
    touch ./strelka/results/variants/somatic.snvs.vcf.gz
    touch ./strelka/results/variants/somatic.snvs.vcf.gz.tbi
    touch ./strelka/results/variants/somatic.indels.vcf.gz
    touch ./strelka/results/variants/somatic.indels.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: 2.9.9
    END_VERSIONS
    """

}
