process MANTA {

    /*

    Run Manta Indel Caller

    */

   label 'process_high'
   tag "Running Manta on ${somatic_meta.somatic_name}"
    
   container "mgibio/manta_somatic-cwl:1.6.0"

   input:
        tuple val(somatic_meta), val(tumor_sequencing_type), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bed), path(bed_index)
        tuple path(reference_fa), path(reference_fai)
    output:
        tuple val(somatic_meta), path("*_manta")

    script:
    def exome_flag = (tumor_sequencing_type == "exome" || tumor_sequencing_type == "exome_ffpe") ? "--exome" : ""
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"

    """
    /usr/bin/manta/bin/configManta.py \
        --normalBam $normal_bam \
        --tumorBam $tumor_bam \
        --referenceFasta $reference_fa \
        ${exome_flag} \
        --callRegions $bed \
        --runDir ./${prefix}_manta/

    ./${prefix}_manta/runWorkflow.py -j $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    """
    mkdir -p ${prefix}_manta/results/variants
    touch ${prefix}_manta/results/variants/candidateSmallIndels.vcf.gz
    """

}
