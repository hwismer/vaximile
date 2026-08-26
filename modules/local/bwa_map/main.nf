process BWA_MAP {

    /*
        Map fastq files using BWA. Outputs a sorted BAM file and its index
        Reads groups are created using metadata information and currently are basically just the same name.
        Creates read group solely based on provided metadata from samplesheet. Any readgroup information
        present in the FASTQs is ignored.
    */

    cpus 8
    memory "32GB"
    container "iarcbioinfo/bwa-mem2-tools:v1.0"
    
    tag "BWA Alignment on ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
        tuple path(reference_fa), path(reference_index)
        path bwa_index

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}.sam")

    script:
    """
    NEW_RG="@RG\\tID:${meta.sample_name}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:${meta.molecule}_${meta.sequencing_type}"

    bwa-mem2 mem -t $task.cpus -R \$NEW_RG $reference_fa $fastq1 $fastq2 > "${meta.sample_name}_${meta.molecule}.sam"

    """
}
