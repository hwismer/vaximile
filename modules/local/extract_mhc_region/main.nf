process EXTRACT_MHC_REGION {
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam), path(bai)
    
    output:
        tuple val(meta), path("${meta.sample_name}_hla_regions.bam") 

    script:
        """
        samtools view --threads $task.cpus -h -b -f 4 $bam > unmapped.bam
        samtools view --threads $task.cpus -h -b $bam chr6:28510120-33480577 > mhc.bam
        samtools merge --threads $task.cpus -o hla_regions.bam mhc.bam unmapped.bam
        samtools collate --threads $task.cpus -o "${meta.sample_name}_hla_regions.bam" hla_regions.bam
        """

}
