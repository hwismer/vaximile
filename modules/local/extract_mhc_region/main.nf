process EXTRACT_MHC_REGION {
    label 'process_medium'
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam), path(bai)
    
    output:
        tuple val(meta), path("*_hla_regions.bam"), emit: bam
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        samtools view --threads $task.cpus -h -b -f 4 $bam > unmapped.bam
        samtools view --threads $task.cpus -h -b $bam chr6:28510120-33480577 > mhc.bam
        samtools merge --threads $task.cpus -o hla_regions.bam mhc.bam unmapped.bam
        samtools collate --threads $task.cpus -o "${prefix}_hla_regions.bam" hla_regions.bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}"
        """
        touch ${prefix}_hla_regions.bam
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: 1.23.1
        END_VERSIONS
        """

}
