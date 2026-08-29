process SAMTOOLS_FLAGSTAT {

    label 'process_high'
    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Samtools flagstat on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
       tuple val(meta), path("*.flagstat"), emit: flagstat
       path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    samtools flagstat --threads $task.cpus $bam > ${prefix}.flagstat
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}.flagstat
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """

}
