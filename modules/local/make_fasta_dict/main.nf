process MAKE_FASTA_DICT {

    // Generate picard fasta.dict file for use with GATK tools

    label 'process_medium'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Creating reference dict for $fasta"

    input:
        tuple path(fasta), path(fai)

    output:
        path("*.dict"), emit: dict
        path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    gatk CreateSequenceDictionary \
        R=${fasta} \
        O=${prefix}.dict
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.dict
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}

//*******************************************************************************************************************
