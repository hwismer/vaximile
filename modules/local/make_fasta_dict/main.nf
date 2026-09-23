process MAKE_FASTA_DICT {

    // Generate picard fasta.dict file for use with GATK tools

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Creating reference dict for $fasta"

    input:
        tuple path(fasta), path(fai)

    output:
        path("*.dict"), emit: dict

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    gatk CreateSequenceDictionary \
        R=${fasta} \
        O=${prefix}.dict
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.dict
    """
}

//*******************************************************************************************************************
