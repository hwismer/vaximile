process MAKE_FASTA_DICT {

    // Generate picard fasta.dict file for use with GATK tools

    cpus 4
    memory "16GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Creating reference dict for $fasta"

    input:
        tuple path(fasta), path(fai)

    output:
        path("${fasta.baseName}.dict"), emit: dict

    script:
    """
    gatk CreateSequenceDictionary \
        R=${fasta} \
        O=${fasta.baseName}.dict
    """
}

//*******************************************************************************************************************
