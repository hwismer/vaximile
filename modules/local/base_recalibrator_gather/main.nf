process BASE_RECALIBRATOR_GATHER {

    /*
    Gathers scattered BaseRecalibrator recal tables from separate intervals and outputs the final table.
    */
    
    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "GatherBQSRReports on ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), val(molecule), path(recal_tables)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_recal_table.table")
    
    script:
    """
    gatk GatherBQSRReports \
        ${recal_tables.collect { "-I ${it}" }.join(' ')} \
        -O "${sample_name}_${molecule}_recal_table.table"
    """

}
