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
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_recal_table.table"), emit: table
        path "versions.yml", topic: versions
    
    script:
    """
    gatk GatherBQSRReports \
        ${recal_tables.collect { "-I ${it}" }.join(' ')} \
        -O "${sample_name}_${molecule}_recal_table.table"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.sample_name}_${meta.molecule}_recal_table.table"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """

}
