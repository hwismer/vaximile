process PULL_VEP_PVAC_PLUGINS {

    storeDir './vaximile_resources/vep_pvac_plugins'

    // Pulls the VEP plugins necessary to run pvactools. Runs locally to ensure internet connection.

    label 'process_single'
    executor "local"

    container "griffithlab/pvactools:6.0.3"

    tag "Pulling Frameshift and Wildtype VEP plugins"

    output:
        path("VEP_plugins"), emit: plugins
        path "versions.yml", topic: versions

    script:
    """
    git clone https://github.com/Ensembl/VEP_plugins.git

    pvacseq install_vep_plugin VEP_plugins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pvacseq --version 2>&1 | tail -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p VEP_plugins
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: 6.0.3
    END_VERSIONS
    """

}
