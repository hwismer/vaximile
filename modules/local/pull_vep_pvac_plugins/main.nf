process PULL_VEP_PVAC_PLUGINS {

    storeDir './vaximile_resources/vep_pvac_plugins'

    // Pulls the VEP plugins necessary to run pvactools. Runs locally to ensure internet connection.

    cpus 1
    memory "4GB"
    executor "local"

    container "griffithlab/pvactools:6.0.3"

    tag "Pulling Frameshift and Wildtype VEP plugins"

    output:
        path("VEP_plugins"), emit: plugins

    script:
    """
    git clone https://github.com/Ensembl/VEP_plugins.git

    pvacseq install_vep_plugin VEP_plugins

    """

}
