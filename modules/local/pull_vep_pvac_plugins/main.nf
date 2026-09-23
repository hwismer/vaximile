process PULL_VEP_PVAC_PLUGINS {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    storeDir workflow.stubRun ? null : './vaximile_resources/vep_pvac_plugins'

    // Pulls the VEP plugins necessary to run pvactools. Runs locally to ensure internet connection.

    label 'process_single'
    executor "local"

    // scratch false, written here rather than left to the config. These run on the local
    // executor - on the node Nextflow was launched from, not as SLURM jobs - so an
    // institutional `scratch = true` would put the download in that node's /tmp, which is
    // sized for nothing like this, while the scratch reservation that makes `scratch`
    // safe elsewhere (`--gres=scratch:...`) only ever applies to SLURM jobs. A directive in
    // the module outranks a generic process scope in config, so this holds whatever
    // config the pipeline is run with.
    scratch false
    container "griffithlab/pvactools:6.0.3"

    tag "Pulling Frameshift and Wildtype VEP plugins"

    output:
        path("VEP_plugins"), emit: plugins
        // No versions.yml here: storeDir only short-circuits when EVERY declared
        // output is already in the store. An absent versions.yml made this process
        // re-run on a populated store and then fail moving its result on top of the
        // copy already there ("unable to remove target: Directory not empty").

    script:
    """
    git clone https://github.com/Ensembl/VEP_plugins.git

    pvacseq install_vep_plugin VEP_plugins

    """

    stub:
    """
    mkdir -p VEP_plugins
    """

}
