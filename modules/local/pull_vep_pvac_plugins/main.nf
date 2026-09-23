process PULL_VEP_PVAC_PLUGINS {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    storeDir workflow.stubRun ? null : './vaximile_resources/vep_pvac_plugins'

    // The two VEP plugins pVACseq depends on - Frameshift and Wildtype, the only ones
    // VEP_ANNOTATE loads - fetched directly from the pVACtools repository. Runs locally for
    // the network access.
    //
    // This used to clone all of Ensembl's VEP_plugins from an unpinned master and then run
    // `pvacseq install_vep_plugin` inside the griffithlab/pvactools image to drop these two
    // files into it. That subcommand is nothing but two file copies, and starting pvacseq to
    // perform them loads the whole CLI and its compiled dependencies, which died with
    // "Illegal instruction (core dumped)" on the launch node's CPU. None of that is needed to
    // obtain two Perl files, and no Ensembl plugin is used.
    //
    // Pinned to the tag the PVACSEQ and PVACFUSE modules run. The image this replaced was
    // 6.0.3 against their 7.0.1; both files are identical across those versions, so no
    // result was affected, but the pin should follow those modules.

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

    tag "Pulling Frameshift and Wildtype VEP plugins"

    output:
        path("VEP_plugins"), emit: plugins

    script:
    def pvactools_version = '7.0.1'
    """
    mkdir -p VEP_plugins
    for plugin in Wildtype Frameshift; do
        curl -fsSL --retry 3 -o "VEP_plugins/\${plugin}.pm" \
            "https://raw.githubusercontent.com/griffithlab/pVACtools/v${pvactools_version}/pvactools/tools/pvacseq/VEP_plugins/\${plugin}.pm"
        # curl -f fails on an HTTP error, but not on an empty body
        test -s "VEP_plugins/\${plugin}.pm"
    done
    """

    stub:
    """
    mkdir -p VEP_plugins
    touch VEP_plugins/Wildtype.pm VEP_plugins/Frameshift.pm
    """

}
