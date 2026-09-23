process PULL_VEP_PVAC_PLUGINS {

    storeDir workflow.stubRun ? null : './vaximile_resources/vep_pvac_plugins'

    // Download the Frameshift and Wildtype VEP plugins pVACseq needs, pinned to the pVACtools version in use.

    label 'process_single'
    executor "local"

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
