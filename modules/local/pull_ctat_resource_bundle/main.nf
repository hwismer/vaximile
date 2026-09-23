process PULL_CTAT_RESOURCE_BUNDLE {

    storeDir workflow.stubRun ? null : './vaximile_resources/ctat_resource_dir'

    // Pulls the hg38 CTAT resource bundle needed for STARfusion and other tools.

    label 'process_single'
    executor "local"

    scratch false
    tag "Pulling CTAT plug-n-play resource bundle"

    output:
        path("./GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir"), emit: ctat_resource_dir

    script:
    """
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
    tar -xzf GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz

    """

    stub:
    """
    mkdir -p ./GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir
    """

}
