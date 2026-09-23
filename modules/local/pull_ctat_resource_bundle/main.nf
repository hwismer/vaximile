process PULL_CTAT_RESOURCE_BUNDLE {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    storeDir workflow.stubRun ? null : './vaximile_resources/ctat_resource_dir'

    // Pulls the hg38 CTAT resource bundle needed for STARfusion and other tools.

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


//***************************************************************************************************************************
// FILE INDEXING OPERATIONS
