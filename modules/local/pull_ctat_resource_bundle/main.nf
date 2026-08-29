process PULL_CTAT_RESOURCE_BUNDLE {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    storeDir workflow.stubRun ? null : './vaximile_resources/ctat_resource_dir'

    // Pulls the hg38 CTAT resource bundle needed for STARfusion and other tools.

    label 'process_single'
    executor "local"

    tag "Pulling CTAT plug-n-play resource bundle"

    output:
        path("./GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir"), emit: ctat_resource_dir
        // No versions.yml here: storeDir only short-circuits when EVERY declared
        // output is already in the store. An absent versions.yml made this process
        // re-run on a populated store and then fail moving its result on top of the
        // copy already there ("unable to remove target: Directory not empty").

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
