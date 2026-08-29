process PULL_ARRIBA_RESOURCES {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    storeDir workflow.stubRun ? null : './vaximile_resources/arriba_hg38'

    // Pulls arribra resources from release 2.5.1. Runs locally incase job nodes don't have internet.

    label 'process_single'
    executor "local"

    tag "Pulling Arriba resources v2.5.1"

    output:
        tuple path("./arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz"), 
            path("./arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz"),
            path("./arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3"), emit: resources
        // No versions.yml here: storeDir only short-circuits when EVERY declared
        // output is already in the store. An absent versions.yml made this process
        // re-run on a populated store and then fail moving its result on top of the
        // copy already there ("unable to remove target: Directory not empty").
    
    script:

        
    """
        wget https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz
        ls
        echo "done"
        tar -xzf arriba_v2.5.1.tar.gz
        ls

    """

    stub:

    """
    mkdir -p ./arriba_v2.5.1/database
    touch ./arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz
    touch ./arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz
    touch ./arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3
    """

}
