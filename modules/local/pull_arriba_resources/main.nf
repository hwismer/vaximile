process PULL_ARRIBA_RESOURCES {

    storeDir './vaximile_resources/arriba_hg38'

    // Pulls arribra resources from release 2.5.1. Runs locally incase job nodes don't have internet.

    label 'process_single'
    executor "local"

    tag "Pulling Arriba resources v2.5.1"

    output:
        tuple path("./arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz"), 
            path("./arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz"),
            path("./arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3"), emit: resources
    
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
