process PULL_ARRIBA_RESOURCES {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    storeDir workflow.stubRun ? null : './vaximile_resources/arriba_hg38'

    // Pulls arribra resources from release 2.5.1. Runs locally incase job nodes don't have internet.

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
