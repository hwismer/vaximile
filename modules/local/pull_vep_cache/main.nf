process PULL_VEP_CACHE {

    /*
        Pull the Ensembl VEP cache when --vep_cache is not supplied.

        Release 115 is pinned to match the VEP the annotation modules actually run
        (ensemblorg/ensembl-vep:release_115.0, bioconda ensembl-vep=115). VEP refuses a
        cache whose version does not match its own, so this is not a free parameter - it
        has to move whenever those modules move.

        homo_sapiens / GRCh38 matches the pipeline's default reference, the GENCODE GRCh38
        primary assembly. A run against another assembly or species has to supply its own
        cache with --vep_cache; VEP will otherwise fail on the assembly mismatch.

        This is the plain Ensembl cache, not the refseq or merged flavour, because the VEP
        modules pass neither --refseq nor --merged.

        About 24 GiB. It unpacks to homo_sapiens/115_GRCh38/, so the directory emitted here
        is the one --dir_cache wants.
    */

    // Not under -stub-run: the stub writes an empty placeholder, and storing that would
    // make a later real run skip the download and annotate against an empty cache.
    storeDir workflow.stubRun ? null : './vaximile_resources/vep_cache'

    // process_long as well as process_single: base.config sets a 4 h default for every
    // process, and config directives override the ones written in a module, so a plain
    // `time` here would be ignored. process_long is the 20 h tier.
    label 'process_single'
    label 'process_long'

    executor "local"

    // scratch false, written here rather than left to the config. These run on the local
    // executor - on the node Nextflow was launched from, not as SLURM jobs - so an
    // institutional `scratch = true` would put the download in that node's /tmp, which is
    // sized for nothing like this, while the scratch reservation that makes `scratch`
    // safe elsewhere (`--gres=scratch:...`) only ever applies to SLURM jobs. A directive in
    // the module outranks a generic process scope in config, so this holds whatever
    // config the pipeline is run with.
    scratch false
    tag "Pulling Ensembl VEP cache release 115 GRCh38"

    output:
        path("vep_cache"), emit: cache

    script:
    """
    mkdir -p vep_cache

    # --continue so a dropped transfer resumes rather than restarting 24 GiB from zero.
    wget --continue --tries=3 \\
        https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz

    tar -xzf homo_sapiens_vep_115_GRCh38.tar.gz -C vep_cache
    rm homo_sapiens_vep_115_GRCh38.tar.gz

    # Fail here rather than letting every VEP task fail hours later on a missing cache.
    test -d vep_cache/homo_sapiens/115_GRCh38
    """

    stub:
    """
    mkdir -p vep_cache/homo_sapiens/115_GRCh38
    """
}
