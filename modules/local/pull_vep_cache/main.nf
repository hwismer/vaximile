process PULL_VEP_CACHE {

    // Download the Ensembl VEP cache (release 115, GRCh38, ~24 GiB) when --vep_cache is not given.
    // The release must match the VEP version the annotation modules run.

    // Not under -stub-run: the stub writes an empty placeholder, and storing that would
    // make a later real run skip the download and annotate against an empty cache.
    storeDir workflow.stubRun ? null : './vaximile_resources/vep_cache'

    // process_long as well as process_single: base.config sets a 4 h default for every
    // process, and config directives override the ones written in a module, so a plain
    // `time` here would be ignored. process_long is the 20 h tier.
    label 'process_single'
    label 'process_long'

    executor "local"

    // scratch false: this runs on the launch node, whose /tmp is too small for the download.
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
