process PULL_VEP_CACHE {

    // Download the Ensembl VEP cache (release 115, GRCh38, ~24 GiB) when --vep_cache is not given.
    // The release must match the VEP version the annotation modules run.

    storeDir workflow.stubRun ? null : './vaximile_resources/vep_cache'

    label 'process_single'
    label 'process_long'

    executor "local"

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
