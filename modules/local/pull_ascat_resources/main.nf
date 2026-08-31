process PULL_ASCAT_RESOURCES {

    // Not under -stub-run: the stub writes empty placeholders, and storing those
    // would make a later real run skip the download and use empty resources.
    //
    // The path carries the chr-prefix flag, unlike the other PULL_* stores. This process
    // is the only one whose output depends on an input: with add_chr_prefix set it rewrites
    // every loci file with sed. storeDir is keyed on the path alone, not on the task hash,
    // so a single shared directory would hand back chr-prefixed loci to a later run that
    // asked for unprefixed ones - silently, and only visible as wrong ASCAT calls.
    storeDir workflow.stubRun ? null : "./vaximile_resources/ascat_hg38_${add_chr_prefix ? 'chr' : 'nochr'}"

    label 'process_single'
    executor "local"

    tag "Pulling ASCAT resources"


    input:
        val(add_chr_prefix)

    output:
        path("G1000_loci_hg38"), emit: loci
        path("G1000_alleles_hg38"), emit: alleles
        path("GC_G1000_hg38"), emit: GC
        path("RT_G1000_hg38"), emit: RT
        // No versions.yml here: storeDir only short-circuits when EVERY declared
        // output is already in the store. An absent versions.yml made this process
        // re-run on a populated store and then fail moving its result on top of the
        // copy already there ("unable to remove target: Directory not empty").

    script:
        
    """
    wget https://zenodo.org/records/14008443/files/G1000_loci_WGS_hg38.zip
    unzip G1000_loci_WGS_hg38.zip -d G1000_loci_hg38


	if ${add_chr_prefix ? 'true' : 'false'}; then
        echo "adding chr prefix to loci files"
        for i in {1..22} X; do
            f=G1000_loci_hg38/G1000_loci_hg38_chr\${i}.txt
            [ -f "\$f" ] || { echo "ERROR: missing \$f -- check filename pattern" >&2; exit 1; }
            sed -i '/^chr/! s/^/chr/' "\$f"
        done
        echo "loci chr1 now starts: \$(head -1 G1000_loci_hg38/G1000_loci_hg38_chr1.txt)"
    fi

    wget https://zenodo.org/records/14008443/files/G1000_alleles_WGS_hg38.zip
    unzip G1000_alleles_WGS_hg38.zip -d G1000_alleles_hg38

    wget https://zenodo.org/records/14008443/files/GC_G1000_WGS_hg38.zip
    unzip GC_G1000_WGS_hg38.zip -d GC_G1000_hg38

    wget https://zenodo.org/records/14008443/files/RT_G1000_WGS_hg38.zip
    unzip RT_G1000_WGS_hg38.zip -d RT_G1000_hg38
    """

    stub:
    """
    mkdir -p G1000_loci_hg38
    mkdir -p G1000_alleles_hg38
    mkdir -p GC_G1000_hg38
    mkdir -p RT_G1000_hg38
    """


}
