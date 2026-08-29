process PULL_ASCAT_RESOURCES {

    label 'process_single'
    executor "local"

    tag "Pulling ASCAT resources"

    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        val(add_chr_prefix)

    output:
        path("G1000_loci_hg38"), emit: loci
        path("G1000_alleles_hg38"), emit: alleles
        path("GC_G1000_hg38"), emit: GC
        path("RT_G1000_hg38"), emit: RT
        path "versions.yml", topic: versions

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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p G1000_loci_hg38
    mkdir -p G1000_alleles_hg38
    mkdir -p GC_G1000_hg38
    mkdir -p RT_G1000_hg38
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """


}
