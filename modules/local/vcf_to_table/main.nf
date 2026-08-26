process VCF_TO_TABLE {

    cpus 2
    memory "8GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Exporint $vcf to tsv table"

    input:
        tuple val(meta), val(file_name), path(vcf), path(vcf_index)

    output:
        tuple val(meta), path("${file_name}.tsv")

    script:
        """
        gatk VariantsToTable \
            -V $vcf \
            -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC -F AF -F set -F FILTER -F CSQ \
            -GF AD -GF DP -GF GT -GF AF \
            -GF RDP -GF RAF -GF RAD -GF RADF -GF RADR -GF TX -GF GX \
            -O "${file_name}.tsv"
        """


}
