process VCF_TO_TABLE {

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Exporting $vcf to tsv table"

    input:
        tuple val(meta), val(file_name), path(vcf), path(vcf_index)

    output:
        tuple val(meta), path("*.tsv"), emit: tsv
        path "versions.yml", topic: versions

    script:
        def prefix = task.ext.prefix ?: "${file_name}"
        """
        gatk VariantsToTable \
            -V $vcf \
            -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC -F AF -F set -F FILTER -F CSQ \
            -GF AD -GF DP -GF GT -GF AF \
            -GF RDP -GF RAF -GF RAD -GF RADF -GF RADR -GF TX -GF GX \
            -O "${prefix}.tsv"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${file_name}"
        """
        touch "${prefix}.tsv"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: 4.6.1.0
        END_VERSIONS
        """


}
