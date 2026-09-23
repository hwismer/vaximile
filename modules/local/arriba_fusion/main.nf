process ARRIBA_FUSION {

    // Arriba is quick after alignment and takes no thread option.
    label 'process_high_memory'
    conda "bioconda::arriba=2.5.1"

    tag "Running Arriba fusion calling on ${meta.sample_name}"

    input:
        tuple val(meta), path(star_bam), path(star_bam_index)
        tuple path(reference_fa), path(reference_index)
        path gtf
        tuple path(arriba_blacklist), path(arriba_known_fusions), path(arriba_protein_domains)

    output:
        tuple val(meta), path("*_arriba_fusions.tsv"), emit: arriba_fusions
        tuple val(meta), path("*_arriba_fusions.discarded.tsv"), emit: discarded_fusions

    script:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        arriba -x $star_bam \
            -g $gtf \
            -a $reference_fa \
            -b $arriba_blacklist \
            -k $arriba_known_fusions \
            -p $arriba_protein_domains \
            -o "${prefix}_arriba_fusions.tsv" \
            -O "${prefix}_arriba_fusions.discarded.tsv"
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
        """
        touch ${prefix}_arriba_fusions.tsv
        touch ${prefix}_arriba_fusions.discarded.tsv
        """

}
