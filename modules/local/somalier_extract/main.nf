process SOMALIER_EXTRACT {


    label 'process_medium'
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier extract on ${meta.sample_name} at ${sites}"

    input:
        tuple val(meta), path(bam), path(bai)
        tuple path(reference_fa), path(fai)
        path(sites)

    output:
        tuple val(meta), path("${meta.sample_name}.somalier"), emit: somalier
        path "versions.yml", topic: versions
        

    script:
    """
    somalier extract \
        -s $sites \
        -f $reference_fa \
        $bam 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(somalier --version 2>&1 | grep -Eo '[0-9]+\.[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_name}.somalier
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: 0.3.2
    END_VERSIONS
    """
}
