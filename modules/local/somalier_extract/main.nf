process SOMALIER_EXTRACT {


    label 'process_medium'
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier extract on ${meta.sample_name} at ${sites}"

    input:
        tuple val(meta), path(bam), path(bai)
        tuple path(reference_fa), path(fai)
        path(sites)

    output:
        tuple val(meta), path("${meta.sample_name}.somalier")
        

    script:
    """
    somalier extract \
        -s $sites \
        -f $reference_fa \
        $bam 
    """

    stub:
    """
    touch ${meta.sample_name}.somalier
    """
}
