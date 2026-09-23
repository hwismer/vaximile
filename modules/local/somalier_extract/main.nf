process SOMALIER_EXTRACT {


    // somalier has no internal threading - its README parallelises across samples,
    // which Nextflow already does. nf-core's somalier modules are process_low.
    label 'process_low'
    conda "bioconda::somalier=0.3.2-0 bioconda::htslib=1.23.1"

    tag "somalier extract on ${meta.sample_name} at ${sites}"

    input:
        tuple val(meta), path(bam), path(bai)
        tuple path(reference_fa), path(fai)
        path(sites)

    output:
        // Glob, not the exact name: `somalier extract` derives the filename from the
        // BAM's SM read-group tag, which need not equal meta.sample_name. Naming it
        // exactly made the process fail with a missing output whenever they differed.
        tuple val(meta), path("*.somalier"), emit: somalier
        

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
