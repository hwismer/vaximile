process MUTECT2_LEARN_READ_ORIENTATION {
    
    cpus 4
    memory "12GB"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Learning read orientation for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(f1r2s)

    output:
        tuple val(somatic_meta), path("${somatic_meta.somatic_name}_orientmodel.tar.gz")

    script:
    
    def f1r2_as_input = f1r2s.collect { f1r2 ->
            "-I ${f1r2}"
        }.join(' ')

    """
    gatk LearnReadOrientationModel \
        $f1r2_as_input \
        -O "${somatic_meta.somatic_name}_orientmodel.tar.gz"
    """
}
