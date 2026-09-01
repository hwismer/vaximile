process MUTECT2_LEARN_READ_ORIENTATION {
    
    label 'process_low_memory'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Learning read orientation for ${somatic_meta.somatic_name}"

    input:
        tuple val(somatic_meta), path(f1r2s)

    output:
        tuple val(somatic_meta), path("*_orientmodel.tar.gz"), emit: tar
        path "versions.yml", topic: versions

    script:

    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"
    def f1r2_as_input = f1r2s.collect { f1r2 ->
            "-I ${f1r2}"
        }.join(' ')

    """
    gatk LearnReadOrientationModel \
        $f1r2_as_input \
        -O "${prefix}_orientmodel.tar.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -Eo 'v[0-9.]+' | head -1 | tr -d 'v')
    END_VERSIONS
    """

    stub:

    def prefix = task.ext.prefix ?: "${somatic_meta.somatic_name}"

    """
    touch ${prefix}_orientmodel.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.1.0
    END_VERSIONS
    """
}
