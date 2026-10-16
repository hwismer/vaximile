process MOSDEPTH {

    // Coverage over the capture regions, for MultiQC.

    label 'process_low'
    conda "bioconda::mosdepth=0.3.10"

    tag "mosdepth on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai), path(capture_bed)

    output:
        tuple val(meta), path("*.mosdepth.global.dist.txt"), path("*.mosdepth.region.dist.txt"), path("*.mosdepth.summary.txt"), emit: dist

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    mosdepth --threads ${task.cpus} --by ${capture_bed} --no-per-base ${prefix} ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}.mosdepth.global.dist.txt ${prefix}.mosdepth.region.dist.txt ${prefix}.mosdepth.summary.txt
    """
}
