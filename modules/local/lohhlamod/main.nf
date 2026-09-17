process LOHHLAMOD {

    label 'process_high'

    // No conda spec: lohhlamod is an R package installed from its own repository, published
    // to neither bioconda nor CRAN, so the container is the only way to get it.
    //
    // Third-party and unofficial - svm-zhang publishes no image, and this one's CI builds
    // `lohhlamod:local` from the checkout. Its build history matches upstream's Dockerfile
    // step for step (rocker/r-ver:4.1.0, the same apt and remotes::install_deps layers, the
    // same symlinks into /usr/local/bin), and `lohhlamod --help` inside it prints exactly
    // the interface R/cli.R declares.
    //
    // Pinned by digest, not by tag: the repository publishes only `latest`, which can be
    // repointed at a different image at any time, and a rebuilt `latest` would silently
    // change what this process runs. The digest also fixes the one thing the build history
    // cannot tell us - `COPY . .` does not record which commit of lohhla-mod was baked in.
    container "hmitzinger/lohhlamod@sha256:8d37ec357a5a5e4f04b6adbd7c69fc73156273527de0132a8f7d5eeb7398a293"

    tag "Detecting HLA LOH in ${somatic_name}"

    input:
        // purityploidy is ASCAT's, and is optional: unset, lohhlamod falls back to a
        // hardcoded ploidy of 2 and purity of 0.5, which is a guess rather than a failure.
        // ASCAT runs per pair and is allowed to fail on this cluster, so pairs without it
        // still get an LOH call rather than silently disappearing from the results.
        tuple val(somatic_name), val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(hla_fasta), path(purityploidy)

    output:
        tuple val(somatic_name), val(meta), path("${somatic_name}_lohhla"), emit: loh_dir
        tuple val(somatic_name), val(meta), path("${somatic_name}_lohhla/${somatic_name}.loh.res.tsv"), emit: loh_res
        // Per-gene intermediates, emitted separately because LOHHLAPLOT reads them: it
        // takes a directory to look for <HLAGene>.rds in, and passing the whole result
        // directory would also drag the filtered BAMs into the plotting task.
        tuple val(somatic_name), val(meta), path("${somatic_name}_lohhla/*.rds"), emit: loh_rds
        path "versions.yml", topic: versions

    script:
    def args = task.ext.args ?: ''
    // ASCAT writes "AberrantCellFraction\tPloidy"; lohhlamod reads a table whose first
    // column is dropped and which must carry TumorPloidy and TumorPurityNGS by name. The
    // reshape is two numbers, so it happens here rather than in a module of its own.
    def tstates = purityploidy ? "--tstates tstates.tsv" : ''
    def build_tstates = purityploidy ? """
    printf 'SampleID\\tTumorPloidy\\tTumorPurityNGS\\n' > tstates.tsv
    awk 'NR == 2 { print "${somatic_name}" "\\t" \$2 "\\t" \$1 }' ${purityploidy} >> tstates.tsv
    """ : ''
    """
    $build_tstates
    lohhlamod --subject ${somatic_name} \\
        --tbam $tumor_bam \\
        --nbam $normal_bam \\
        --hlaref $hla_fasta \\
        $tstates \\
        --threads $task.cpus \\
        $args \\
        --outdir ${somatic_name}_lohhla

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${somatic_name}_lohhla
    touch ${somatic_name}_lohhla/${somatic_name}.loh.res.tsv
    touch ${somatic_name}_lohhla/hla_a.rds ${somatic_name}_lohhla/hla_b.rds ${somatic_name}_lohhla/hla_c.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unknown: unknown
    END_VERSIONS
    """
}
