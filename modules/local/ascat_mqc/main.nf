process ASCAT_MQC {

    label 'process_single'

    conda "conda-forge::pillow=12.3.0 python=3.12"

    tag "Preparing ${somatic_name} ASCAT results for MultiQC"

    input:
        tuple val(somatic_name), val(meta), path(plots), path(purityploidy), path(metrics)
        path montage_script

    output:
        tuple val(somatic_name), val(meta), path("*_ascat_mqc.png"), emit: png
        tuple val(somatic_name), val(meta), path("*_ascatmetrics.tsv"), emit: tsv
        path "versions.yml", topic: versions

    script:
    // Flat layout rather than the grid the LOH sheet uses: ASCAT's plots have no two-axis
    // structure to lay out against, they are just a set. Two columns because the genome
    // profiles are wide; the script keeps each panel's own aspect ratio, so the square
    // sunrise plot sits in a row of its own height rather than being stretched to match.
    """
    python3 $montage_script ${plots} \\
        --out ${somatic_name}_ascat_mqc.png \\
        --layout flat \\
        --cols 2 \\
        --strip-prefix ${somatic_name}

    python3 - <<'METRICS' > ${somatic_name}_ascatmetrics.tsv
import csv
import pathlib
import sys

# Purity and ploidy come from one file and the QC metrics from another, both written by
# ASCAT as a header and a single row. They are merged into one row here so the report has
# one ASCAT line per pair rather than two tables to read against each other.
#
# The columns of ascat.metrics() are whatever the installed ASCAT version returns, so they
# are passed through as found rather than named here: a fixed list would silently drop
# columns a newer ASCAT adds.
def read_row(path):
    with open(path) as handle:
        rows = list(csv.reader(handle, delimiter="\\t"))
    if len(rows) < 2:
        return [], []
    return rows[0], rows[1]

pp_head, pp_row = read_row("${purityploidy}")
qc_head, qc_row = read_row("${metrics}")

# Pair first: ASCAT keys neither file by sample, and MultiQC reads the first column as the
# row key, so without this every pair would land on the same row and overwrite it.
head = ["Pair"] + pp_head + qc_head
row = ["${somatic_name}"] + pp_row + qc_row

writer = csv.writer(sys.stdout, delimiter="\\t", lineterminator="\\n")
writer.writerow(head)
writer.writerow(row)
METRICS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pillow: \$(python3 -c "import PIL; print(PIL.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${somatic_name}_ascat_mqc.png
    printf 'Pair\\tAberrantCellFraction\\tPloidy\\n${somatic_name}\\t0.5\\t2.0\\n' > ${somatic_name}_ascatmetrics.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pillow: 12.3.0
    END_VERSIONS
    """
}
