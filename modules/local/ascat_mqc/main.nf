process ASCAT_MQC {

    label 'process_single'

    conda "conda-forge::pillow=12.3.0 python=3.12"

    tag "Preparing ${somatic_name} ASCAT results for MultiQC"

    input:
        tuple val(somatic_name), val(meta), path(plots), path(purityploidy), path(metrics)

    output:
        tuple val(somatic_name), val(meta), path("*_ascat_mqc.png"), emit: png
        tuple val(somatic_name), val(meta), path("*_ascatmetrics.tsv"), emit: tsv

    script:
    // Tile ASCAT's plots into one sheet for MultiQC, keeping each panel's aspect ratio.
    def plot_list = [plots].flatten().collect { "'${it}'" }.join(", ")
    """
    python3 - <<'PY'
import pathlib
import sys

from PIL import Image, ImageDraw, ImageFont

# Tile ASCAT's plots into one sheet per pair, two across in name order. Panels keep their aspect
# ratio: ASCAT's square sunrise plot sits beside genome-wide profiles three times as wide as tall.
PAIR = "${somatic_name}"
PANELS = [${plot_list}]
COLS = 2
CELL_W, LABEL_H, PAD = 900, 34, 8

panels = sorted((pathlib.Path(p) for p in PANELS if pathlib.Path(p).exists()), key=lambda p: p.stem)
if not panels:
    sys.exit("No panels to tile")
cells = [(i % COLS, i // COLS, path) for i, path in enumerate(panels)]
n_cols, n_rows = min(COLS, len(panels)), (len(panels) + COLS - 1) // COLS


def label(path):
    stem = path.stem[len(PAIR):] if path.stem.startswith(PAIR) else path.stem
    return stem.strip("._- ").replace(".", " ") or path.stem


# Row heights are per row, so a wide row does not sit in the whitespace of the square sunrise plot.
scaled_h = {p: round(CELL_W * Image.open(p).height / Image.open(p).width) for p in panels}
row_h = {}
for _col, row, path in cells:
    row_h[row] = max(row_h.get(row, 0), scaled_h[path])
row_y, offset = {}, PAD
for row in range(n_rows):
    row_y[row] = offset
    offset += row_h[row] + LABEL_H + PAD

font = ImageFont.load_default(size=22)
sheet = Image.new("RGB", (PAD + n_cols * (CELL_W + PAD), offset), "white")
draw = ImageDraw.Draw(sheet)
for col, row, path in cells:
    x, y = PAD + col * (CELL_W + PAD), row_y[row]
    draw.text((x, y + 6), label(path), fill="black", font=font)
    panel = Image.open(path).convert("RGB")
    scale = min(CELL_W / panel.width, row_h[row] / panel.height)
    panel = panel.resize((round(panel.width * scale), round(panel.height * scale)), Image.LANCZOS)
    sheet.paste(panel, (x + (CELL_W - panel.width) // 2, y + LABEL_H))
sheet.save(f"{PAIR}_ascat_mqc.png", optimize=True)
PY

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

    """

    stub:
    """
    touch ${somatic_name}_ascat_mqc.png
    printf 'Pair\\tAberrantCellFraction\\tPloidy\\n${somatic_name}\\t0.5\\t2.0\\n' > ${somatic_name}_ascatmetrics.tsv
    """
}
