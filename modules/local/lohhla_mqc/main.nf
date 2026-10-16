process LOHHLA_MQC {

    // Loads LOHHLAmod results into MultiQC

    label 'process_single'

    conda "conda-forge::poppler=26.09.0 conda-forge::pillow=12.3.0 python=3.12"

    tag "Preparing ${somatic_name} HLA LOH results for MultiQC"

    input:
        tuple val(somatic_name), val(meta), path(plot_dir), path(loh_res)

    output:
        tuple val(somatic_name), val(meta), path("*_hla_loh_mqc.png"), emit: png
        tuple val(somatic_name), val(meta), path("*_lohres.tsv"), emit: tsv

    script:
    // Rasterise the LOH PDFs and tile them into one sheet per pair for MultiQC (genes across, plots down).
    """
    mkdir -p panels
    for pdf in ${plot_dir}/*.pdf; do
        [ -e "\$pdf" ] || continue
        pdftoppm -png -r 100 -singlefile "\$pdf" "panels/\$(basename "\$pdf" .pdf)"
    done

    python3 - <<'PY'
import pathlib
import sys

from PIL import Image, ImageDraw, ImageFont

# Tile the LOH panels into one sheet per pair: HLA genes across, plots down. Panel names are
# <gene>.<plot>.png, e.g. hla_a.logR.png. Panels keep their aspect ratio.
OUT = "${somatic_name}_hla_loh_mqc.png"
ROW_ORDER = ["t_dp", "n_dp", "tn_dp", "logR", "baf"]
CELL_W, LABEL_H, PAD = 900, 34, 8

panels = sorted(pathlib.Path("panels").glob("*.png"))
if not panels:
    sys.exit("No panels to tile")
table = {}
for path in panels:
    gene, _, plot = path.stem.partition(".")
    table.setdefault(gene, {})[plot] = path
genes = sorted(table)
present = {plot for entry in table.values() for plot in entry}
plots = [p for p in ROW_ORDER if p in present] + sorted(present - set(ROW_ORDER))
cells = [(c, r, table[gene].get(plot)) for r, plot in enumerate(plots) for c, gene in enumerate(genes)]

# Row heights are per row, so a short row does not sit in the whitespace of a tall one.
scaled_h = {p: round(CELL_W * Image.open(p).height / Image.open(p).width) for p in panels}
row_h = {}
for _col, row, path in cells:
    if path is not None:
        row_h[row] = max(row_h.get(row, 0), scaled_h[path])
row_y, offset = {}, PAD
for row in range(len(plots)):
    row_y[row] = offset
    offset += row_h.get(row, 0) + LABEL_H + PAD

font = ImageFont.load_default(size=22)
sheet = Image.new("RGB", (PAD + len(genes) * (CELL_W + PAD), offset), "white")
draw = ImageDraw.Draw(sheet)
for col, row, path in cells:
    if path is None:
        continue
    x, y = PAD + col * (CELL_W + PAD), row_y[row]
    draw.text((x, y + 6), path.stem.strip("._- ").replace(".", " "), fill="black", font=font)
    panel = Image.open(path).convert("RGB")
    scale = min(CELL_W / panel.width, row_h[row] / panel.height)
    panel = panel.resize((round(panel.width * scale), round(panel.height * scale)), Image.LANCZOS)
    sheet.paste(panel, (x + (CELL_W - panel.width) // 2, y + LABEL_H))
sheet.save(OUT, optimize=True)
PY

    # The result table is keyed by HLA gene alone - hla_a, hla_b, hla_c - which is unique
    # within a pair and not across them. MultiQC reads the first column as the row key and
    # merges every file it matches into one table, so without the pair in the key a second
    # pair overwrites the first rather than appending to it.
    awk -F'\t' -v OFS='\t' -v pair="${somatic_name}" '
        NR == 1 { \$1 = "Pair_HLAGene"; print; next }
        { \$1 = pair "_" \$1; print }
    ' ${loh_res} > ${somatic_name}_lohres.tsv

    """

    stub:
    """
    touch ${somatic_name}_hla_loh_mqc.png
    touch ${somatic_name}_lohres.tsv
    """
}
