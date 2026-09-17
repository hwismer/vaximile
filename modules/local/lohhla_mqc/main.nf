process LOHHLA_MQC {

    label 'process_single'

    conda "conda-forge::poppler=26.09.0 conda-forge::pillow=12.3.0 python=3.12"

    tag "Preparing ${somatic_name} HLA LOH results for MultiQC"

    input:
        tuple val(somatic_name), val(meta), path(plot_dir), path(loh_res)

    output:
        // MultiQC picks up any *_mqc.png as a custom-content image and titles the section
        // from the file name, so the pair is encoded there: it is the only place that
        // survives into the report, and two pairs sharing a name would collide.
        tuple val(somatic_name), val(meta), path("*_hla_loh_mqc.png"), emit: png
        tuple val(somatic_name), val(meta), path("*_lohres.tsv"), emit: tsv
        path "versions.yml", topic: versions

    script:
    // lohhlaplot writes PDFs, which MultiQC cannot embed - it takes png, jpg and jpeg only -
    // so they are rasterised here rather than lost from the report, and tiled into one sheet
    // per pair. Fifteen separate images, three genes by five plot types, is fifteen report
    // sections to scroll past for a single pair; as one sheet it is one section, and the
    // three genes sit side by side, which is the comparison worth making.
    //
    // Genes across and plot types down rather than the other way round: five columns
    // squeezes each panel to a fifth of the report width, where three leaves them legible.
    // The PDFs stay the published artefact for anything needing full resolution.
    """
    mkdir -p panels
    for pdf in ${plot_dir}/*.pdf; do
        [ -e "\$pdf" ] || continue
        pdftoppm -png -r 100 -singlefile "\$pdf" "panels/\$(basename "\$pdf" .pdf)"
    done

    cat > montage.py <<'MONTAGE'
import pathlib
import sys

from PIL import Image, ImageDraw, ImageFont

pair, plot_dir, out = sys.argv[1], pathlib.Path(sys.argv[2]), sys.argv[3]

# lohhlaplot's order: the two coverage tracks, their ratio, then the derived signals.
# Anything it grows later still gets tiled, after these.
KINDS = ["t_dp", "n_dp", "tn_dp", "logR", "baf"]

panels = {}
for png in plot_dir.glob("*.png"):
    gene, _, kind = png.stem.partition(".")
    panels.setdefault(gene, {})[kind] = png
if not panels:
    sys.exit(f"No rasterised panels found for {pair}")

genes = sorted(panels)
present = {k for gene in panels.values() for k in gene}
kinds = [k for k in KINDS if k in present] + sorted(present - set(KINDS))

CELL_W, LABEL_H, PAD = 900, 34, 8
probe = Image.open(next(iter(panels[genes[0]].values())))
CELL_H = round(CELL_W * probe.height / probe.width)
font = ImageFont.load_default(size=22)

sheet = Image.new(
    "RGB",
    (PAD + len(genes) * (CELL_W + PAD), PAD + len(kinds) * (CELL_H + LABEL_H + PAD)),
    "white",
)
draw = ImageDraw.Draw(sheet)

for row, kind in enumerate(kinds):
    for col, gene in enumerate(genes):
        x = PAD + col * (CELL_W + PAD)
        y = PAD + row * (CELL_H + LABEL_H + PAD)
        draw.text((x, y + 6), f"{gene}  {kind}", fill="black", font=font)
        src = panels[gene].get(kind)
        if src is None:
            continue
        panel = Image.open(src).convert("RGB").resize((CELL_W, CELL_H), Image.LANCZOS)
        sheet.paste(panel, (x, y + LABEL_H))

sheet.save(out, optimize=True)
MONTAGE

    python3 montage.py ${somatic_name} panels ${somatic_name}_hla_loh_mqc.png

    # The result table is keyed by HLA gene alone - hla_a, hla_b, hla_c - which is unique
    # within a pair and not across them. MultiQC reads the first column as the row key and
    # merges every file it matches into one table, so without the pair in the key a second
    # pair overwrites the first rather than appending to it.
    awk -F'\t' -v OFS='\t' -v pair="${somatic_name}" '
        NR == 1 { \$1 = "Pair_HLAGene"; print; next }
        { \$1 = pair "_" \$1; print }
    ' ${loh_res} > ${somatic_name}_lohres.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pdftoppm: \$(pdftoppm -v 2>&1 | head -1 | sed 's/pdftoppm version //')
        pillow: \$(python3 -c "import PIL; print(PIL.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${somatic_name}_hla_loh_mqc.png
    touch ${somatic_name}_lohres.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pdftoppm: 26.09.0
    END_VERSIONS
    """
}
