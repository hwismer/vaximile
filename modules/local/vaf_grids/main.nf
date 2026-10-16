process VAF_GRIDS {

    // A patient's VAF histograms and scatter plots, one panel per pair or normal, two across.

    label 'process_low'
    conda "python=3.10 pandas=2.1 matplotlib=3.8"

    tag "VAF plots for ${patient}"

    input:
        tuple val(patient), path(tables)

    output:
        tuple val(patient), path("*_mqc.png"), emit: png, optional: true

    script:
    def table_list = [tables].flatten().collect { "'${it}'" }.join(", ")
    """
    python3 - <<'PY'
import csv
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATIENT = "${patient}"
TABLES = [${table_list}]

VAF_BIN = 0.05
COLS = 2
# Plots from the panel tables <pair or sample>_<plot>.tsv, written by SOMATIC_QC and GERMLINE_QC.
PLOTS = {"somatic_vaf": "histogram", "somatic_normal_vaf": "scatter", "somatic_rna_vaf": "scatter",
         "germline_vaf": "histogram"}


# One panel per pair or sample, at most two across, adding rows as needed. panels: (name, xlab, ylab, points).
def plot_grid(path, panels, plot):
    histogram = PLOTS[plot] == "histogram"
    cols = min(len(panels), COLS)
    rows = math.ceil(len(panels) / cols)
    width, height = (5.5, 3.6) if histogram else (5, 5.2)
    fig, axes = plt.subplots(rows, cols, figsize=(width * cols, height * rows), squeeze=False)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")
    for ax, (name, xlab, ylab, points) in zip(axes.flat, panels):
        xs = [x for x, _ in points]
        if histogram:
            ax.hist(xs, bins=[i * VAF_BIN for i in range(round(1 / VAF_BIN) + 1)], color="#4C72B0", edgecolor="white")
            ax.set_xlim(0, 1)
        else:
            if plot == "somatic_rna_vaf":
                ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=0.8)
            ax.scatter(xs, [y for _, y in points], s=14, alpha=0.6, color="#4C72B0", edgecolors="none")
            # Axes run just past 0 and 1 so points at the edges are not cut off.
            ax.set_xlim(-0.03, 1.03)
            ax.set_ylim(-0.03, 1.03)
            ax.set_aspect("equal")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(f"{name} (n = {len(points):,})", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


panels = {plot: {} for plot in PLOTS}
for path in TABLES:
    plot = next(p for p in PLOTS if path.endswith(f"_{p}.tsv"))
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\\t"))
    if rows:
        name = rows[0]["panel"]
        points = [(float(r["x"]), float(r["y"]) if r["y"] else None) for r in rows if r["x"]]
        panels[plot][name] = (name, rows[0]["xlab"], rows[0]["ylab"], points)
for plot, by_name in panels.items():
    if by_name:
        plot_grid(f"{PATIENT}_{plot}_mqc.png", [by_name[n] for n in sorted(by_name)], plot)
PY
    """

    stub:
    """
    touch ${patient}_somatic_vaf_mqc.png ${patient}_somatic_normal_vaf_mqc.png ${patient}_somatic_rna_vaf_mqc.png ${patient}_germline_vaf_mqc.png
    """
}
