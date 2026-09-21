#!/usr/bin/env python3
"""Tile a set of PNG panels into one sheet for the MultiQC report.

MultiQC gives every image its own report section, so a pair with fifteen LOH
panels or eight ASCAT plots buries everything else. Tiling them into one sheet
per pair keeps the report navigable, and puts panels that are read against each
other - the three HLA genes, the before and after GC correction plots - side by
side.

Two layouts:

  --layout grid   panel names are <col>.<row>.png, e.g. hla_a.logR.png, and the
                  grid is the cross product. Column and row order follow
                  --row-order where given, with anything unlisted appended.
  --layout flat   panels have no such structure, so they are laid out in reading
                  order across --cols columns.

Panels keep their aspect ratio: ASCAT emits a square sunrise plot alongside wide
genome profiles, and forcing those into one cell shape distorts them.
"""

import argparse
import pathlib
import sys

from PIL import Image, ImageDraw, ImageFont

CELL_W = 900
LABEL_H = 34
PAD = 8


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panels", nargs="+", type=pathlib.Path)
    parser.add_argument("--out", required=True)
    parser.add_argument("--layout", choices=["grid", "flat"], default="flat")
    parser.add_argument("--cols", type=int, default=2, help="flat layout only")
    parser.add_argument(
        "--row-order",
        default="",
        help="grid layout only: comma-separated row order, unlisted rows appended",
    )
    parser.add_argument(
        "--strip-prefix",
        default="",
        help="dropped from the start of every label, e.g. the pair name",
    )
    return parser.parse_args()


def label_for(path, strip_prefix):
    stem = path.stem
    if strip_prefix and stem.startswith(strip_prefix):
        stem = stem[len(strip_prefix) :]
    return stem.strip("._- ").replace(".", " ") or path.stem


def grid_cells(panels, row_order):
    """(column, row) placement from <col>.<row>.png names."""
    table = {}
    for path in panels:
        col, _, row = path.stem.partition(".")
        table.setdefault(col, {})[row] = path
    cols = sorted(table)
    present = {row for entry in table.values() for row in entry}
    preferred = [r for r in row_order if r in present]
    rows = preferred + sorted(present - set(preferred))
    return [
        (c, r, table[col].get(row))
        for r, row in enumerate(rows)
        for c, col in enumerate(cols)
    ], len(cols), len(rows)


def flat_cells(panels, cols):
    ordered = sorted(panels, key=lambda p: p.stem)
    placed = [(i % cols, i // cols, path) for i, path in enumerate(ordered)]
    rows = (len(ordered) + cols - 1) // cols
    return placed, min(cols, len(ordered)), rows


def main():
    args = parse_args()
    panels = [p for p in args.panels if p.exists()]
    if not panels:
        sys.exit("No panels to tile")

    row_order = [r for r in args.row_order.split(",") if r]
    if args.layout == "grid":
        cells, n_cols, n_rows = grid_cells(panels, row_order)
    else:
        cells, n_cols, n_rows = flat_cells(panels, args.cols)

    # Row heights are per row, not one height for the whole sheet. ASCAT emits a square
    # sunrise plot beside genome-wide profiles three times as wide as they are tall; a
    # single height sized to the tallest leaves every wide panel floating in whitespace.
    scaled_h = {
        p: round(CELL_W * Image.open(p).height / Image.open(p).width) for p in panels
    }
    row_h = {}
    for _col, row, path in cells:
        if path is not None:
            row_h[row] = max(row_h.get(row, 0), scaled_h[path])
    row_y, offset = {}, PAD
    for row in range(n_rows):
        row_y[row] = offset
        offset += row_h.get(row, 0) + LABEL_H + PAD

    font = ImageFont.load_default(size=22)
    sheet = Image.new("RGB", (PAD + n_cols * (CELL_W + PAD), offset), "white")
    draw = ImageDraw.Draw(sheet)

    for col, row, path in cells:
        if path is None:
            continue
        x = PAD + col * (CELL_W + PAD)
        y = row_y[row]
        draw.text((x, y + 6), label_for(path, args.strip_prefix), fill="black", font=font)
        panel = Image.open(path).convert("RGB")
        scale = min(CELL_W / panel.width, row_h[row] / panel.height)
        panel = panel.resize(
            (round(panel.width * scale), round(panel.height * scale)), Image.LANCZOS
        )
        sheet.paste(panel, (x + (CELL_W - panel.width) // 2, y + LABEL_H))

    sheet.save(args.out, optimize=True)


if __name__ == "__main__":
    main()
