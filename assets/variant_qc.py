#!/usr/bin/env python3
"""Variant QC plots and metrics for the MultiQC report.

Writes MultiQC custom-content files (*_mqc.json). Files with the same plot id from
different pairs or samples are merged by MultiQC into one plot or table.

  somatic  - one tumour/normal pair: VAFs, caller overlap, filtering funnel,
             substitution spectrum, TMB, contamination, and the pVACseq tiers
  germline - one normal: VAF histogram, het/hom and Ti/Tv, caller overlap
"""

import argparse
import csv
import gzip
import itertools
import json
import statistics
from collections import Counter, namedtuple

Record = namedtuple("Record", "chrom pos ref alt info calls")

SOMATIC = {"id": "somatic", "name": "Somatic variants",
           "description": "Somatic calls for each tumour/normal pair, from the three callers through to the pVACseq candidates"}
GERMLINE = {"id": "germline", "name": "Germline variants",
            "description": "Germline calls for each normal sample"}
NEOANTIGENS = {"id": "neoantigens", "name": "Neoantigens",
               "description": "pVACseq MHC class I results for each tumour/normal pair"}

SOMATIC_CALLERS = {"mutect": "Mutect2", "deepsomatic": "DeepSomatic", "strelka": "Strelka"}
GERMLINE_CALLERS = {"deepvariant": "DeepVariant", "haplotypecaller": "HaplotypeCaller", "strelka": "Strelka"}

NONSYNONYMOUS = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
                 "inframe_insertion", "inframe_deletion", "protein_altering_variant"}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
SPECTRUM = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
TIERS = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "ProbPos",
         "PoorTranscript", "Anchor", "Subclonal", "LowExpr", "NoExpr", "Poor"]
VAF_BIN = 0.02


# ---------------------------------------------------------------------------
# Reading

def read_vcf(path):
    """Return the sample names, the VEP CSQ field names (or None), and the records."""
    opener = gzip.open if path.endswith(".gz") else open
    samples, csq_fields, records = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = line.split("Format: ")[1].split('"')[0].split("|")
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = fields[9:]
                continue
            fmt = fields[8].split(":") if len(fields) > 8 else []
            calls = {s: dict(zip(fmt, v.split(":"))) for s, v in zip(samples, fields[9:])}
            records.append(Record(fields[0], int(fields[1]), fields[3], fields[4], fields[7], calls))
    return samples, csq_fields, records


def key(record):
    return (record.chrom, record.pos, record.ref, record.alt)


def ad_vaf(call):
    """VAF from AD (ref,alt), or None when AD is missing or has no reads."""
    try:
        ref, alt = (int(x) for x in call["AD"].split(",")[:2])
    except (KeyError, ValueError):
        return None
    return alt / (ref + alt) if ref + alt else None


def field_float(call, name):
    try:
        return float(call[name].split(",")[0])
    except (KeyError, ValueError):
        return None


def dna_vaf(call):
    """Tumour DNA VAF: the bam-readcount AF if present, otherwise from AD."""
    vaf = field_float(call, "AF")
    return vaf if vaf is not None else ad_vaf(call)


def info_value(info, name):
    for item in info.split(";"):
        if item.startswith(name + "="):
            return item[len(name) + 1:]
    return None


def is_snv(record):
    return len(record.ref) == 1 and len(record.alt) == 1


def bed_size(path):
    """Total size of a BED file in base pairs, with overlapping intervals merged."""
    opener = gzip.open if path.endswith(".gz") else open
    intervals = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            chrom, start, end = line.split("\t")[:3]
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    total = 0
    for spans in intervals.values():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start > cur_end:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
            else:
                cur_end = max(cur_end, end)
        total += cur_end - cur_start
    return total


# ---------------------------------------------------------------------------
# Metrics

def histogram(values):
    """Counts of values in [0, 1] per bin, keyed by the bin's lower edge."""
    n = round(1 / VAF_BIN)
    counts = Counter(min(int(v / VAF_BIN), n - 1) for v in values)
    return {f"{i * VAF_BIN:.2f}": counts.get(i, 0) for i in range(n)}


def ti_tv(records):
    snvs = [(r.ref, r.alt) for r in records if is_snv(r)]
    ti = sum(1 for s in snvs if s in TRANSITIONS)
    tv = len(snvs) - ti
    return round(ti / tv, 2) if tv else None


def spectrum(records):
    counts = Counter()
    for r in records:
        if is_snv(r):
            ref, alt = r.ref, r.alt
            if ref in "AG":
                ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
            counts[f"{ref}>{alt}"] += 1
    return {c: counts.get(c, 0) for c in SPECTRUM}


def caller_overlap(caller_files, names):
    """Count variants by the combination of callers that called them."""
    sets = {}
    for path in caller_files:
        caller = next(c for c in names if path.endswith(f"_{c}_variants.vcf.gz"))
        sets[caller] = {key(r) for r in read_vcf(path)[2]}
    order = [c for c in names if c in sets]
    counts = Counter(tuple(c for c in order if k in sets[c]) for k in set().union(*sets.values()))
    data = {}
    for size in range(1, len(order) + 1):
        for combo in itertools.combinations(order, size):
            if size == 1:
                label = f"{names[combo[0]]} only"
            elif size == len(order):
                label = "All three"
            else:
                label = " + ".join(names[c] for c in combo)
            data[label] = counts.get(combo, 0)
    return data, {names[c]: len(sets[c]) for c in order}


def genotype(call):
    alleles = [a for a in call.get("GT", "").replace("|", "/").split("/") if a not in ("", ".")]
    if len(alleles) < 2 or set(alleles) == {"0"}:
        return None
    return "het" if len(set(alleles)) > 1 else "hom"


# ---------------------------------------------------------------------------
# Writing

def plot_vaf_histogram(path, vafs, title, xlabel):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.hist(vafs, bins=[i * VAF_BIN for i in range(round(1 / VAF_BIN) + 1)], color="#4C72B0", edgecolor="white")
    ax.set_xlim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Variants")
    ax.set_title(f"{title} (n = {len(vafs):,})")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_mqc(path, plot_id, section, description, parent, plot_type, data, pconfig=None, headers=None):
    if not data:
        return
    content = {
        "id": plot_id,
        "section_name": section,
        "description": description,
        "parent_id": parent["id"],
        "parent_name": parent["name"],
        "parent_description": parent["description"],
        "plot_type": plot_type,
        "pconfig": {"id": plot_id, **(pconfig or {})},
        "data": data,
    }
    if headers:
        content["headers"] = headers
    with open(path, "w") as fh:
        json.dump(content, fh, indent=1)


def vaf_scatter_points(pairs):
    return [{"x": round(x, 3), "y": round(y, 3)} for x, y in pairs]


# ---------------------------------------------------------------------------
# Modes

def somatic(args):
    pair = args.pair
    _, _, consensus = read_vcf(args.consensus)
    _, csq_fields, final = read_vcf(args.final)

    # Tumour VAF and tumour vs normal VAF, on the consensus calls before the population filter.
    tumour_vafs, normal_pairs = [], []
    for r in consensus:
        t = ad_vaf(r.calls.get(args.tumor, {}))
        n = ad_vaf(r.calls.get(args.normal, {}))
        if t is not None:
            tumour_vafs.append(t)
            if n is not None:
                normal_pairs.append((t, n))
    # A matplotlib image rather than a MultiQC plot, since MultiQC only draws horizontal bars.
    # The MultiQC step gives each pair's image its own section.
    plot_vaf_histogram(f"{pair}_somatic_vaf_mqc.png", tumour_vafs, f"{pair}: tumour VAF", "Tumour VAF")
    write_mqc(f"{pair}_somatic_normal_vaf_mqc.json", "somatic_normal_vaf", "VAF - tumour vs normal",
              "Tumour against normal DNA VAF of each consensus somatic call, before the population filter. True somatic "
              "calls lie on the x-axis. Calls with alt reads in the normal suggest germline leakage, tumour cells in the "
              "normal, or a normal that is not from the same person.",
              SOMATIC, "scatter", {pair: vaf_scatter_points(normal_pairs)},
              {"title": "Somatic: tumour vs normal VAF", "xlab": "Tumour VAF", "ylab": "Normal VAF",
               "xmin": 0, "xmax": 1, "ymin": 0, "ymax": 1})

    # RNA vs DNA VAF of the final calls.
    rna_pairs = []
    for r in final:
        call = r.calls.get(args.tumor, {})
        d, rna = dna_vaf(call), field_float(call, "RAF")
        if d is not None and rna is not None and (field_float(call, "RDP") or 0) > 0:
            rna_pairs.append((d, rna))
    # One plot per pair, so pairs don't overlap; the id makes each its own section.
    write_mqc(f"{pair}_somatic_rna_vaf_mqc.json", f"somatic_rna_vaf_{pair}", f"VAF - RNA vs DNA - {pair}",
              "Tumour DNA against RNA VAF of the final somatic calls that have RNA coverage. Mutations near the diagonal "
              "are expressed in proportion to their DNA fraction; those on the x-axis are not expressed.",
              SOMATIC, "scatter", {pair: vaf_scatter_points(rna_pairs)},
              {"title": f"Somatic: RNA vs DNA VAF - {pair}", "xlab": "Tumour DNA VAF", "ylab": "Tumour RNA VAF",
               "xmin": 0, "xmax": 1, "ymin": 0, "ymax": 1})

    # Caller overlap, and the funnel from caller calls to pVACseq candidates.
    overlap, per_caller = caller_overlap(args.callers, SOMATIC_CALLERS)
    write_mqc(f"{pair}_somatic_overlap_mqc.json", "somatic_overlap", "Calls - caller overlap",
              "PASS calls from each somatic caller, by which callers made them. Only calls made by at least two callers "
              "are kept in the consensus.",
              SOMATIC, "bargraph", {pair: overlap},
              {"title": "Somatic: caller overlap", "ylab": "Variants", "cpswitch_counts_label": "Variants"})

    report = []
    if args.pvacseq:
        with open(args.pvacseq) as fh:
            report = list(csv.DictReader(fh, delimiter="\t"))
    funnel = {f"{name} PASS": count for name, count in per_caller.items()}
    funnel["Consensus"] = len(consensus)
    funnel["After population filter"] = len(final)
    if args.pvacseq:
        funnel["With pVACseq peptides"] = len(report)
        funnel["pVACseq Pass"] = sum(1 for row in report if row["Tier"] == "Pass")
    write_mqc(f"{pair}_somatic_funnel_mqc.json", "somatic_funnel", "Calls - filtering funnel",
              "Number of somatic variants left at each step, from each caller's PASS calls to pVACseq's Pass tier.",
              SOMATIC, "table", {pair: funnel},
              {"title": "Somatic: filtering funnel", "col1_header": "Pair"},
              {c: {"title": c, "format": "{:,.0f}", "scale": False} for c in funnel})

    write_mqc(f"{pair}_somatic_spectrum_mqc.json", "somatic_spectrum", "Calls - substitution spectrum",
              "SNV substitution classes of the final somatic calls. An excess of C>T can indicate FFPE damage and an "
              "excess of C>A oxidative damage.",
              SOMATIC, "bargraph", {pair: spectrum(final)},
              {"title": "Somatic: substitution spectrum", "ylab": "SNVs", "cpswitch_counts_label": "SNVs",
               "cpswitch_c_active": False})

    # Summary metrics.
    consequence = csq_fields.index("Consequence") if csq_fields and "Consequence" in csq_fields else None
    nonsynonymous = 0
    if consequence is not None:
        for r in final:
            csq = info_value(r.info, "CSQ") or ""
            terms = {t for ann in csq.split(",") for t in ann.split("|")[consequence].split("&")}
            nonsynonymous += bool(terms & NONSYNONYMOUS)
    capture_mb = bed_size(args.bed) / 1e6
    contamination = None
    with open(args.contamination) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        if rows:
            contamination = float(rows[0]["contamination"]) * 100
    final_vafs = [v for v in (dna_vaf(r.calls.get(args.tumor, {})) for r in final) if v is not None]
    metrics = {
        "Variants": len(final),
        "SNVs": sum(1 for r in final if is_snv(r)),
        "Indels": sum(1 for r in final if len(r.ref) != len(r.alt)),
        "Median VAF": round(statistics.median(final_vafs), 3) if final_vafs else None,
        "TMB": round(nonsynonymous / capture_mb, 2) if capture_mb else None,
        "Ti/Tv": ti_tv(final),
        "Contamination": round(contamination, 2) if contamination is not None else None,
        "Normal VAF >= 5%": round(100 * sum(1 for _, n in normal_pairs if n >= 0.05) / len(normal_pairs), 1)
                            if normal_pairs else None,
    }
    headers = {
        "Variants": {"title": "Variants", "description": "Final somatic calls, after the population filter", "format": "{:,.0f}"},
        "SNVs": {"title": "SNVs", "format": "{:,.0f}"},
        "Indels": {"title": "Indels", "format": "{:,.0f}"},
        "Median VAF": {"title": "Median VAF", "description": "Median tumour DNA VAF of the final calls", "format": "{:,.3f}", "max": 1},
        "TMB": {"title": "TMB", "description": f"Non-synonymous coding variants per Mb of capture region ({capture_mb:.1f} Mb)", "format": "{:,.2f}"},
        "Ti/Tv": {"title": "Ti/Tv", "description": "Transition/transversion ratio of the final SNVs", "format": "{:,.2f}"},
        "Contamination": {"title": "Contamination", "description": "Cross-sample contamination of the tumour, from GATK CalculateContamination", "suffix": "%", "format": "{:,.2f}"},
        "Normal VAF >= 5%": {"title": "Normal VAF >= 5%", "description": "Consensus calls with at least 5% alt reads in the normal; should be near 0", "suffix": "%", "format": "{:,.1f}", "max": 100},
    }
    write_mqc(f"{pair}_somatic_metrics_mqc.json", "somatic_metrics", "Calls - summary metrics",
              "Summary of the final somatic calls for each tumour/normal pair.",
              SOMATIC, "table", {pair: metrics}, {"title": "Somatic: metrics", "col1_header": "Pair"}, headers)

    # pVACseq tiers and Pass candidates per HLA allele.
    if report:
        tiers = Counter(row["Tier"] for row in report)
        write_mqc(f"{pair}_neoantigen_tiers_mqc.json", "neoantigen_tiers", "Tiers",
                  "Number of variants in each pVACseq tier (MHC class I aggregated report).",
                  NEOANTIGENS, "bargraph", {pair: {t: tiers.get(t, 0) for t in TIERS}},
                  {"title": "Neoantigens: pVACseq tiers", "ylab": "Variants", "cpswitch_counts_label": "Variants"})
        alleles = Counter(row["Allele"] for row in report if row["Tier"] == "Pass")
        write_mqc(f"{pair}_neoantigen_alleles_mqc.json", "neoantigen_alleles", "Pass candidates per HLA allele",
                  "Number of Pass-tier candidates whose best peptide binds each HLA allele. Compare with the HLA LOH "
                  "results: candidates on a lost allele will not be presented by the tumour.",
                  NEOANTIGENS, "bargraph", {pair: dict(sorted(alleles.items()))},
                  {"title": "Neoantigens: Pass candidates per allele", "ylab": "Candidates", "cpswitch_counts_label": "Candidates"})


def germline(args):
    sample = args.sample
    _, _, records = read_vcf(args.consensus)
    vafs, het_vafs, gts = [], [], Counter()
    for r in records:
        call = r.calls.get(args.sample_column, {})
        gt, vaf = genotype(call), ad_vaf(call)
        gts[gt] += 1
        if vaf is not None:
            vafs.append(vaf)
            if gt == "het":
                het_vafs.append(vaf)
    write_mqc(f"{sample}_germline_vaf_mqc.json", "germline_vaf", "VAF",
              "VAF of the consensus germline calls. Heterozygous variants should peak at 0.5 and homozygous ones at 1. "
              "A heterozygous peak that is broad or shifted away from 0.5 suggests contamination or a mixed sample.",
              GERMLINE, "linegraph", {sample: histogram(vafs)},
              {"title": "Germline: VAF", "xlab": "VAF", "ylab": "Variants", "xmin": 0, "xmax": 1, "ymin": 0})

    overlap, _ = caller_overlap(args.callers, GERMLINE_CALLERS)
    write_mqc(f"{sample}_germline_overlap_mqc.json", "germline_overlap", "Calls - caller overlap",
              "PASS calls from each germline caller, by which callers made them. Only calls made by at least two callers "
              "are kept in the consensus.",
              GERMLINE, "bargraph", {sample: overlap},
              {"title": "Germline: caller overlap", "ylab": "Variants", "cpswitch_counts_label": "Variants"})

    metrics = {
        "Variants": len(records),
        "SNVs": sum(1 for r in records if is_snv(r)),
        "Indels": sum(1 for r in records if len(r.ref) != len(r.alt)),
        "Het": gts["het"],
        "Hom alt": gts["hom"],
        "Het/Hom": round(gts["het"] / gts["hom"], 2) if gts["hom"] else None,
        "Ti/Tv": ti_tv(records),
        "Median het VAF": round(statistics.median(het_vafs), 3) if het_vafs else None,
    }
    headers = {c: {"title": c, "format": "{:,.0f}"} for c in ("Variants", "SNVs", "Indels", "Het", "Hom alt")}
    headers["Het/Hom"] = {"title": "Het/Hom", "description": "Heterozygous to homozygous-alt ratio; typically 1.3-2", "format": "{:,.2f}"}
    headers["Ti/Tv"] = {"title": "Ti/Tv", "description": "Transition/transversion ratio; about 2.0-2.1 for genomes and 2.8-3.0 for exomes", "format": "{:,.2f}"}
    headers["Median het VAF"] = {"title": "Median het VAF", "description": "Median VAF of heterozygous calls; should be near 0.5", "format": "{:,.3f}", "max": 1}
    write_mqc(f"{sample}_germline_metrics_mqc.json", "germline_metrics", "Calls - summary metrics",
              "Summary of the consensus germline calls for each normal.",
              GERMLINE, "table", {sample: metrics}, {"title": "Germline: metrics", "col1_header": "Sample"}, headers)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    modes = parser.add_subparsers(dest="mode", required=True)

    s = modes.add_parser("somatic")
    s.add_argument("--pair", required=True)
    s.add_argument("--tumor", required=True, help="Tumour sample name in the VCFs")
    s.add_argument("--normal", required=True, help="Normal sample name in the VCFs")
    s.add_argument("--callers", nargs="+", required=True, help="Normalised PASS VCF of each caller")
    s.add_argument("--consensus", required=True, help="Consensus VCF, before the population filter")
    s.add_argument("--final", required=True, help="Final annotated VCF, after the population filter")
    s.add_argument("--contamination", required=True, help="GATK CalculateContamination table")
    s.add_argument("--bed", required=True, help="Capture region BED")
    s.add_argument("--pvacseq", help="pVACseq MHC class I aggregated report")

    g = modes.add_parser("germline")
    g.add_argument("--sample", required=True, help="Name to show in the report")
    g.add_argument("--sample-column", required=True, help="Sample name in the VCFs")
    g.add_argument("--callers", nargs="+", required=True, help="Normalised PASS VCF of each caller")
    g.add_argument("--consensus", required=True, help="Consensus germline VCF")

    args = parser.parse_args()
    {"somatic": somatic, "germline": germline}[args.mode](args)


if __name__ == "__main__":
    main()
