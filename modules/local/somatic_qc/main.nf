process SOMATIC_QC {

    // Filtering metrics, substitution types, TMB, contamination and pVACseq tiers for one pair,
    // and the points for the patient's VAF plots, which VAF_GRIDS draws.

    label 'process_low'
    conda "python=3.10 pandas=2.1"

    tag "Somatic QC for ${meta.somatic_name}"

    input:
        tuple val(meta), path(caller_vcfs), path(consensus_vcf, stageAs: 'consensus/*'),
            path(final_vcf, stageAs: 'final/*'), path(contamination), path(capture_bed), path(pvacseq_report)

    output:
        tuple val(meta), path("*_mqc.json"), emit: mqc
        tuple val(meta), path("*_vaf.tsv"), emit: vaf_tables

    script:
    def callers = [caller_vcfs].flatten().collect { "'${it}'" }.join(", ")
    def pvacseq = pvacseq_report ? "'${pvacseq_report}'" : "None"
    """
    python3 - <<'PY'
import csv
import gzip
import json
import statistics
from collections import Counter, namedtuple

PAIR = "${meta.somatic_name}"
TUMOR = "${meta.tumor_meta.sample_name}"
NORMAL = "${meta.normal_meta.sample_name}"
CALLER_VCFS = [${callers}]
CONSENSUS = "${consensus_vcf}"
FINAL = "${final_vcf}"
CONTAMINATION = "${contamination}"
CAPTURE_BED = "${capture_bed}"
PVACSEQ = ${pvacseq}

SOMATIC = {"id": "somatic", "name": "Somatic",
           "description": "Somatic calls for each tumor/normal pair, from the three callers through to the pVACseq candidates"}
NEOANTIGENS = {"id": "neoantigens", "name": "Neoantigens",
               "description": "pVACseq MHC class I results for each tumor/normal pair"}
CALLER_NAMES = {"mutect": "Mutect2", "deepsomatic": "DeepSomatic", "strelka": "Strelka"}
NONSYNONYMOUS = {"missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
                 "inframe_insertion", "inframe_deletion", "protein_altering_variant"}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
SPECTRUM = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
TIERS = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "ProbPos",
         "PoorTranscript", "Anchor", "Subclonal", "LowExpr", "NoExpr", "Poor"]

Record = namedtuple("Record", "chrom pos ref alt info calls")


# Sample names, VEP CSQ field names (or None), and records.
def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    samples, csq_fields, records = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = line.split("Format: ")[1].split('"')[0].split("|")
                continue
            fields = line.rstrip("\\n").split("\\t")
            if line.startswith("#"):
                samples = fields[9:]
                continue
            fmt = fields[8].split(":") if len(fields) > 8 else []
            calls = {s: dict(zip(fmt, v.split(":"))) for s, v in zip(samples, fields[9:])}
            records.append(Record(fields[0], int(fields[1]), fields[3], fields[4], fields[7], calls))
    return samples, csq_fields, records


# VAF from AD (ref,alt), or None when AD is missing or has no reads.
def ad_vaf(call):
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


# Tumor DNA VAF: the bam-readcount AF if present, otherwise from AD.
def dna_vaf(call):
    vaf = field_float(call, "AF")
    return vaf if vaf is not None else ad_vaf(call)


def info_value(info, name):
    for item in info.split(";"):
        if item.startswith(name + "="):
            return item[len(name) + 1:]
    return None


def is_snv(record):
    return len(record.ref) == 1 and len(record.alt) == 1


# Total size of a BED file in base pairs, with overlapping intervals merged.
def bed_size(path):
    opener = gzip.open if path.endswith(".gz") else open
    intervals = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            chrom, start, end = line.split("\\t")[:3]
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


# Points for one panel of a VAF plot, drawn by VAF_GRIDS. A row with no point keeps an empty panel.
def write_panel(path, xlab, ylab, points):
    with open(path, "w") as fh:
        fh.write("panel\\txlab\\tylab\\tx\\ty\\n")
        fh.writelines(f"{PAIR}\\t{xlab}\\t{ylab}\\t{x:.4f}\\t{'' if y is None else f'{y:.4f}'}\\n" for x, y in points)
        if not points:
            fh.write(f"{PAIR}\\t{xlab}\\t{ylab}\\t\\t\\n")


# A MultiQC custom-content file. Files with the same plot id from different pairs are merged into one plot or table.
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


_, _, consensus = read_vcf(CONSENSUS)
_, csq_fields, final = read_vcf(FINAL)

# VAF plot points. Tumor VAF histogram: the final calls, after every filter.
final_vafs = [v for v in (dna_vaf(r.calls.get(TUMOR, {})) for r in final) if v is not None]
write_panel(f"{PAIR}_somatic_vaf.tsv", f"{TUMOR} DNA VAF", "Variants", [(v, None) for v in final_vafs])

# Tumor vs normal VAF: the consensus calls before the population filter, which removes the
# germline leakage this plot is meant to show.
normal_pairs = []
for r in consensus:
    t = ad_vaf(r.calls.get(TUMOR, {}))
    n = ad_vaf(r.calls.get(NORMAL, {}))
    if t is not None and n is not None:
        normal_pairs.append((t, n))
write_panel(f"{PAIR}_somatic_normal_vaf.tsv", f"{TUMOR} DNA VAF", f"{NORMAL} DNA VAF", normal_pairs)

# RNA vs DNA VAF of the final calls that have RNA coverage.
rna_pairs = []
for r in final:
    call = r.calls.get(TUMOR, {})
    d, rna = dna_vaf(call), field_float(call, "RAF")
    if d is not None and rna is not None and (field_float(call, "RDP") or 0) > 0:
        rna_pairs.append((d, rna))
write_panel(f"{PAIR}_somatic_rna_vaf.tsv", f"{TUMOR} DNA VAF", f"{TUMOR} RNA VAF", rna_pairs)

# Filtering metrics: variants left at each step, from each caller's PASS calls to pVACseq's Pass tier.
report = []
if PVACSEQ:
    with open(PVACSEQ) as fh:
        report = list(csv.DictReader(fh, delimiter="\\t"))
funnel = {}
for caller, name in CALLER_NAMES.items():
    for path in CALLER_VCFS:
        if path.endswith(f"_{caller}_variants.vcf.gz"):
            funnel[f"{name} PASS"] = len(read_vcf(path)[2])
funnel["Consensus"] = len(consensus)
funnel["After population filter"] = len(final)
if PVACSEQ:
    funnel["With pVACseq peptides"] = len(report)
    funnel["pVACseq Pass"] = sum(1 for row in report if row["Tier"] == "Pass")
write_mqc(f"{PAIR}_somatic_funnel_mqc.json", "somatic_funnel", "Filtering Metrics",
          "Number of somatic variants left at each step, from each caller's PASS calls to pVACseq's Pass tier.",
          SOMATIC, "table", {PAIR: funnel},
          {"title": "Somatic: filtering metrics", "col1_header": "Pair"},
          {c: {"title": c, "format": "{:,.0f}", "scale": False} for c in funnel})

write_mqc(f"{PAIR}_somatic_spectrum_mqc.json", "somatic_spectrum", "Substitution Types",
          "SNV substitution classes of the final somatic calls. An excess of C>T can indicate FFPE damage and an "
          "excess of C>A oxidative damage.",
          SOMATIC, "bargraph", {PAIR: spectrum(final)},
          {"title": "Somatic: substitution types", "ylab": "SNVs", "cpswitch_counts_label": "SNVs",
           "cpswitch_c_active": False})

# Summary metrics.
consequence = csq_fields.index("Consequence") if csq_fields and "Consequence" in csq_fields else None
nonsynonymous = 0
if consequence is not None:
    for r in final:
        csq = info_value(r.info, "CSQ") or ""
        terms = {t for ann in csq.split(",") for t in ann.split("|")[consequence].split("&")}
        nonsynonymous += bool(terms & NONSYNONYMOUS)
capture_mb = bed_size(CAPTURE_BED) / 1e6
contamination = None
with open(CONTAMINATION) as fh:
    rows = list(csv.DictReader(fh, delimiter="\\t"))
    if rows:
        contamination = float(rows[0]["contamination"]) * 100
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
    "Median VAF": {"title": "Median VAF", "description": "Median tumor DNA VAF of the final calls", "format": "{:,.3f}", "max": 1},
    "TMB": {"title": "TMB", "description": f"Non-synonymous coding variants per Mb of capture region ({capture_mb:.1f} Mb)", "format": "{:,.2f}"},
    "Ti/Tv": {"title": "Ti/Tv", "description": "Transition/transversion ratio of the final SNVs", "format": "{:,.2f}"},
    "Contamination": {"title": "Contamination", "description": "Cross-sample contamination of the tumor, from GATK CalculateContamination", "suffix": "%", "format": "{:,.2f}"},
    "Normal VAF >= 5%": {"title": "Normal VAF >= 5%", "description": "Consensus calls with at least 5% alt reads in the normal; should be near 0", "suffix": "%", "format": "{:,.1f}", "max": 100},
}
write_mqc(f"{PAIR}_somatic_metrics_mqc.json", "somatic_metrics", "Summary Metrics",
          "Summary of the final somatic calls for each tumor/normal pair.",
          SOMATIC, "table", {PAIR: metrics}, {"title": "Somatic: metrics", "col1_header": "Pair"}, headers)

# pVACseq tiers and Pass candidates per HLA allele.
if report:
    tiers = Counter(row["Tier"] for row in report)
    write_mqc(f"{PAIR}_neoantigen_tiers_mqc.json", "neoantigen_tiers", "Tiers",
              "Number of variants in each pVACseq tier (MHC class I aggregated report).",
              NEOANTIGENS, "bargraph", {PAIR: {t: tiers.get(t, 0) for t in TIERS}},
              {"title": "Neoantigens: pVACseq tiers", "ylab": "Variants", "cpswitch_counts_label": "Variants"})
    alleles = Counter(row["Allele"] for row in report if row["Tier"] == "Pass")
    write_mqc(f"{PAIR}_neoantigen_alleles_mqc.json", "neoantigen_alleles", "Pass candidates per HLA allele",
              "Number of Pass-tier candidates whose best peptide binds each HLA allele. Compare with the HLA LOH "
              "results: candidates on a lost allele will not be presented by the tumor.",
              NEOANTIGENS, "bargraph", {PAIR: dict(sorted(alleles.items()))},
              {"title": "Neoantigens: Pass candidates per allele", "ylab": "Candidates", "cpswitch_counts_label": "Candidates"})

# pVACseq's ML classifier: its call ("Accept (0.71)"), present only when the classifier ran.
if report and "ML Prediction (score)" in report[0]:
    calls = Counter(row["ML Prediction (score)"].split(" ")[0] for row in report)
    write_mqc(f"{PAIR}_neoantigen_ml_mqc.json", "neoantigen_ml", "ML predictions",
              "pVACseq's machine-learning call for each variant: Accept, Review or Reject. NA means the "
              "classifier could not score the variant.",
              NEOANTIGENS, "bargraph", {PAIR: {c: calls.get(c, 0) for c in ["Accept", "Review", "Reject", "NA"]}},
              {"title": "Neoantigens: ML predictions", "ylab": "Variants", "cpswitch_counts_label": "Variants"})
PY
    """

    stub:
    """
    touch ${meta.somatic_name}_somatic_metrics_mqc.json ${meta.somatic_name}_somatic_vaf.tsv ${meta.somatic_name}_somatic_normal_vaf.tsv ${meta.somatic_name}_somatic_rna_vaf.tsv
    """
}
