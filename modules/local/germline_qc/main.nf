process GERMLINE_QC {

    // Caller overlap, het/hom and Ti/Tv for one normal, and the points for its VAF histogram,
    // which VAF_GRIDS draws.

    label 'process_low'
    conda "python=3.10 pandas=2.1"

    tag "Germline QC for ${meta.sample_name}"

    input:
        tuple val(meta), path(caller_vcfs), path(consensus_vcf, stageAs: 'consensus/*')

    output:
        tuple val(meta), path("*_mqc.json"), emit: mqc
        tuple val(meta), path("*_vaf.tsv"), emit: vaf_tables

    script:
    def callers = [caller_vcfs].flatten().collect { "'${it}'" }.join(", ")
    """
    python3 - <<'PY'
import gzip
import itertools
import json
import statistics
from collections import Counter, namedtuple

SAMPLE = "${meta.sample_name}_${meta.molecule}"
SAMPLE_COLUMN = "${meta.sample_name}"
CALLER_VCFS = [${callers}]
CONSENSUS = "${consensus_vcf}"

GERMLINE = {"id": "germline", "name": "Germline", "description": "Germline calls for each normal sample"}
CALLER_NAMES = {"deepvariant": "DeepVariant", "haplotypecaller": "HaplotypeCaller", "strelka": "Strelka"}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

Record = namedtuple("Record", "chrom pos ref alt calls")


def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    samples, records = [], []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if line.startswith("#"):
                samples = fields[9:]
                continue
            fmt = fields[8].split(":") if len(fields) > 8 else []
            calls = {s: dict(zip(fmt, v.split(":"))) for s, v in zip(samples, fields[9:])}
            records.append(Record(fields[0], int(fields[1]), fields[3], fields[4], calls))
    return records


# VAF from AD (ref,alt), or None when AD is missing or has no reads.
def ad_vaf(call):
    try:
        ref, alt = (int(x) for x in call["AD"].split(",")[:2])
    except (KeyError, ValueError):
        return None
    return alt / (ref + alt) if ref + alt else None


def is_snv(record):
    return len(record.ref) == 1 and len(record.alt) == 1


def ti_tv(records):
    snvs = [(r.ref, r.alt) for r in records if is_snv(r)]
    ti = sum(1 for s in snvs if s in TRANSITIONS)
    tv = len(snvs) - ti
    return round(ti / tv, 2) if tv else None


def genotype(call):
    alleles = [a for a in call.get("GT", "").replace("|", "/").split("/") if a not in ("", ".")]
    if len(alleles) < 2 or set(alleles) == {"0"}:
        return None
    return "het" if len(set(alleles)) > 1 else "hom"


# A MultiQC custom-content file. Files with the same plot id from different samples are merged into one plot or table.
def write_mqc(path, plot_id, section, description, plot_type, data, pconfig=None, headers=None):
    content = {
        "id": plot_id,
        "section_name": section,
        "description": description,
        "parent_id": GERMLINE["id"],
        "parent_name": GERMLINE["name"],
        "parent_description": GERMLINE["description"],
        "plot_type": plot_type,
        "pconfig": {"id": plot_id, **(pconfig or {})},
        "data": data,
    }
    if headers:
        content["headers"] = headers
    with open(path, "w") as fh:
        json.dump(content, fh, indent=1)


records = read_vcf(CONSENSUS)
vafs, het_vafs, gts = [], [], Counter()
for r in records:
    call = r.calls.get(SAMPLE_COLUMN, {})
    gt, vaf = genotype(call), ad_vaf(call)
    gts[gt] += 1
    if vaf is not None:
        vafs.append(vaf)
        if gt == "het":
            het_vafs.append(vaf)

# VAF histogram points, drawn by VAF_GRIDS with any other normals of the patient.
with open(f"{SAMPLE}_germline_vaf.tsv", "w") as fh:
    fh.write("panel\\txlab\\tylab\\tx\\ty\\n")
    fh.writelines(f"{SAMPLE}\\t{SAMPLE_COLUMN} VAF\\tVariants\\t{v:.4f}\\t\\n" for v in vafs)
    if not vafs:
        fh.write(f"{SAMPLE}\\t{SAMPLE_COLUMN} VAF\\tVariants\\t\\t\\n")

# Caller overlap: PASS calls by the combination of callers that made them.
sets = {}
for caller in CALLER_NAMES:
    for path in CALLER_VCFS:
        if path.endswith(f"_{caller}_variants.vcf.gz"):
            sets[caller] = {(r.chrom, r.pos, r.ref, r.alt) for r in read_vcf(path)}
combos = Counter(tuple(c for c in sets if k in sets[c]) for k in set().union(*sets.values()))
overlap = {}
for size in range(1, len(sets) + 1):
    for combo in itertools.combinations(sets, size):
        if size == 1:
            label = f"{CALLER_NAMES[combo[0]]} only"
        elif size == len(sets):
            label = "All three"
        else:
            label = " + ".join(CALLER_NAMES[c] for c in combo)
        overlap[label] = combos.get(combo, 0)
write_mqc(f"{SAMPLE}_germline_overlap_mqc.json", "germline_overlap", "Caller Overlap",
          "PASS calls from each germline caller, by which callers made them. Only calls made by at least two callers "
          "are kept in the consensus.",
          "bargraph", {SAMPLE: overlap},
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
write_mqc(f"{SAMPLE}_germline_metrics_mqc.json", "germline_metrics", "Summary Metrics",
          "Summary of the consensus germline calls for each normal.",
          "table", {SAMPLE: metrics}, {"title": "Germline: metrics", "col1_header": "Sample"}, headers)
PY
    """

    stub:
    """
    touch ${meta.sample_name}_${meta.molecule}_germline_metrics_mqc.json ${meta.sample_name}_${meta.molecule}_germline_vaf.tsv
    """
}
