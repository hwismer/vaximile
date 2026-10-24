process HLA_CALLS_PVAC {

    // HLA alleles for pVACtools, and a MultiQC table of them: class I from OptiType, class II from
    // HLA-HD. The table also notes class I genes where HLA-HD disagrees. Calls are trimmed to two
    // fields and null alleles dropped.

    label 'process_low'

    conda "python=3.10"

    tag "HLA alleles for ${meta.sample_name}"

    input:
        tuple val(meta), path(optitype_result), path(optitype_pdf), path(hlahd_result)

    output:
        tuple val(meta), path("*_hla_calls.csv"), emit: pvac_calls
        tuple val(meta), path("*_hla_alleles_mqc.json"), emit: mqc

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    def library = "${meta.sample_name}_${meta.molecule}"
    """
    #!/usr/bin/env python3

    import csv
    import json
    import re

    CLASS_I = ["A", "B", "C"]
    CLASS_II = ["DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"]

    # Two-field names without "HLA-": "HLA-A*02:01:01" -> "A*02:01". Untyped and null (N) alleles are dropped.
    def two_field(*names):
        matches = [(re.match(r"(?:HLA-)?(\\w+\\*\\d+:\\d+)", n.strip()), n.strip()) for n in names]
        return sorted(m.group(1) for m, n in matches if m and not n.endswith("N"))

    with open("${optitype_result}") as fh:
        row = next(csv.DictReader(fh, delimiter="\\t"))
    optitype = {g: two_field(row[f"{g}1"], row[f"{g}2"]) for g in CLASS_I}

    # HLA-HD: gene, allele 1, allele 2 ("-" when homozygous).
    hlahd = {}
    with open("${hlahd_result}") as fh:
        for line in fh:
            gene, first, second = (line.rstrip("\\n").split("\\t") + ["", ""])[:3]
            hlahd[gene] = two_field(first, first if second.strip() == "-" else second)

    notes = [f"{g}: HLA-HD {', '.join(hlahd.get(g, [])) or 'none'}"
             for g in CLASS_I if optitype[g] != hlahd.get(g, [])]

    # pVACtools names class I alleles with "HLA-" and class II without, and pairs the DQ and DP chains itself.
    class_i = sorted({"HLA-" + a for g in CLASS_I for a in optitype[g]})
    class_ii = sorted({a for g in CLASS_II for a in hlahd.get(g, [])})
    if not class_i:
        raise SystemExit("No class I alleles typed by OptiType")
    with open("${prefix}_hla_calls.csv", "w") as fh:
        fh.write(",".join(class_i + class_ii) + "\\n")

    row = {g: ", ".join(optitype[g]) for g in CLASS_I}
    row["class_ii"] = ", ".join(class_ii)
    row["notes"] = "; ".join(notes) or "HLA-HD agrees"
    content = {
        "id": "hla_pvacseq_alleles",
        "section_name": "pVACseq alleles",
        "description": "Alleles passed to pVACseq, at two-field resolution with null alleles removed: class I from OptiType, class II from HLA-HD. Notes lists class I genes where HLA-HD's calls differ.",
        "parent_id": "hla",
        "parent_name": "HLA",
        "plot_type": "table",
        "pconfig": {"id": "hla_pvacseq_alleles", "title": "HLA: pVACseq alleles"},
        "headers": {
            "A": {"title": "HLA-A (OptiType)"},
            "B": {"title": "HLA-B (OptiType)"},
            "C": {"title": "HLA-C (OptiType)"},
            "class_ii": {"title": "Class II (HLA-HD)"},
            "notes": {"title": "Notes"},
        },
        "data": {"${library}": row},
    }
    with open("${prefix}_hla_alleles_mqc.json", "w") as fh:
        json.dump(content, fh, indent=1)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch "${prefix}_hla_calls.csv" "${prefix}_hla_alleles_mqc.json"
    """
}
