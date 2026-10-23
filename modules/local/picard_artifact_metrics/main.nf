process PICARD_ARTIFACT_METRICS {

    // Library-level DNA damage, for MultiQC: 8-oxoguanine oxidation (G>T, read as C>A on the
    // other strand) from CollectOxoGMetrics, and the per-substitution error rates from
    // CollectSequencingArtifactMetrics, summarised into one table row.

    label 'process_low'
    conda "bioconda::gatk4=4.6.1.0"
    container "broadinstitute/gatk:4.6.1.0"

    tag "Artifact metrics on ${meta.sample_name}"

    input:
        tuple val(meta), path(bam), path(bai)
        tuple path(reference_fa), path(reference_index)
        path reference_dict

    output:
        tuple val(meta), path("*.oxog_metrics.txt"), path("*_artifacts_mqc.json"), emit: mqc

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    gatk CollectOxoGMetrics \
        -I $bam \
        -R $reference_fa \
        -O ${prefix}.oxog_metrics.txt

    gatk CollectSequencingArtifactMetrics \
        -I $bam \
        -R $reference_fa \
        -O ${prefix} \
        --FILE_EXTENSION .txt

    python3 - <<'PY'
import json

# Picard metrics files: comment lines, then a tab-separated table, then a blank line.
def read_metrics(path):
    with open(path) as fh:
        lines = [l.rstrip("\\n") for l in fh if l.strip() and not l.startswith("#")]
    header = lines[0].split("\\t")
    return [dict(zip(header, l.split("\\t"))) for l in lines[1:] if "\\t" in l]

def qscores(path):
    return {r["REF_BASE"] + ">" + r["ALT_BASE"]: float(r["TOTAL_QSCORE"]) for r in read_metrics(path)}

pre_adapter = qscores("${prefix}.pre_adapter_summary_metrics.txt")
bait_bias = qscores("${prefix}.bait_bias_summary_metrics.txt")
worst_pre = min(pre_adapter, key=pre_adapter.get)
worst_bait = min(bait_bias, key=bait_bias.get)

row = {
    "oxog_q": pre_adapter.get("G>T"),
    "deamination_q": pre_adapter.get("C>T"),
    "worst_pre_adapter_q": pre_adapter[worst_pre],
    "worst_pre_adapter": worst_pre,
    "worst_bait_bias_q": bait_bias[worst_bait],
    "worst_bait_bias": worst_bait,
}
q_column = {"format": "{:,.1f}", "scale": "RdYlGn", "min": 0, "max": 60}
content = {
    "id": "sequencing_artifacts",
    "section_name": "Sequencing Artifacts",
    "description": "Phred-scaled error rate of each base substitution from Picard CollectSequencingArtifactMetrics: lower means more artefactual substitutions. Pre-adapter errors arise before adapter ligation, such as 8-oxoguanine oxidation (G>T, seen as C>A on the other strand) and cytosine deamination (C>T, typical of FFPE); bait-bias errors come from hybrid capture.",
    "plot_type": "table",
    "pconfig": {"id": "sequencing_artifacts", "title": "Sequencing Artifacts"},
    "headers": {
        "oxog_q": {"title": "OxoG Q (G>T)", "description": "Pre-adapter Q for G>T, the 8-oxoguanine oxidation artifact", **q_column},
        "deamination_q": {"title": "Deamination Q (C>T)", "description": "Pre-adapter Q for C>T, the cytosine deamination artifact", **q_column},
        "worst_pre_adapter_q": {"title": "Worst pre-adapter Q", "description": "Lowest pre-adapter Q across all substitutions", **q_column},
        "worst_pre_adapter": {"title": "Worst pre-adapter", "description": "Substitution with the lowest pre-adapter Q"},
        "worst_bait_bias_q": {"title": "Worst bait-bias Q", "description": "Lowest bait-bias Q across all substitutions", **q_column},
        "worst_bait_bias": {"title": "Worst bait-bias", "description": "Substitution with the lowest bait-bias Q"},
    },
    "data": {"${prefix}": row},
}
with open("${prefix}_artifacts_mqc.json", "w") as fh:
    json.dump(content, fh, indent=1)
PY
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}_${meta.molecule}"
    """
    touch ${prefix}.oxog_metrics.txt ${prefix}_artifacts_mqc.json
    """
}
