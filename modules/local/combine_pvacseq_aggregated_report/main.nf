process COMBINE_PVACSEQ_AGGREGATED_REPORT {

    label 'process_low'
    conda "python=3.10 pandas=2.1"
    
    tag "Combing pVACseq reports for ${patient}"

    input:
        tuple val(patient), path(reports), path(vcfs)

    output:
        tuple val(patient), path("*_pvacseq_reports.tsv"), emit: tsv

    script:

    def prefix = task.ext.prefix ?: "${patient}"
    def report_list = [reports].flatten().collect { "'${it}'" }.join(", ")
    def vcf_list = [vcfs].flatten().collect { "'${it}'" }.join(", ")

    """
    #!/usr/bin/env python3
    
    import gzip
    import pandas as pd

    def read_ad(vcf):
        # AD for every sample in the VCF, keyed the way pVACseq builds its IDs:
        # chrom-start-stop-ref-alt, with start = POS - 1 for SNVs/MNVs and POS for indels.
        ad = {}
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                fields = line.rstrip("\\n").split("\\t")
                if line.startswith("#"):
                    samples = fields[9:]
                    continue
                chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                start = pos - 1 if len(ref) == len(alt) else pos
                fmt = fields[8].split(":")
                ad[(chrom, start, ref, alt)] = {
                    sample: dict(zip(fmt, values.split(":"))).get("AD")
                    for sample, values in zip(samples, fields[9:])
                }
        return samples, ad

    def report_key(variant_id):
        chrom, start, _stop, ref, alt = variant_id.rsplit("-", 4)
        return (chrom, int(start), ref, alt)

    # MERGING REPORTS, WITH EACH SAMPLE'S AD FROM THE PAIR'S SOMATIC VCF
    dfs = []
    for report, vcf in zip([${report_list}], [${vcf_list}]):
        report_df = pd.read_csv(report, sep = "\\t")
        report_df["sample"] = report.split(".")[0]
        samples, ad = read_ad(vcf)
        keys = report_df["ID"].map(report_key)
        for sample in samples:
            report_df[f"{sample} AD"] = keys.map(lambda k: ad.get(k, {}).get(sample))
        dfs.append(report_df)
    df = pd.concat(dfs)
            
    # SORTING
    tier_order = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "PoorTranscript", 
        "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
    
    df["Tier"] = pd.Categorical(df["Tier"], categories=tier_order, ordered=True)
    
    df["sum_rank"] = (df["Allele Expr"].rank(method="min", ascending=False) + df["%ile MT"].rank(method="min") + df["IC50 MT"].rank(method="min"))
    
    df = df.sort_values(["Tier", "sum_rank", "IC50 MT", "Gene", "AA Change"],kind="mergesort").drop(columns="sum_rank").reset_index(drop=True)

    # pVACseq's ML call, when the classifier ran, next to the tier.
    if "ML Prediction (score)" in df.columns:
        ml = df.pop("ML Prediction (score)")
        df.insert(df.columns.get_loc("Tier") + 1, "ML Prediction (score)", ml)
    df.to_csv("${prefix}_pvacseq_reports.tsv", sep = "\\t", index = False)
    """

    stub:

    def prefix = task.ext.prefix ?: "${patient}"

    """
    touch ${prefix}_pvacseq_reports.tsv
    """

}
