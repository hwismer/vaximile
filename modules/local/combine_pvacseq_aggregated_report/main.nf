process COMBINE_PVACSEQ_AGGREGATED_REPORT {

    label 'process_low'
    conda "python=3.10 pandas=2.1"
    
    tag "Combing pVACseq reports for ${patient}"

    input:
        tuple val(patient), path(reports)

    output:
        tuple val(patient), path("*_pvacseq_reports.tsv"), emit: tsv
        path "versions.yml", topic: versions

    script:

    def prefix = task.ext.prefix ?: "${patient}"
    def file_list = reports.collect { "'${it.name}'" }.join(", ")

    """
    #!/usr/bin/env python3
    
    import pandas as pd
    
    files = [${file_list}]
    
    # CREATING MERGED DATAFRAME
    df = None
    for f in files:
        if df is None:
            df = pd.read_csv(f, sep = "\\t")
            df["sample"] = f.split(".")[0]
        else:
            new_df = pd.read_csv(f, sep = "\\t")
            new_df["sample"] = f.split(".")[0]
            df = pd.concat([df,new_df])
            
    # SORTING
    tier_order = ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "PoorTranscript", 
        "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
    
    df["Tier"] = pd.Categorical(df["Tier"], categories=tier_order, ordered=True)
    
    df["sum_rank"] = (df["Allele Expr"].rank(method="min") + df["%ile MT"].rank(method="min") + df["IC50 MT"].rank(method="min"))
    
    df = df.sort_values(["Tier", "sum_rank", "IC50 MT", "Gene", "AA Change"],kind="mergesort").drop(columns="sum_rank").reset_index()
    df.to_csv("${prefix}_pvacseq_reports.tsv", sep = "\\t", index = False)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """

    stub:

    def prefix = task.ext.prefix ?: "${patient}"

    """
    touch ${prefix}_pvacseq_reports.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """

}
