process HLAHD_TO_TSV {

    label 'process_low'

    conda "python=3.10 pandas=2.1"

    tag "Converting HLAHD to TSV for  ${meta.sample_name}"

    input:
        tuple val(meta), val(sample_name), path(hlahd_result)

    output:
		tuple val(meta), path("*_hlahd.tsv"), emit: hlahd_tsv

    script:
    def prefix = task.ext.prefix ?: "${sample_name}"
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    
    test = pd.read_csv("${hlahd_result}", sep = "\\t", header = None, names=range(10))
    test = test.rename(columns = {0:"locus"})
    test.insert(0, "sample", "${sample_name}")
    test["calls"] = test.iloc[:, 2:].apply(
        lambda r: " - ".join([x for x in r if pd.notna(x) and x != "-" and x != "Not typed"]),
            axis=1
            )

    wide = (
        test.pivot(index="sample", columns="locus", values="calls")
              .reset_index()
              )
    wide.to_csv("${prefix}_hlahd.tsv", sep = "\\t", index = False)
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_name}"
    """
    touch ${prefix}_hlahd.tsv
    """
}
