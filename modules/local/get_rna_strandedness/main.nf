process GET_RNA_STRANDEDNESS {

    label 'process_single'

    conda "python=3.10 pandas=2.1"

    tag "Predicting RNA strandedness on ${meta.sample_name}"
    
    input:
        tuple val(meta), path(salmon_quant)

    output:
        tuple val(meta), path("*_strandedness.txt"), emit: strand_txt

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    #!/usr/bin/env python3

    import json
    import sys
    
    with open("${salmon_quant}/lib_format_counts.json", "r") as f:
        data = json.load(f)

    expected_format = data.get("expected_format")
    print(expected_format)
    if expected_format[1] == "U":
        strandedness = "XS"
    elif expected_format[1] == "S":
        if expected_format[2] == "R":
            strandedness = "RF"
        elif expected_format[2] == "F":
            strandedness = "FR"
    else:
        print("Unable to parse strandedness. Check salmon output")
        sys.exit(1)
        
    with open("${prefix}_strandedness.txt", "w") as f:
        f.write(strandedness + "\\n")

    """


}
