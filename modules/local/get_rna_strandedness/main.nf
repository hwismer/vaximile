process GET_RNA_STRANDEDNESS {

    label 'process_single'

    conda "python=3.10 pandas=2.1"

    tag "Predicting RNA strandedness on ${meta.sample_name}"
    
    input:
        tuple val(meta), path(salmon_quant)

    output:
        tuple val(meta), path("*_strandedness.txt"), emit: strand_txt
        path "versions.yml", topic: versions

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
            # Without this branch `strandedness` stays unbound and the script dies with a
            # NameError instead of the intended message.
            print("Unable to parse strandedness. Check salmon output")
            sys.exit(1)
    else:
        print("Unable to parse strandedness. Check salmon output")
        sys.exit(1)
        
    with open("${prefix}_strandedness.txt", "w") as f:
        f.write(strandedness + "\\n")

    # versions.yml must be written by this interpreter: the script block runs under
    # python, so a bash heredoc here would be a syntax error.
    with open("versions.yml", "w") as _vf:
        _vf.write('"${task.process}":\\n    python: 3.10\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch ${prefix}_strandedness.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """


}
