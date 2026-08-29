process HLA_CALLS_PVAC {

    label 'process_low'

    conda "python=3.10 pandas=2.1"

    tag "HLA consensus calls for ${meta.sample_name}"

    input:
        tuple val(meta), path(optitype_result), path(optitype_pdf), path(hlahd_result)

    output:
		tuple val(meta), path("*_hla_calls.csv"), emit: pvac_calls
		path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import csv
    
    class_i = {"A":{}, "B":{}, "C":{}}
    class_ii = {"DRB1":{}, "DQA1":{}, "DQB1":{}, "DPA1":{}, "DPB1":{}}
    
    opti = pd.read_csv("${optitype_result}", sep = "\\t")
    
    # Optitype Parsing
    for col in ["A1", "A2", "B1", "B2", "C1", "C2"]:
        short_allele = opti[col][0]
        for allele in class_i:
            if short_allele.startswith(allele):
                long_allele = "HLA-" + short_allele
                if long_allele in class_i[allele]:
                    class_i[allele][long_allele] += 1
                else:
                    class_i[allele][long_allele] = 1
                    
    with open("${hlahd_result}") as f:
        for line in f:
            line = line.rstrip("\\n")
            line = line.split("\\t")[:3]

            allele = line[0]
            type1 = line[1]
            type2 = line[2]
            
            if type1 != "Not typed":
                type1 = type1.split("*")[0] + "*" + type1.split("*")[1][:5]
                if type2 == "-":
                    type2 = type1
                else:
                    type2 = type2.split("*")[0] + "*" + type2.split("*")[1][:5]
                    
                if allele in class_i:
                    if type1 in class_i[allele]:
                        class_i[allele][type1] += 1
                    else:
                        class_i[allele][type1] = 1
                        
                    if type2 in class_i[allele]:
                        class_i[allele][type2] += 1
                    else:
                        class_i[allele][type2] = 1
                
                elif allele in class_ii:
                    if type1 in class_ii[allele]:
                        class_ii[allele][type1] += 1
                    else:
                        class_ii[allele][type1] = 1
                        
                    if type2 in class_ii[allele]:
                        class_ii[allele][type2] += 1
                    else:
                        class_ii[allele][type2] = 1
                        
    alleles = []
    
    for allele in class_i:
        top2 = sorted(class_i[allele].items(), key=lambda x: x[1], reverse=True)[:2]
        allele_set = {k for k, v in top2}
        for a in allele_set:
            alleles.append(a)
            
    for allele in class_ii:
        top2 = sorted(class_ii[allele].items(), key=lambda x: x[1], reverse=True)[:2]
        allele_set = {k for k, v in top2}
        for a in allele_set:
            alleles.append(a)
            
    with open("${prefix}_hla_calls.csv", "w", newline="",encoding="utf-8") as f:
        writer = csv.writer(f,lineterminator="\\n")
        writer.writerow(alleles)

    # versions.yml must be written by this interpreter: the script block runs under
    # python, so a bash heredoc here would be a syntax error.
    with open("versions.yml", "w") as _vf:
        _vf.write('"${task.process}":\\n    python: 3.10\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    """
    touch "${prefix}_hla_calls.csv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """
}
