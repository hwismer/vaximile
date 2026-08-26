process HLA_CALLS_PVAC {

    cpus 2
    memory "4GB"

    conda "python=3.10 pandas=2.1"

    tag "HLA consensus calls for ${meta.sample_name}"

    input:
        tuple val(meta), path(optitype_result), path(optitype_pdf), path(hlahd_result)

    output:
		tuple val(meta), path("${meta.sample_name}_hla_calls.csv"), emit: pvac_calls

    script:
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
            
    with open("${meta.sample_name}_hla_calls.csv", "w", newline="",encoding="utf-8") as f:
        writer = csv.writer(f,lineterminator="\\n")
        writer.writerow(alleles)

    """
}
