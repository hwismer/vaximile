
process HLA_CALLS_PVAC {

    cpus 1
    memory "4GB"

    conda "python=3.10 pandas=2.1"

    input:
        tuple val(meta), path(optitype_result), path(optitype_pdf), path(hlahd_result)

    output:
		tuple val(meta), path("${meta.sample_name}_hla_calls.csv")

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


process HLAHD_HLA_CALLS {

    /*

    Parse the hla-hd final_result file to a .csv compatible with pvac (all hla alleles on single line)

    */

    cpus 1
    memory "2GB"

    conda "python=3.10 pandas=2.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/pvac_input", mode:"copy"

    input:
        tuple val(meta), path(hla_result)

    output:
        tuple val(meta), path("${meta.sample_name}_hla_calls.csv")

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import csv
    
    typed_alleles = []
    class_i_skips = ["HLA-E","HLA-F","HLA-G","HLA-H","HLA-J","HLA-K","HLA-L","HLA-V"]

    for line in open("${hla_result}", "r"):
        split_allele_line = line.strip().split("\t")
        for whole_allele in split_allele_line[1:]:
            if whole_allele != "-" and whole_allele != "Not typed":
                hla = whole_allele.split("*")[0]
                if hla not in class_i_skips:
                    hla_type = whole_allele.split("*")[1]
                    hla_type_pvac_res = hla_type[:5]
                    final_allele = hla + "*" + hla_type_pvac_res
                    typed_alleles.append(final_allele)

    with open("${meta.sample_name}_hla_calls.csv", "w", newline="",encoding="utf-8") as f:
        writer = csv.writer(f,lineterminator="\\n")
        writer.writerow(typed_alleles)

    """
        

}



process OPTITYPE_HLA_CALLS {

    /*

        Use the optitype tsv file to create pvac-compatible .csv (all hla alleles listed on single line)

    */
    
    cpus 1
    memory "4GB"

    conda "python=3.10 pandas=2.1"

    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/", mode: "copy"

    input:
        tuple val(meta), path(optitype_result_tsv)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_hla_pvacinput.csv")

    script:
        """
        #!/usr/bin/env python3

        import pandas as pd

        df = pd.read_csv("$optitype_result_tsv", sep = "\t")
        
        alleles = []
        for allele in ["A1", "A2", "B1", "B2", "C1", "C2"]:
            alleles.append("HLA-" + df[allele].iloc[0])
        allele_csv = ",".join(alleles)

        with open("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_hla_pvacinput.csv", "w") as f:
            print(allele_csv, file=f)

        """

}

process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    container "fred2/optitype:latest"
    cpus 16

    memory "128GB"

    input:
        tuple val(meta), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("optitype_out/*_result.tsv"), path("optitype_out/*_coverage_plot.pdf")


    script:

        def molecule_flag = meta.molecule.toLowerCase()
        """
        cat << EOF > OptiType.ini
        [mapping]
        razers3=/usr/local/bin/razers3
        threads=${task.cpus}
        [ilp]
        solver=cbc
        threads=${task.cpus}
        [behavior]
        deletebam=true
        unpaired_weight=0
        use_discordant=false
        EOF
        
        which python
        which OptiTypePipeline.py
        
        pwd
        ls -lh

        python /usr/local/bin/OptiType/OptiTypePipeline.py \
            -i $fastq1 $fastq2 \
            --$molecule_flag \
            -c OptiType.ini \
            --prefix "${meta.sample_name}_${meta.molecule}" \
            --outdir optitype_out
        """
}


process HLAHD {

    /*

    Run HLA-HD to perform HLA typing on patient fastqs

    This has given me some trouble before with temp directories and issues where bowtie never initializes.
    This should be fixed by unzipping fastqs at the beginning of the script.
    One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.


    */
    
    cpus 16
    memory "32GB"

    container "griffithlab/hlahd:1.0"
    
    //publishDir "${params.outdir}/${meta.somatic_sample}/HLA/hlahd/${meta.sample_name}_hlahd", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
    
    output:
        tuple val(meta), path("./${meta.sample_name}/result/${meta.sample_name}_final.result.txt"), emit: final_hla_calls
        tuple val(meta), path("./${meta.sample_name}/result/"), emit: result_dir

    script:
        """
        ulimit -n 1024
        
        gzip -dc $fastq1 > fastq_r1.fastq
        gzip -dc $fastq2 > fastq_r2.fastq

        /opt/hlahd/bin/hlahd.1.6.1.sh \
            -f /opt/hlahd/freq_data \
            -t $task.cpus \
            fastq_r1.fastq fastq_r2.fastq \
            /opt/hlahd/HLA_gene.split.txt \
            /opt/hlahd/dictionary \
            "${meta.sample_name}" \
            ./
        """

}
