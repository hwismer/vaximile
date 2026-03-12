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

    hlahd = pd.read_csv("$hla_result", sep = "\t", names = ["HLA", "Allele 1", "Allele 2"], nrows=21)
    hlahd = hlahd[(hlahd["Allele 1"] != "Not typed") & (hlahd["Allele 2"] != "Not typed")]

    allele_2_new = []
    allele_1_new = []
    for allele_1, allele_2 in zip(hlahd["Allele 1"], hlahd["Allele 2"]):
        allele_1_split = allele_1.split(":")
        allele_1 = allele_1_split[0] + ":" + allele_1_split[1]
        allele_1_new.append(allele_1)
                        
        if allele_2 == "-":
            allele_2_new.append(allele_1)
        else:
            allele_2_split = allele_2.split(":")
            allele_2 = allele_2_split[0] + ":" + allele_2_split[1]
            allele_2_new.append(allele_2)

    hlahd["Allele 1"] = allele_1_new
    hlahd["Allele 2"] = allele_2_new
    #hlahd = hlahd[hlahd["HLA"].isin(["A","B","C","DRB1","DQA1","DQB1"])]

    alleles = set()
    for allele, allele_1,allele_2 in zip(hlahd["HLA"],hlahd["Allele 1"], hlahd["Allele 2"]):
            if allele != "A" and allele != "B" and allele != "C":
                alleles.add(allele_1.split("-")[1])
                alleles.add(allele_2.split("-")[1])
            else:
                alleles.add(allele_1)
                alleles.add(allele_2)
                                                                    
    alleles = list(alleles)

    with open("${meta.sample_name}_hla_calls.csv", "w", newline="",encoding="utf-8") as f:
        writer = csv.writer(f,lineterminator="\\n")
        writer.writerow(alleles)

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
process POSTPROCESS_OPTITYPE {

    /*

    Parse optitype output folder to get just the tsv and pdf output

    */

    cpus 1
    memory "4GB"

    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/optitype/${meta.sample_name}_optitype", mode: "copy"

    input:
        tuple val(meta), path(optitype_output_dir)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype.tsv"), emit: result_tsv
        tuple val(meta), path("${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype_coverage.pdf"), emit: coverage_plot

    script:
        """
        RESULT_FILE=\$(find ${optitype_output_dir}/* -name "*_result.tsv" | head -n 1)    
        mv "\$RESULT_FILE" "${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype.tsv"
        
        COVERAGE_FILE=\$(find ${optitype_output_dir}/* -name "*_coverage_plot.pdf" | head -n 1)    
        mv "\$COVERAGE_FILE" "${meta.sample_name}_${meta.sample_type}_${meta.molecule}_optitype_coverage.pdf"

        """
}


process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    container "fred2/optitype:release-v1.3.1"

    cpus 8

    memory "150GB"

    input:
        tuple val(meta), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_optitype")

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
        
        python /usr/local/bin/OptiType/OptiTypePipeline.py \
            -i $fastq1 $fastq2 \
            --$molecule_flag \
            -c OptiType.ini \
            --outdir "${meta.sample_name}_${meta.molecule}_optitype"
        """
}


process HLAHD {

    /*

    Run HLA-HD to perform HLA typing on patient fastqs

    This has given me some trouble before with temp directories and issues where bowtie never initializes.
    This should be fixed by unzipping fastqs at the beginning of the script.
    One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.


    */
    
    cpus 6
    memory "24GB"

    container "griffithlab/hlahd:1.0"
    
    publishDir "${params.outdir}/${meta.somatic_sample}/HLA/hlahd/${meta.sample_name}_hlahd", mode: "copy"

    input:
        tuple val(meta), path(fastq1), path(fastq2)
    
    output:
        tuple val(meta), path("./${meta.sample_name}/result/${meta.sample_name}_final.result.txt"), emit: hla_calls
        path("./${meta.sample_name}/result/")

    script:
        """
        ulimit -n 1024
        echo "\$TMPDIR"
        #mkdir -p tmp
        #mkdir -p /tmp/

        #export TMPDIR=/tmp/

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
