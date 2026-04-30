
process HLA_COMBINE_FASTQS {

    cpus 1
    memory "8GB"

    input:
        tuple val(sample1_meta), path(sample1_fastq1), path(sample1_fastq2)
        tuple val(sample2_meta), path(sample2_fastq1), path(sample2_fastq2)

    output:
        tuple val(sample1_meta),
            path("merged_${sample1_meta.sample_name}_${sample2_meta.sample_name}_R1.fastq.gz"),
            path("merged_${sample1_meta.sample_name}_${sample2_meta.sample_name}_R1.fastq.gz")


    script:
    """
    cat $sample1_fastq1 $sample2_fastq1 > "merged_${sample1_meta.sample_name}_${sample2_meta.sample_name}_R1.fastq.gz"
    cat $sample1_fastq2 $sample2_fastq2 > "merged_${sample1_meta.sample_name}_${sample2_meta.sample_name}_R2.fastq.gz"

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
