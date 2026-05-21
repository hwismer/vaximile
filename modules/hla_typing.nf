process MHC_REGION_FASTQS {
    
    cpus 4
    memory "16GB"

    conda "bioconda::samtools=1.23.1 bioconda::bedtools=2.31.1 bioconda::htslib=1.23.1"
    
    tag "Extracting MHC regions from ${meta.sample_name}"

    input:
    	tuple val(meta), path(bam), path(bai)

    output:
    	tuple val(meta), path("${meta.sample_name}_R1.fastq"), path("${meta.sample_name}_R2.fastq")

    script:
    def mhc_region = 'chr6:28510120-33480577'

    """
    set -euo pipefail

    # Build BED of non-primary contigs
    samtools idxstats $bam \
      | awk '
          \$1 != "*" &&
          \$1 != "chr6" &&
          \$1 !~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ &&
          \$1 != "MT"
        ' \
      | awk '{print \$1 "\\t0\\t" \$2}' \
      > nonprimary.bed

    # chr6 MHC interval
    samtools view \
        -b \
        -F 0x904 \
        ${bam} \
        ${mhc_region} \
        > chr6_mhc.bam

    # unmapped reads
    samtools view \
        -b \
        -f 4 \
        -F 0x904 \
        ${bam} \
        > unmapped.bam

    # reads on non-primary contigs
    samtools view \
        -b \
        -F 0x904 \
        -L nonprimary.bed \
        ${bam} \
        > nonprimary.bam

    # merge selected reads
    samtools merge \
        -f \
        merged.bam \
        chr6_mhc.bam \
        unmapped.bam \
        nonprimary.bam

    # BAM -> paired FASTQ
    samtools collate -Ou merged.bam \
      | samtools fastq \
            -1 ${meta.sample_name}_R1.fastq \
            -2 ${meta.sample_name}_R2.fastq \
			-0 /dev/null \
			-s /dev/null \
            -
    """
}



process HLAHD_TO_TSV {

    cpus 2
    memory "4GB"

    conda "python=3.10 pandas=2.1"

    tag "Converting HLAHD to TSV for  ${meta.sample_name}"

    input:
        tuple val(meta), path(hlahd_result)

    output:
		tuple val(meta), path("${meta.sample_name}_hlahd.tsv"), emit: hlahd_tsv

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    
    test = pd.read_csv("${hlahd_result}", sep = "\\t", header = None)
    test = test.rename(columns = {0:"locus"})
    test.insert(0, "sample", "${meta.sample_name}")
    test["calls"] = test.iloc[:, 2:].apply(
        lambda r: " - ".join([x for x in r if pd.notna(x) and x != "-" and x != "Not typed"]),
            axis=1
            )

    wide = (
        test.pivot(index="sample", columns="locus", values="calls")
              .reset_index()
              )
    wide.to_csv("${meta.sample_name}_hlahd.tsv", sep = "\t", index = False)
    """
}




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



process OPTITYPE {

    /*
        HLA type patient's fastqs using OPTITYPE
        This process can be fairly memory intensive. 
        One potential improvement here would be to map BAMs first, then extract reads from HLA region and run.

    */
    
    container "fred2/optitype:latest"
    cpus 8
    memory "32GB"

    tag "Optitype calls for ${meta.sample_name}"

    input:
        tuple val(meta), path(fastq1), path(fastq2)

    output:
        tuple val(meta), path("optitype_out/*_result.tsv"), path("optitype_out/*_coverage_plot.pdf"), emit: hla_calls

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
    
    cpus 8
    memory "32GB"
    container "griffithlab/hlahd:1.0"

    tag "HLA-HD on ${meta.sample_name}"
    
    input:
        tuple val(meta), path(fastq1), path(fastq2)
    
    output:
        tuple val(meta), path("./${meta.sample_name}/result/${meta.sample_name}_final.result.txt"), emit: final_hla_calls
        tuple val(meta), path("./${meta.sample_name}/result/"), emit: result_dir

    script:
        """
        /opt/hlahd/bin/hlahd.1.6.1.sh \
            -f /opt/hlahd/freq_data \
            -t $task.cpus \
            $fastq1 $fastq2 \
            /opt/hlahd/HLA_gene.split.txt \
            /opt/hlahd/dictionary \
            "${meta.sample_name}" \
            ./
        """

}


process EXTRACT_MHC_REGION {
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam), path(bai)
    
    output:
        tuple val(meta), path("${meta.sample_name}_hla_regions.bam") 

    script:
        """
        samtools view --threads $task.cpus -h -b -f 4 $bam > unmapped.bam
        samtools view --threads $task.cpus -h -b $bam chr6:28510120-33480577 > mhc.bam
        samtools merge --threads $task.cpus -o hla_regions.bam mhc.bam unmapped.bam
        samtools collate --threads $task.cpus -o "${meta.sample_name}_hla_regions.bam" hla_regions.bam
        """

}


process BAM_TO_FASTQ {

    /*

    */
    
    cpus 4
    memory "16GB"
    conda "bioconda::samtools=1.23.1 bioconda::htslib=1.23.1"

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path("${meta.sample_name}_${meta.molecule}_R1.fastq"), path("${meta.sample_name}_${meta.molecule}_R2.fastq")


    script:
        """
        samtools fastq --threads $task.cpus -1 ${meta.sample_name}_${meta.molecule}_R1.fastq -2 ${meta.sample_name}_${meta.molecule}_R2.fastq -n $bam
        """

}


