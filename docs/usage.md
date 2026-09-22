# vaximile: Usage

vaximile predicts tumour neoantigens from paired tumour/normal bulk DNA and matched tumour
RNA sequencing. See the [README](../README.md) for what the pipeline runs.

## Samplesheet

`--samplesheet` takes a CSV with one row per FastQ pair. DNA and RNA rows for the same
tumour go in the same sheet.

```csv title="samplesheet.csv"
patient,somatic_name,sample_name,sample_type,sequencing_type,sex,capture_kit,fastqr1,fastqr2
PatientX,PatientX_T1_N1,Tumor1,Tumor,exome,XX,twist_2,t1_r1.fq.gz,t1_r2.fq.gz
PatientX,PatientX_T1_N1,Normal1,Normal,exome,XX,twist_2,n1_r1.fq.gz,n1_r2.fq.gz
PatientX,PatientX_T1_N1,Tumor1,Tumor,rna,XX,,t1_rna_r1.fq.gz,t1_rna_r2.fq.gz
```

| Column | Description |
| ------ | ----------- |
| `patient` | Groups samples for the MultiQC report, the somalier relatedness check and the combined pVACseq report. |
| `somatic_name` | Identifies the tumour/normal pair. Rows sharing a value are called together, so this is what pairs a tumour with its normal and its RNA. |
| `sample_name` | Name of the individual library. |
| `sample_type` | `Tumor` or `Normal`. |
| `sequencing_type` | `exome`, `genome`, `exome_FFPE`, `genome_FFPE` or `rna`. Anything else drops the row at the DNA/RNA branch. |
| `sex` | `XX` or `XY`, or blank. Used only by ASCAT. |
| `capture_kit` | Must match a `kit` in the capture kit sheet. Blank for RNA rows. |
| `fastqr1`, `fastqr2` | Gzipped FastQs, ending `.fq.gz` or `.fastq.gz`. |

### A normal shared between pairs

List it once per pair, with a different `somatic_name` each time:

```csv
PatientX,PatientX_T1_N1,Normal1,Normal,exome,XX,twist_2,n1_r1.fq.gz,n1_r2.fq.gz
PatientX,PatientX_T2_N1,Normal1,Normal,exome,XX,twist_2,n1_r1.fq.gz,n1_r2.fq.gz
```

Rows sharing `patient` + `sample_name` + molecule are one library and are processed
**once** - one fastp, alignment, duplicate marking, BQSR, germline call, HLA typing and
somalier extraction. Only the per-pair steps run again per `somatic_name`.

## Capture kits

`--capture_kits` maps kit names to BED targets. Intervals are scattered per kit, so every
kit named in the samplesheet needs a row.

```csv title="capture_kits.csv"
kit,bed
agilent_v7,/beds/AGV7_GRCh38_chr.bed
twist_2,/beds/TwistExome_GRCh38_chr.bed
```

## Running Example

```bash
nextflow run . -profile singularity,conda \
    --samplesheet ./samplesheet.csv \
    --capture_kits ./capture_kits.csv \
    --outdir ./vaximile_out \
    -resume
```

## Parameters

### Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--samplesheet` | **Required.** Path to comma-separated file describing the patients, samples and FastQ files. | — |
| `--capture_kits` | **Required.** Path to comma-separated file mapping capture kit names to their BED target files. | — |
| `--vep_cache` | Path to a local Ensembl VEP cache directory. Optional; downloaded if omitted. | — |
| `--outdir` | **Required.** The output directory where the results will be saved. | `./vaximile_out` |

### Reference genome options

Reference genome, annotation and transcriptome files.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--reference_fa` | Reference genome FASTA. Gzipped input is decompressed and indexed automatically. | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz` |
| `--reference_includes_chr_prefix` | Set to false when the reference contigs are named without a 'chr' prefix. | `true` |
| `--gtf` | Gene annotation GTF used for RNA alignment, quantification and fusion calling. | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz` |
| `--transcriptome_reference` | Transcriptome FASTA used to build the salmon index. | `https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz` |
| `--human_ref_peptides` | Reference proteome peptide FASTA used by pVACtools for reference proteome similarity. | `https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz` |

### Alignment index options

Prebuilt indices. Any left unset are built from the reference during the run.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--bwa_index` | Directory containing a prebuilt minibwa index. | — |
| `--star_index` | Directory containing a prebuilt STAR 2.7.10 index. | — |
| `--salmon_index` | Directory containing a prebuilt salmon index. | — |

### Variant calling options

Interval scattering behaviour for the GATK-based callers.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--interval_padding` | Bases of padding added around each calling interval. | `100` |
| `--scatter_count` | Number of interval shards to scatter variant calling across. | `30` |

### HLA typing options

Reference data for MHC region extraction and HLA allele calling.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--hla_fasta` | HLA class I reference FASTA. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/master/data/abc_complete.fasta` |
| `--hla_fasta_fai` | Index for the HLA reference FASTA. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/master/data/abc_complete.fasta.fai` |
| `--hla_kmers` | Unique HLA k-mer file used by the typing tools. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/master/data/abc_v14.uniq` |
| `--hla_freq` | HLA allele population frequency table. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/master/data/HLA_FREQ.txt` |

### GATK resource bundle options

Known-sites VCFs used for BQSR, contamination estimation and somatic filtering. Each must have a matching .tbi.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--common_germline` | Common germline biallelic SNP sites for contamination estimation. | `gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz` |
| `--known_sites_dbsnp` | dbSNP known sites VCF for BQSR. | `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz` |
| `--known_sites_1000g_snps` | 1000 Genomes high-confidence SNPs for BQSR. | `gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz` |
| `--known_indels` | Known indels VCF for BQSR. | `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz` |
| `--gnomad` | gnomAD allele-frequency-only VCF for Mutect2. | `gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz` |
| `--pon` | Panel of normals VCF for Mutect2. | `gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz` |
| `--hapmap` | HapMap VCF for germline variant recalibration. | `gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz` |
| `--mills` | Mills gold standard indels VCF. | `gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` |
| `--somalier_sites` | Sites VCF used by somalier for relatedness checks. | `https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz` |

### Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--publish_dir_mode` | Method used to save pipeline results to output directory. | `copy` |
| `--version` | Display version and exit. | `false` |
| `--help` | Display help text. | `false` |

## Reproducibility

Pin a release so the same code runs each time:

```bash
nextflow run hwismer/vaximile -r 1.0.0 -profile conda --samplesheet ...
```
