# vaximile

<img src="assets/images/vaximile_logo.png" alt="vaximile" align="right" width="128" height="128">

**Tumour neoantigen discovery from paired tumour/normal bulk DNA and tumour RNA sequencing.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.0-23aa62.svg)](https://www.nextflow.io/)

vaximile takes a patient's tumour and matched normal exome or genome libraries plus tumour
RNA-seq and produces ranked neoantigen predictions, along with the evidence behind them:
somatic and germline variants, HLA genotypes and HLA loss of heterozygosity, allele-specific
copy number, expression, and fusions.

## Steps

**DNA preprocessing** — [fastp](https://github.com/OpenGene/fastp) trimming,
[minibwa](https://github.com/lh3/minibwa) alignment, duplicate marking and GATK base
recalibration, then `samtools` QC and [somalier](https://github.com/brentp/somalier)
relatedness checks to catch sample swaps.

**Somatic variant calling** — [Mutect2](https://gatk.broadinstitute.org),
[Strelka2](https://github.com/Illumina/strelka) with
[Manta](https://github.com/Illumina/manta), and
[DeepSomatic](https://github.com/google/deepsomatic). Each callset is filtered to `PASS`
and normalised, then combined on an **n−1 consensus**: a variant is kept when at least two
of the three callers report it.

**Germline variant calling** — GATK `HaplotypeCaller` with CNN scoring,
[Strelka2](https://github.com/Illumina/strelka) in germline mode, and
[DeepVariant](https://github.com/google/deepvariant), combined on the same 2-of-3 rule. The
consensus callset also drives proximal variant phasing for pVACseq.

**HLA typing** — [OptiType](https://github.com/FRED-2/OptiType) (class I) and
[HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) (class I and II) per library, whose
calls go to pVACtools.

**CNV** — [ASCAT](https://github.com/VanLoo-lab/ascat) allele-specific copy number, which
also supplies the tumour purity and ploidy the LOH step uses.

**HLA LOH** — [mhcflow](https://github.com/svm-zhang/mhcflow) types each normal against a
sample-specific HLA reference, each tumour is realigned against the reference inferred for *its
own normal*, and [lohhla-mod](https://github.com/svm-zhang/lohhla-mod) calls loss over the
pair. mhcflow needs Novoalign, which you supply with `--novoalign_dir`; without it this step
is skipped.

**RNA alignment and quantification** — [STAR](https://github.com/alexdobin/STAR) alignment
and [salmon](https://combine-lab.github.io/salmon/) quantification at transcript and gene
level.

**RNA fusion calling** — [Arriba](https://github.com/suhrig/arriba) and
[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion).

**Neoantigen prediction** — [Ensembl VEP](https://www.ensembl.org/vep) annotation with DNA
and RNA coverage and expression, then [pVACseq](https://pvactools.readthedocs.io) on
somatic variants and pVACfuse on fusions.

**Report** — one [MultiQC](https://multiqc.info/) report per patient, including the HLA LOH
and ASCAT results, somatic and germline variant QC (VAFs, caller overlap, filtering funnel,
substitution spectrum, TMB), pVACseq tier summaries, capture coverage from
[mosdepth](https://github.com/brentp/mosdepth), and library damage (8-oxoG oxidation,
deamination) from Picard.

## Parameters

### Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--samplesheet` | **Required.** Path to comma-separated file describing the patients, samples and FastQ files. | — |
| `--capture_kits` | **Required.** Path to comma-separated file mapping capture kit names to their BED target files. | — |
| `--outdir` | **Required.** The output directory where the results will be saved. | `./vaximile_out` |

### Reference genome options

Reference genome, annotation and transcriptome files.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--reference_fa` | Reference genome FASTA. Gzipped input is decompressed and indexed automatically. | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz` |
| `--reference_includes_chr_prefix` | Set to false when the reference contigs are named without a 'chr' prefix. | `true` |
| `--gtf` | Gene annotation GTF used for RNA alignment, quantification and fusion calling. | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz` |
| `--transcriptome_reference` | Transcriptome FASTA used to build the salmon index. | `https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz` |
| `--human_ref_peptides` | Reference proteome peptide FASTA used by pVACtools for reference proteome similarity. | `https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz` |

### Alignment index options

Prebuilt indices. Any left unset are built from the reference during the run.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--bwa_index` | Directory containing a prebuilt minibwa index. When not given, the index is built and saved to `vaximile_resources/bwa_index/` in the launch directory. | — |
| `--star_index` | Directory containing a prebuilt STAR index. When not given, the index is built and saved to `vaximile_resources/star_index/STARGenomeDir/` in the launch directory. | — |
| `--salmon_index` | Directory containing a prebuilt salmon index. When not given, the index is built and saved under `vaximile_resources/salmon_index/` in the launch directory. | — |

### Downloaded resource options

Resources the pipeline downloads when not given, once, into `./vaximile_resources/` in the launch directory, reusing them on later runs. Each default is the download it makes; giving a parameter skips that download. See [Saved resources](#saved-resources) for where each one is saved.

| Parameter | Description | Default |
| --------- | ----------- | ------- |
| `--vep_cache` | Path to a local Ensembl VEP cache directory. Optional; downloaded if omitted. | Downloads `https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz` |
| `--vep_plugins` | Directory containing the pVACtools Frameshift.pm and Wildtype.pm VEP plugins. Downloaded if omitted. | Downloads both plugins from `https://github.com/griffithlab/pVACtools/tree/v7.1.5/pvactools/tools/pvacseq/VEP_plugins` |
| `--ctat_resource_lib` | CTAT genome library for STAR-Fusion: the `ctat_genome_lib_build_dir` directory of a plug-n-play bundle. Downloaded if omitted. | Downloads `https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz` |
| `--arriba_blacklist` | Arriba blacklist. Downloaded with the Arriba v2.5.1 release if omitted. | `blacklist_hg38_GRCh38_v2.5.1.tsv.gz` from `https://github.com/suhrig/arriba/releases/download/v2.5.1/arriba_v2.5.1.tar.gz` |
| `--arriba_known_fusions` | Arriba known fusions. Downloaded with the Arriba v2.5.1 release if omitted. | `known_fusions_hg38_GRCh38_v2.5.1.tsv.gz` from the same release |
| `--arriba_protein_domains` | Arriba protein domains (GFF3). Downloaded with the Arriba v2.5.1 release if omitted. | `protein_domains_hg38_GRCh38_v2.5.1.gff3` from the same release |
| `--ascat_resources` | Directory of ASCAT 1000 Genomes resources: `G1000_loci_hg38/`, `G1000_alleles_hg38/`, `GC_G1000_hg38/` and `RT_G1000_hg38/`, as the pipeline downloads them. Downloaded if omitted. | Downloads the four 1000G WGS hg38 archives from `https://zenodo.org/records/14008443`, adding a `chr` prefix to the loci when `--reference_includes_chr_prefix` is set |

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
| `--hla_fasta` | HLA class I reference FASTA. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/9326ac2ef9b112b9bc5f25222854296097b58be5/data/abc_complete.fasta` |
| `--hla_fasta_fai` | Index for the HLA reference FASTA. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/9326ac2ef9b112b9bc5f25222854296097b58be5/data/abc_complete.fasta.fai` |
| `--hla_kmers` | Unique HLA k-mer file used by the typing tools. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/9326ac2ef9b112b9bc5f25222854296097b58be5/data/abc_v14.uniq` |
| `--hla_freq` | HLA allele population frequency table. | `https://raw.githubusercontent.com/jason-weirather/hla-polysolver/9326ac2ef9b112b9bc5f25222854296097b58be5/data/HLA_FREQ.txt` |
| `--novoalign_dir` | Directory containing the `novoalign` and `novoindex` binaries, needed for HLA LOH. **Without it, HLA LOH is skipped**; HLA typing with OptiType and HLA-HD still runs. Novoalign is proprietary and no longer on bioconda, so download it from [Novocraft](https://www.novocraft.com/support/download/). Novoalign V4 also needs its `novoalign.lic` in the same directory; V3 runs without one. | — |

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

### Samplesheet 
One sample per line. A complete somatic grouping consists of a tumor sample, normal sample, and a tumor RNA sample.

- **`patient`** — identifies a shared patient field that can contain multiple somatic groups. The pipeline concludes with a patient-level MultiQc reports with all patient-specific sample metrics.

- **`somatic_name`** — identifies the somatic comparison being done. Each somatic_name should have a single tumor dna sample, normal dna sample, and tumor RNAseq. In the case that multiple tumor samples share one normal sample, the normal sample will still need its own line in the samplesheet.

- **`sample_name`** — name of the individual sample.

- **`sample_type`** — either Tumor or Normal.

- **`sequencing_type`** — one of: exome, genome, rna, exome_FFPE, genome_FFPE.

- **`sex`** — either XX, XY, or left blank. Used only when running ASCAT.

- **`capture_kit`** — name of the capture kit used. Must match a kit name in the capture_kits.csv. Left blank by RNA samples.

- **`fastqr{1.2}`** — R1 and R2 fastqs for the sample

#### Example:

`samplesheet.csv`

| patient  | somatic_name   | sample_name | sample_type | sequencing_type | sex | capture_kit | fastqr1         | fastqr2         |
| -------- | -------------- | ----------- | ----------- | --------------- | --- | ----------- | --------------- | --------------- |
| PatientX | PatientX_T1_N1 | Tumor1      | Tumor       | exome           | XX  | twist_2     | t1_r1.fq.gz     | t1_r2.fq.gz     |
| PatientX | PatientX_T1_N1 | Normal1     | Normal      | exome           | XX  | twist_2     | n1_r1.fq.gz     | n1_r2.fq.gz     |
| PatientX | PatientX_T1_N1 | Tumor1      | Tumor       | rna             | XX  |             | t1_rna_r1.fq.gz | t1_rna_r2.fq.gz |

`capture_kits.csv`

| kit     | bed                             |
| ------- | ------------------------------- |
| twist_2 | /beds/TwistExome_GRCh38_chr.bed |

```bash
nextflow run . -profile singularity,conda --samplesheet ./samplesheet.csv --capture_kits ./capture_kits.csv --outdir ./vaximile_out -resume
```
Some modules are container-only while others are conda-only. To execute a full run without problems, both will need to be enabled.


## Output

Publishing is driven by the `output {}` block in `vaximile.nf`. Results for each patient sit
under `<outdir>/<patient>/`, split into tumor/normal pairs and single libraries.

### Directory structure

```text
<outdir>/
└── <patient>/
    ├── <patient>_report.html          MultiQC report
    ├── <patient>_pvacseq_reports.tsv  combined MHC class I report across the patient's pairs
    ├── somatic_samples/               one directory per tumor/normal pair
    │   └── <somatic_name>/
    │       ├── somatic/               somatic VCF (+ index) and its tabular form
    │       ├── neoantigen_prediction/ pVACseq and pVACfuse output directories
    │       ├── hla_loh/               lohhlamod results and per-gene plots
    │       └── CNV/ascat/             ASCAT results
    └── samples/                       one directory per library
        └── <sample_name>/
            ├── alignment/             duplicate-marked BAM, recalibrated BAM, STAR BAM (+ indexes)
            ├── germline/              germline VCF (+ index) and its tabular form, normals only
            ├── hla_typing/
            │   ├── optitype/          class I calls and coverage plot
            │   ├── hlahd/             class I and II calls
            │   └── <sample_name>_hla_calls.csv  combined calls, formatted for pVACtools
            └── salmon/                gene-level abundance, RNA libraries only
```

A library shared between pairs - a normal used as the control for several tumors - is
processed once, so its alignments, germline calls and HLA typing appear once under
`samples/` rather than beneath each pair.

### Saved resources

Downloaded resources and the alignment indexes the pipeline builds are saved in the launch
directory, not the output directory, so later runs from the same directory reuse them. To
use them from another directory, give the parameter shown.

```text
vaximile_resources/
├── vep_cache/vep_cache/                --vep_cache
├── vep_pvac_plugins/VEP_plugins/       --vep_plugins
├── ctat_resource_dir/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/
│   └── ctat_genome_lib_build_dir/      --ctat_resource_lib
├── arriba_hg38/arriba_v2.5.1/database/ --arriba_blacklist, --arriba_known_fusions, --arriba_protein_domains
├── ascat_hg38_chr/                     --ascat_resources (ascat_hg38_nochr without the chr prefix)
├── bwa_index/                          --bwa_index
├── star_index/STARGenomeDir/           --star_index
└── salmon_index/<name>_salmon_index/   --salmon_index
```

### Results by stage

#### Quality control

`<patient>_report.html` is the MultiQC report for the patient. It aggregates fastp reports
for every DNA and RNA library, samtools and mosdepth coverage, GATK base recalibration
tables, STAR and salmon summaries, VEP summaries, and somalier relatedness. Its custom
sections cover the somatic calls (filtering metrics, substitution types, TMB, contamination,
and VAF plots), the germline calls (caller overlap, het/hom, Ti/Tv, and VAF), HLA typing and
LOH, ASCAT copy number, the pVACseq tiers, and sequencing artifacts from Picard
CollectOxoGMetrics and CollectSequencingArtifactMetrics. It ends with how the run was
launched and every parameter it used, including where downloaded or built resources went.

The somalier tables in that report are the check for sample swaps - confirm the tumor and
normal of each pair relate to each other before trusting downstream calls.

#### Variant calling

`somatic/` holds the merged, filtered and annotated somatic callset. Three callers
contribute: Mutect2, Strelka2 (with Manta for indel candidates) and DeepSomatic. Each is
filtered to `PASS` and normalised, then combined with GATK3 `CombineVariants` under
`--minimumN 2` - an **n-1 consensus**, so a variant is kept when at least two of the three
callers report it - and filtered to variants below 0.1% frequency in every population. The
`.tsv` alongside it is the same content flattened for spreadsheet use.

`germline/` holds the consensus germline callset for the normal sample, built the same way
from three callers: GATK `HaplotypeCaller` (CNN-scored and tranche-filtered), Strelka2 in
germline mode, and DeepVariant, combined under the same 2-of-3 rule. This VCF is also the
phasing input for proximal variant detection in pVACseq.

#### HLA typing

`hla_typing/` holds class I alleles from OptiType and class I and II alleles from HLA-HD,
both run per library, and the two combined into one allele list for pVACtools. The
normal's list is what feeds pVACtools, since HLA type is germline.

#### HLA loss of heterozygosity

`hla_loh/` holds the lohhlamod result table for the pair - per-allele copy number, the
mismatch logR p-value and BAF - alongside per-gene coverage, logR and BAF plots. The tumor
is realigned against the HLA reference mhcflow inferred for its own normal, so the two BAMs
sit on one subject-specific reference.

#### Copy number

`CNV/ascat/` publishes the whole ASCAT result set for a pair - segments, CNV calls,
purity/ploidy estimates, BAF and LogR tables, run metrics and plots. Its purity and ploidy
are what the HLA LOH analysis uses.

#### Neoantigen prediction

`neoantigen_prediction/` contains the pVACseq output for somatic variants and the pVACfuse
output for RNA fusions called by Arriba and STAR-Fusion. `<patient>_pvacseq_reports.tsv` is
the MHC class I aggregated report combined across the patient's tumor/normal pairs, with the
allelic depth of each sample from the pair's somatic VCF.

#### Expression

`salmon/` publishes `quant.genes.sf`, the gene-level TPM table salmon aggregates from
`--geneMap`. The transcript-level `quant.sf` and the rest of the salmon run directory are
used internally - for VCF expression annotation, strandedness prediction and MultiQC -
rather than published.

#### Alignments

`alignment/` holds three BAMs per library, each with its index, distinguished by suffix:
`_markdup.bam` after duplicate marking, `_bqsr.bam` after base recalibration, and
`_STAR_sorted.bam` for RNA. Mutect2 and HaplotypeCaller read the recalibrated BAMs;
Strelka, Manta, DeepSomatic and DeepVariant read the duplicate-marked ones.

### Not published by default

Per-caller intermediate VCFs, RNA fusion tables and the transcript-level salmon output
remain in the work directory. Add them to the `publish:` block in `vaximile.nf` and give
them a `path {}` in the `output {}` block to change that.

## Credits

vaximile is developed by Harrison Wismer at UCSF.

## Citations

Tool references are collected in [CITATIONS.md](CITATIONS.md).
