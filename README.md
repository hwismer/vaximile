# vaximile

<img src="docs/images/vaximile_logo.png" alt="vaximile" align="right" width="128" height="128">

**Tumor neoantigen discovery using paired tumour-normal bulk DNA and RNA sequencing.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

vaximile takes a patient's tumour and matched normal exome or genome libraries plus tumour
RNA-seq, and produces ranked neoantigen predictions along with the evidence behind them:
somatic and germline variant calls, HLA class I and II genotypes, allele-specific copy
number, RNA expression, and RNA fusion calls.

## Pipeline summary

1. **Read QC** — [`fastp`](https://github.com/OpenGene/fastp) on all DNA and RNA libraries
2. **DNA alignment and pre-processing** — [`minibwa`](https://github.com/lh3/minibwa),
   then GATK best-practices duplicate marking and base recalibration
3. **Alignment QC** — `samtools flagstat`/`coverage`/`idxstats`, and
   [`somalier`](https://github.com/brentp/somalier) relatedness checks to catch sample swaps
4. **HLA typing** — [`OptiType`](https://github.com/FRED-2/OptiType) (class I) and
   [`HLA-HD`](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) (class I and II), run per library
   and on a merged per-pair BAM
5. **Somatic variant calling** — three callers on the tumour/normal pair:
   [`Mutect2`](https://gatk.broadinstitute.org),
   [`Strelka2`](https://github.com/Illumina/strelka) with
   [`Manta`](https://github.com/Illumina/manta), and
   [`DeepSomatic`](https://github.com/google/deepsomatic). Each callset is filtered to
   `PASS`, normalised, then combined into an **n−1 consensus**: a variant is retained if at
   least two of the three callers report it
6. **Germline variant calling** — three callers on the normal sample: GATK
   `HaplotypeCaller` with CNN variant scoring,
   [`Strelka2`](https://github.com/Illumina/strelka) in germline mode, and
   [`DeepVariant`](https://github.com/google/deepvariant). Combined into an **n−1
   consensus** on the same 2-of-3 rule. The consensus germline callset is also what drives
   proximal variant phasing for pVACseq
7. **Copy number** — [`ASCAT`](https://github.com/VanLoo-lab/ascat) (vendored nf-core module)
8. **RNA processing** — [`STAR`](https://github.com/alexdobin/STAR),
   [`kallisto`](https://pachterlab.github.io/kallisto/) and
   [`salmon`](https://combine-lab.github.io/salmon/) quantification
9. **Fusion calling** — [`Arriba`](https://github.com/suhrig/arriba) and
   [`STAR-Fusion`](https://github.com/STAR-Fusion/STAR-Fusion)
10. **Annotation** — [`Ensembl VEP`](https://www.ensembl.org/vep) plus DNA/RNA coverage and
    expression annotation of the somatic VCF
11. **Neoantigen prediction** — [`pVACseq`](https://pvactools.readthedocs.io) on somatic
    variants and `pVACfuse` on fusions
12. **Reporting** — per-patient [`MultiQC`](https://multiqc.info/) report

## Usage

> [!NOTE]
> A local Ensembl VEP cache is required. See [docs/usage.md](docs/usage.md#prerequisites).

Prepare a samplesheet describing each library:

```csv title="samplesheet.csv"
patient,somatic_name,sample_name,sample_type,sequencing_type,sex,capture_kit,fastqr1,fastqr2
PatientX,PatientX_Tumor1_Normal1,Tumor1,Tumor,exome,XX,twist_2,t1_r1.fq.gz,t1_r2.fq.gz
PatientX,PatientX_Tumor1_Normal1,Normal1,Normal,exome,XX,twist_2,n1_r1.fq.gz,n1_r2.fq.gz
PatientX,PatientX_Tumor1_Normal1,Tumor1,Tumor,rna,XX,,t1_rna_r1.fq.gz,t1_rna_r2.fq.gz
```

and a sheet mapping capture kits to their targets:

```csv title="capture_kits.csv"
kit,bed
twist_2,/beds/TwistExome_GRCh38_chr.bed
```

Then run:

```bash
nextflow run . -profile conda --samplesheet ./samplesheet.csv --capture_kits ./capture_kits.csv --vep_cache ./vep/vep_data/ --outdir ./vaximile_out -resume
```

Full documentation: [usage](docs/usage.md) · [output](docs/output.md)

## Repository layout

This pipeline follows the nf-core directory layout:

```text
main.nf                              entry workflow, publish and output blocks
nextflow.config                      parameter defaults, profiles, manifest
nextflow_schema.json                 parameter documentation and validation
conf/
  base.config                        resource defaults and process_* labels
  modules.config                     per-module tool arguments
  test.config, test_full.config      test profiles
  ucsf_krummellab.config             UCSF SLURM/Apptainer settings
workflows/vaximile/main.nf           the VAXIMILE workflow
subworkflows/local/<name>/main.nf    17 subworkflows
modules/local/<name>/main.nf         94 single-process modules
modules/nf-core/                     vendored nf-core modules
assets/                              input schemas, example inputs, test fixtures
docs/                                usage, output, migration notes
```

The reorganisation into this layout kept every process body unchanged. What that leaves
outstanding for `nf-core lint` is recorded in
[docs/nf-core-migration.md](docs/nf-core-migration.md).

## Credits

vaximile is developed by Harrison Wismer at UCSF.

## Citations

Tool references are collected in [CITATIONS.md](CITATIONS.md).
