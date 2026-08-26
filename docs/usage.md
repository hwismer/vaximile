# vaximile: Usage

## Introduction

vaximile predicts tumour neoantigens from paired tumour/normal bulk DNA sequencing and
matched bulk RNA sequencing. A run produces somatic and germline variant calls, HLA class I
and II genotypes, RNA fusion calls, and pVACseq/pVACfuse neoantigen predictions, plus a
per-patient MultiQC report.

## Prerequisites

### VEP cache (required)

Download the Ensembl VEP cache matching your reference genome and pass the directory with
`--vep_cache`. There is no default; the pipeline fails at startup if the directory is
missing.

<https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/>

### Prebuilt indices (optional)

Any index left unset is built from `--reference_fa` / `--transcriptome_reference` during the
run, which adds substantial wall time to a first run.

| Parameter          | Tool                  |
| ------------------ | --------------------- |
| `--bwa_index`      | bwa-mem2              |
| `--star_index`     | STAR 2.7.10           |
| `--kallisto_index` | kallisto              |
| `--salmon_index`   | salmon                |

## Samplesheet input

`--samplesheet` takes a comma-separated file with the header shown below. One row per
FastQ pair. DNA and RNA rows for the same tumour go in the same sheet.

```csv title="samplesheet.csv"
patient,somatic_name,sample_name,sample_type,sequencing_type,sex,capture_kit,fastqr1,fastqr2
PatientX,PatientX_Tumor1_Normal1,Tumor1,Tumor,exome,XX,twist_2,t1_r1.fq.gz,t1_r2.fq.gz
PatientX,PatientX_Tumor1_Normal1,Normal1,Normal,exome,XX,twist_2,n1_r1.fq.gz,n1_r2.fq.gz
PatientX,PatientX_Tumor1_Normal1,Tumor1,Tumor,rna,XX,,t1_rna_r1.fq.gz,t1_rna_r2.fq.gz
```

| Column            | Description                                                                                                                                  |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `patient`         | Patient identifier. Groups samples for the MultiQC report, somalier relatedness check and the combined pVACseq report.                        |
| `somatic_name`    | Identifies the tumour/normal pair. Rows sharing a value are called together, so this is what pairs a tumour with its normal and its RNA.      |
| `sample_name`     | Name of the individual library.                                                                                                              |
| `sample_type`     | `Tumor` or `Normal`.                                                                                                                         |
| `sequencing_type` | `exome`, `genome`, `exome_FFPE`, `genome_FFPE` or `rna`. Anything else leaves `meta.molecule` null and the row is dropped at the DNA/RNA branch. |
| `sex`             | `XX` or `XY`. Passed to ASCAT for copy number calling.                                                                                       |
| `capture_kit`     | Must match a `kit` value in the capture kit sheet. Leave blank for RNA rows.                                                                  |
| `fastqr1`         | Gzipped FastQ, read 1. Must end `.fq.gz` or `.fastq.gz`.                                                                                      |
| `fastqr2`         | Gzipped FastQ, read 2. Must end `.fq.gz` or `.fastq.gz`.                                                                                      |

## Capture kit input

`--capture_kits` maps kit names to their BED target files. Intervals are scattered per kit,
so every kit named in the samplesheet needs a row here.

```csv title="capture_kits.csv"
kit,bed
agilent_v7,/beds/AGV7_GRCh38_chr.bed
twist_2,/beds/TwistExome_GRCh38_chr.bed
```

## Running the pipeline

```bash
nextflow run . \
    -profile conda \
    --samplesheet ./samplesheet.csv \
    --capture_kits ./capture_kits.csv \
    --vep_cache ./vep/vep_data/ \
    --outdir ./vaximile_out \
    -resume
```

Parameters can be supplied as a file instead, which is easier to version:

```bash
nextflow run . -profile conda -params-file assets/params_example.json -resume
```

### Reference contig naming

`--reference_includes_chr_prefix` must match your reference. It selects the correct ASCAT
resource files and the MHC region coordinates used for HLA typing. Setting it wrongly
produces empty extractions rather than an error.

## Profiles

Use `-profile` to select a software provisioning method. Multiple profiles are
comma-separated and later entries override earlier ones.

- `conda` — the only fully supported option today. Every local module declares a conda spec.
- `docker`, `singularity`, `apptainer` — the vendored nf-core `ascat` module carries a
  container, but the local modules do not, so these profiles do not yet run the whole
  pipeline on their own.
- `test` — minimal settings for a smoke test; see `conf/test.config`.

### Institutional configuration

`conf/ucsf_krummellab.config` holds the SLURM/Apptainer settings for the UCSF cluster:

```bash
nextflow run . -c conf/ucsf_krummellab.config -params-file params.json -resume
```

Use `-c` (merge), not the `-C` (replace) this config was previously used with. Parameter
defaults now live in `nextflow.config`, which `-C` would skip entirely.

## Resource requests

Most local modules set `cpus` and `memory` inside the process body. `conf/base.config`
provides defaults and the `process_single`/`process_low`/`process_medium`/`process_high`
labels used by vendored nf-core modules. To retune a module without editing it, add a
`withName:` block to a custom config passed with `-c`.

## Reproducibility

Pin a release with `-r` so the same code runs each time:

```bash
nextflow run hwismer/vaximile -r 1.0.0 -profile conda --samplesheet ...
```
