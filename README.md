# vaximile

<img src="docs/images/vaximile_logo.png" alt="vaximile" align="right" width="128" height="128">

**Tumour neoantigen discovery from paired tumour/normal bulk DNA and tumour RNA sequencing.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.0-23aa62.svg)](https://www.nextflow.io/)

vaximile takes a patient's tumour and matched normal exome or genome libraries plus tumour
RNA-seq and produces ranked neoantigen predictions, along with the evidence behind them:
somatic and germline variants, HLA genotypes and HLA loss of heterozygosity, allele-specific
copy number, expression, and fusions.

## What it runs

**DNA** — `fastp`, [minibwa](https://github.com/lh3/minibwa), duplicate marking and GATK
base recalibration, then `samtools` QC and [somalier](https://github.com/brentp/somalier)
relatedness checks to catch sample swaps.

**Somatic variants** — [Mutect2](https://gatk.broadinstitute.org),
[Strelka2](https://github.com/Illumina/strelka) with
[Manta](https://github.com/Illumina/manta), and
[DeepSomatic](https://github.com/google/deepsomatic). Each callset is filtered to `PASS`
and normalised, then combined on an **n−1 consensus**: a variant is kept when at least two
of the three callers report it.

**Germline variants** — GATK `HaplotypeCaller` with CNN scoring,
[Strelka2](https://github.com/Illumina/strelka) in germline mode, and
[DeepVariant](https://github.com/google/deepvariant), combined on the same 2-of-3 rule. The
consensus callset also drives proximal variant phasing for pVACseq.

**HLA** — [OptiType](https://github.com/FRED-2/OptiType) (class I) and
[HLA-HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) (class I and II) per library, and
[mhcflow](https://github.com/svm-zhang/mhcflow) typing that feeds
[lohhla-mod](https://github.com/svm-zhang/lohhla-mod) for HLA loss of heterozygosity, each
tumour measured against its own normal.

**Copy number** — [ASCAT](https://github.com/VanLoo-lab/ascat), which also supplies the
purity and ploidy the LOH analysis uses.

**RNA** — [STAR](https://github.com/alexdobin/STAR) and
[salmon](https://combine-lab.github.io/salmon/) quantification, with
[Arriba](https://github.com/suhrig/arriba) and
[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) for fusions.

**Neoantigens** — [Ensembl VEP](https://www.ensembl.org/vep) annotation with DNA and RNA
coverage and expression, then [pVACseq](https://pvactools.readthedocs.io) on somatic
variants and pVACfuse on fusions.

**Report** — one [MultiQC](https://multiqc.info/) report per patient, including the HLA LOH
and ASCAT results.

## Quick start

Reference data downloads itself on the first run - the VEP cache, ASCAT and Arriba
resources, the CTAT bundle - into `./vaximile_resources/`, and later runs reuse it from
there. The VEP cache is the big one at ~24 GiB, so pass `--vep_cache` to point at one you
already have. Indices are built from the reference unless you pass `--bwa_index`,
`--star_index` or `--salmon_index`. See [docs/usage.md](docs/usage.md#prerequisites).

One row per library. A normal shared by two tumours is listed once per pair; RNA rows leave
`capture_kit` empty:

```csv title="samplesheet.csv"
patient,somatic_name,sample_name,sample_type,sequencing_type,sex,capture_kit,fastqr1,fastqr2
PatientX,PatientX_T1_N1,Tumor1,Tumor,exome,XX,twist_2,t1_r1.fq.gz,t1_r2.fq.gz
PatientX,PatientX_T1_N1,Normal1,Normal,exome,XX,twist_2,n1_r1.fq.gz,n1_r2.fq.gz
PatientX,PatientX_T1_N1,Tumor1,Tumor,rna,XX,,t1_rna_r1.fq.gz,t1_rna_r2.fq.gz
```

```csv title="capture_kits.csv"
kit,bed
twist_2,/beds/TwistExome_GRCh38_chr.bed
```

```bash
nextflow run . -profile conda --samplesheet ./samplesheet.csv --capture_kits ./capture_kits.csv --outdir ./vaximile_out -resume
```

Some modules are container-only, so a real run wants a container engine enabled alongside
conda. `conf/ucsf_krummellab.config` does that for the UCSF cluster:

```bash
nextflow run . -c conf/ucsf_krummellab.config -params-file params.json -resume
```

Before a real run, `nextflow run . -profile test -stub-run` exercises the whole DAG in
seconds with no data, no containers and no conda.

## Layout

```text
main.nf                  entry workflow, publishing
nextflow.config          parameters, profiles
nextflow_schema.json     parameter validation and --help
conf/                    resources, per-module arguments, cluster and test profiles
workflows/vaximile/      the pipeline body
subworkflows/local/      grouped stages
modules/local/           one process each
assets/                  input schemas, tiling script, test fixtures
docs/                    usage and output
```

## Documentation

[usage](docs/usage.md) · [output](docs/output.md) · [citations](CITATIONS.md)

## Credits

vaximile is developed by Harrison Wismer at UCSF.
