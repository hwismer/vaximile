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

## Fast checks before a real run

```bash
nextflow run . -profile test -stub-run
```

Every module has a `stub:` block, so this executes the **entire DAG** - all 132 tasks -
in about 10 seconds, offline, with no containers, no conda and no data. Stubs only create
the files each process declares as output, so what it verifies is the wiring, not the
science:

- every output declaration actually resolves (a glob that matches nothing fails here)
- tuple arities line up between each process and its call sites
- **how many times each process runs**

That last one matters. A reference channel built as a queue instead of a value channel is
consumed by the first task, so an aligner silently processes one sample and skips the
rest - which no amount of linting, `nextflow inspect` or `-preview` will reveal, because
none of them execute tasks. Check counts against your samplesheet:

```bash
grep -oE 'Submitted process > [A-Za-z0-9_:]+' .nextflow.log | sed 's/.*://' | sort | uniq -c | sort -rn
```

The `test` profile points every reference at an empty placeholder under
`assets/test/refs/`, purely so the stub run stays offline - Nextflow stages inputs even
under `-stub-run`, and the real defaults are remote `https://` and `gs://` URIs.

ASCAT is skipped via `ext.when = false` in `conf/test.config`: it is the one vendored
nf-core module, and its stub still runs `Rscript -e "library(ASCAT)"` to capture a
version, which needs its container.

## Software versions

Every module writes a `versions.yml`, collected into
`<outdir>/pipeline_info/software_versions.yml`, so a set of results records the tool
versions that produced it. Modules whose tool has a usable version flag query it at run
time; the rest report the version pinned in their own `conda`/`container` directive.

Two things that file will show you, both worth acting on:

- `optitype: 1.3.1` - pinned to `fred2/optitype:release-v1.3.1`, that repo's newest tag.
  It was previously `:latest`, a floating tag last pushed in 2018.

  Upgrading to the current OptiType (1.5.0, on bioconda) is a **migration, not a repin**:
  1.5.0 is a CLI rewrite. The entry point became a click group (`optitype run` rather than
  `OptiTypePipeline.py`), and the `OptiType.ini` config this module writes is replaced by
  `--solver` / `--mapper` / `--razers3` / `--threads` flags. Output file naming has not been
  verified against the module's `optitype_out/*_result.tsv` glob either. Worth doing, but it
  needs a real HLA-typing run to validate.
- `unknown: unknown` - the four modules that declare neither `conda` nor `container`
  (plus two whose tool could not be identified) rely on whatever is on the host `PATH`.

## Profiles

Use `-profile` to select a software provisioning method. Multiple profiles are
comma-separated and later entries override earlier ones.

Software provisioning is **mixed**, and no single profile yet covers the whole pipeline.
Of the 94 local modules:

| Provisioning declared     | Modules |
| ------------------------- | ------- |
| Both `conda` and `container` | 39   |
| `conda` only              | 31      |
| `container` only          | 20      |
| Neither                   | 4       |

`-profile conda` now resolves software for 70 of 94 modules. The 20 container-only ones
still need a container engine, so a run today wants **both** conda and a container engine
enabled, which is what `conf/ucsf_krummellab.config` does.

Every conda spec was checked against bioconda with `conda search` and pins the same
version the container provides, with two exceptions noted below.

The four modules with neither — `combine_fastqs`, `prepare_fasta`,
`pull_arriba_resources`, `pull_ctat_resource_bundle` — rely on tools on the host `PATH`.

### Why the remaining 20 are still container-only

| Modules | Blocker |
| ------- | ------- |
| 4 × GATK3 (`merge_*_vcfs`, `phase_vcf_combine_variants`, `phase_vcf_rbphasing`) | bioconda's `gatk` 3.6 is a wrapper that needs the licensed jar registered manually; it cannot install unattended |
| 4 × VAtools | bioconda only has 6.0.1; the container pins 5.2.0, a major-version gap |
| 2 × STAR | bioconda only has 2.7.11b; the container pins 2.7.10a. STAR indices are version-sensitive, so switching would invalidate a prebuilt `--star_index` |
| 2 × Strelka, 1 × Manta | scripts call `configure*Workflow.py` by absolute container path, and Manta's bioconda floor (1.28) is far above the pinned 1.6.0 |
| `deepvariant` | bioconda has an exact 1.10.0, but the script hardcodes `/opt/deepvariant/bin/run_deepvariant`. Convertible with a one-line script change |
| `deepsomatic`, `hlahd` | not packaged in bioconda (HLA-HD is licence-restricted) |
| `optitype` | needs a `config.ini` and an ILP solver that the container supplies |
| `bamreadcount` | runs `bam_readcount_helper.py`, a script that exists only in the CWL image |
| `pull_vep_pvac_plugins` | bioconda has no pVACtools 6.0.3 |
| `phase_vcf_sort_vcf` | runs `java -jar /usr/picard/picard.jar`; bioconda `picard` 3.4.0 exists but needs the call rewritten to `picard SortVcf` |

### Reusing a prebuilt bwa-mem2 index

`--bwa_index` takes a directory. When it is set the pipeline skips `CREATE_BWA_INDEX`;
when it is not, the index is built once and published to `./resources/bwa/`.

The directory must contain the five files `bwa-mem2 index` produces, **named after the
prepared reference FASTA**:

```
<reference_fa stem>_prc.fa.0123
<reference_fa stem>_prc.fa.amb
<reference_fa stem>_prc.fa.ann
<reference_fa stem>_prc.fa.bwt.2bit.64
<reference_fa stem>_prc.fa.pac
```

The naming is not incidental. `BWA_MAP` passes the reference FASTA to bwa-mem2 as the
index prefix, so bwa-mem2 looks for `<reference_fa>.0123` and friends. An index built from
a FASTA with a different filename is unusable even if it is otherwise perfectly valid.
`PREPARE_FASTA` decompresses and renames the reference to `<stem>_prc.fa`, which is why
that suffix appears.

The pipeline validates this at launch and names any missing file, rather than letting
bwa-mem2 fail per-sample once alignment starts. The easiest way to get a valid directory
is to run once without `--bwa_index` and reuse `./resources/bwa/`.

Note for indices published before this was fixed: `CREATE_BWA_INDEX` used to publish with
Nextflow's default `publishDir` mode, which is **symlink**, so `./resources/bwa/` held
links into `work/`. Once `work/` was cleaned those links dangled - they still list in the
directory but no longer resolve, which is why validation could report a file as missing
while showing it in the same message. The module now copies. If you have such a directory,
rebuild the index by running once without `--bwa_index`.

### GATK versions are pinned deliberately — do not unify them

Three different GATK generations are in use, and the split is load-bearing. Each conda
spec pins exactly the version its container provided, and the modules carry comments
saying so.

| Modules | GATK | Why it cannot move |
| ------- | ---- | ------------------ |
| 19 modules (Mutect2, BQSR, pileups, interval/VCF utilities) | `gatk4=4.6.1.0` | current baseline |
| `haplotype_caller_scatter`, `haplotype_caller_cnn_score_variants`, `haplotype_caller_filter_variants` | `gatk4=4.3.0.0` | the CNN germline chain. `CNNScoreVariants` was deprecated in favour of `NVScoreVariants` and is not in current GATK4, so bumping these to 4.6.1.0 breaks the chain |
| `merge_germline_vcfs`, `merge_somatic_vcfs`, `phase_vcf_combine_variants`, `phase_vcf_rbphasing` | GATK3 3.6, container only | `CombineVariants` and `ReadBackedPhasing` were dropped in GATK4 and have no equivalent. bioconda's `gatk` 3.x is a wrapper needing the licensed jar registered by hand, so no conda spec is possible |

The practical rule: a GATK module's version is part of its behaviour. `HaplotypeCaller`
defaults, `FilterMutectCalls` filters and the CNN tranche models all differ between
releases, so treat any version change as a change to results and revalidate rather than
assuming it is a maintenance bump.

### Two conda pins that are not exact

- **`bwa-mem2=2.2.1`** — the container tag (`iarcbioinfo/bwa-mem2-tools:v1.0`) names the
  toolset, not bwa-mem2, so the bundled version is unstated. 2.2.1 is the long-standing
  stable release and is the most likely match, but it is inferred rather than verified.
  Confirm before trusting conda and container to produce identical alignments.
- **`ensembl-vep=115`** — the container is `release_115.0`; bioconda publishes this as
  `115` (plus patches `115.1`, `115.2`). Same VEP release, and it matches a release-115
  cache.

- `test` — minimal settings for a smoke test; see `conf/test.config`.

### Execution reports and the `ps` requirement

`timeline`, `report` and `trace` are **disabled** in `nextflow.config`, deliberately.

Enabling any of them makes Nextflow inject a guard into every task wrapper that runs
`command -v ps || exit 1`. It is a hard failure, not a warning — the task dies before the
tool runs, leaving an empty `.command.out` and a single line in `.command.err`:

```
Command 'ps' required by nextflow to collect task metrics cannot be found
```

Several images this pipeline uses have no `procps` (`alexdobin/star`,
`google/deepvariant`, `google/deepsomatic`, `staphb/bcftools`), and a container cannot see
the host's `ps`. `dag` does not inject the guard and stays enabled.

If you are running a configuration where every process does have `ps`, request the reports
per-run rather than re-enabling them globally:

```bash
nextflow run . -profile conda -with-report -with-timeline -with-trace
```

This is worth doing when you can, since `execution_report.html` is the only practical way
to right-size the resource tiers in `conf/base.config`.

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
