# vaximile: Usage

vaximile predicts tumour neoantigens from paired tumour/normal bulk DNA and matched tumour
RNA sequencing. See the [README](../README.md) for what the pipeline runs.

## Prerequisites

Nothing has to be downloaded by hand. The VEP cache, ASCAT and Arriba resources and the
CTAT bundle are fetched on the first run into `./vaximile_resources/` and reused after
that, and any index left unset is built from the reference.

Two things are worth supplying if you already have them, because both are slow to produce:

| Parameter | Notes |
| --------- | ----- |
| `--vep_cache` | Skips a ~24 GiB download. Point at the directory *containing* `homo_sapiens/`. Must be release **115** - VEP rejects a cache whose version differs from its own - and the plain cache, not `refseq` or `merged`. Required if you are not running human GRCh38, since the automatic download is fixed to `homo_sapiens_vep_115_GRCh38`. |
| `--bwa_index` | minibwa. See [reusing an index](#reusing-a-prebuilt-minibwa-index) - the filenames matter. |
| `--star_index` | STAR. Built with `star=2.7.11b`, and STAR indices are version-specific, so an index built by an older STAR is not reusable. |
| `--salmon_index` | salmon `1.11.4`. |

The `PULL_*` processes run on the local executor, so the head node needs outbound network
access.

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

Parameters can come from a file instead, which is easier to version:

```bash
nextflow run . -profile conda -params-file assets/params_example.json -resume
```

## Before a real run

```bash
nextflow run . -profile test -stub-run
```

Every module has a `stub:` block, so this runs the entire DAG in seconds - offline, no
containers, no conda, no data. It verifies wiring, not science: that every output
declaration resolves, that tuple arities match at each call site, and **how many times each
process runs**.

That last one matters. A reference channel built as a queue instead of a value channel is
consumed by the first task, so an aligner silently processes one sample and skips the rest,
which no amount of linting or `-preview` will reveal. Check the counts against your
samplesheet:

```bash
grep -oE 'Submitted process > [A-Za-z0-9_:]+' .nextflow.log | sed 's/.*://' | sort | uniq -c | sort -rn
```

The `test` profile points every reference at an empty placeholder under `assets/test/refs/`
purely to keep this offline, since Nextflow stages inputs even under `-stub-run`. ASCAT is
skipped there via `ext.when = false`, because its stub runs `Rscript` to capture a version
and so needs its container.

## Reusing a prebuilt minibwa index

`--bwa_index` takes a directory holding the two files `minibwa index` produces, **named
after the prepared reference FASTA**:

```
<reference_fa stem>_prc.fa.l2b
<reference_fa stem>_prc.fa.mbw
```

The naming is not incidental: `BWA_MAP` passes the reference FASTA to minibwa as the index
prefix, so an index built from a differently named FASTA is unusable even if it is
otherwise valid. `PREPARE_FASTA` renames the reference to `<stem>_prc.fa`, which is where
that suffix comes from. The pipeline validates this at launch and names any missing file.

The simplest way to get a valid directory is to run once without `--bwa_index` and reuse
`./resources/bwa/`.

**bwa-mem2 indices are not reusable.** They are a different set of five files, and pointing
`--bwa_index` at one fails at launch saying so. Because minibwa alignments are not
identical to bwa-mem2's, do not mix BAMs from before and after the switch within a cohort.



## Software versions

Every module writes a `versions.yml`, collected into
`<outdir>/pipeline_info/software_versions.yml`, so results record the tools that produced
them. `unknown: unknown` marks the few modules that declare neither `conda` nor
`container` and rely on the host `PATH`.

GATK versions are pinned per module on purpose and should not be unified: the CNN germline
chain needs 4.3.0.0, `CombineVariants` and `ReadBackedPhasing` exist only in GATK3, and
everything else is on 4.6.1.0. Treat any GATK change as a change to results.

