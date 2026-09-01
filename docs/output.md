# vaximile: Output

## Introduction

This document describes the output produced by the pipeline. Publishing is driven by the
`output {}` block in `main.nf`, which nests results under
`<outdir>/<patient>/<somatic_name>/`.

## Directory structure

```text
<outdir>/
├── pipeline_info/
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   └── pipeline_dag.html
└── <patient>/
    ├── multiqc/
    │   └── <patient>_report.html
    ├── pvactools_report/
    │   └── combined MHC class I aggregated report across the patient's pairs
    ├── <somatic_name>/          per tumour/normal pair
    │   ├── variants/     somatic VCF (+ index) and its tabular form
    │   ├── ascat/        CNV segments, purity/ploidy, BAF/LogR, metrics and plots
    │   ├── hla/          pVACtools-formatted allele list for the merged pair BAM
    │   └── pvactools/    pVACseq and pVACfuse prediction directories
    └── <sample_name>/           per library
        ├── alignment/    duplicate-marked BAM, recalibrated BAM, STAR BAM (+ indexes)
        ├── germline/     germline VCF (+ index) and its tabular form, normals only
        ├── hla/
        │   ├── optitype/ class I calls and coverage plot
        │   ├── hlahd/    class I and II calls
        │   └── pVACtools-formatted allele list
        └── salmon/       gene-level abundance, RNA samples only
```

Outputs split by whether they belong to a pair or to a single library. A library shared
between pairs - a normal used as the control for several tumours - is processed once, so
its germline calls and per-sample HLA typing sit under the sample rather than being
duplicated under each pair.

## Results by stage

### Quality control

`multiqc/<patient>_report.html` aggregates, per patient: fastp reports for every DNA and
RNA library, GATK base recalibration tables, samtools flagstat/coverage/idxstats, STAR
alignment logs, salmon quantification summaries, OptiType and HLA-HD calls, VEP summaries
for both germline and somatic annotation, and somalier relatedness tables.

The somalier pairs/samples tables in that report are the check for sample swaps — confirm
the tumour and normal of each pair relate to each other before trusting downstream calls.

### Variant calling

`variants/` holds the merged, filtered and annotated somatic callset. Three callers
contribute: Mutect2, Strelka2 (with Manta for indel candidates) and DeepSomatic. Each is
filtered to `PASS` and normalised, then combined with GATK3 `CombineVariants` under
`--minimumN 2` - an **n-1 consensus**, so a variant is kept when at least two of the three
callers report it. The published VCF is therefore a consensus product rather than any
single caller's raw output. The `.tsv` alongside it is the same content flattened for
spreadsheet use.

`germline/` holds the consensus germline callset for the normal sample, built the same way
from three callers: GATK `HaplotypeCaller` (CNN-scored and tranche-filtered), Strelka2 in
germline mode, and DeepVariant, combined under the same 2-of-3 rule. This VCF is also the
phasing input for proximal variant detection in pVACseq.

### HLA typing

Class I alleles come from OptiType, class I and II from HLA-HD. Both run on each individual
library and on a per-`somatic_name` merged BAM; the merged calls are what feed pVACtools,
since they use the most reads available for the patient.

### Neoantigen prediction

`pvactools/` contains the pVACseq output for somatic variants and the pVACfuse output for
RNA fusions called by Arriba and STAR-Fusion. `<patient>/pvactools_report/` is the
aggregated MHC class I report combined across the patient's tumour/normal pairs.

### Expression

`salmon/` publishes `quant.genes.sf`, the gene-level TPM table salmon aggregates from
`--geneMap`. The transcript-level `quant.sf` and the rest of the salmon run directory are
produced but consumed internally - for VCF expression annotation, strandedness prediction
and MultiQC - rather than published.

### Alignments

`alignment/` holds three BAMs per library, each with its index, distinguished by suffix:
`_markdup.bam` after duplicate marking, `_bqsr.bam` after base recalibration, and
`_STAR_sorted.bam` for RNA. The somatic callers read the recalibrated BAMs; Strelka and
Manta read the duplicate-marked ones.

### Copy number

`ascat/` publishes the whole ASCAT result set for a pair - segments, CNV calls,
purity/ploidy estimates, BAF and LogR tables, run metrics and plots.

## Not published by default

ASCAT copy number output, per-caller intermediate VCFs, aligned BAMs, RNA fusion tables and
the built alignment indices remain in the work directory. Add them to the `publish:` block
in `main.nf` and give them a `path {}` in the `output {}` block to change that.
