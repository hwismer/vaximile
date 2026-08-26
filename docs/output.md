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
    └── <somatic_name>/
        ├── variants/     somatic VCF (+ index) and its tabular form
        ├── germline/     germline VCF (+ index) from the normal sample
        ├── hla/
        │   ├── optitype/ class I calls and coverage plot
        │   ├── hlahd/    class I and II calls
        │   └── pVACtools-formatted allele list
        ├── pvactools/    pVACseq and pVACfuse prediction directories
        └── kallisto/     gene-level abundance for the RNA sample
```

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
contribute: Mutect2, Strelka2 (with Manta for indel candidates) and DeepSomatic. Calls are
merged and postprocessed before annotation, so the published VCF is a consensus product
rather than any single caller's raw output. The `.tsv` alongside it is the same content
flattened for spreadsheet use.

`germline/` holds HaplotypeCaller output for the normal sample, CNN-scored and filtered.
This VCF is also used as the phasing input for proximal variant detection in pVACseq.

### HLA typing

Class I alleles come from OptiType, class I and II from HLA-HD. Both run on each individual
library and on a per-`somatic_name` merged BAM; the merged calls are what feed pVACtools,
since they use the most reads available for the patient.

### Neoantigen prediction

`pvactools/` contains the pVACseq output for somatic variants and the pVACfuse output for
RNA fusions called by Arriba and STAR-Fusion. `<patient>/pvactools_report/` is the
aggregated MHC class I report combined across the patient's tumour/normal pairs.

### Expression

`kallisto/` publishes gene-level abundance. Transcript-level abundance and salmon
quantification are produced but consumed internally for VCF expression annotation rather
than published.

## Not published by default

ASCAT copy number output, per-caller intermediate VCFs, aligned BAMs, RNA fusion tables and
the built alignment indices remain in the work directory. Add them to the `publish:` block
in `main.nf` and give them a `path {}` in the `output {}` block to change that.
