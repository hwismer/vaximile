# vaximile: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

### Changed

- Reorganised the repository into the nf-core directory layout: `main.nf` entry point,
  `workflows/vaximile/`, one directory per module under `modules/local/`, one directory per
  subworkflow under `subworkflows/local/`, and parameter defaults moved to
  `nextflow.config` and `nextflow_schema.json`. Process bodies are unchanged. See
  [docs/nf-core-migration.md](docs/nf-core-migration.md) for what remains.

- Processes that use metadata inside their script now declare it as explicit `val()`
  inputs instead of reaching into the meta map. 30 modules and 12 caller files changed:
  each module takes named arguments (`sample_name`, `molecule`, `sequencing_type`,
  `somatic_name`, `tumor_sample_name`, `normal_sample_name`, `tumor_sequencing_type`) and
  the workflow projects those fields at the call site. Nested access such as
  `somatic_meta.tumor_meta.sample_name` is gone from every script body.

  The meta map stays as element 0 of each tuple, so no `join`, `groupTuple`, `branch` or
  `combine(by:)` key changed - the projection maps are always chained *after* the keyed
  operation. Verified mechanically: every process input tuple arity matches the tuple
  built at its call site, no keyed operation differs from the previous commit, and
  `nextflow run . -profile test -preview` builds the same DAG as before.

  `tag` directives and `ext.prefix` defaults still read the map. Both are cosmetic, and
  narrowing the change to script bodies kept the diff reviewable.
- Added verified conda specs to 39 container-only modules, so `-profile conda` now covers
  70 of 94 local modules instead of 31. Containers are retained; every spec was checked
  against bioconda with `conda search` and pins the container's version (exceptions:
  `bwa-mem2=2.2.1`, inferred because the container tag does not state a version, and
  `ensembl-vep=115`, which is how bioconda publishes release 115.0). The 20 still
  container-only have real blockers — GATK3's licensed jar, unpackaged DeepSomatic and
  HLA-HD, major version gaps, and absolute container paths in scripts — tabulated in
  docs/usage.md.
- Converted all 94 local modules to nf-core process conventions. Each now carries a
  `process_*` label instead of its own `cpus`/`memory` (tiers in `conf/base.config`), 67
  support `task.ext.prefix`, and 19 read `task.ext.args` from `conf/modules.config`.
  Output filenames are unchanged. 27 modules cannot use `ext.prefix` without risking a
  glob that captures their own inputs; those are enumerated in
  [docs/nf-core-migration.md](docs/nf-core-migration.md).
- Resource requests now come from eight tiers. Each module was placed in the smallest tier
  meeting or exceeding its previous CPU and memory request, so 79 of 94 request somewhat
  more than before and none requests less.

### Fixed

- Every parameter is now printed at the start of a run, grouped by the sections in
  `nextflow_schema.json`, with `*` marking values that differ from the schema default.
  nf-schema's `paramsSummaryLog()` prints only non-default values, which hid every
  reference URI and GATK resource VCF a run depends on - exactly what you want recorded
  alongside a set of results.
- Added `default` to the 20 parameters in `nextflow_schema.json` that have one in
  `nextflow.config`. Previously only 5 declared a default, so `--help` under-reported and
  nothing could tell an overridden value from a default one.
- `--bwa_index` now works. Previously the main workflow bypassed the `BWA_INDEX`
  subworkflow entirely when the parameter was set, emitting one channel item per index
  file instead of a single item holding the whole set - so `BWA_MAP` did not receive a
  usable index. Both paths now go through `BWA_INDEX` and emit the same shape, and a
  supplied directory is validated at launch: it must contain the five files `bwa-mem2
  index` writes, named after the prepared reference (`<stem>_prc.fa.<ext>`), because
  `BWA_MAP` passes that FASTA to bwa-mem2 as the index prefix. Missing files are listed
  by name rather than surfacing as a per-sample bwa-mem2 failure mid-run.
- `APPLY_BQSR_GATHER` computed `sorted_bams` and then passed the unsorted `bams` to
  `GatherBamFiles`, which concatenates without re-sorting. Shard order out of
  `groupTuple` is not guaranteed, so the merged BAM could be mis-ordered. It now uses the
  sorted list, matching `MUTECT2_GATHER_VCFS` and `HAPLOTYPE_CALLER_GATHER_VCFS`.
- `PIPELINE_COMPLETION`'s `workflow.onComplete` handler threw
  `NullPointerException: Cannot get property 'success' on null object` on every run,
  because `workflow` resolves to null inside a closure invoked from a named workflow body.
  The metadata object is now captured before the closure.
- `HLAHD` declared its directory output as the glob `*/result/`. A glob ending in `/`
  never matches, so the task failed with "Missing output file(s) `*/result/`" even though
  HLA-HD had run and written its results. Introduced when the explicit
  `./${meta.sample_name}/result/` path became a glob during the ext.prefix conversion. Now
  `path("*/result", type: 'dir')`. No other module has a trailing-slash glob.
- `SALMON_QUANT` passed `--libType` after `-1/-2`, which salmon rejects outright
  ("The (--libType/-l) option must precede the input files"). Introduced when
  `--libType A` and `--validateMappings` — which sat on opposite sides of the read
  files — were collapsed into a single `$args` placed after them. `$args` now precedes
  `-1/-2`. The other 18 modules using `ext.args` were audited for the same reordering;
  10 also had non-contiguous flags but all are order-insensitive option-only CLIs
  (VEP, GATK, bcftools, DeepVariant/DeepSomatic) or place `$args` after every
  positional (pVACseq, pVACfuse).
- Disabled `timeline`, `report` and `trace` in `nextflow.config`. Enabling them makes
  Nextflow inject `command -v ps || exit 1` into every task wrapper, which killed tasks
  running in containers without `procps` (STAR, DeepVariant, DeepSomatic, bcftools) before
  their tool ran. The symptom is an empty `.command.out` and only the "Command 'ps'
  required by nextflow" line in `.command.err`. Request the reports per-run with
  `-with-report`/`-with-timeline`/`-with-trace` instead. `dag` is unaffected and stays on.
- `nextflow.config` interpolated `${manifest.name}` inside the `validation.help.command`
  string, where `manifest` is not in scope. This failed config parsing outright, so every
  `nextflow` invocation from the repository root aborted before compiling anything.

- `mhcflow` referenced its conda environment as `./envs/mhcflow.yml`, which resolved
  against the launch directory. The file now lives beside the module and is referenced with
  `${moduleDir}`.
