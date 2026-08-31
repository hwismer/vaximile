# vaximile: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

### Changed

- Replaced `GATK MarkDuplicatesSpark` with the nf-core `SAMTOOLS_SORMADUP` module
  (samtools 1.24: cat | collate | fixmate | sort | markdup). The module is vendored
  unmodified at `modules/nf-core/samtools/sormadup/` and recorded in `modules.json`; a local
  `BAM_MARKDUPLICATES` subworkflow adapts it to this pipeline's conventions.

  **Duplicate marking results will change.** samtools markdup and MarkDuplicatesSpark use
  different algorithms, so duplicate flags, and therefore depth and variant calls, will
  shift. `-S` is set so supplementary alignments of a duplicate template are also flagged,
  which is what the Spark tool did. Optical duplicate tagging (`-d`) is not set, as it needs
  a flowcell-specific pixel distance.

  Output names (`<sample>_<molecule>_markdup.bam`) and the downstream
  `tuple(meta, bam, bai)` shape are unchanged, so no consuming module was touched.

  Three things the swap required:

  - **A `.bai`, not a `.csi`.** The nf-core module's `--write-index` produces a `.csi`;
    every downstream module here declares `path(bai)`, and Strelka and Manta read `.bai`
    specifically. A local `INDEX_BAM` step makes the `.bai`, which keeps the nf-core module
    unpatched and updatable.
  - **Meta bridging.** The module tags and names from `meta.id`, which this pipeline's meta
    map lacks. `BAM_MARKDUPLICATES` adds an `id` for the call and strips it from the output
    again - the meta map is the join key in `DNA_ALIGN_AND_PREPROC` and the grouping key for
    HLA typing and somatic pairs, so an extra key would have silently broken those joins.
  - **Mixed `versions` topic shapes.** Local modules emit a `versions.yml` path; unmodified
    nf-core modules emit a `(process, tool, version)` tuple from `eval()`. `collectFile` read
    the process name in those tuples as a filename, dropping the version from
    `software_versions.yml` and writing a process-shaped junk file into `pipeline_info/`.
    The collector in `main.nf` now renders tuples to YAML and reads the `.yml` files to text
    so both shapes merge. This also fixes the same latent bug for the vendored ASCAT module,
    which does not run under the test profile.

  MultiQC now reports a duplicate rate, which it previously could not: MarkDuplicatesSpark
  was not run with `--metrics-file`, so no duplicate metrics existed at all.

  Resourcing drops from `process_very_high` (16 CPU / 96 GB) plus a
  `--gres=scratch:600G` reservation to the module's `process_medium` (4 CPU / 32 GB) with no
  scratch request. Raise it for WGS if sorting proves slow; note that samtools sort writes
  its temporary files into the task work directory rather than a cluster scratch mount.

- Switched DNA alignment from `bwa-mem2` 2.2.1 to [`minibwa`](https://github.com/lh3/minibwa)
  0.7 (`quay.io/biocontainers/minibwa:0.7--h118bc1c_0`). `BWA_MAP` now runs `minibwa map`
  and `CREATE_BWA_INDEX` runs `minibwa index`; the `-R` read-group string, threading and
  SAM-to-stdout behaviour are unchanged, so no downstream module was touched.

  **Alignments are not identical to bwa-mem2.** Minibwa uses a different algorithm - bwa-mem
  seeding with minimap2 chaining and base alignment - so variant calls will shift slightly.
  Do not mix BAMs from the two aligners within a cohort; re-align rather than resume.

  The index is now two files (`.l2b`, `.mbw`) instead of five, so an existing
  `./resources/bwa/` is **not** reusable. `--bwa_index` rejects a leftover bwa-mem2 index at
  launch with a message naming the aligner change and telling you to rebuild, rather than the
  generic missing-file report - the directory is a valid index, just for the wrong aligner.
  Rebuild by running once without `--bwa_index`. Index construction needs
  ~18x the genome size in RAM (~56 GB for GRCh38), less than bwa-mem2's, so
  `process_high_memory` is unchanged.

  Note that minibwa does not support alternate contigs. The default reference is
  `GRCh38.primary_assembly`, which has none; a reference carrying alts must not be used.

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

- Removed `versions.yml` from the three `storeDir` modules (`PULL_VEP_PVAC_PLUGINS`,
  `PULL_ARRIBA_RESOURCES`, `PULL_CTAT_RESOURCE_BUNDLE`). `storeDir` only short-circuits
  when *every* declared output is already in the store, so adding a `versions.yml` the
  store had never held made these re-download on each run and then fail moving the result
  on top of the copy already there - `mv: inter-device move failed ... unable to remove
  target: Directory not empty`. These processes fetch reference data; their tool version
  was not meaningful anyway.
- `CREATE_BWA_INDEX` published with Nextflow's default `publishDir` mode, which is
  **symlink**, so `./resources/bwa/` held links into `work/`. After `work/` was cleaned
  those dangled, and passing that directory to `--bwa_index` failed validation with a
  contradictory message: files reported missing while being listed as present in the same
  error. `exists()` follows symlinks and returns false for a broken one, whereas `list()`
  still shows the name. The module now copies, and validation reports broken symlinks
  separately from absent files with the actual remedy. The other three index builders
  already used `mode: "copy"`.
- Migrated OptiType to a pinned 1.5.0, with both a conda spec and the matching
  biocontainer, replacing the floating `fred2/optitype:latest` (last pushed 2018). This
  required rewriting the invocation: 1.5.0's entry point is the click group `optitype run`
  rather than `OptiTypePipeline.py`, `-i` is `multiple=True` so each read file needs its
  own flag, and the `OptiType.ini` file the module wrote is replaced by
  `--solver`/`--threads`/`--ilp-threads`. Output naming is unchanged, so the module's
  `optitype_out/*_result.tsv` globs still match. Also removed `which`/`pwd`/`ls` debug
  probes: under `set -e`, `which OptiTypePipeline.py` exits non-zero and aborted the task.
- Fixed the `ensembl-vep` version command in five modules. It ended `grep -Eo '[0-9.]+$'`,
  and Groovy consumed the unescaped `$'` during interpolation, so bash received an
  unterminated quote and the task died with "unexpected EOF while looking for matching `''".
  The `$` anchor is gone. Every other version command was audited for the same hazard.
- Fixed a bash `versions.yml` heredoc being appended to five modules whose script runs
  under `python3`/`Rscript` (`get_rna_strandedness`, `hla_calls_pvac`, `hlahd_to_tsv`,
  `combine_pvacseq_aggregated_report`, `kallisto_tximport`). They now write the file from
  their own interpreter. Introduced by the versions change and invisible to the stub run,
  whose stub blocks are bash.
- `STAR_FUSION` declared `fastq1`/`fastq2` inputs and never passed them, and never set
  `--CPU`, so its 16-CPU label did nothing. Both fixed. Fusion calls may change, since
  STAR-Fusion now has read-level evidence it previously lacked.
- `get_rna_strandedness` left `strandedness` unbound when salmon reported a stranded
  library whose orientation was neither `R` nor `F`, dying with a `NameError` instead of
  the intended message.
- `SOMALIER_EXTRACT` declared its output as `${meta.sample_name}.somalier`, but the tool
  names the file from the BAM's `SM` read-group tag; it now globs `*.somalier`.
- `ADD_VCF_GT_FIELD` copied its input's `.tbi` onto a freshly written *uncompressed* VCF,
  which cannot have a tabix index. The bogus index is gone, along with `FILTER_VCF`'s
  index input, which was staged and never read.
- `star_align` and `kallisto_quant` emitted `path("*")` catch-alls that globbed the whole
  work directory, publishing staged inputs as results. Nothing consumed them.
- `POSTPROCESS_STRELKA` passed `--threads` twice to one `bcftools concat` and indexed an
  intermediate it then discarded.
- The three GATK3 modules hard-coded `-Xmx16g`, decoupled from their label's memory; the
  heap is now derived from `task.memory`.
- The three `PULL_*` modules with `storeDir` no longer populate the real resource store
  under `-stub-run`, where they write empty placeholders that a later real run would have
  reused instead of downloading.
- Unified the split bcftools pin (five modules on 1.23, two on 1.23.1) onto 1.23.1.
- `MHCFLOW` had a leftover debug `ls` and wrote to a different name than it declared;
  `DEEPSOMATIC` was missing a space before a line continuation; `SAMTOOLS_COVERAGE`'s conda
  spec pulled bedtools and htslib it never used. Plus assorted stale comments and typos.
- Every module now reports its tool version. Each writes a `versions.yml` to Nextflow's
  `versions` topic channel, which `main.nf` collects into
  `<outdir>/pipeline_info/software_versions.yml`. `topic:` rather than `emit:` means no
  per-process wiring: 94 modules would otherwise each need threading through their
  subworkflow. 62 modules query the tool directly; the other 32 report the version pinned
  in their own `conda`/`container` directive, since their tool has no usable version flag.
  Stub blocks always report the pinned literal, so a stub run stays offline.
- Gave every module's outputs `emit:` names (48 had none). This was a prerequisite:
  adding `versions.yml` makes every process multi-output, so a call site consuming the
  result bare would break. 62 call sites now select an emit name explicitly.
- Added a `stub:` block to all 94 local modules, and made `-profile test` self-contained,
  so `nextflow run . -profile test -stub-run` exercises the whole DAG - 132 tasks - in
  ~10 seconds offline with no containers, conda or data. This is the first check in this
  repo that actually executes tasks, and therefore the first that can catch a wrong output
  declaration, a tuple-arity mismatch, or a process running the wrong number of times.
- Fixed capture-kit BED paths being passed as bare strings. A relative path in
  `capture_kits.csv` failed with "Not a valid path value" once a task was submitted; they
  are now `file(..., checkIfExists: true)`, so a wrong path is a startup error.
- Fixed `MHCFLOW` declaring `output: tuple path(meta), ...` where `meta` is the metadata
  map, not a file. Nextflow looked for a file named by the map's string form and failed
  with "Missing output file(s) [somatic_name:..., patient:...]". Found by the first stub run.
- Fixed `--bwa_index` causing `BWA_MAP` to align only one sample. The prebuilt index was
  emitted with `channel.fromList`, a one-item *queue* channel, which the first `BWA_MAP`
  task consumes - so every other sample was silently skipped. It is now `channel.value`,
  which is read without being consumed. The auto branch was unaffected because a process
  output that emits exactly once is treated as a value channel.
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
