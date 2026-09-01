# vaximile: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

### Changed

- `APPLY_BQSR_GATHER` moves from `process_medium` to `process_very_high`, 4 CPU to 16, since
  `-@` now feeds `samtools merge` and this process absorbed the compression work the removed
  `SORT_BAM` used to do. The k-way merge itself is serial - the threads go to BAM
  (de)compression - so the gain tapers well before 16, and the tier reserves 96 GB that a
  streaming merge does not need.

- `APPLY_BQSR_GATHER` now merges its shards with `samtools merge` instead of concatenating
  them with `gatk GatherBamFiles`, and emits an indexed BAM. `SORT_BAM` is gone with no
  replacement.

  Concatenation cannot produce a sorted BAM here. `APPLY_BQSR_SCATTER` passes `-L <shard>`
  and GATK emits every read *overlapping* the interval, not only those starting inside it, so
  a read spanning a boundary is written near the start of the later shard while beginning
  before the end of the earlier one. Ordering the shards correctly is necessary but not
  sufficient - the seams are still out of order, by less than a read length, which is why the
  header still claimed `SO:coordinate` and only `samtools index` objected. That was
  previously absorbed by re-sorting the entire gathered BAM; a k-way merge of
  already-sorted shards gets there in one streaming pass instead, so the sort is eliminated
  rather than hidden.

  `-c` is required, not cosmetic. Every shard carries the same `@RG`, and without it merge
  treats those as colliding IDs and renames all but one - verified locally, `ID:S1` became
  `ID:S1-7A2E5CD9` and the reads split across both - which would have broken every
  downstream GATK step. `-p` does the same for `@PG`. `--write-index` with the `##idx##`
  output syntax produces the `.bai` in the same pass, so no separate index step is needed
  either.

  Verified against real samtools on two shards reproducing the reported overlap: concatenation
  fails with the same `Unsorted positions on sequence #1: 52263835 followed by 52263738`,
  while the merge command as written in the module emits positions in order, keeps a single
  `@RG ID:S1` with all reads tagged to it, reports `SO:coordinate`, re-indexes cleanly and
  preserves the read count. Stub run is 132 tasks, two fewer than before, with `SORT_BAM`
  removed and nothing added, and `GET_PILEUP_SUMMARIES` still receiving bam plus bai.

- Reverted `SORT_BAM` -> `INDEX_BAM` on the BQSR path. Replacing the sort with a plain
  `samtools index` was wrong: the gathered BQSR BAM is **not** fully coordinate-sorted, and
  a real run failed with

  ```
  [E::hts_idx_push] Unsorted positions on sequence #1: 52263835 followed by 52263738
  samtools index: failed to create index for "..._bqsr.bam"
  ```

  `APPLY_BQSR_SCATTER` passes `-L <shard>`, and GATK emits every read *overlapping* the
  interval, not only those starting inside it. A read spanning a shard boundary is therefore
  written near the start of the later shard while beginning before the end of the earlier
  one, so concatenating the shards puts positions slightly out of order at each boundary.
  The offsets are under one read length, which is why the header still reads `SO:coordinate`
  and only indexing notices. `APPLY_BQSR_GATHER` ordering its shards by filename is necessary
  but not sufficient.

  The guard held: `samtools index` refused rather than writing a corrupt index, so this
  surfaced as a task failure instead of bad downstream calls. `INDEX_BAM` is still used after
  duplicate marking, where SAMTOOLS_SORMADUP does its own `samtools sort` and the output
  really is sorted.

- Bumped `nf-schema` from 2.1.1 to 2.8.0 in **both** places that pin it. Every run printed

  ```
  WARN: Unrecognized config option 'validation.help.enabled'
  WARN: Unrecognized config option 'validation.help.command'
  ```

  and the warnings were not cosmetic - `nextflow run . --help` failed with
  `Unable to create help message: Specified param 'true' does not exist in JSON schema`.
  2.1.1 introduced the `validation.help` scope but predates Nextflow 26, which does not
  recognise the way that version registers it. Under 2.8.0 the warnings are gone and the
  help message renders, including the custom `validation.help.command`.

  `conf/ucsf_krummellab.config` carries its own `plugins` block, because a merged one
  replaces rather than extends the block in `nextflow.config`. It was still pinned to 2.1.1,
  and being the more specific config it won - so anyone running with it kept both warnings
  and the broken help even on a checkout whose `nextflow.config` said 2.8.0. The two pins
  have to move together.

  Verified: with `-c conf/ucsf_krummellab.config` the warnings are gone and the help message
  renders, where before the fix that exact command reproduced both warnings and the
  `Specified param 'true'` error. The stub run is 134 tasks with a DAG identical to the
  previous commit, and parameter validation still rejects an unknown `--notaparam`.

  `--help` also now exits cleanly. nf-schema's `HelpObserver` prints the message in
  `onFlowCreate` and then calls `session.cancel()`, which stops tasks from being submitted
  but does not stop the entry workflow's body from running. It carried on to build channels
  from parameters a help request never supplies, so `--help` printed the help and then either
  failed parameter validation or parked on the `output {}` block - the same
  `DataflowVariable.get()` hang already documented for `-preview` in docs/nf-core-migration.md.
  The entry workflow now returns immediately when help was requested, after the observer has
  already printed it.

  The guard tests `params.containsKey('helpFull')` rather than reading `params.helpFull`:
  `helpFull` and `showHidden` belong to the plugin and are not declared here, so reading one
  that was not passed emits "Access to undefined parameter" and lets the guard fall through.

  Verified: `--help`, `--helpFull`, `--help <param>` and `--help` under
  `-c conf/ucsf_krummellab.config` all exit 0 with no errors and no warnings, stable over
  three repetitions; per-parameter help returns the shorter targeted message. A normal stub
  run is unaffected at 134 tasks, and an unknown `--notaparam` is still rejected.

- Cleaned out dead and redundant work found by auditing the DAG.

  Removed outright:

  - `COMBINE_FASTQS` - imported in `workflows/vaximile/main.nf` but never invoked. The
    module is gone with it.
  - `--intervals_file` - declared in `nextflow.config`, the schema, `assets/params_example.json`
    and `conf/test.config`, but never read. `SPLIT_INTERVALS` takes its intervals from
    `--capture_kits`. Leaving it in `params_example.json` would have failed schema validation
    for anyone copying that file once the parameter was dropped.
  - STAR's `--quantMode GeneCounts` and the `star_gene_quant` emit - the `ReadsPerGene.out.tab`
    it produced was consumed by nothing, and salmon now covers gene-level counts.
  - `PULL_ASCAT_RESOURCES`'s conda spec, which declared samtools only to feed the versions.yml
    that had to be dropped for `storeDir`.
  - The 13 `.unique{...}` calls in the MultiQC channel mappings. They existed to absorb the
    duplicate sample rows that library deduplication now prevents, so they had become a way
    of hiding duplication rather than avoiding it.

  `SORT_BAM` is replaced by `INDEX_BAM`. It ran a full `samtools sort` over every BQSR'd BAM,
  but `APPLY_BQSR_GATHER` already sorts its shards by filename so `GatherBamFiles`
  concatenates in coordinate order, and `APPLY_BQSR_SCATTER` uses `-L` with no padding
  (`SPLIT_INTERVALS` is called with 0). The gathered BAM is therefore already sorted and only
  needed an index. `samtools index` refuses an unsorted BAM, so if that reasoning is ever
  wrong this fails loudly instead of writing a bad index.

  The tumour/normal pairing that was hand-rolled three times - branch, join, rebuild the same
  five-field somatic_meta - is now `pair_tumor_normal()`. It is arity-agnostic, so the BAM
  sites (meta, bam, bai) and the pileup site (meta, pileup) share it, and about 45 lines
  collapse to three calls.

  `STAR_ALIGN` moves from `process_very_high` to `process_max`, 16 to 32 threads;
  `--runThreadN` already took `task.cpus`. STAR scales sub-linearly past ~16 threads, so
  expect well under 2x, and on a busy cluster the wider reservation may cost more in queue
  time than it saves.

  Verified by stub run: 134 tasks, and a per-process diff against the previous commit shows
  exactly `SORT_BAM` replaced by `INDEX_BQSR_BAM` and nothing else moved - the somatic callers
  still receive the same number of pairs through the new helper.

- Dropped kallisto; salmon is now the only RNA quantifier. `CREATE_KALLISTO_INDEX`,
  `KALLISTO_QUANT` and `KALLISTO_TXIMPORT` are removed, along with `--kallisto_index`. The
  two ran side by side producing the same transcript and gene abundances.

  Gene-level TPM comes from salmon's own `--geneMap`, which takes the GTF and writes
  `quant.genes.sf` next to `quant.sf`, rather than from a separate tximport step. That
  removes a whole process and its R dependency stack (tximport, rtracklayer, dplyr, readr)
  from the conda environments. It is copied out as `<sample>_<molecule>.gene_tpm.tsv`,
  since salmon names it identically for every sample, and published under `salmon/` where
  `kallisto/` used to be.

  The VCF expression annotators were adapted to salmon's column names. vcf-expression-annotator
  has no salmon parser, so both now go through `custom`: transcript-level reads `quant.sf`
  with `-i Name -e TPM` in place of the `kallisto` format, and gene-level reads
  `quant.genes.sf` with `-i Name` where tximport's table had called that column `ENSEMBLID`.
  Both keep `--ignore-ensembl-id-version`, since the Ensembl cDNA FASTA carries versioned IDs.

  **The gene symbol column is gone.** tximport's table was ENSEMBLID/TPM/Gene; salmon's is
  Name/Length/EffectiveLength/TPM/NumReads, with no symbol. Nothing consumed that column -
  the annotator matches on Ensembl IDs - but it was useful in the published file. Recovering
  it means a gene_id-to-gene_name map from the GTF, which is an awk step rather than a
  reason to keep R.

  Verified by stub run: 134 tasks, and a per-process diff against the previous commit shows
  exactly the three kallisto processes gone and nothing else changed. `quant.sf` stages into
  the transcript annotator and the gene TPM table into the gene annotator, and the table
  publishes to `<patient>/<sample>/salmon/<sample>_<molecule>.gene_tpm.tsv`. The salmon run
  itself has not been executed, so the `--geneMap` output and the annotators' parsing of it
  are unverified against real data.

- Fixed sample-level outputs publishing into a literal `null` directory. Deduplicating
  libraries replaced the scalar `somatic_name` on sample-level metas with a `somatic_names`
  list, but nine publish path closures in `main.nf` still read `meta.somatic_name`, so
  germline VCFs, per-sample HLA calls and the RNA gene abundance landed in
  `<outdir>/<patient>/null/...`. Stub runs did not catch it: the paths resolve and the run
  succeeds, it is only the directory name that is wrong.

  `publish_scope()` now picks the pair for pair-level outputs and the sample for
  sample-level ones. This changes where those files land - germline calls move from
  `<patient>/<somatic_name>/germline/` to `<patient>/<sample_name>/germline/`, and per-sample
  HLA likewise - which is the point: a library shared between pairs is processed once and no
  longer belongs to exactly one of them. Pair-level outputs (somatic VCFs, pVACtools reports,
  the merged-BAM HLA call) are unchanged.

- `PULL_ASCAT_RESOURCES` now uses `storeDir`, so its four Zenodo downloads are fetched once
  and reused rather than re-downloaded on every run. It was the only `PULL_*` module without
  one; `PULL_VEP_CACHE`, `PULL_ARRIBA_RESOURCES`, `PULL_CTAT_RESOURCE_BUNDLE` and
  `PULL_VEP_PVAC_PLUGINS` already had it.

  Its store path carries the chr-prefix flag - `./vaximile_resources/ascat_hg38_chr` or
  `ascat_hg38_nochr`. This is the only `PULL_*` process whose output depends on an input:
  with `--reference_includes_chr_prefix` it rewrites every loci file with sed. storeDir is
  keyed on the path alone rather than on the task hash, so a single shared directory would
  return chr-prefixed loci to a later run that asked for unprefixed ones - silently, and
  visible only as wrong ASCAT calls.

  Its `versions.yml` output was removed, as on the other four. storeDir only short-circuits
  when every declared output is already in the store, and a `topic:`-routed versions.yml
  never lands there, so the process would re-run on a populated store and then fail moving
  its results on top of the existing copies. ASCAT resources therefore no longer contribute
  a samtools version to `software_versions.yml`, matching the other `PULL_*` modules.

  Verified by stub run (137 tasks, unchanged) and by evaluating the directive line directly
  with both flag values, since `workflow.stubRun` short-circuits it during a stub run: the
  two settings write to separate stores as intended.

- `--vep_cache` is now optional. When it is not given, `PULL_VEP_CACHE` downloads the
  release-115 human GRCh38 cache (~24 GiB) from Ensembl's FTP into
  `./vaximile_resources/vep_cache`, and later runs reuse it from the `storeDir`. Passing
  `--vep_cache` still skips the download and uses the given directory, and it was removed
  from the schema's required list.

  The release is pinned to 115 rather than exposed as a parameter: VEP rejects a cache whose
  version differs from its own, so it has to move with the pinned `ensembl-vep` 115 the
  annotation modules run. Species and assembly are likewise fixed to homo_sapiens/GRCh38 to
  match the default reference - anything else has to supply `--vep_cache`. The plain cache is
  used rather than refseq or merged, because the VEP modules pass neither flag.

  The process carries both `process_single` and `process_long`. A `time` directive written
  in the module would have been silently overridden by the 4 h default `conf/base.config`
  sets for every process, since config directives beat module ones; `process_long` is the
  20 h tier, and the two labels apply cumulatively.

  Verified by stub run on both branches: with `--vep_cache` supplied the DAG is unchanged at
  137 tasks and PULL_VEP_CACHE does not run; without it the run is 138 tasks, the cache
  directory is staged into every VEP task as `vep_cache/homo_sapiens/`, and the process
  resolves to 20 h / 1 CPU once the test profile's 1 h resourceLimits clamp is lifted. The
  download URL, its ~23.5 GiB size and the `homo_sapiens/115_GRCh38/` layout inside the
  tarball were each checked against the live FTP; the tarball itself was not downloaded.

- Samples shared between tumour/normal pairs are now processed once instead of once per
  pair. The samplesheet carries one row per (library, pair), so a normal used as the control
  for two tumours is listed twice with different `somatic_name`s. Because `somatic_name` was
  part of the meta map, and the meta map is a channel item's identity, those rows were
  distinct items and every sample-level process ran twice on byte-identical data.

  Measured on a two-tumour/one-normal samplesheet, these all ran twice for the shared normal:
  `FASTP`, `BWA_MAP`, `SAMTOOLS_SORMADUP`, `INDEX_BAM`, base recalibration and BQSR (scatter
  and gather), `SORT_BAM`, `GET_PILEUP_SUMMARIES`, `SOMALIER_EXTRACT`, the samtools QC trio,
  the whole HaplotypeCaller chain plus `DEEPVARIANT` and `STRELKA_GERMLINE`, and HLA typing
  (`MHC_REGION_FASTQS`, `OPTITYPE`, `HLAHD`). With N tumours on one normal it was N times,
  and the normal's germline VCF was written N times to the same published name.

  `dedupe_libraries()` collapses rows sharing patient + sample_name + molecule at the entry
  point, replacing the scalar `somatic_name` with a `somatic_names` list. `expand_pairs()`
  and `fan_out_pairs()` restore the per-pair view at the seven tumour/normal join sites, so
  nothing downstream of those joins changed shape - roughly 140 of the 155 `somatic_name`
  references are past the join and untouched.

  Two forms are needed because the joins are not all alike: most key on the name as a string,
  but `PVAC_VCF_PHASING` joins `somatic_meta.normal_meta` against the germline VCF channel and
  `somatic_meta.tumor_meta` against the tumour BAMs - whole meta maps as keys, which only
  match if the channel side carries the scalar `somatic_name` too.

  Rows treated as one library must agree. If rows sharing a `sample_name` list different
  FastQs, or disagree on `sample_type`, `sequencing_type`, `capture_kit` or `sex`, the run
  fails at launch naming the conflict rather than silently processing whichever row sorted
  first. Distinct data needs a distinct `sample_name`.

  Verified two ways. On the existing single-pair test samplesheet the DAG is unchanged -
  137 tasks, and a per-process count diff against the previous commit is empty, so nothing
  shifted for the non-shared case. On a two-tumour/one-normal sheet the total falls from 232
  to 191 tasks (18%), every sample-level step for the shared normal drops from 2 to 1, and
  the pair-level steps stay at 2 with each pair staging the correct BAMs: Tumor1+Normal1 and
  Tumor2+Normal1 against the same singly-produced normal. Both validators were triggered and
  produce the intended message.

  This invalidates existing `-resume` caches, since every DNA task rehashes.

- Fixed `SAMTOOLS_SORMADUP` failing on every sample with `[main_cat] ERROR: input is not
  BAM or CRAM`. The module's pipeline opens with `samtools cat`, which reads BAM and CRAM
  only, and `BWA_MAP` was handing it the SAM minibwa writes. MarkDuplicatesSpark accepted
  SAM, so nothing had needed a BAM before. Stub runs could not catch this: a stub `touch`es
  its outputs, so no tool ever saw the file's format.

  `BWA_MAP` now pipes `minibwa map` into `samtools view -1` and emits an unsorted BAM.
  Converting in the pipe rather than in a following process means the SAM is never written
  at all, which for WGS removes a multi-hundred-GB intermediate - so this is now faster than
  the bwa-mem2 pipeline it replaced, not merely fixed. `-1` (fast compression) is used
  because SORMADUP reads the file once.

  `BWA_MAP` consequently has no container: it needs minibwa and samtools in one task, and
  every published minibwa image carries minibwa alone. It is conda-only, falling back to the
  host `PATH` under the container profiles as several modules already do. See docs/usage.md.

  Verified against real samtools rather than stubs: the failure reproduces exactly on a SAM,
  the new BAM feeds the full cat/collate/fixmate/sort/markdup chain to exit 0, the @RG line
  survives it, output is SO:coordinate, a planted duplicate pair is flagged (2 of 6 reads),
  and `samtools index` produces the .bai INDEX_BAM expects.

- Gave the two per-sample alignment steps more cores. `BWA_MAP` moves from
  `process_high` to `process_very_high` (8 -> 16 CPU, 48 -> 96 GB), matching `STAR_ALIGN`
  and the tier's stated purpose. minibwa's `-t` scales well, so this is close to a
  straight speedup.

  `SAMTOOLS_SORMADUP` gets the same 16 CPU / 96 GB, but set by name in `conf/base.config`
  rather than by label: it is a vendored nf-core module kept byte-identical to upstream so
  `nf-core modules update` keeps working, its `process_medium` label lives in that file,
  and a label cannot be added from config.

  Expect less from the extra cores here than from the aligner's. The module pipes five
  samtools stages and passes `--threads` to four of them, so the nominal thread request is
  4x this number, and only `samtools sort` really parallelises - the middle stages hand off
  uncompressed BAM and stay mostly idle. If sorting is still the bottleneck, the next lever
  is `ext.args4 = '-m 2G'`, which cuts temp-file spilling, rather than more cores.

  Both now request 96 GB where they previously asked for 48 and 32. Neither tool needs that
  much; it comes bundled with the tier. On a busy cluster the larger reservation may cost
  more in queue time than it saves in runtime.

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
