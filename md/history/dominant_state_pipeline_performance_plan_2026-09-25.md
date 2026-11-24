# Historical pipeline performance proposal (2026-09-25)

Recorded: 2026-09-25

This is the original investigation and proposal. For implemented cache changes,
measured benchmarks, and current remaining work, see the
[2026-10-09 performance report](../dominant_state_pipeline_performance_2026-10-09.md).

## Purpose and method

This session ran the full 55-step `--scope all` pipeline four separate
times to verify code changes (the calibration-bug fix, the RDG-branch
fix, and now the gene-set revision). Every run took 30-45+ minutes. This
note is a real, timing-based investigation of *why*, done by reading
`dominant_pipeline_last_run_summary.csv`/live log timestamps from the
in-progress verification run rather than guessing, plus reading the
actual scoring-loop code for the steps that timing showed were heaviest.
It ends with a concrete, priority-ordered plan.

## What the real numbers say, so far (from the in-progress verification run)

| step | wall time | what it does |
|---|---|---|
| fst_pages | 4m 13s | scan `coverage_index.fst` pages for every marker gene's selected transcript across 3,857 samples |
| clean_cds | 4m 47s | compute clean-CDS usage + all 9 (now: fewer, see below) dominant-state scores per sample |
| score_audit | 8s | percentile-rank recalibration audit |
| state_activity_metadata | 3s | calibrated within-state activity ranks |
| auxiliary_state_screen | 5m 19s | screen ~15 candidate auxiliary states against metadata |
| adjusted_enrichment | 2m 2s | per-field/interaction adjusted enrichment + Fisher tests |
| control_aware | 10s | case-vs-control effect sizes |
| **broad_metadata** | **14m 21s** | **n-wise metadata interaction scan (see below) -- the single heaviest step measured so far** |
| glmnet | (running) | elastic-net regression per state + 12x stability-selection resamples |
| ~46 more RDG/hierarchical/LOO/consensus steps | not yet reached | structurally likely to include the next-heaviest costs (see "LOO-style steps" below) |

Already, 8 steps into 55, roughly 31 minutes have elapsed. `broad_metadata`
alone is worse than every other step measured combined.

## Finding 1: `broad_metadata`'s 4-way interaction scan is the single biggest, most reducible cost found so far

`dominant_state_broad_metadata_model.R` builds every 2-, 3-, and 4-way
combination of 14 metadata fields (`max_interaction_order <- 4L`,
`nwise_priority_fields` has 14 entries): `choose(14,2) + choose(14,3) +
choose(14,4) = 91 + 364 + 1001 = 1456` combinations, each scored with its
own full pass over the sample table (`score_group_terms()`), printed
live as "Scoring n-wise term i/1456". This step alone took 14m21s in this
run.

Checking the *previous* run's output file
(`dominant_state_broad_metadata_nwise_terms.csv`, 152MB, 434,862 rows)
directly:

- **1,001 of the 1,456 scored combinations (68.7%) are 4-way terms.**
- **280,179 of the 434,862 output rows (64.4%) are 4-way terms** -- so
  4-way terms are not being disproportionately filtered out by the
  `min_nwise_n=8`/`min_nwise_out=30` thresholds; they are passing at
  roughly the same rate as lower orders, and dominate both the compute
  and the output size.
- A representative 4-way row: `Cancer_type=osteosarcoma | CONDITION=KD |
  TISSUE=bone | AUTHOR=Jang` -- a single-author, single-cancer-type,
  single-condition stratum. Rows like this are exactly the kind of
  narrow, non-independently-replicated stratum
  `dominant_state_open_analysis_2026-09-24.md`'s Finding 6 already warned
  about trusting at face value (there, for the uORF on/off catalog) --
  4-way metadata interactions in a compendium this heterogeneous are
  structurally unlikely to have enough independent studies per cell to
  mean much, even when they clear the raw sample-count filter.

**Recommendation**: cap `max_interaction_order` at 3 by default (drop the
1,001 four-way combinations, keeping the 455 lower-order ones), or --
better -- switch to apriori-style pruning: only expand a k-way
combination into its (k+1)-way children if the k-way term already cleared
a loose significance/effect-size bar. This is the standard trick used to
make combinatorial interaction mining tractable, and directly targets the
combinations most likely to be real (built from ingredients that already
showed something) rather than scoring the full combinatorial space
uniformly. Expected effect: since 4-way terms are the dominant cost
component of the dominant-cost step, this alone should cut a large
fraction of this step's ~14 minutes, and of its consumers' data volume.
Keep the full order-4 scan available behind an opt-in flag for anyone who
specifically wants the deep combinatorial sweep.

## Finding 2: none of the heavy loops use the server's 48 cores

`nproc` on this server reports **48 cores**. Every heavy loop found in
the codebase this session is a plain serial `lapply()`/`for` loop:

- `broad_metadata`'s 1,456-term n-wise scan (`rbindlist(lapply(seq_along(nwise_combos), ...))`)
- `dominant_state_auxiliary_module_screen.R`'s per-axis metadata
  enrichment loops
- `dominant_rdg_hierarchical_effect_model.R`'s nested
  `for (i in seq_len(nrow(eligible))) { for (study_id in studies) {...} }`
  leave-one-study-out refitting loop -- structurally the most expensive
  *pattern* in the codebase, since it is O(genes x studies): with ~82
  genes and 50-130 studies per design family (matching
  `dominant_rdg_mechanistic_motif_model_summary.csv`'s own `n_genes`/
  `n_studies` columns), this is thousands of individual model refits
- `dominant_rdg_model_agreement_loo.R`'s equivalent LOO structure
- `dominant_rdg_mechanistic_motif_model.R`'s block-CV loop (already
  itself a leave-one-block-out refit inside a larger leave-one-gene-out
  scheme)

Every one of these iterations is embarrassingly parallel (each term/gene/
study refit is fully independent of the others, reads shared input data,
and writes to its own output row). None of them are wrapped in
`parallel::mclapply()`/`BiocParallel::bplapply()` despite both being
already-available, already-used-elsewhere-in-this-ecosystem tools
(`massiveNGSpipe` already uses `BiocParallel` for exactly this kind of
thing). The ~270-330% CPU utilization observed for this R process across
this whole session (a handful of cores, out of 48, mostly from
`data.table`'s own internal multithreading, not from these R-level loops)
confirms these loops are running effectively single-threaded on a machine
built for much higher parallelism.

**Recommendation**: wrap the five loop sites above in
`parallel::mclapply(..., mc.cores = N)` with a conservative, configurable
`N` (e.g. 8-16, not all 48 -- this is a shared server, and the existing
memory-usage caution from earlier this session applies here too). This is
a low-risk change (each iteration is already independent and
side-effect-free besides its own return value) with a potentially large
payoff: an 8-16x reduction in wall time for exactly the steps timing
already shows are the heaviest, without touching the underlying
statistics.

## Finding 3: cache invalidation is coarser than it needs to be, and its own log message doesn't say why

Reading `scripts/dominant_pipeline_cache.R` directly: every step's cache
fingerprint includes not just its own script and declared dependencies,
but a fixed, shared list of four "extra source files"
(`dominant_state_module_definitions.R`, `dominant_coverage_cache.R`,
`dominant_internal_orf_signal_helpers.R`, `dominant_htmlwidgets.R`),
passed identically to *every one of the 55 steps* regardless of whether
that step actually uses them. **Editing any one of these four files
invalidates the cache for all 55 steps**, not just the ones that
genuinely depend on it. This session's own verification runs are a real
example: editing `dominant_state_module_definitions.R` for the gene-set
changes correctly needed to invalidate `fst_pages`/`clean_cds` onward, but
also invalidated late-pipeline bookkeeping steps
(`rdg_review_feedback`, `rdg_validation_labels`, etc.) that have no real
dependency on the gene panel's content.

Separately, `pipeline_cache_status()`'s failure reason is always the same
generic string -- `"inputs, source files, or cache version changed"` --
regardless of which specific dependency or source file actually changed
(confirmed directly in `dominant_pipeline_cache.R`: this string is
returned unconditionally whenever `!identical(fingerprint,
manifest$fingerprint)`, with no diff of *what* differs). This made it
impossible, all session, to tell from the log alone whether a given
step's rebuild was expected (a real upstream change) or worth
investigating as a bug.

**Recommendation**: (1) scope `extra_source_files` per step to only the
helper files that step's own script actually `source()`s, rather than a
blanket list for all 55; (2) when a fingerprint mismatch is found, diff
the old and new fingerprint's `sources`/`dependencies` data frames and
report which specific path's mtime/size/md5 changed in the log message,
turning "inputs, source files, or cache version changed" into something
like "changed: scripts/dominant_state_module_definitions.R (md5)". Low
risk, pure observability improvement, and would have saved real
diagnostic time this session.

## Finding 4: the fast scopes already exist -- use them

`pipeline_scopes()` already defines `core` (`manual_translons ->
fst_pages -> clean_cds -> score_audit -> state_activity_metadata`) and
`metadata` (adds `auxiliary_state_screen -> adjusted_enrichment ->
control_aware -> broad_metadata -> glmnet -> latent`) as named,
narrower scopes, in addition to the much larger `rdg` scope and the
default `all`. **This session's own verification runs used `--scope all`
every time**, including for changes (like the gene-set edits) whose
correctness could have been confirmed with `--scope metadata` alone,
skipping the ~40 RDG/hierarchical/consensus/therapeutic-prediction steps
entirely for that specific check.

**Recommendation**: this is a workflow fix, not a code fix, and the
cheapest win in this whole plan. Use `--scope core` for anything that only
touches gene selection/normalization/basic scoring, `--scope metadata`
for anything touching module definitions or metadata-driven enrichment,
and reserve `--scope all` for final, complete verification before
treating a change as done. Document this explicitly in
`dominant_cell_state_scientific_note.md` or a `HANDOFF.md`-style file so
it isn't only discovered by reading `pipeline_scopes()`'s source.

## Priority-ordered plan

| # | change | expected impact | risk | effort |
|---|---|---|---|---|
| 1 | Use `--scope core`/`--scope metadata` instead of `--scope all` for targeted verification | Immediate, large, for the common case | None (no code change) | Trivial |
| 2 | Cap `broad_metadata`'s `max_interaction_order` at 3, or add apriori pruning for order 4 | Cuts the single heaviest measured step (~14m+) by a large fraction | Low (order-4 terms are mostly narrow/underpowered anyway; keep as opt-in) | Small |
| 3 | Parallelize the 5 identified serial embarrassingly-parallel loops with `mclapply`/`BiocParallel` (conservative core count) | Likely the single largest aggregate win, since it applies to the heaviest steps (`broad_metadata`, `auxiliary_state_screen`, and all LOO-style RDG/hierarchical steps) | Low-medium (needs a sanity check that results are bit-identical to the serial version once parallelized) | Medium |
| 4 | Scope `extra_source_files` per-step instead of blanket-including 4 files in every step's fingerprint | Prevents unnecessary full-pipeline invalidation from unrelated helper-file edits | Low | Small |
| 5 | Make cache-miss reasons specific (diff fingerprints, name the changed file) | Pure observability; saves diagnostic time on every future session | None | Small |

None of these have been implemented yet -- this is a plan, matching the
same pending-review status as the gene-set revision plan. #1 can be
adopted immediately with no code change at all; #2-#5 would need their
own small, verified patches the same way the gene-set changes were
(implement, parse-check, run, compare before/after numbers).
