Dominant state / fatigue pipeline code review

Scope: `dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R`, the
root compatibility launcher, and the analysis scripts it orchestrates, with
special attention to FST-heavy steps.

Findings

1. High impact: pipeline steps were always treated as dirty.

   The entry point reran every selected script even when all outputs and inputs
   were unchanged. This made `clean_cds` and `uorf_regulation` expensive on
   routine reruns because both call `coverageByTranscriptFST()` across target
   transcripts.

Fix implemented: step-level cache manifests in
`dominant_cell_states/.pipeline_cache/`. A step now skips when its declared
outputs exist and the cache fingerprint still matches.

2. High impact: no manual invalidation flag existed for internal model changes.

   Output files alone cannot tell whether a script was intentionally changed or
   whether the biological model should be recomputed. The cache key now includes
   the step script md5, the shared module-definition script md5, selected input
   file signatures, upstream output signatures, and a user-provided cache
   version.

Use `--cache-version some-new-name` or set
`DOMINANT_PIPELINE_CACHE_VERSION=some-new-name` when the model logic changed
outside the declared scripts or when you want a reproducible full refresh
without deleting outputs.

Bootstrap existing outputs with `--seed-cache` or
`run_dominant_cell_state_pipeline("all", seed_cache = TRUE)`. This creates
manifests for already completed CSV/RDS outputs without running the heavy
scripts.

3. Medium impact: FST readiness was not part of rerun logic.

   The readiness step now fingerprints the required/additional page lists and
   any page paths named in those files. If a missing FST page is later
   downloaded, the readiness cache invalidates automatically.

4. Medium impact: repeated FST reads can still happen inside one cold run.

   Fix implemented: `dominant_coverage_cache.R` stores one coverage RDS per
   `(gene, tx_id, display-region, coverage-index signature)`. Both `clean_cds`
   and `uorf_regulation` now use this shared accessor, so a forced cold run
   should not re-read the same transcript coverage twice unless the display
   window or coverage index changes.

5. Medium impact: dependency declarations are now explicit but conservative.

   Each step has a declared output set and dependency set in
   `run_dominant_cell_state_pipeline.R`. This is much easier to reason about than
   implicit script side effects, but when a new output or input is added to a
   model script, the maps should be updated in the same commit.

6. Low/medium impact: uORF regulation repeatedly loaded transcript regions.

   Fix implemented: `dominant_uorf_regulation_inference.R` now preloads leader
   ranges once beside the existing CDS/mRNA preload and passes those objects into
   `measure_uorf_features()`. This removes one `loadRegion()` round trip per
   marker transcript.

How to run

```sh
Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope all
Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope all --seed-cache
Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope all --force
Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope all --force-step clean_cds,uorf_regulation
Rscript dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R --scope all --cache-version model-v2
```

Internal checks

```sh
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); source("dominant_cell_states/tests/testthat/test-pipeline-cache.R")'
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); source("dominant_cell_states/tests/testthat/test-coverage-cache.R")'
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); source("dominant_cell_states/tests/integration/timing_pipeline_cache.R")'
```

Current optimization boundary

The cache now answers the main operational problem: do not rerun an analysis
that is already complete and unchanged, and do not re-read the same transcript
coverage across the two FST consumers. The next boundary is page-level sharing:
`coverageByTranscriptFST()` may still read the same underlying page for nearby
transcripts when their transcript windows differ. If cold runs remain too slow,
the next optimization should group target transcripts by FST page and read each
page once before slicing transcript windows.

Sanity pass, 2026-06-18

1. Result inventory and invariant checks passed.

   Added `dominant_cell_states/tests/integration/result_sanity_review.R`. It reads
   every generated CSV, inventories PNG/PDF/HTML outputs, checks the main atlas
   and validation dossier tables, and writes:

   - `dominant_cell_states/results/test_reports/result_sanity_report.md`
   - `dominant_cell_states/results/test_reports/result_sanity_csv_inventory.csv`
   - `dominant_cell_states/results/test_reports/result_sanity_figure_inventory.csv`
   - `dominant_cell_states/results/test_reports/result_sanity_key_checks.csv`
   - `dominant_cell_states/results/test_reports/result_sanity_top_primary_candidates.csv`

   Current result: 445 CSV files readable, no CSV read errors, one expected empty
   counterfactual table, 3458 figure/HTML files, no empty figures, no private
   htmlwidget dependency folders, 1680 atlas cards, 184 primary validation rows,
   and 3857 clean-CDS runs.

2. Shared HTML widget dependencies are now content-aware.

   `dominant_htmlwidgets.R` no longer deletes and recopies an existing shared
   widget dependency directory when the source dependency has identical content.
   This matters for figure-heavy reruns where each plotly HTML save previously
   rewrote the same `plotly`, `htmlwidgets`, `jquery`, and related folders.
   Added `tests/testthat/test-htmlwidgets-shared-libs.R` to cover unchanged and
   changed dependency-copy behavior.

3. FST page readiness is explicit after download.

   `dominant_cell_state_required_fst_pages.R` now has
   `DOMINANT_CHECK_LOCAL_FSTS=auto` behavior: it avoids probing absent mounts,
   but fills `local_file_exists` once the local page directory exists. The
   current `required_fst_pages.csv` was refreshed from the existing page table
   and reports 252/252 required pages present.

   Note: a forced full recomputation of the `fst_pages` audit is still slower
   than expected because it re-evaluates transcript/page mappings for the full
   290-gene dominant-state target set. That is not a duplicate-gene bug; the RDG
   atlas covers 120 measurable RDG genes, while the FST target list covers all
   dominant-state marker/prior genes. If this step becomes a frequent bottleneck,
   add a fast metadata-only local-readiness refresh path or cache transcript
   display-region/page mappings separately.

4. Validation dossier plot readability improved.

   `rdg_validation_dossier_priority.png` now includes a visible review-tier
   legend with human-readable labels. The previous plot used color without an
   explanatory legend. The `rdg_validation_dossiers` step was rerun and completed
   successfully in about 4 seconds.

5. Verification run.

   The following internal checks passed after the patch:

   - `test_pipeline_cache.R`
   - `test_coverage_cache.R`
   - `test_htmlwidgets_shared_libs.R`
   - `test_internal_orf_signal_metrics.R`
   - `test_pipeline_invariants.R`
