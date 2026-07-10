# Dominant Cell States

Standalone analysis workspace for dominant-state and ribosome decision graph
(RDG) discovery. This directory is being split out of the RiboCrypt source tree
so the analysis code, Shiny browser, notes, and tests can live in their own
GitHub repository.

## Current Status

The project is a review-grade discovery system, not a final claim engine. It
combines dominant-state scores, clean-CDS quantification, RDG branch allocation,
model-agreement layers, evidence-gap triage, active-review batches, browser
links, feedback labels, and a Shiny review app.

The latest algorithm focus is merged-replicate translon shifts: find condition
contrasts where clean-CDS allocation changes relative to leader or overlapping
uORF allocation. The current first-stop table is:

```text
results/dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv
```

That table collapses raw context rows into one row per gene/transcript/design
and shift direction, with recurrence class, pooling benefit/loss class, output
allocation context, and next review action. Current results collapse 95
reviewable context rows to 60 candidate-direction rows. No candidate is
currently claim-ready; many strong viral rows first need protocol/start-site or
pooling audit.

## Important Files To Read First

- `AGENTS.md`: local operating instructions for Codex or another coding agent.
- `md/dominant_cell_state_scientific_note.md`: scientific state, interpretation
  rules, current result counts, and next-model priorities.
- `md/dominant_rdg_atlas_browser_tutorial.md`: how to use the Shiny browser and
  review surfaces.
- `NEWS.md`: chronological implementation history.
- `scripts/run_dominant_cell_state_pipeline.R`: pipeline registry and entry
  point.
- `tests/integration/test_pipeline_invariants.R`: broad invariant checks for the
  pipeline, app startup, and key output tables.

## Standalone Migration

Target standalone path:

```text
/home/roler/Desktop/forks/dominant_cell_states
```

Current source path before migration:

```text
/home/roler/Desktop/forks/RiboCrypt/dominant_cell_states
```

Recommended GitHub setup:

```bash
cd /home/roler/Desktop/forks
git clone git@github.com:YOUR_USER/dominant_cell_states.git dominant_cell_states

rsync -a --info=progress2 \
  --exclude='.git/' \
  /home/roler/Desktop/forks/RiboCrypt/dominant_cell_states/ \
  /home/roler/Desktop/forks/dominant_cell_states/
```

This checkout is now structured as an R package plus a local analysis-results
tree:

- `R/`: package helper code for path discovery, result paths, RiboCrypt loading,
  pipeline startup, and `run_dcs_app()`.
- `tests/testthat/`: fast package tests that do not require Bio_data, FST pages,
  or a mounted S drive.
- `scripts/`: executable analysis and pipeline scripts.
- `inst/shiny/rdg_atlas_browser/`: packaged Shiny app files.
- `manual_translons/`: curated manual translon text input only.
- `results/`: ignored local result tree. This holds generated tables, FST page
  manifests, app-derived curated CSVs under `results/curated_inputs/`, caches,
  static coverage assets, and other artifacts.

Do not store FST pages, ORFik annotations, BAMs, reference indexes, or other
Bio_data artifacts in this repo. Those are external runtime inputs expected to
come from the user's Bio_data/S-drive setup.

Generated outputs have been moved under `results/`. Ignored compatibility
symlinks from earlier cleanup attempts should not be kept at the repo root; if
they appear, remove them with `dcs_remove_legacy_result_links()`.

Commit only source-like files, not generated outputs:

```bash
cd /home/roler/Desktop/forks/dominant_cell_states
git add .gitignore DESCRIPTION NAMESPACE R tests AGENTS.md README.md NEWS.md md
git add scripts inst tests man
git add manual_translons/translons_to_be_added_to_manual.txt
git add -n .
```

Inspect the dry run before committing. The directory is about 4.7 GB locally,
and normal GitHub should not receive the generated tables, caches, or static
coverage artifacts.

## RiboCrypt Dependency

This analysis still depends on local RiboCrypt code and should load RiboCrypt
with `devtools::load_all()`, not an installed package. From the old nested
layout, the pipeline discovers the RiboCrypt repo root automatically. From the
standalone layout, set a local RiboCrypt path explicitly for scripts or agent
work:

```bash
export RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt
```

The pipeline entry point now separates the RiboCrypt package root from the
analysis directory. From a standalone checkout it should discover the current
analysis directory and load RiboCrypt from `RIBOCRYPT_REPO`.

## Typical Commands

Always run R commands with `LC_ALL` unset. For pipeline runs, source the
pipeline runner and let it load local RiboCrypt after it has configured the
analysis runtime environment:

```bash
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt \
env -u LC_ALL R --vanilla -q -e 'source("scripts/run_dominant_cell_state_pipeline.R"); run_dominant_cell_state_pipeline(steps = "merged_replicate_translon_shift", dry_run = TRUE, rdg_plot_mode = "tables")'
```

For table/model iterations, prefer table mode:

```r
run_dominant_cell_state_pipeline(
  steps = "merged_replicate_translon_shift",
  force = TRUE,
  rdg_plot_mode = "tables"
)
```

For broad verification:

```bash
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt \
env -u LC_ALL R --vanilla -q -e 'source("scripts/run_dominant_cell_state_pipeline.R"); setup_runtime_env(find_analysis_dir()); devtools::load_all(Sys.getenv("RIBOCRYPT_REPO")); source("tests/integration/test_pipeline_invariants.R")'
```

For package-structure checks that do not touch Bio_data:

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::test()'
```

For `R CMD build`/`R CMD check`, use a clean checkout without local `results/`
or set `TMPDIR` to a filesystem with enough free space. `.Rbuildignore`
excludes results from the final package, but R still stages the local source
tree in a temporary build directory first.

If real libraries or references are needed and the S drive is mounted, the
usual ORFik experiment is:

```r
df <- ORFik::read.experiment("human_all_merged_l50")
```

## Current Review Priorities

1. Use `merged_replicate_translon_shift_candidate_queue.csv` for the
   clean-CDS/uORF redistribution review goal.
2. Start with `next_review_action`, `pooling_benefit_class`,
   `pooling_loss_class`, and `missing_evidence_flags`, not the score alone.
3. Treat protocol/start-site bias, direction conflicts, missing static
   coverage, missing RDG annotation, and single-run dominance as blockers for
   biological claims.
4. Use the Shiny browser Workbench and Inspection Assets for manual review.
5. Feed completed labels back through review-feedback and validation-label
   bridge steps before recalibrating motif/context filters.
