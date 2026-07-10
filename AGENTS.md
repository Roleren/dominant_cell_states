# Local Agent Notes

This directory is being split out of
`/home/roler/Desktop/forks/RiboCrypt/dominant_cell_states` into a standalone
repo at `/home/roler/Desktop/forks/dominant_cell_states`.

## First Steps

1. Read `README.md`, `md/dominant_cell_state_scientific_note.md`,
   `md/dominant_rdg_atlas_browser_tutorial.md`, and `NEWS.md` before making
   changes.
2. Catch up on local changes and generated outputs before editing.
3. Always inspect related code and output schemas before making assumptions.
   If needed, inspect ORFik/Bioconductor-level behavior before patching.
4. Prefer correctness over speed. If unsure, make a focused test.

## RiboCrypt Loading

This analysis depends on local RiboCrypt code. Always use
`devtools::load_all()` and do not rely on an installed RiboCrypt version.

When working from the standalone repo, use:

```bash
export RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt
env -u LC_ALL R --vanilla -q -e 'devtools::load_all(Sys.getenv("RIBOCRYPT_REPO")); ...'
```

For pipeline runs, prefer sourcing the pipeline runner first and letting
`run_dominant_cell_state_pipeline()` load RiboCrypt after runtime setup:

```bash
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt \
env -u LC_ALL R --vanilla -q -e 'source("scripts/run_dominant_cell_state_pipeline.R"); run_dominant_cell_state_pipeline(steps = "merged_replicate_translon_shift", dry_run = TRUE, rdg_plot_mode = "tables")'
```

When working from the old nested RiboCrypt repo, this is also acceptable:

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); ...'
```

Use `env -u LC_ALL R --vanilla -q -e 'devtools::load_all(...); ...'` to avoid
`LC_ALL='C.UTF-8'` warnings from `testthat`/`withr`.

## Data Access

FST pages, ORFik annotations, BAMs, reference indexes, and other Bio_data
artifacts are external runtime inputs. Do not copy them into this repository or
try to package them. The repo should contain code, docs, tests, and small
curated inputs only.

If real libraries or references are needed:

```r
df <- ORFik::read.experiment("human_all_merged_l50")
```

If the S drive is not mounted:

```r
if (!dir.exists("/media/roler/S")) {
  message("Mounting the S drive!")
  system("udisksctl mount -b /dev/sda2")
}
```

## Git And Output Policy

The local analysis folder is large and contains generated tables, caches, and
static coverage assets. Do not blindly `git add .`.

Commit source-like files first:

- `README.md`, `AGENTS.md`, `.gitignore`, `NEWS.md`, `DESCRIPTION`,
  `NAMESPACE`, `md/`
- `R/`
- `scripts/`
- `inst/shiny/`
- `tests/testthat/`
- `tests/`
- manual translon text input files

Generated output directories, caches, huge CSVs, FST page manifests, static
coverage assets, app-derived curated CSVs under `results/curated_inputs/`, and
readlength/BAM inputs live under ignored `results/` or external Bio_data
locations. Keep them local unless the user explicitly asks for an artifact
strategy. Do not keep root-level compatibility symlinks; remove any that
appear.

## Current Algorithmic State

The current main discovery target is merged-replicate clean-CDS/uORF
redistribution. The first table to inspect is:

```text
dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv
```

The canonical local path is now:

```text
results/dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv
```

Interpret it as a claim gate:

- `next_review_action` says the next job.
- `pooling_benefit_class` says what pooling gained.
- `pooling_loss_class` says what pooling may have hidden.
- `missing_evidence_flags` says why a candidate is not claim-ready.

Current state: 95 raw review rows collapse to 60 candidate-direction rows. No
candidate is currently claim-ready. Viral `PPP1R15A`, `DDIT3`, `ATF4`, `TFE3`,
`IFIH1`, and `FBXL4` are important review leads, but most currently require
protocol/start-site/pooling audit before interpretation.

## Verification

For a targeted rerun:

```r
source("scripts/run_dominant_cell_state_pipeline.R")
run_dominant_cell_state_pipeline(
  steps = "merged_replicate_translon_shift",
  force = TRUE,
  rdg_plot_mode = "tables"
)
```

For broad verification:

```r
setup_runtime_env(find_analysis_dir())
devtools::load_all(Sys.getenv("RIBOCRYPT_REPO"))
source("tests/integration/test_pipeline_invariants.R")
```

For fast package-structure tests:

```r
devtools::test()
```

The invariant suite parses scripts, checks Shiny app startup, verifies pipeline
registry/output assumptions, and checks key model output tables.
