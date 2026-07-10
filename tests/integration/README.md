# Integration Checks

Fast package tests live in `tests/testthat/` and run with:

```sh
env -u LC_ALL R --vanilla -q -e 'devtools::test()'
```

The scripts in this directory are longer local verification and diagnostic
checks for the current generated result tree. Run them from the standalone repo
with `RIBOCRYPT_REPO` pointing at the local RiboCrypt checkout:

```sh
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt env -u LC_ALL R --vanilla -q -e 'source("tests/integration/test_pipeline_invariants.R")'
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt env -u LC_ALL R --vanilla -q -e 'source("tests/integration/timing_pipeline_cache.R")'
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt env -u LC_ALL R --vanilla -q -e 'source("tests/integration/result_sanity_review.R")'
```

`test_pipeline_invariants.R` checks the current analysis outputs and catches
human-review failure modes: unsorted candidate ranks, missing atlas priority
scores, broken grouped-branch denominators, missing manual translon output, and
missing count-aware branch-allocation outputs.

`timing_pipeline_cache.R` does not run the heavy analysis. It times cache-status
checks for the current pipeline outputs and writes
`results/test_reports/pipeline_cache_status.csv`.

`result_sanity_review.R` reads generated CSVs, inventories figures/HTML, checks
the main atlas/review tables, and writes compact Markdown/CSV sanity reports
under `results/test_reports/` without failing on expected empty products.
