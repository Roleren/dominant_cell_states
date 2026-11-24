# ORFik coverage-function benchmarks and improvements

Benchmarked on 2026-10-09 using the local ORFik checkout at
`~/livemount/forks/ORFik` (starting commit `45fda30`, version 1.29.15), local
RiboCrypt loaded with `devtools::load_all()`, and real dominant-cell-state
coverage inputs. Raw timings, profiles, input manifests, original function
snapshots, validation logs, and the package diff are under
`results/orfik_performance_2026-10-09/`. Generated artifacts remain ignored.

## What changed

`coverageByTranscriptFST()` now makes one chromosome-aware page/exon interval
join for the input batch. It shares reads of overlapping or nearby slices from
the same page within that call, then extracts independent output tables in the
original transcript/exon order. Numeric offsets stay numeric, eliminating the
old transpose-to-character and format/parse conversions. Requested column order
and the existing per-page minus-strand reversal are preserved.

Coalesced reads are limited to 8,192 rows and eight times the number of requested
rows in a block. An individual exon longer than 8,192 rows retains its original
read size. This reduces repeated decoding while bounding additional transient
data. It deliberately reads some gap positions and discards them; the output
contains only the original requested positions. The function retains no cache
between calls, so new indexes/pages are visible immediately. The pipeline's
existing index-timestamp cache invalidation remains applicable.

`coverage_to_dt()` now expands positive Rle runs directly into count, position,
and gene vectors. It no longer constructs a separate one-base IRanges entry
for every retained nucleotide. Integer/fractional values, factor levels,
positions, frames, fractions, and empty-profile warnings are preserved.
`coveragePerTiling()` and its downstream consumers inherit this improvement
when they request `drop.zero.dt = TRUE`. Its default dense output path remains
unchanged.

The batched chromosome key also prevents cross-chromosome page matches for a
multi-chromosome list element. Empty list elements now return empty tables.
These edge cases have explicit tests; production gene inputs each use one
chromosome and reproduce their existing coverage exactly.

The changes are installed in `/usr/local/lib/R/site-library/ORFik` as 1.29.15.
Fresh R sessions use them; existing sessions require an ORFik reload or restart.
No version bump or scientific-output regeneration was needed.

## Inputs and measurement design

| Gene | Transcript | Exon slices | Display nucleotides | Available libraries |
| --- | --- | ---: | ---: | ---: |
| ATF4 | ENST00000674920 | 3 | 1,338 | 3,857 |
| FMR1 | ENST00000370475 | 17 | 2,160 | 3,857 |
| BOD1L1 | ENST00000040738 | 26 | 9,294 | 3,857 |

These are the actual leader/CDS display coordinates stored by this pipeline,
not whole-transcript lengths. FMR1 uses the `.75` transcript. The index is
`~/livemount/Bio_data/processed_data/all_samples-Homo_sapiens/collection_tables_indexed/coverage_index.fst`,
with the unchanged 2026-02-04 timestamp recorded in `real_inputs.csv`.

Each transcript was tested with four actual library columns and all 3,857
columns. Real library-summed profiles supply the table-conversion inputs.
Additional tests use the first 30 distinct cached genes in sorted filename
order (manifest in `source_verification/panel_inputs.csv`) and seeded synthetic
100-transcript, 5,000-nucleotide profiles with 1% and 50% nonzero positions.
Variable-length dense profiles were also tested as an alternative optimization.

Candidates were first evaluated in memory against saved unmodified functions;
source edits followed demonstrated speedups and output-equivalence checks.
A fresh-session pass then benchmarked the edited package against those original
functions. The installed original implementations were separately confirmed to
match the original checkout by comparing their deparsed bodies before installing
the changes.

Timing order alternates original/candidate within each pair. Each case has
5–9 paired repetitions, with multiple calls per repetition for small operations.
Every case warms both implementations and verifies output equivalence before
timing. Garbage collection and Rprof runs happen outside the measured calls.
Reported numbers are medians per call; raw repetitions and min/max values are
retained. Final summaries distinguish the prototype pass from the edited-source
pass. These are warm operating-system-cache reads, **not cold-disk timings**.

The shared server has 48 logical CPUs (two Xeon Gold 6126 processors) and about
502 GiB RAM. data.table uses four threads; the existing FST reader policy of five
threads is preserved. Other server workloads and benchmark processes were
active, so tiny changes are not treated as convincing improvements. Rprof uses
the platform's 10 ms sampling interval and supplies qualitative attribution,
not microsecond measurements.

## Final timings

<!-- BEGIN GENERATED TIMINGS -->
Edited-source results; seconds per call. The four-library batch uses the
same four actual columns for both implementations.

| Function | Input | Before (s) | After (s) | Speedup | Paired repetitions |
| --- | --- | ---: | ---: | ---: | ---: |
| coverageByTranscriptFST | ATF4_4_libraries | 0.02200 | 0.01467 | 1.50× | 7 |
| coverageByTranscriptFST | ATF4_all_libraries | 0.39500 | 0.16000 | 2.47× | 5 |
| coverage_to_dt | ATF4_summed_profile | 0.01380 | 0.00980 | 1.41× | 7 |
| coverageByTranscriptFST | FMR1_4_libraries | 0.03533 | 0.02300 | 1.54× | 7 |
| coverageByTranscriptFST | FMR1_all_libraries | 2.18800 | 1.34900 | 1.62× | 5 |
| coverage_to_dt | FMR1_summed_profile | 0.01430 | 0.01010 | 1.42× | 7 |
| coverageByTranscriptFST | BOD1L1_4_libraries | 0.05033 | 0.03100 | 1.62× | 7 |
| coverageByTranscriptFST | BOD1L1_all_libraries | 3.58800 | 2.59900 | 1.38× | 5 |
| coverage_to_dt | BOD1L1_summed_profile | 0.01830 | 0.01315 | 1.39× | 7 |
| coverageByTranscriptFST | 30_gene_panel_4_libraries | 0.65400 | 0.14700 | 4.45× | 7 |
| coverage_to_dt | 100x5000_density_0.01 | 0.01900 | 0.01400 | 1.36× | 7 |
| coverage_to_dt | 100x5000_density_0.5 | 0.14900 | 0.12400 | 1.20× | 7 |

`final_paired_summary.csv` additionally retains observed ranges, the number
of faster paired trials, and a paired bootstrap interval for the median
speedup (5,000 resamples, fixed seed). These intervals describe this shared
server run; they do not capture uncertainty across machines or cold disks.
<!-- END GENERATED TIMINGS -->

## Why the FST reader improved

| Input | Original page reads | Improved page reads | Original rows decoded into tables | Improved rows decoded into tables | Largest improved read |
| --- | ---: | ---: | ---: | ---: | ---: |
| ATF4 | 3 | 1 | 1,338 | 2,041 | 2,041 |
| FMR1 | 17 | 8 | 2,160 | 9,294 | 4,442 |
| BOD1L1 | 26 | 14 | 9,294 | 21,414 | 7,981 |

`fst_read_audit.csv` records these counts and independently compares the returned
values with the saved production tracks. FMR1's original all-library profile
spends about 98% of sampled time in FST retrieval. Avoiding repeated decoding
there matters much more than only accelerating the interval join.

The largest combined BOD1L1 read is 7,981 rows. Across all 3,857 libraries its
values occupy approximately 117 MiB for integer storage or 235 MiB for double
storage, before table/container overhead.
This is a transient data estimate, not a measured process-memory peak. High
worker counts must budget for these buffers and the requested output tables.

## Validation

- ORFik's coverage suite: 341 passing assertions, no failures/warnings/skips.
  Subsequent focused edge-case suite: 104 passing assertions.
- RiboCrypt coverage-download/QC suites: 29 passing assertions.
- Full dominant-cell-state package suite: 228 passing assertions.
- Candidate validation: 800 seeded combinations of sparse/dense profiles,
  integer/fractional counts, gene-name policies, and frame columns, plus
  multi-page/strand/column FST fixtures. Final checks use zero numerical tolerance.
- Edited-source and installed-package checks compare ATF4, FMR1, and BOD1L1
  directly with their existing production tracks. Every timed case also checks
  the original and improved outputs, including attributes.
- Regression tests cover page boundaries, unsorted/duplicate exons, both
  strands, independently mutable result tables, moved page paths, invalid
  inputs, updated page contents, empty elements, chromosome matching, coverage
  positions, gene factors, fractions, frames, and empty-profile warnings.

Production analysis tables and reference/coverage inputs were not rewritten.
The new regression tests live in ORFik's
`tests/testthat/test_coverage_hotpaths.R`; the implementation changes are in
`R/grl_helpers.R` and `R/coverage_helpers.R`.

## Other functions and rejected changes

`coverageByTranscriptC()`, `coverageByTranscriptSum()`, `pmapToTranscriptF()`,
and the complete `coveragePerTiling()` calls on covRle and weighted GRanges were
profiled using the same real coordinates, with three and 300 transcript entries.
Full measurements are in `screen_timings.csv` and the `screen_*.out.csv` profiles.

Applying the existing chromosome-Views extraction path to smaller
`coverageByTranscriptC()` batches improved medians by only about 2% at three
transcripts and 300; at 6,000 transcripts it was effectively unchanged.
That experiment was rejected. Its paired measurements are saved under
`coverage_extraction/`.

Replacing dense-output data.table grouping with direct sequence construction
improved the two variable-window stress cases by approximately 1.5% and 0%.
That experiment was also rejected. The initial measurements are retained under
`coverage_tables/`; the production dense path stays unchanged.

## Remaining expensive work

All-library FST retrieval remains the largest cost in these transcript reads.
The bounded strategy trades additional gap rows for fewer reads; highly isolated
exons have less opportunity to benefit. A page-format change or selective
compressed-block reuse would require FST-level investigation and a separate
memory/I/O benchmark.

For weighted-read `coveragePerTiling()`, constructing genomic coverage and
dispatching/extracting S4 range objects remain substantial costs. Reusing covRle
within a caller is already beneficial, but avoiding a full library import or
global coverage build requires a caller-level design and validation against
real regional sums. The large cold manual-candidate activity stage is not solved
by these FST/table changes.

Warm pipeline runs often bypass these functions through valid per-transcript
caches. The measured function speedups therefore do not imply the same speedup
for the full pipeline, or for independent model/bootstrap stages. No full cold
pipeline runtime claim is made here.

## Reproduction

The existing artifact directory contains the original pre-edit function
snapshot. To rerun the final comparison from this repository:

```bash
env -u LC_ALL R --vanilla -q -e 'source("scripts/benchmarks/benchmark_orfik_source_verification.R")'
env -u LC_ALL R --vanilla -q -e 'source("scripts/benchmarks/validate_orfik_hotpath_candidates.R")'
env -u LC_ALL R --vanilla -q -e 'source("scripts/benchmarks/audit_orfik_fst_reads.R")'
env -u LC_ALL R --vanilla -q -e 'source("scripts/benchmarks/profile_orfik_consumers.R")'
```

For a new before/after experiment, set `ORFIK_BENCHMARK_CAPTURE_BASELINE=1`
while the checkout is still unmodified, then run
`scripts/benchmarks/benchmark_orfik_hotpaths.R`. The scripts refuse to silently capture an
edited package as a missing baseline. `ORFIK_BENCHMARK_DIR` selects an isolated
artifact directory. Keep the original `baseline_functions.qs` when comparing
subsequent edits. Real-input runs require the project's existing transcript
coverage caches and external megafst pages; the ORFik regression tests use only
small synthetic fixtures.
