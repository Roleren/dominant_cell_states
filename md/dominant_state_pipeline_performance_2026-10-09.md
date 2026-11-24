# Pipeline performance and incremental caching

Date: 2026-10-09. Measurements use the real human compendium on the shared
analysis server, with local RiboCrypt loaded through `devtools::load_all()`.
Generated logs, timings, validation snapshots, and isolated benchmark caches
are saved under `results/performance_2026-10-09/`.

## Baseline and plan

The completed October 8 FMR1 coordinate correction took 3,216 seconds across
50 core/RDG steps. Its manual-translon step alone took 1,389 seconds (43% of
step time), followed by clean CDS (210), candidate discovery (174), FST page
selection (166), RDG graphs (148), graph geometry (113), motif modeling (107),
branch covariates (104), and uORF regulation (92). State scores were unchanged
by that coordinate correction. The benchmark preserves that actual run summary;
it is a historical baseline, not a new controlled cold run.

Work prioritized avoiding repeated work without changing statistical methods:
cache fixed discovered-manual-candidate activity separately; cache per-gene
selection and measurements; correct dependency signatures and invalidation;
retire successive candidate-ranking exports; benchmark real computations and
verify output equivalence. Statistical layers that answer different questions
remain separate.

## Implemented behavior

- Discovered SNAI2/MMP2 manual-candidate activity is cached independently of the
  fixed-coordinate queue. Its key includes candidate/mRNA ranges, resolved
  pshifted coverage paths and timestamps, unique-mapper policy, scoring code,
  and ORFik version. Fixed user-coordinate edits do not rescore these six ORFs.
- Transcript selection, clean-CDS vectors, uORF feature tables, and RDG
  sequence annotations have
  per-gene keys. Inputs include relevant ranges/translons, canonical selection
  and exon-QC results where applicable, library sizes, index signature, and
  calculation/helper code. A changed range invalidates that gene; a changed
  library collection invalidates every affected gene. Unrelated plotting or
  modeling code is excluded from calculation signatures.
- Raw transcript coverage already had a per-transcript cache. It is retained
  and prevalidated against the unchanged megafst index as authorized by the
  user. The index timestamp/size/path is the collection update marker. The
  user rebuilds the index when replacing megafst and does not independently
  replace single chromosome pages; there is therefore no per-page scan.
- The runner now tracks the megafst index, experiment CSV, total library-size
  vector, annotations, merged coverage, and applicable CDS count-table inputs.
  External inputs use timestamp/size signatures. The CDS count-table fallback is currently disabled in clean-CDS scoring;
  its input is watched conservatively so enabling it cannot reuse stale results.
  A count-table update invalidates that step even if the index is unchanged. Missing-to-present input changes also invalidate caches.
- Project tables use memoized content hashes (up to 256 MiB) so rewriting an
  identical table does not cascade into recomputation. Larger tables and external binary inputs use
  file signatures. Output edits invalidate the producing step. Cache misses
  identify changed source/dependency files.
- Helper dependencies follow actual local `source()` calls transitively,
  including filenames assigned to variables. Comments and file-existence
  probes are excluded; these previously caused unrelated scripts to enter
  every step's dependency graph. The
  HTML-widget helper remains shared because the runner sources it for every
  step. The cache audit uses the same helper discovery and output-mode rules.
- State-only consumers no longer depend on unused per-gene diagnostics/QC.
  Clean-CDS scoring no longer depends on the FST page-discovery summary
  (which contains a generation timestamp); manual-range consumers track the
  ranges RDS rather than unused candidate manifests. uORF inference tracks
  its diagnostic input rather than the unused full gene-expression matrix.
  Models using the whole gene-expression matrix still legitimately depend on
  that matrix. Global scores, joins, multiple-testing corrections, and cohort
  models are not treated as independent per-gene computations.
- Table-mode cache checks exclude suppressed figures. A conditional IFIH1 lung
  plot is no longer a universally required output.
- Candidate ranking now has one table, `candidate_rank` and `candidate_score`.
  Intermediate tables/copies/rank aliases were removed, preserving additive
  score weights and tie ordering. All active consumers use the current schema.
  [Ranking history](history/candidate_ranking_history.md) records the retired
  exports. Dated server-review narrative and the superseded performance plan were moved to explicit history notes. Two unused legacy workflows were archived under `md/history/legacy_scripts/`.

## Benchmarks and validation

The incremental benchmark calls the real uORF measurement function using an
isolated cache. Existing raw coverage caches remain available. These are cold
**derived-measurement** timings, not cold disk/megafst timings.

| Scenario | Genes requested | Genes computed | Seconds |
| --- | ---: | ---: | ---: |
| Initial panel, withholding FMR1 | 69 | 69 | 71.452 |
| Add FMR1 | 70 | 1 | 4.331 |
| Repeat unchanged panel | 70 | 0 | 3.466 |
| Edit one FMR1 range | 70 | 1 | 3.645 |
| Rebuilt index timestamp | 70 | 70 | 71.401 |

Every unchanged-input measurement is checked against the production per-gene
result. Range and index changes are simulated in memory; production annotations,
coverage, and curated inputs are untouched. The executable benchmark is
`scripts/benchmarks/benchmark_dominant_incremental_measurements.R`; exact timings are in
`incremental_gene_benchmark.csv` under the artifact directory.

A cold discovered-manual-candidate activity pass took 2,190.776 seconds on this
server, versus 4.540 seconds for a forced repeat with the activity cache and
0.494 seconds for the normal outer-cache skip. The candidate CSV remained
byte-identical. The cold pass reached about 53 GiB RSS, making memory use a
material constraint on parallelizing this computation. The historical October 8
manual step took 1,389 seconds; shared-server contention prevents treating the
cold-time difference as a regression measurement.

The corrected 50-step core/RDG rerun completed with 47 successful executions
and three unchanged steps skipped, totaling 2,134.574 seconds of step time
(35.6 minutes). This pass includes the mapping correction and a refreshed
derived cache, so it is not a controlled cache-only comparison with October 8.
Two fully warm runs skipped all 50 steps in 7.035 and 8.178 seconds of wall time, including
runner setup and cache checks within the existing R process. The final dependency
refresh also completed without errors. Producer comparisons below retain raw
coverage caches and measure cold versus warm derived objects, not cold disk I/O. The FST rebuild used the final guarded vectorized mapper;
its earlier run was 320.997 seconds before that function signature refreshed.

| Producer | Derived-cache rebuild (seconds) | Forced cached repeat (seconds) |
| --- | ---: | ---: |
| Discovered manual activity | 2190.776 | 4.540 |
| FST transcript/page selection | 258.754 | 35.851 |
| Clean-CDS scoring | 421.832 | 45.390 |
| uORF regulation | 175.622 | 69.903 |
| RDG graphs with cached annotation | 165.268 | 94.125 |

Forced repeats still load annotations, categorize translons, assemble cohort
outputs, and write tables; the outer cache skips these costs completely. Invalid
or nonfinite measurement results are deliberately not cached. These are sequential
real-data measurements on the shared server, not a dedicated-machine benchmark.

Initial cache-only comparisons reproduced the old measurements and
ranking exactly; a subsequent RDG annotation check exposed a pre-existing
coordinate-mapping problem, so those old results are not used as the sole oracle.

Local ORFik native mapping uses `counter <= txSize` while indexing an exon
array of length `txSize`. An unmappable input can exhaust the search and read
past that array. Repeated annotations of the same EEF2/EIF4B features produced
negative two-base sentinel coordinates or million-base genomic coordinates in
a transcript-coordinate table. The project now checks that every query exon
lies within its paired transcript before invoking the native mapper, then checks
mapped coordinate bounds. Unmappable features are excluded deterministically;
valid pairs use the vectorized mapping path. Focused tests cover both strands,
spliced ranges, unrelated transcript pairs, zero-width and absent mappings, and
batch/scalar equivalence. Final cached measurements are compared to a separate
uncached run using these same guards. Corrections to malformed old features are
reported explicitly rather than claimed as unchanged biology.

The correction removes 32 malformed non-CDS features, all with zero observed
counts. Eighteen genes consequently lack eligible ORFs and leave the ORF-specific
analysis panel: 88 genes/329 features become 70 genes/279 features. Those genes
remain in the main 129-gene clean-CDS analysis. All retained feature measurements,
including every BOD1L1 and FMR1 measurement, are identical to their prior values.
The main state-score, clean-CDS expression, gene-diagnostic, and manual-candidate
CSVs are byte-identical. Corrected cached feature measurements also match the
independent uncached run byte-for-byte. Current candidate scores reproduce the
unchanged formula to within 9e-15; ranks can change because eligible inputs and
population summaries changed. See `mapping_correction_summary.json`,
`mapping_removed_features.csv`, and `output_equivalence.csv` in the benchmark
directory. Table-mode validation preserves existing static figures; figures for
affected features should be regenerated before using them for interpretation.

Helpers also load into each step's environment, and computation signatures
inspect function calls without evaluating unrelated global variables. Signatures containing
exon-QC `data.table`s compare serialized values, avoiding false misses from
transient self-reference pointers. Biologically empty uORF results are cached
too; genuine coverage failures raise errors and cannot be cached.

The package suite passes 149 assertions with no failures, warnings, or skips,
including same-size index/count-table timestamp updates, output edits, identical
project-table rewrites, transitive helper dependencies, and one-gene changes.
The broad invariant suite passes script/app parsing, Shiny startup, and registry
checks, then stops at a pre-existing missing server artifact:
`results/dominant_rdg_sdrive_readiness/sdrive_experiment_status.csv`.
This suite is not reported as passing.

## Thread benchmark

The same design-aware score-semantics step ran in isolated output directories
with unchanged real inputs. All eight scientific CSVs were byte-identical to
production in each pass.

| data.table threads | Seconds |
| ---: | ---: |
| 1 | 90.384 |
| 4 | 51.391 |
| 24 (previous server default) | 95.630 |

Four threads were 1.86 times faster than 24 in this comparison. The runner now
defaults to at most four threads, preserving an already smaller thread setting.
`DOMINANT_DT_THREADS` accepts a positive integer override. This is a sequential
single-pass comparison on a shared 48-CPU server; load varied substantially, so
it is directional evidence rather than a dedicated-machine scaling study.
`data_table_thread_benchmark.csv` and `thread_benchmark.R` preserve the test.

## Remaining expensive work and next priorities

A genuinely cold manual-activity computation still imports the large all-merged
pshifted library and computes regional coverage repeatedly. Its cache removes
this expense from fixed-coordinate edits, but a changed merged library must
recompute it. A future optimization should use indexed regional access or
restrict reads to candidate/flanking regions, with exact score equivalence
against the current ORFik implementation.

Even warm per-gene caches leave compendium annotation loading, translon categorization,
joins, cohort aggregation, and output serialization. Large feature and metadata
CSV tables still have to be read by consumers. Binary intermediate tables and
shared process-level input loading could reduce this, but require coordinated
schema/reader changes rather than an unmeasured format substitution.

The historical September metadata benchmark identified a four-way interaction
scan as expensive. This change does not cap its interaction order, prune tests,
or change statistical coverage. The metadata scope was not included in the
October 8 baseline. Benchmark that scope separately before deciding whether to
parallelize or change the search space. Similarly, hierarchical leave-one-out,
transport, graph-geometry, and motif models are real remaining computations;
parallel execution needs deterministic-output and memory-budget validation on
this shared server. Candidate discovery still runs matched-control bootstraps
and applies population-wide score normalization and multiplicity corrections.
Caching those per gene would require careful separation of reusable fits from
global adjustments; the current bootstrap RNG stream also depends on traversal
order. A stable per-gene seed scheme would need an explicit statistical
validation, rather than silently changing existing resampling results.

The installed server ORFik lacks `cdsExonCoverageQC`; the existing clean-CDS
code warns and emits empty exon QC on this server. This is a pre-existing
scientific-validation limitation, not repaired by caching. Production output
comparisons must explicitly check that caching has not changed this behavior.

## Reproduction and cache controls

Run the normal producer first, then:

```bash
RIBOCRYPT_REPO=/home/rstudio/livemount/forks/RiboCrypt \
env -u LC_ALL R --vanilla -q -f scripts/benchmarks/benchmark_dominant_incremental_measurements.R
```

`force = TRUE` reruns a pipeline step while retaining valid measurement caches.
`DOMINANT_MEASUREMENT_CACHE=off` bypasses the new derived-measurement caches;
existing raw transcript-coverage caches still follow the index signature.
A new megafst index causes all raw transcript keys and corresponding derived
keys to miss. Rerun from `fst_pages`/`clean_cds` after rebuilding coverage, with
updated experiment metadata and total library sizes. Normal dependency checks
then refresh downstream results. The benchmark directory includes the exact
R scripts for the forced cold/warm passes and final validation.

## Serialization follow-up

The subsequent full qs2 migration and verification are documented in [qs2_migration_2026-10-09.md](qs2_migration_2026-10-09.md). All five forced cache producers succeeded, and all 50 core/RDG steps were skipped in the final 6.05-second warm run. Serialization improvements do not remove the remaining inference and table/model construction costs.
