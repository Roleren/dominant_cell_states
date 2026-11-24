# ORFik accessors and transcript mapping: 2026-10-10

A subsequent same-day refactor moves wrapper preparation and output factories
into C++ and supersedes the remaining-cost discussion below. See the
[ORFik binding benchmark note](../../ORFik/benchmark/transcript_mapping/README.md).
Its baseline is the optimized implementation described in this report.

This audit follows the FST/coverage-table audit of 2026-10-09. It targets the
local ORFik source used by RiboCrypt and dominant_cell_states, particularly
compressed GRangesList accessors and paired transcript mapping. External
annotations and coverage caches were not changed.

## Measurements and reproducibility

The retained baseline consists of the original R functions in `baseline.qs`
and independent copies of both original native translation units under
`results/orfik_mapping_2026-10-10/baseline_native/`. The native copies are
compiled under different symbols: saving an old `.Call` wrapper alone would
otherwise benchmark the newly loaded native code on both sides.

Final comparisons use optimized `-O2` builds for both native implementations.
`R CMD INSTALL --preclean` builds the edited package, and benchmark loading
uses `devtools::load_all(..., compile = FALSE)` to preserve that build.
Earlier development measurements mixed sourceCpp's optimized build with
load_all's `-O0` debug build; those ratios are superseded. The debug comparison
is retained under `debug_build_audit/`, not used for speedup claims.

Run from the project root:

```bash
env -u LC_ALL R --vanilla -q -e 'source("scripts/benchmark_orfik_mapping_source.R")'
```

The runner loads local ORFik and RiboCrypt, checks exact outputs including
attributes before timing, alternates implementation order, and records five
paired repetitions with collection outside timed regions. It saves raw
seconds, medians, ranges, and session information under
`results/orfik_mapping_2026-10-10/source_verification/`. These local ignored
artifacts are required to reproduce the original baseline comparison.

## Accessors

The correct name slot for a compressed GRangesList is
`grl@partitioning@NAMES`. A separate 20,000-iteration expression benchmark
measured `names(grl)` at 2.035 microseconds versus 0.144 microseconds for the
slot, approximately **14.1-fold**, with no measured allocation for either.
This dispatch saving is small in absolute terms. ORFik already overrides
several accessors with slot-based implementations, so replacing every call
would add maintenance cost without a substantial pipeline gain.

The useful wrapper shortcut is obtaining the unlisted strand vector directly
from `grl@unlistData@strand`, instead of constructing and unlisting a strand
RleList. On 1,000 transcripts with 17 exons, the isolated operation fell from
about 1.54 ms to 0.13 ms (11.8-fold). The paired mapping wrapper now uses this
and obtains partition widths directly from partition end offsets.

For `strandPerGroup(..., FALSE)` and `seqnamesPerGroup(..., FALSE)`, expanding
an entire character vector wastes work on very exon-rich lists. Indexing the
Rle first helps there but slows smaller inputs. The implementation switches
only when the unlisted list has more than 65,536 elements and more than 32
exons per group on average. Small inputs retain the original path. The 33-exon crossover series measured
0.63 to 0.59 ms at 2,048 groups, 1.42 to 0.65 ms at 4,096 groups, and
3.06 to 0.78 ms at 8,192 groups. These synthetic accessor lists use one
chromosome and strand. Tests cover both paths, names, and empty groups. Direct-slot variants of groupings
and widthPerGroup were not adopted because measured gains were negligible.

## Native mapping changes

The C++ mapping code was not optimal. It copied input integer vectors into
std::vectors and allocated intermediate output buffers. To-transcript searches
were linear even for long ordered transcripts; the loop could search beyond
the paired transcript or beyond the available exon array. From-transcript
searches were unbounded and did not validate query coordinates.

The implementation now borrows integer inputs, holds direct integer pointers
through the hot loops, writes final integer outputs directly, checks dimensions and coordinate limits, and bounds searches to the
paired transcript. To-transcript mapping uses a linear search for at most eight
exons, binary search for longer ordered groups, and a bounded linear fallback
for irregular exon order. From-transcript mapping builds cumulative exon
widths and allocates only required output rows when empty exons are removed.
Long query loops check for user interrupts. Prefix arithmetic checks integer
limits and avoids intermediate coordinate overflow.

Correctness fixes also cover empty inputs, paired chromosome mismatches, and
scalar GRanges/IRanges references with three or more rows. Named native index
vectors preserve the original plain output schema. Invalid from-transcript
coordinates now produce an explicit error; unmappable to-transcript queries
retain the existing start=0/end=-1 sentinel. The public signatures are unchanged.

## Final matched-build timings

Medians in milliseconds per call. Speedup is original/edited; below 1 means
slower. Small differences in wrapper calls overlap the recorded timing ranges
and should not be interpreted as statistically established gains. Native
to-transcript inputs contain 100,000 queries to the last exon of one reference;
native from-transcript inputs contain 10,000 paired references and retain only
the terminal mapped exon. Both are search stress tests. Group accessors use
`keep.names = FALSE`.

| Input | Original ms | Edited ms | Speedup |
|---|---:|---:|---:|
| `to_real_exons_3` | 19.000 | 17.800 | 1.07× |
| `from_real_3_dropTRUE` | 18.000 | 19.000 | 0.95× |
| `to_real_exons_300` | 24.000 | 23.000 | 1.04× |
| `from_real_300_dropTRUE` | 19.000 | 19.000 | 1.00× |
| `to_real_exons_6000` | 95.000 | 93.000 | 1.02× |
| `from_real_6000_dropTRUE` | 27.000 | 27.000 | 1.00× |
| `to_native_100k_3exons_+` | 1.233 | 0.800 | 1.54× |
| `to_native_100k_3exons_-` | 1.267 | 0.800 | 1.58× |
| `from_native_10k_3exons_dropTRUE` | 0.233 | 0.267 | 0.88× |
| `to_native_100k_20exons_+` | 3.900 | 1.067 | 3.66× |
| `to_native_100k_20exons_-` | 3.933 | 1.033 | 3.81× |
| `from_native_10k_20exons_dropTRUE` | 1.667 | 1.000 | 1.67× |
| `to_native_100k_200exons_+` | 35.900 | 1.500 | 23.93× |
| `to_native_100k_200exons_-` | 30.067 | 1.467 | 20.50× |
| `from_native_10k_200exons_dropTRUE` | 38.333 | 6.333 | 6.05× |
| `strandPerGroup_1000x17exons_FALSE` | 0.140 | 0.140 | 1.00× |
| `seqnamesPerGroup_1000x17exons_FALSE` | 0.140 | 0.140 | 1.00× |
| `strandPerGroup_5000x200exons_FALSE` | 9.570 | 0.630 | 15.19× |
| `seqnamesPerGroup_5000x200exons_FALSE` | 9.680 | 0.630 | 15.37× |

Short from-transcript native calls cost approximately 33 microseconds more per
10,000 mappings (14% in this case). This retains boundary and coordinate
validation. The full from-transcript wrapper was broadly unchanged. Full
to-transcript wrapper medians improved 2–6%; native search gains do not imply
an equivalent full-pipeline improvement.

## Verification

- Existing ORFik genomic-range tests: 119 assertions passed.
- Combined targeted coverage, genomic-range, and mapping suites: 671 assertions
  passed before twelve additional accessor checks; the final mapping suite
  alone passed all 220 assertions against the installed optimized package.
- Six hundred randomized independent base-enumeration comparisons passed for
  both directions, both strands, and both removeEmpty settings.
- The same candidate kernels passed those comparisons under UndefinedBehaviorSanitizer
  with halt-on-error enabled. This is not an AddressSanitizer run.
- dominant_cell_states package tests: 228 assertions passed, including checked
  transcript mapping integration.

## Remaining costs and scope

Native last-exon workloads deliberately stress exon searching and are not a
prediction of whole-pipeline speedup. Full wrapper cases use actual cached
ATF4 (3 exons), FMR1 (17 exons), and BOD1L1 (26 exons) transcript structures,
repeated to 3, 300, and 6,000 transcripts. Much of wrapper time remains in
S4 construction, subsetting, validation, strand splitting and result assembly.
An initial checked implementation made short native calls slower. A separately
benchmarked direct-pointer loop recovered that cost before the final build.
Any remaining regressions remain visible in the final table.

A further useful design change would let repeated mapping operations reuse a
validated transcript catalogue and cumulative-width index. In particular,
from-transcript wrappers currently duplicate reference exons per query before
calling native code. Avoiding that requires an explicit indexed API and
careful strand/order semantics, rather than changing an accessor. No global
annotation cache or implicit cache invalidation policy was added here.
