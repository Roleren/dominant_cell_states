# Methods

## Study design

This study is a secondary analysis of public ribosome-profiling and matched
RNA-expression data. No new biological samples were generated. Candidate
discovery and focused validation were separated conceptually and in code. The
discovery layer screened grouped leader-uORF versus clean-CDS allocation
contrasts in a human Ribo-seq compendium. A candidate was not treated as
validated from its discovery score. Focused validation reconstructed exact
case and control contexts, collapsed sequencing runs to biological units,
formed a non-overlapping union of annotated uORFs, quantified a disjoint
clean-CDS region, and applied count, frame, shape, coordinate and
position-deletion gates. Raw length-preserving reprocessing of the tRNA-Glu
and matched tRNA-Arg studies, the total-RNA audit, feature-wide specificity
ranking and raw UV reprocessing were performed as separate falsification
analyses.

The public studies used in the focused analyses were:

| Perturbation or state | Accession | Cell line | Biological units |
|---|---|---|---|
| tRNA-Glu(UUC) overexpression | GSE77347; PRJNA310107 | MDA-MB-231 | 2 case, 2 control; three technical lanes per unit |
| tRNA-Arg(CCG) overexpression | GSE77317; PRJNA310096 | MDA-MB-231 | 2 case, 2 control |
| UV irradiation | GSE141459; PRJNA593553 | HeLa | 2 UV, 2 untreated |
| ABCE1 loss | GSE144163; PRJNA602917 | HCT116 | 2 targeting, 2 non-targeting control |
| Original MDA-LM2 state | PRJNA310015 | MDA-MB-231 derivatives | 2 MDA-LM2, 2 parental |
| CN34-LM1a state | PRJNA309995 | CN34 derivatives | 2 CN34-LM1a, 2 parental |
| Later MDA-LM2 state | GSE186639; PRJNA774895 | MDA-MB-231 derivatives | 2 MDA-LM2, 2 parental |

Sample identities, treatments and replicate labels were recovered from the
public study records and the local extended metadata table. The UV compendium
metadata had incorrectly labelled the cell context as HEK293/kidney; the GEO
record and raw sample accessions establish that the profiled samples were HeLa
cells. The public accession, rather than the erroneous compendium label, was
used throughout the focused analysis.

## Software and reference data

Processed compendium coverage was accessed through a local development checkout
of RiboCrypt version 1.15.2 and ORFik version 1.29.11. RiboCrypt was loaded with
`devtools::load_all()` rather than from an installed binary. Genomic and
transcript operations used Bioconductor classes and ORFik [12]. Analyses were
run with the R 4.4.0 development snapshot dated 18 March 2024 (revision 86148)
under Ubuntu 24.04.2 with `data.table` 1.18.99,
`GenomicRanges` 1.56.2, `GenomicFeatures` 1.55.4, `GenomicAlignments` 1.39.4
and `Rsamtools` 2.19.4. The raw tRNA and UV workflows used Skewer 0.2.2 [13],
STAR 2.7.4a [14] and samtools 1.19.2. The final submission will include a
machine-readable session-information file and immutable code release; the
versions above record the current working analysis and must be refreshed after
the final rerun.

Human reads were analysed against GRCh38 with the local Ensembl release 101
reference and transcript database. The focal transcript was
ENST00000339995, corresponding to RefSeq NM_001418. External FST coverage
pages, annotations, BAM files and genome indexes are runtime inputs and are not
copied into the source repository.

## EIF4G2 feature definition

Transcript positions are one-based. The focal transcript has a 308-nt leader,
uORF1 at positions 23–73, uORF2 at 128–220 and a main CDS at 309–3032. uORF1
is 51 nt long, begins with ATG and has a moderate Kozak context. uORF2 is 93 nt
long, begins with TTG and has a weak Kozak context. The main CDS is 2,724 nt and
begins with the known GUG initiation codon. The two uORFs do not overlap one
another or the main CDS.

For exact-union validation, all leader and overlapping uORF annotations were
converted to transcript positions and unioned before counting, preventing
double-counting of overlapping features. Any uORF-union positions were removed
from the annotated clean CDS. For EIF4G2, the uORF union was 144 nt and the
unique clean CDS was 2,724 nt. Individual-uORF results were quantified
separately; the mechanistically resolved claims were restricted to uORF1
because uORF2 had weaker phase support.

## Recovery of biological units

Exact case and control runs were recovered using study, author, cell line,
tissue, condition, perturbation, fraction, cancer type, sex and timepoint
fields, excluding the selected case–control axis and explicitly ignored
technical axes. Biological-unit identifiers were assigned conservatively. An
explicit replicate label was used first; otherwise the fallback order was SRA
Experiment, BioSample, sample name and Run. This order was required for the
tRNA-Glu study because the SRA submission assigned separate Experiment and
BioSample identifiers to sequencing lanes from the same labelled biological
replicate. Coverage was summed across runs within a biological unit before any
replicate-level result was calculated.

The tRNA-Glu experiment therefore contained two, not six, biological units per
arm. Lane-level values and lane subsets were used solely to diagnose technical
confounding and were never used as independent observations in statistical
tests or sample-size reporting.

## Ribosome-occupancy metrics

For each biological unit, counts were summed over the union-uORF and clean-CDS
positions. The primary within-transcript quantity was

`clean-CDS allocation = clean-CDS counts / (clean-CDS counts + union-uORF counts)`.

The contrast was case allocation minus control allocation after aggregating
counts within each arm. A negative contrast therefore denotes relative
redistribution from the clean CDS towards the uORF union. This is an occupancy
allocation; it is not interpreted as an initiation probability, protein
abundance or peptide-presentation rate.

Counts per million (CPM) were calculated separately for the uORF union and
clean CDS using the total coverage signal recorded for each sample. Group
ratios were calculated after averaging biological-unit CPM values, with a 0.5
pseudocount. Density ratios divided counts by feature length. Replicate-range
separation was the distance between the nearest case and control allocation
ranges in the observed direction; a positive value denotes complete
separation.

Reading phase within a uORF was assigned relative to its annotated first
nucleotide. When uORFs overlapped, phase was evaluated only at positions for
which all covering uORF annotations agreed. Phase-0 fraction was phase-0 counts
divided by the total phase-informative uORF counts. Frame-usage values for the
whole library and dominant read lengths were carried from the compendium QC
metadata as independent technical diagnostics.

## Exact-union validation gates

The focused validator required at least two biological units in each arm; at
least 30 union-uORF and 100 unique-clean-CDS counts in each aggregated arm; and
phase-0 fractions of at least 0.40 in both arms, with at least half of each
unit's uORF counts at phase-informative positions. The whole-library frame gate
required a mean frame-0 difference no larger than 15 percentage points and a
mean frame-0 value of at least 45% in both arms.

The broad-shape gate required no biological unit to place more than 35% of its
uORF counts at one position or more than 65% at its three largest positions,
and required at least 10% of uORF positions to be non-zero. Effects were
recomputed after shifting both uORF and clean-CDS coordinates by each integer
from −2 to +2 nt. A direction gate required the same sign at every shift and a
minimum absolute allocation contrast of 0.03.

Because a genuine short ORF may contain a translated pause, a failed shape gate
did not automatically erase an otherwise broad effect. Every uORF position was
deleted in turn, and the three highest pooled uORF positions were also deleted
together. A shape-flagged result could be retained as
`C_shape_flag_leaveout_robust` only when the biological-unit ranges were
separated, all other gates passed, every single-position deletion preserved the
direction with an absolute contrast of at least 0.03, and deletion of the
pooled top three positions also preserved the direction and magnitude of at
least 0.03. The tRNA-Glu EIF4G2 result met this class rather than being relabelled
as an ordinary shape pass.

## tRNA-Glu lane and read-length audit

For each of the 12 tRNA-Glu/control runs, we calculated the union-uORF and
clean-CDS counts, allocation, CPM, uORF phase and carried whole-library frame
and dominant-length metrics. Five subsets were declared from the observed
lane design: all lanes; all three lanes from the fully condition-matched first
replicate; lane 3 from both replicates, for which every sample was
30-nt-dominant; lanes 1 and 2 from replicate 1, for which every sample was
28-nt-dominant; and the partially confounded lanes 1 and 2 from replicate 2.
The subset audit asked whether the direction required the partial condition–
length linkage. It did not create additional biological replication.

## Raw tRNA preprocessing and per-length P-site calibration

The 12 tRNA-Glu/control and four tRNA-Arg/control ribosome-footprint FASTQs
were downloaded from ENA and checked against the repository-declared MD5
values and with `gzip -t` before processing. Skewer was run in tail mode
against `AGATCGGAAGAGCACACGTCT`, retaining 15--40-nt inserts. Reads were
aligned uniquely to the same GRCh38 STAR index used for the UV analysis, with
sorted BAM output, `--outFilterMultimapNmax 1`,
`--outFilterMismatchNoverLmax 0.1` and removal of non-canonical unannotated
intron motifs. Indexed BAMs were required to pass `samtools quickcheck` and
`idxstats`.

P-site offsets were inferred independently for each run and read length using
the 4,670 non-overlapping CDS-start windows described below. The physical-
offset and independent phase gates were identical to the raw UV analysis,
except that the accepted length range was 15--40 nt. Target P sites were then
mapped to EIF4G2 ENST00000339995. The study-wide category admitted a length
only when it passed in every run of a perturbation. Because the tRNA-Glu study
contains three deliberately distinct technical lanes, a complementary
lane-matched category admitted a length only when it passed in all four case
and control runs from the same lane. Every tRNA-Glu lane was required to
contribute a passing length before lanes were summed into the two true
biological units per arm. The all-length, 20--23-nt and 28--32-nt lane-matched
strata each had their own complete-lane requirement. The tRNA-Arg libraries
were treated as one lane with two units per arm.

The decisive raw tRNA gate required minimum counts and phase-0 dominance for
the exact 51-nt uORF1 in both arms, in addition to the union-uORF allocation
requirements. The uORF1-specific and union-uORF allocation contrasts had to
have the same direction, absolute uORF1-specific magnitude of at least 0.03,
and non-overlapping biological-unit ranges for both metrics. Passing categories
also required an absolute union contrast of at least 0.05, coordinate-shift
robustness and single- and greedy-three-position deletion robustness. These
criteria were fixed before the lane-3 case and tRNA-Arg target profiles were
examined.

## Matched tRNA specificity and codon analysis

To evaluate post-selection specificity, we quantified individual,
non-overlapping leader uORFs of at least five codons in both the tRNA-Glu and
tRNA-Arg studies. For each feature, uORF allocation was individual-uORF counts
divided by the sum of individual-uORF and clean-CDS counts. The tRNA-specific
effect was the tRNA-Glu allocation contrast minus the matched tRNA-Arg
contrast. Features were ranked at minimum uORF-count thresholds of 5, 10 and 20
in every arm of both experiments.

The outlier score subtracted the eligible-background median from the EIF4G2
specificity effect and divided by the unscaled median absolute deviation. It
was used only as a robust descriptive distance. It is not a parametric z-score.
The empirical top-rank fraction was 1 divided by the number of eligible
features and was not treated as an inferential P value.

Codons were counted in each uORF and clean CDS. Relative glutamate-codon
fraction was the difference between the fraction of GAA/GAG codons in the uORF
and clean CDS. The analogous CGG contrast was used for the tRNA-Arg control.
Associations with allocation effects were evaluated with two-sided Spearman
rank correlations in the strict 20-count set (`n = 18` features). No multiple-
testing adjustment was applied because the four reported correlations were
predefined mechanistic falsification tests rather than a discovery scan; all
four tests and exact P values are reported.

## Matched total-RNA and isoform audit

Source-provided Cufflinks FPKM tables for parallel ribosome-footprint and total-
transcript libraries were retrieved from GSE77347. EIF4G2 rows were identified
by RefSeq tracking identifier. For each isoform and assay, FPKM was averaged
across the two control and two tRNA-Glu biological replicates. We calculated
the log2 tRNA-Glu/control ratio, the fraction of summed EIF4G2 isoform-level
FPKM attributable to each transcript, and an RPF/total-transcript translation-
efficiency proxy. NM_001418 was prespecified as the focal isoform because it
maps to ENST00000339995. These processed short-read estimates were used to test
whether transcript abundance or major isoform replacement could plausibly
explain the 13.7-percentage-point occupancy-allocation shift; they were not used
as definitive isoform-resolved measurements.

## Raw UV preprocessing

Raw single-end FASTQ files for SRR11569080–SRR11569083 were downloaded from SRA
and verified with `gzip -t`. The source library contained four random 5′
nucleotides and six random nucleotides immediately before the constant 3′
adapter `CACTCGGGCACCAAGGA`. Skewer was run in tail mode against
`NNNNNNCACTCGGGCACCAAGGA`, initially retaining lengths 19–38 nt; the first
four nucleotides were then removed from sequence and quality strings, yielding
15–34-nt biological inserts. Trimmed FASTQs were gzip-validated.

Reads were aligned uniquely to the GRCh38 STAR index with sorted BAM output,
`--outFilterMultimapNmax 1`, `--outFilterMismatchNoverLmax 0.1` and removal of
non-canonical unannotated intron motifs. BAMs were indexed and passed
`samtools quickcheck`. Secondary, supplementary and unmapped alignments were
excluded in R.

## Per-run, per-length UV P-site calibration

P-site offsets were detected independently for every run and read length using
ORFik CDS-start metagene profiles. Calibration transcripts had leaders of at
least 40 nt and CDSs of at least 300 nt. A 220-nt window spanning 40 nt upstream
and 180 nt downstream of the first CDS nucleotide was formed for each eligible
transcript, and overlapping genomic windows were removed. The final calibration
set contained 4,670 non-overlapping CDS-start windows.

For each run, `ORFik::detectRibosomeShifts` was applied to 15–34-nt reads with
at least 100 reads per accepted length and at least 25 reads near the
translation-initiation site. Detector output was accepted only when the ORFik
start offset was non-positive and the implied P-site displacement was smaller
than the footprint. One impossible 16-nt, −18-nt result in control replicate 2
was rejected before EIF4G2 mapping.

An independent calibration gate then shifted five-prime positions by the
inferred offset and counted phase over the first 150 CDS nucleotides. A length
passed within a run when at least 100 CDS-start reads were available, phase-0
fraction was at least 0.40 and exceeded the larger alternative phase by at
least 0.05. Decisive shared lengths had to pass in all four runs. They were 17,
19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31 and 33 nt.

Calibrated P-sites were mapped to ENST00000339995 and assigned to uORF1,
uORF2, clean CDS or the remaining transcript. Results were summarized for each
exact length, all calibrated reads, all shared well-phased lengths, shared
20–23-nt short footprints and shared 28–32-nt canonical footprints. The three
shared summary categories were required to have exactly two units per arm,
minimum feature counts, non-overlapping replicate ranges, dominant uORF phase
0 in both arms, coordinate robustness and single- and greedy-three-position
deletion robustness. Position deletion was greedy for the raw UV analysis: at
each step, the position that minimized the absolute remaining allocation
contrast was removed, and this was repeated three times.

## Statistical analysis

Biological-unit sample sizes are stated for every contrast. No technical lane
was counted as a biological replicate. With two biological units per arm, the
focused contrast analyses were intentionally effect-size and falsification
based; no small-sample parametric P value is reported for the primary tRNA,
UV or ABCE1 allocation contrasts. Exact unit values, aggregate counts, range
separation, direction under coordinate perturbation and worst-case position
deletions are reported instead. Error bars in future figures must be defined in
their legends and all individual points must remain visible.

The only inferential tests in the completed focused analysis are the two-sided
Spearman correlations in the matched codon analysis. No observations were
excluded after examining an outcome except through the explicit offset,
feature-count and phase gates described above. The later MDA-LM2 dataset was
retained as a non-replication and classified as a frame-imbalance risk rather
than omitted.

## Reproducibility and code review

Exact-union helpers have unit tests for biological-unit collapse, ignored
technical axes, unambiguous phase assignment and feature metrics. Raw-length
helpers have unit tests for length strata, physical offset validation,
calibration phase gates and worst-position deletion. At the current checkpoint,
the package test suite passes 66 tests, and the broader integration suite passes
script parsing, application startup, pipeline-registry, output-schema and model
invariants. Source files are checked for whitespace errors and the UV shell
workflow passes `bash -n`.

## Data availability

All primary sequencing data analysed here are public under GSE77347/
PRJNA310107, GSE77317/PRJNA310096, GSE141459/PRJNA593553,
GSE144163/PRJNA602917, PRJNA310015, PRJNA309995 and GSE186639/PRJNA774895.
The final paper will deposit the minimum machine-readable evidence tables,
figure-source data, exact sample manifest, checksums and frozen environment
metadata in a DOI-minting repository. **Repository DOI: [required before
submission].** Large public FASTQs and reference indexes will not be
redistributed.

## Code availability

The working source is in the standalone `dominant_cell_states` repository and
depends on a local development checkout of RiboCrypt. The final paper will
release the exact source commit, tagged RiboCrypt dependency, tests, command
lines and a DOI-minted archive sufficient to rebuild every reported table and
figure. **Public repository URL, commit and archive DOI: [required before
submission].** The current code must be shared with editors and referees on
request even before the archival release.

## Use of generative AI

OpenAI Codex was used during exploratory code review, analysis implementation,
quality-control planning and preparation of an initial manuscript draft. The
authors inspected the underlying code, output tables, primary literature and
all resulting text, and remain responsible for the analyses, interpretation
and final manuscript. This statement is a working disclosure and should be
updated to match the exact tools and the journal policy at submission.

## Methods references

Methods references 12–14 are listed with the main references. Add formal
dataset citations and any final environment/container references when the
archive is deposited.
