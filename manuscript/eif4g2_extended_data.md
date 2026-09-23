# Extended Data plan and draft legends

Nature typically permits up to ten Extended Data display items. These items
should carry the exhaustive falsification record so the four main figures can
remain readable. Extended Data are peer-reviewed evidence, not a place to hide
contradictory results.

## Extended Data Figure 1 | Discovery queue and validation gates

Show the full progression from grouped branch contrasts to 95 inspectable rows,
60 collapsed candidate directions and the focused exact-union gate. Include
the distribution of missing-evidence flags and validation classes. Make clear
that the focal candidate was selected after review and that discovery and
validation are not statistically independent.

## Extended Data Figure 2 | EIF4G2 transcript geometry and phase assignment

Show ENST00000339995 exon structure, transcript-coordinate conversion, uORF1,
uORF2, clean CDS, unambiguous phase positions, Kozak contexts, nucleotide
sequence and peptide sequence. Include a diagram of union construction and the
disjoint clean-CDS denominator.

## Extended Data Figure 3 | Complete raw tRNA reconstruction and lane audit

Show P-site offset, external CDS-start phase, target counts and exact-uORF1
phase for all 12 tRNA-Glu/control and four tRNA-Arg/control raw libraries and
every 15--40-nt class. Include the study-wide and per-lane shared-length sets,
all-length/short/canonical lane-coverage gates, run accessions and biological-
unit membership. Beside this, retain the original compendium lane audit with
dominant lengths, whole-library frame usage, CPM and allocation, explicitly
bracketing technical splits. Agreement between the independent raw rebuild and
old tracks is a reconstruction check, not additional replication.

## Extended Data Figure 4 | Feature-wide specificity at all count thresholds

Show the complete non-overlapping leader-uORF background at minimum counts of
5, 10 and 20 in every arm. Include tRNA-Glu effect, tRNA-Arg effect,
Glu-minus-Arg specificity, ranks, background medians and unscaled MADs. Do not
label the robust distances as z-scores.

## Extended Data Figure 5 | Total-RNA and EIF4G2 isoform audit

Show all source-reported EIF4G2 RefSeq isoforms in paired RPF and total-RNA
libraries, their abundance, isoform fractions and RPF/total-RNA proxy. Highlight
NM_001418 without suppressing low-abundance isoforms.

## Extended Data Figure 6 | Complete UV P-site calibration

Show inferred offsets, physical-offset gate, CDS-start read depth, phase-0
fraction and phase margin for each run and read length. Include all failed
lengths and the excluded 16-nt, −18-nt detector result.

## Extended Data Figure 7 | UV exact-length and position-deletion audit

Show all exact 15–34-nt EIF4G2 feature counts, allocations, target phases,
coordinate shifts, every single-position deletion and the order of greedy
three-position deletion. Low-count lengths should remain visible but be
distinguished from decisive shared-phase categories.

## Extended Data Figure 8 | ABCE1 recurrence and topology

Show the two ABCE1-targeting and two non-targeting units, uORF1/uORF2/clean-CDS
counts, allocation, phase and positional topology. Contrast increased uORF1
stop-end occupancy after ABCE1 loss with the different UV profile. This is
recurrence of allocation, not proof of one shared mechanism.

## Extended Data Figure 9 | Metastatic support and non-replication

Show the original MDA-LM2, CN34-LM1a and later MDA-LM2 datasets together. The
later dataset must show its +0.039 opposite effect, 21-nt versus 26-nt
condition-linked length difference and 17.28-percentage-point frame-0
imbalance. The legend should state that a general metastatic-state claim was
rejected.

## Extended Data Figure 10 | CDK1 adjacent branch and architectural ambiguity

Show the tRNA-Glu-specific CDK1 union shift beside EIF4G2, its six overlapping
near-cognate/other uORFs, position-17 ambiguity and tRNA-Arg deletion failure.
Use this only as pathway-level support. The current combined asset is
`results/dominant_rdg_trna_eif4g2_cdk1_axis_packet/trna_eif4g2_cdk1_axis.pdf`.

## Extended Data tables or source-data files

Large, machine-readable tables should be archived rather than converted into
unreadable figure panels:

1. Complete sample/run/biological-unit manifest.
2. Exact contrast metrics and validation gates.
3. Biological-unit and individual-uORF metrics.
4. Coordinate and position-deletion values.
5. Feature-wide tRNA specificity and codon counts.
6. Raw tRNA offset/calibration/run-length, lane-coverage and specificity tables.
7. UV offset/calibration/run-length tables.
8. Software versions, checksums and commands.

Each archived table should include a data dictionary and stable identifier.
