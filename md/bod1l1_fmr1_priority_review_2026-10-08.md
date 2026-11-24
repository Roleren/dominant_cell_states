# BOD1L1 and FMR1 priority-panel review

An [expert-facing ribosome-profiling report](bod1l1_fmr1_ribosome_profiling_report_2026-10-08.md)
describes the observations, biological hypotheses, and measurement limitations
without requiring familiarity with the pipeline. Its sample-level audit found
that the BOD1L1 G3BP comparison pools arsenite and water-control libraries in
balanced proportions and uses simplified knockout/rescue labels; the reported
effect is not treatment-specific.

The user corrected the FMR1 overlapping start from 143 to **110** after the
first completed run below. On ENST00000370475, tx 110–322 is **CTG–TAA**
(213 bases). The source queue has been corrected and the superseding core/RDG
run completed all 50 steps successfully. The earlier GGC-boundary concern applies to the old
143–322 input, not the corrected overlapping feature. The leader inputs remain
18–98 (CTC–TGA) and 40–57 (GTG–TGA).
An independent in-frame scan confirms no internal stop in any of these three
ranges: their first in-frame stops occur at 98, 57, and 322, respectively.
The position-18 CTC initiation boundary remains a separate review caveat.

Both genes are included as analysis-only priorities in
`curated_inputs/curated_additional_gene_priorities.csv`. The loader unions this
tracked input with any legacy priorities under `results/curated_inputs/`.
Neither gene is assigned to a scored dominant-state module.

BOD1L1 (BOD1L in older literature) protects stalled replication forks and
limits damaging resection ([Higgs et al., 2015](https://pubmed.ncbi.nlm.nih.gov/26166705/)).
This motivates replication-stress analysis, but does not establish that its
ribosome occupancy is a consistent positive proliferation or negative
senescence marker across the human compendium.

FMR1 encodes FMRP, an RNA-binding translation regulator
([Darnell et al., 2011](https://pubmed.ncbi.nlm.nih.gov/21784246/)).
A newer neuronal study distinguishes RNA-granule regulation from changes in
ribosome stall positions
([neuronal RNA-granule study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13384498/)).
These mechanisms motivate gene-level analysis without treating FMR1 abundance
as a directional ISR or mTOR/TOP activity marker.

Module membership would require reproducible, directionally coherent evidence
in matched conditions and a check against cell-type and study confounding.
The measurement remains ribosome occupancy/clean-CDS allocation, not protein
production or pathway activity by itself.

The user subsequently added BOD1L1 tx 47–307 as an overlapping uORF and
FMR1 tx 18–98 and 40–57 as leader uORFs plus tx 143–322 as an overlapping
uORF to the manual text queue. These map against the selected transcripts
ENST00000040738 and ENST00000370475, respectively. The initial run without
these additions was stopped and superseded by a run starting with manual
validation. The generated manual manifest and downstream feature tables are
the authority for category acceptance and clean-CDS masking.

The refreshed manual manifest accepted all four entries with matching
positional categories. Generated manual files now live exclusively under
`results/manual_translons/`; the source directory retains the editable queue.

To process the full existing queue, manual transcript lookup now falls back
to unique canonical annotations for genes absent from the measured panel;
it retains measured choices and rejects ambiguous fallback mappings. Verified
fallbacks on this server are IFIT2/ENST00000371826,
ZNFX1/ENST00000396105, PELO/ENST00000274311, and HK2/ENST00000290573.

An independent sequence/category check on the selected transcripts found:

| Gene | Transcript coordinates | Category | First triplet | Stop codon |
|---|---|---|---|---|
| BOD1L1 | 47–307 | overlapping uORF | TTG | TGA |
| FMR1 | 18–98 | leader uORF | CTC | TGA |
| FMR1 | 40–57 | leader uORF | GTG | TGA |
| FMR1 | 143–322 | overlapping uORF | GGC | TAA |

The CTC/GGC boundaries need clarification as intentional noncanonical/RAN
candidates versus transcript-coordinate mismatch. Their positional categories
and occupancy measurements can be computed as supplied, but category acceptance
does not validate an initiation site. The raw coordinate check is saved locally
under `results/targeted_bod1l1_fmr1_2026-10-08/`.

## Corrected overlapping-start rerun (start 110)

The corrected run finished on 2026-10-08 at 15:37 UTC, with all 50 core/RDG
steps marked `ok`. The selected manual manifest and RDG annotation both
contain FMR1 tx 110–322; no FMR1 feature with start 143 remains in the current
annotation. Each of the 3,857 samples has measurements for the corrected
213-base overlapping feature, alongside the unchanged 81- and 18-base
leader features. The coupling table now assigns the CTG start to
`near_cognate`; this sequence classification is not itself evidence of
initiation. The coupling result remains `mixed_or_weak_coupling`.

Across samples, overlapping-feature raw counts increase from 193,460 for the
old 180-base range to 207,594 for the corrected 213-base range. These totals
are not directly comparable expression estimates without the length change.
FMR1's clean-CDS mask remains 61 bases, leaving 1,838 measured clean-CDS
bases. All 3,857 state-score rows remain exactly equal to the preserved
pre-addition table. The panel and atlas counts remain 129 genes and 1,148
cards across 82 genes.

FMR1's viral-infection and starvation consensus labels remain
`B_replicated_robust_consensus`, and DNA damage/senescence remains
`C_multi_model_branch_consensus`. In the corrected starvation result,
leave-one-study-out is `consensus_loo_support_robust` and grouped transport
omission is `consensus_transport_direction_robust`, superseding the earlier
transport-fragile result below. DNA damage/senescence remains not evaluable
in those two robustness layers.

The merged-replicate queue still contains one FMR1 DNA-damage/senescence
candidate (PRJNA593553), with `clean_CDS_down_uORF_up` direction and
`measurement_limited_candidate` status. Its clean-CDS and overlapping-uORF
allocation deltas are approximately -0.350 and +0.299. The minimum
two-sided dominant-uORF count is now 29 rather than 25, but the count gate
and the other missing-evidence flags still prevent claim readiness; there
are zero claim-candidate rows. The next action remains
`rebuild_or_review_translon_measurement`. BOD1L1 remains
`G_general_atlas_review`, with no merged-replicate candidate-queue entry.

The completed run summary, sequence/in-frame-stop QC, filtered gene tables,
and exact final assertions are retained in
`results/targeted_fmr1_start110_2026-10-08/`. Manual cache validation also
passes. Package source code was unchanged in this correction; the broader
integration limitations documented below remain outstanding.

## First completed manual-aware rerun (overlapping start 143; superseded)

The superseding run completed all 50 core/RDG steps with status `ok` on
2026-10-08. The measured panel contains 129 genes; the atlas contains
1,148 gene/design cards for 82 genes, including 14 cards each for BOD1L1 and
FMR1. All 3,857 sample state-score rows are exactly equal to the preserved
pre-addition table. No scored module membership was changed.

BOD1L1's additional overlapping feature excludes 169 CDS bases, leaving
8,987 clean CDS bases. FMR1's overlapping feature excludes 61 CDS bases,
leaving 1,838 clean CDS bases. FMR1 now has two leader uORFs and one
overlapping uORF in the RDG feature annotation. These are measurements of
the supplied ranges, subject to the start-coordinate caveat above.

The final consensus tables classify BOD1L1's best result as
`G_general_atlas_review`. FMR1 has `B_replicated_robust_consensus` labels
for viral infection and nutrient starvation, and
`C_multi_model_branch_consensus` for DNA damage/senescence. These labels
are review priorities, not claim-readiness decisions. In particular,
starvation is direction-robust under leave-one-study-out but fragile under
grouped transport omission; DNA damage/senescence is not evaluable in those
robustness layers. The viral auxiliary-state audit also warns that branch
support is weak in that audit.

The merged-replicate claim gate contains one FMR1 candidate, in
DNA damage/senescence (PRJNA593553): `clean_CDS_down_uORF_up`. Its
classification is `measurement_limited_candidate`, with zero claim-candidate
rows. The next action is `rebuild_or_review_translon_measurement`. Missing
evidence includes two-sided uORF count support, browser labeling,
robustness/partial-pooling support, and static coverage/backing-data assets.
BOD1L1 has no entry in this candidate queue.

Filtered per-gene tables, the run summary, coordinate QC, and exact target
assertions are retained locally under
`results/targeted_bod1l1_fmr1_2026-10-08/`. All 74 package checks passed.
The complete integration invariant suite stopped because the desktop
S-drive/static-inspection readiness tables are absent on this server.
An additional diagnostic run excluding only that prerequisite block passed
parsing, Shiny startup, registry, manual-manifest, and candidate-rank checks,
then stopped on the auxiliary antigen-presentation candidate assertion.
Those auxiliary score outputs were not rebuilt in this run; this is an
outstanding broader validation limitation, not a passing full-suite result.
