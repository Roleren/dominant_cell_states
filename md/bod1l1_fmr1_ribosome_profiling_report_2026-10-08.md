# FMR1 and BOD1L1: condition-dependent ribosome occupancy and upstream ORF regulation

**Expert report — 8 October 2026**

FMR1 shows the stronger evidence for condition-dependent redistribution of
ribosome footprints between its upstream overlapping ORF and the main coding
region. Nutrient starvation provides the best replicated signal; UV exposure
produces a substantially larger shift in one study. BOD1L1 shows weaker,
context-specific responses, including an upstream-region increase in a G3BP
perturbation experiment. These findings nominate regulatory hypotheses. They
do not yet demonstrate initiation at the proposed upstream starts,
uORF-mediated repression, or changes in protein synthesis rates.

## Measurement and scope

The analysis comprises 3,857 human ribosome-profiling samples from heterogeneous
studies. Gene-specific comparisons use matched case/control metadata where
available, with additional pooling of same-condition replicates for count-based
allocation estimates. Individual libraries and comparisons are not equivalent
to independent biological studies. Neither gene contributes to the
dominant-state scores used elsewhere in the project: the results discussed
here come from direct measurements of their transcript regions.

Two readouts must be distinguished:

- **Clean-CDS occupancy:** library- and length-normalized footprint density in
  the main CDS after excluding regions overlapping the specified upstream ORFs.
  Reported percentages are model-derived changes from matched contrasts, some
  using shrinkage or pseudocounts, rather than literal ratios of raw pooled
  counts. Without paired RNA-seq adjustment they are not translation-efficiency
  estimates.
- **Regional allocation:** the fraction of measured feature counts assigned to
  CDS, leader ORFs, overlapping ORFs, and other annotated features. Changes are
  given in percentage points. They measure composition, not protein-production
  flux, and differ from the separate density-based relative-usage summaries.

Matched effects are summarized with random-effects models, while several
count/composition models compare regional allocation. Omitting individual
studies or study groups assesses robustness. Ordinary meta-analytic intervals
can appear supportive while more conservative study-held-out calibration
includes zero. Agreement among allocation models is useful sensitivity
evidence, but the models reuse the same reads and are not independent
replication.

## Transcript models and putative upstream ORFs

Coordinates are 1-based, inclusive transcript coordinates on the selected
reference isoforms.

| Gene / transcript | Region | Proposed start | Length | Overlap with annotated main CDS |
|---|---|---|---:|---:|
| FMR1 / ENST00000370475 | 18–98, leader ORF | CTC | 81 nt | 0 |
| FMR1 / ENST00000370475 | 40–57, leader ORF | GTG | 18 nt | 0 |
| FMR1 / ENST00000370475 | **110–322, overlapping ORF** | **CTG** | **213 nt** | **61 nt** |
| BOD1L1 / ENST00000040738 | 16–36, leader ORF | CTG | 21 nt | 0 |
| BOD1L1 / ENST00000040738 | **47–307, overlapping ORF** | **TTG** | **261 nt** | **169 nt** |

The corrected FMR1 overlapping interval is CTG–TAA and has no internal
in-frame stop. Its two leader intervals also terminate at their first in-frame
stops. The position-18 CTC boundary remains unresolved as an initiation site.
BOD1L1 additionally has an annotated N-terminal-extension feature, measured
only in its unique non-CDS extension.

The overlap masks leave 1,838 clean-CDS bases for FMR1 and 8,987 for BOD1L1.
Reference-isoform selection does not establish that every sample expresses only
that isoform. Overlapping-ORF counts sum coverage across the complete nominated
interval; they are not restricted to the proposed reading frame. Main-CDS
ribosomes can therefore contribute to counts in the shared sequence. This is
particularly important for BOD1L1, where approximately 65% of the overlapping
ORF lies within the main CDS; the corresponding proportion is approximately
29% for FMR1. The nested FMR1 leader intervals also share sequence and are
measured separately. Allocation denominators therefore represent measured
features, not mutually exclusive molecular translation outcomes.

## FMR1: starvation-associated redistribution

Across nutrient-starvation comparisons, pooled clean-CDS occupancy changes
little: the atlas estimate is approximately +1.9%, with an ordinary interval
spanning approximately −5.1% to +9.4%. This is not evidence for a consistent
change in main-CDS occupancy across starvation conditions.

The allocation result is more informative. The joint count-based summary
shifts approximately **3.0 percentage points away from clean CDS** and
**2.8 points toward the overlapping ORF**. Other study-aware allocation
estimators give the same direction, with different effect magnitudes. Three
studies pass the strict grouped/count-support gates. Support persists when
individual studies are omitted; grouped omissions retain the direction,
although direction retention does not imply that statistical support persists
in every omitted-group fit.

A branched-chain-amino-acid deprivation comparison, labelled NoBCAA versus WT
in PRJNA472636, is a prominent individual example. It should be treated as a
sample-level review target; starvation conditions are not interchangeable, and
the pooled result spans heterogeneous nutrient perturbations.

**Biological interpretation:** the data are compatible with stress-dependent
changes in leader translation, start-site selection, or relative ribosome
loading. They are also compatible with altered elongation or residence time in
the upstream region. The stable overall CDS-density estimate alongside
changing allocation makes regional regulation a more defensible working
hypothesis than generalized FMRP suppression.

Published work provides a relevant precedent: translation of the native FMR1
CGG-containing leader can form an inhibitory overlapping uORF and regulate
FMRP synthesis, including in unexpanded-repeat contexts.
[Rodriguez et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32066985/).
However, the present CTG110 interval has not been demonstrated to correspond
to that characterized ORF. Published CGG-associated initiation sites include
ACG/GUG starts; sequence alignment, reading frame, and sample repeat length
must establish any correspondence. Cap-dependent scanning and
structured-leader effects are plausible mechanisms, rather than conclusions
from these regional counts.
[Kearse et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4854189/).

## FMR1: a large UV-associated shift in one study

The UV versus control comparison in PRJNA593553 shows approximately
**45% lower clean-CDS footprint density**, together with a count-allocation
change of approximately **−35 percentage points for CDS** and
**+30 points for the overlapping ORF**.

This is the strongest individual redistribution lead. Multiple allocation
estimators agree on its direction. Nevertheless, it comes from one study,
with two runs per side in the pooled comparison, and the minimum
overlapping-ORF count across its two sides is only 29. Independent-study
robustness cannot be evaluated. The narrow ordinary within-study interval
should not be read as precise evidence for a general UV response; conservative
calibration weakens interval support. The current candidate assessment remains
measurement-limited, with no validated coverage/browser review.

The source is a **HeLa ribosome-collision study**, with untreated and UV
libraries, rather than a demonstrated senescence experiment.
[GSE141459 / PRJNA593553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141459).
Its biological setting makes ribotoxic stress, altered elongation, and
collision-linked signaling especially relevant alternatives to increased
upstream initiation. The source study connects ribosome collisions with
ZAK/SAPK and GCN2-mediated stress signaling.
[Source study record](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA593553).

The current evidence is therefore a large, biologically plausible **regional
footprint redistribution requiring validation**, not proof that CTG initiation
represses FMRP synthesis.

## FMR1: infection responses depend on cellular context

The pooled viral-infection estimate is approximately **25% higher clean-CDS
occupancy** across 13 represented studies, 11 with strict-supported CDS
comparisons. The ordinary pooled interval is approximately +11% to +42%, but
heterogeneity is substantial (I² approximately 0.82), and conservative
calibration weakens interval support.

The upstream-allocation response is not uniform. Context-stratified modeling
suggests that A549 comparisons lose overlapping-ORF allocation while the
remaining contexts gain it. Estimated changes are approximately −7.1 versus
+1.9 percentage points. This contrast is a reason to examine individual
cell/virus/protocol combinations, not to assign a universal infection-induced
uORF switch.

Changes in transcript abundance, leader isoform use, viral host-shutoff effects,
and initiation or elongation could all contribute. The analysis does not
establish increased FMRP protein abundance, repeat-associated translation, or a
change in FMRP target regulation.

## BOD1L1: weaker general signals and one upstream-region lead

BOD1L1's overlapping-region and clean-CDS measurements show mixed or weak
coupling across the compendium. No BOD1L1 comparison enters the final
merged-replicate redistribution candidate queue.

| Context | Observed clean-CDS response | Interpretation |
|---|---:|---|
| ISR/ER/translation-stress comparisons | Approximately +7% | Small pooled increase without a strong accompanying upstream-ORF switch; conservative interval calibration includes zero. |
| UV comparison in PRJNA593553 | Approximately −28% | One-study decrease; branch counts are insufficient to establish a reproducible allocation shift. |
| Heat/osmotic/recovery comparisons | Approximately −34% | Heterogeneous comparisons from one study; exploratory rather than a general recovery response. |

The most interesting upstream example is a U2OS comparison in PRJNA1071172,
labelled G3BP1 KO versus WT in the analysis metadata. Overlapping-region
footprint density increases approximately **2.1-fold**, while clean-CDS density
increases only approximately **7%**. The disproportionate upstream response
is compatible with local leader regulation, but the original comparison is
more complex than its simplified label suggests.

Inspection of the eight contributing sample titles shows four KO-labelled
libraries versus four G3BP1-labelled libraries. **Each arm contains two NaAs
libraries and two water-control libraries.** The reported effect thus pools
stress-treated and untreated samples in balanced proportions; it is not a
stress-specific contrast. Balance avoids a simple difference in treatment
proportions between arms but does not resolve a genotype-by-treatment
interaction. The source study also uses G3BP1/2 knockout and rescue/isoform
designs, so these labels should not be interpreted as a verified isolated
G3BP1 knockout versus parental wild type.

G3BP perturbations can affect RNA condensation and stress-associated gene
expression, making the context biologically relevant without establishing
direct G3BP regulation of BOD1L1.
[Liboy-Lugo et al., 2024](https://pubmed.ncbi.nlm.nih.gov/39356796/).
The appropriate follow-up is to separate the water and arsenite comparisons
and establish the exact knockout/rescue backgrounds before interpreting this
as a stress-granule-dependent BOD1L1 response.

BOD1L1, called BOD1L in earlier literature, protects stressed replication
forks by restraining deleterious resection and supporting RAD51-associated
fork protection.
[Higgs et al., 2015](https://pubmed.ncbi.nlm.nih.gov/26166705/).
That function makes stress-responsive expression biologically plausible,
but it does not predict the direction of its immediate ribosome-occupancy
response. The UV decrease cannot establish loss of fork-protection activity,
and the G3BP-associated upstream increase cannot establish a causal connection
between stress granules and genome maintenance.

## Evidence needed to distinguish the regulatory hypotheses

The most informative next analysis is **replicate-level, nucleotide-resolved
coverage with footprint-length and frame stratification**, prioritizing FMR1
starvation/UV and the BOD1L1 G3BP comparison. Leader-only footprints should be
separated from shared CDS sequence, and distributed in-frame coverage
distinguished from isolated peaks, pausing, or non-ribosomal protection. The
UV source makes length-dependent footprint states particularly important.

Paired RNA-seq would distinguish transcript abundance from altered ribosome
engagement. TSS/isoform analysis and sample-specific FMR1 repeat information
would test whether the reference leader is appropriate. The CTC18 boundary
and exact G3BP knockout/rescue contrast require resolution before stronger
mechanistic statements. Shared and nested feature positions should also be
handled explicitly when constructing allocation denominators.

If those checks support translation of the proposed regions, initiation-site
profiling and leader reporters could test start-site use and its effect on
downstream synthesis. Start substitutions must be interpreted alongside
possible changes in GC-rich RNA structure; protein-level measurements are
required to connect footprint redistribution to FMRP or BOD1L1 output.

**Present interpretation:** FMR1 is a credible candidate for stress-dependent
cis-regulation involving its leader/overlapping ORF, with starvation providing
the strongest replicated allocation evidence. BOD1L1 remains a context-specific
exploratory lead. Neither gene currently has a validated causal uORF mechanism
from this analysis.

## Local data provenance

Numerical observations use the completed 8 October 2026 analysis with the
corrected FMR1 start at 110. Filtered outputs provide traceability independently
of pipeline terminology:

- [Matched-condition atlas estimates](../results/targeted_fmr1_start110_2026-10-08/dominant_rdg_atlas_cards.csv)
- [Final redistribution candidate table](../results/targeted_fmr1_start110_2026-10-08/merged_replicate_translon_shift_candidate_queue.csv)
- [Upstream-region/CDS coupling summaries](../results/targeted_fmr1_start110_2026-10-08/dominant_ouorf_cds_coupling_gene_summary.csv)
- [Feature annotation and coordinates](../results/targeted_fmr1_start110_2026-10-08/rdg_feature_annotation.csv)
- [Completed run summary](../results/targeted_fmr1_start110_2026-10-08/run_summary.csv)

The G3BP comparison contains KO/NaAs runs SRR27792535–SRR27792536,
KO/water runs SRR27792541–SRR27792542, G3BP1/NaAs runs
SRR27792533–SRR27792534, and G3BP1/water runs SRR27792539–SRR27792540.
Treatment pooling was verified from the external sample metadata's
`sample_title` field; this audit changes interpretation, not the existing
numerical estimate. The analysis and its simplified sample labels require
further stratification before treatment-specific effects can be reported.

These results are exploratory reanalyses, not gene-specific conclusions
reported by the original source studies.
