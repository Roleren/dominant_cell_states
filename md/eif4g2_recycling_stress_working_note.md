# EIF4G2 uORF1 as a candidate tRNA-specific branch and ribosome-stress sensor

Status: working note, updated 2026-08-18. This document records the current evidence
and open falsification tests. It is not yet a biological claim.

## Current verdict after the 2026-08-18 falsification pass

The strongest surviving observation is now narrower and more interesting than
the initial broad stress/metastasis story:

> EIF4G2 uORF1 is the strongest tRNA-Glu-specific leader-branch outlier among
> strict count-supported, non-overlapping leader uORFs shared with the matched
> tRNA-Arg experiment. The effect is biological-replicate separated, phased,
> robust to coordinate and peak deletion, and occurs without a corresponding
> change in the dominant EIF4G2 mRNA isoform.

This is not yet a claim that the two GAG codons directly cause the response.
A transcriptome-wide codon-dose test is null, so direct codon causality must be
established by reporter mutation. It is also not a claim that the switch is a
general metastatic state: a later independent MDA-LM2 dataset does not
replicate the original direction and has a severe condition-linked footprint
length/frame imbalance.

The biological consequence is unusually concrete. A 2025 Nature
Communications paper independently demonstrated that the same EIF4G2 uORF
produces an HLA-presented peptide and used an EIF4G2-leader reporter plus
uORF-start mutation to show strong induction under mitotic arrest:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC12391379/>. Thus the new candidate
intersection is tRNA-pool-specific control of an already experimentally proven,
immunogenic EIF4G2 uORF—not discovery that the uORF exists or is translated.

## Initial stress candidate, now a supporting branch

The initial hypothesis was that the annotated 51-nt ATG uORF1 in the
EIF4G2/DAP5 5' leader is a reproducible sensor of ribosome collision/recycling
stress. Two unrelated perturbations produce a large, replicate-separated,
phase-coherent redistribution toward this uORF:

1. UV irradiation in HeLa cells (GSE141459 / PRJNA593553).
2. ABCE1 knockout in HCT116 cells (GSE144163 / PRJNA602917).

The two perturbations have different downstream behavior. UV increases uORF1
occupancy while depressing the main clean CDS, whereas ABCE1 loss increases
uORF1 occupancy while retaining or slightly increasing clean-CDS abundance.
This argues against a generic low-coverage leader artifact and suggests that
uORF1 can report distinct collision/recycling states rather than acting as a
simple binary translational repressor.

The raw UV read-length gate now passes. Independently calibrated, shared
well-phased canonical 28--32-nt footprints and UV-associated 20--23-nt
footprints reproduce the switch separately, with complete replicate-range
separation and worst-case position-deletion robustness. Thus the UV recurrence
is no longer plausibly explained by pooling short and canonical footprints or
by applying one incorrect P-site offset across lengths. It remains supporting
recurrence rather than proof of the tRNA mechanism.

## Stronger branch discovered after the initial note: tRNA-Glu and metastasis

The compendium condition label `glu` was initially interpreted as glucose
deprivation. Raw study provenance shows that this was wrong: PRJNA310107 /
GSE77347 is tRNA-Glu(UUC) overexpression in poorly metastatic MDA-MB-231
cells. This correction exposes a more specific and potentially more important
mechanism.

The source SRA submission split each biological replicate into three sequencing
lanes, assigned separate Experiment and BioSample accessions, and stored the
lane name in `FRACTION`. The validation code was extended to prefer explicit
replicate labels over SRA Experiment IDs and to allow lane-like context fields
to be collapsed. All six case and six control lane-runs therefore reduce to two
true biological units in each arm.

Lane-collapsed tRNA-Glu(UUC) overexpression result:

- 2 case and 2 control biological units, each assembled from three lanes.
- Union-uORF counts: 1118 tRNA-Glu versus 2033 control.
- Unique clean-CDS counts: 2653 tRNA-Glu versus 10720 control.
- Clean-CDS allocation delta: -0.1371.
- uORF CPM log2 ratio: +0.7475.
- clean-CDS CPM log2 ratio: -0.3877.
- Biological-unit clean-CDS allocations: 0.6997 and 0.7065 in tRNA-Glu,
  versus 0.8259 and 0.8515 in controls; the ranges are completely separated.
- uORF phase-0 fraction: 0.811 tRNA-Glu versus 0.748 control.
- Coordinate-shift minimum absolute delta across -2 to +2 nt: 0.1201.

The read-length/frame QC contains both a sequencing-lane batch and one partial
condition linkage, so it must be described exactly. Lane 3 is 30-nt-dominant in
both conditions and both biological replicates. Lanes 1 and 2 are 28-nt-dominant
in both arms of replicate 1 and in control replicate 2, but are 30-nt-dominant
in tRNA-Glu replicate 2. Thus read length is partially condition-linked within
replicate 2, although frame-0 usage stays near 65% in every one of the 12 runs
and the aggregate frame-0 difference is only -0.082 percentage points.

The newly emitted technical-run audit provides the more decisive bound. Every
one of the six tRNA-Glu lanes has lower clean-CDS allocation than every one of
the six control lanes: case range 0.6667--0.7205 versus control range
0.8224--0.8625, a run-level separation of 0.1019. This includes the fully
length-matched replicate-1 lanes and the 30-nt lane-3 comparisons in both
replicates. Every individual run also has strong uORF phase-0 support
(0.708--0.838). Moreover, every case lane has higher total-coverage-normalized
uORF CPM than every control (204.5--252.4 versus 112.7--143.3), while every
case lane has lower clean-CDS CPM (502.4--578.0 versus 654.4--721.8). The
partial replicate-2 length linkage is therefore real and is not hidden, but it
cannot explain an effect present across all matched lane subsets and stable to
a +/-2-nt coordinate audit. Raw length-preserving reprocessing would remain
the strongest possible final audit.

A predeclared lane-subset summary makes that bound numerical. Clean-CDS
allocation deltas are -0.1262 using all three lanes of fully lane-matched
replicate 1, -0.1350 using only the 30-nt lane 3 from both replicates, and
-0.1503 using only the 28-nt lanes 1/2 from replicate 1. The partially
confounded replicate-2 lanes 1/2 give a similar -0.1330, not a uniquely larger
effect. All five subsets retain positive lane-range separation and phase-0
dominance in both arms. The reproducible summary is generated by
`scripts/dominant_rdg_eif4g2_trna_lane_subset_audit.R` and written under
`results/dominant_rdg_eif4g2_trna_lane_subset_audit/`. These technical subsets
are an artifact audit, not extra biological replicates or inferential n.

The ordinary shape gate fails because one biological unit has 36.9% of its
union-uORF signal at one coordinate, just above the conservative 35% ceiling.
This is not being hidden or relabeled as an ordinary pass. A formal
leave-one-position-out and pooled-top-three sensitivity audit now assigns
`C_shape_flag_leaveout_robust`: the shape warning remains, but the branch
redistribution survives deletion of the concentrated positions. The recurring
coordinate is transcript position 35, phase 0 of uORF1, corresponding to the
fifth translated codon in the current P-site convention. It is a shared,
well-phased within-uORF pause, not a uORF boundary artifact.

Importantly, the allocation effect does not disappear when that peak is
removed: the aggregate clean-CDS allocation delta remains -0.0888 in the
worst single-position deletion and -0.0513 after removing the three positions
with the highest pooled case-plus-control occupancy (35, 32, and 38). Thus the
peak strengthens the signal but is not its sole cause.

### Matched tRNA-Arg negative control

The same Cell paper and parental cell background contain a separately sequenced
tRNA-Arg(CCG) overexpression experiment (PRJNA310096 / GSE77317). It is a
high-value specificity control because both tRNA-Glu(UUC) and tRNA-Arg(CCG)
were shown by the source paper to promote metastasis.

- 2 case and 2 control biological units.
- Clean-CDS allocation delta: +0.0184, essentially null and in the opposite
  direction to tRNA-Glu.
- Biological-unit allocations: 0.8634 and 0.8416 in tRNA-Arg, versus 0.8325
  and 0.8399 in controls.
- Coordinate-shift minimum absolute delta: 0.0184, below the 0.03 gate.
- uORF phase remains strong, but there is no leader-to-CDS redistribution.

This argues against a generic consequence of tRNA overexpression, stable-line
generation, MDA culture, cycloheximide treatment, the Goodarzi laboratory
protocol, or low EIF4G2 coverage.

### Feature-wide specificity rank

A separate audit compares individual non-overlapping leader uORFs between the
tRNA-Glu and matched tRNA-Arg perturbations. To avoid picking EIF4G2 after the
fact, the audit ranks all features that meet the same minimum count support in
both perturbations and have at least five codons.

For EIF4G2 uORF1:

- tRNA-Glu uORF-allocation delta: +0.13618;
- tRNA-Arg uORF-allocation delta: -0.02260;
- tRNA-Glu-specific delta (Glu minus Arg): +0.15878;
- Glu-effect rank: 2;
- Glu-specificity rank: 1.

The rank is insensitive to three predeclared count thresholds: EIF4G2 remains
rank 1 of 35, 27, and 18 eligible features when every arm is required to have
at least 5, 10, or 20 uORF counts, respectively. The corresponding empirical
top-rank fractions are 0.0286, 0.0370, and 0.0556. Relative to the median and
unscaled median absolute deviation of each eligible background, the EIF4G2
Glu-minus-Arg effect is 13.27, 15.99, and 14.84 raw-MAD units above the median
at those thresholds. These are descriptive outlier scores, not parametric
z-scores or p-values. This is a useful outlier screen, not a substitute for a
biological-replicate-level genome-wide model. The exact EIF4G2 validator
remains the claim gate.

The reproducible audit is
`scripts/dominant_rdg_trna_codon_branch_audit.R`; outputs and the outlier plot
are under `results/dominant_rdg_trna_codon_branch_audit/`.

### Adjacent CDK1 branch: pathway support, not a cleaner replacement

A broader summed-annotation screen identified CDK1 above EIF4G2 by raw
Glu-minus-Arg branch magnitude. That nomination statistic double-counts
overlapping uORFs and is not itself evidence; it was used only to trigger an
exact-union alternative-candidate check. For tRNA-Glu, CDK1 has 705 union-uORF
and 508 unique-CDS counts versus 1451 and 2309 in control, respectively; the
clean-CDS allocation delta is -0.1953. Both two-versus-two biological-unit
ranges are separated by 0.1879, union phase-0 fractions are 0.768 versus 0.717,
and the coordinate-shift minimum absolute effect is 0.1951. Removing the
largest single position leaves -0.0955 and removing the pooled top three leaves
-0.0908. The matched tRNA-Arg contrast is small (+0.0253 clean-CDS allocation),
not replicate-range separated, poorly phased in the case arm, and fails peak
robustness.

This supports a broader, tRNA-Glu-specific cell-cycle/noncanonical-initiation
branch rather than an EIF4G2-only counting accident. It is not the primary
case because the CDK1 leader annotation contains six mutually overlapping
near-cognate/other candidate uORFs. Its largest shared peak is transcript
position 17, simultaneously phase 0 inside uORF1 and the CTG start of nested
uORF2; that architecture is intrinsically harder to assign than the exact,
non-overlapping ATG EIF4G2 uORF1. Local exploratory outputs are under
`results/dominant_rdg_cdk1_trna_glu_validation/` and
`results/dominant_rdg_cdk1_trna_arg_validation/`.

The dual hit is biologically coherent: the 2025 immunogenic-uORF paper
specifically highlights EIF4G2/DAP5 as a mitotic translation factor for CDK1.
The public data do not establish whether the two leader shifts form a causal
EIF4G2-to-CDK1 chain, a shared upstream response, or two parallel effects, so
the observation is a focused follow-up axis rather than a mechanistic claim.
A combined biological-unit and peak-deletion figure plus machine-readable
tables are generated by
`scripts/dominant_rdg_trna_eif4g2_cdk1_axis_packet.R` under
`results/dominant_rdg_trna_eif4g2_cdk1_axis_packet/`; the figure labels CDK1's
architecture as ambiguous and visibly retains the unstable tRNA-Arg deletion
control rather than suppressing it.

There is also an important 2025 novelty boundary: ALKBH8 loss, which perturbs
U34 modification on a set of tRNAs including tRNA-Glu(UUC), was reported to
reduce CDK1 protein while broadly altering A-ending-codon elongation
(<https://www.nature.com/articles/s41467-025-64144-0.pdf>). Thus a generic
tRNA-modification/CDK1-expression connection is not new. That study neither
isolates tRNA-Glu from the other ALKBH8 substrates nor reports the CDK1 leader
uORF redistribution or the matched EIF4G2 branch. It also emphasizes GAA
decoding, whereas EIF4G2 uORF1 contains GAG, further cautioning against a simple
direct-codon story.

This check exposed a reproducibility bug in manual-context resolution:
explicitly ignored technical axes such as `FRACTION` were still included when
requiring one context identity. The validator now excludes both the selected
case/control axis and every explicitly ignored axis; the focused regression
test passes.

### Matched total-RNA audit

GSE77347 contains total-transcript libraries prepared simultaneously with the
RPF libraries. The source-provided Cufflinks estimates independently address
whether the branch shift merely reflects a change in EIF4G2 transcript amount
or dominant isoform.

For NM_001418, the RefSeq transcript corresponding to ENST00000339995:

- mean total-transcript FPKM: 200.041 control versus 194.895 tRNA-Glu;
- total-transcript log2 ratio: -0.0376;
- fraction of summed EIF4G2 isoform-level total-RNA FPKM: 0.7467 control
  versus 0.7213 tRNA-Glu (delta -0.0255);
- mean RPF FPKM: 72.764 control versus 63.271 tRNA-Glu;
- RPF log2 ratio: -0.2017;
- fraction of summed EIF4G2 isoform-level RPF FPKM: 0.4968 control versus
  0.4685 tRNA-Glu (delta -0.0283);
- RPF/total-transcript translation-efficiency proxy log2 ratio: -0.1656.

The target mRNA isoform is therefore individually stable and remains dominant
in total RNA. Its composition shifts by only 2.5--2.8 percentage points in the
source isoform estimates, much smaller than the 13.7-point clean-CDS allocation
shift. This does not make short-read isoform assignment perfect, but it argues
strongly against isoform replacement as the whole explanation. The retrieval
and summary are reproducible with
`scripts/dominant_rdg_eif4g2_goodarzi_rna_audit.R`; outputs are under
`results/dominant_rdg_eif4g2_goodarzi_rna_audit/`. GEO explicitly documents
parallel RPF and total-RNA preparation:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77347>.

### Codon-dose hypothesis did not pass

EIF4G2 uORF1 contains two GAG codons, but a transcriptome-wide test does not
support a general relationship between uORF-versus-CDS Glu-codon enrichment
and the tRNA-Glu branch response. In the strict 18-feature matched analysis:

- tRNA-Glu delta versus relative Glu-codon fraction: Spearman rho 0.129,
  p = 0.610;
- tRNA-Glu-specific delta versus relative Glu-codon fraction: rho 0.117,
  p = 0.645;
- the matched tRNA-Arg negative-control correlation is also nonsignificant.

This falsifies the simple genome-wide codon-dose model. The EIF4G2 response is
locus-specific in these data. Its GAG codons remain a reporter-level hypothesis,
not an inferred mechanism.

### Natural metastatic derivatives from two patients

The same source study profiled two independently selected highly metastatic
derivatives and their isogenic poorly metastatic parental lines. The paper
measured tRNA-Glu(UUC) and tRNA-Arg(CCG) as among the most highly upregulated
tRNAs in both derivatives and confirmed the changes by independent assays.

MDA-LM2 versus MDA-parental:

- 2 versus 2 biological replicates.
- Union-uORF counts: 441 versus 120.
- Unique clean-CDS counts: 1610 versus 1876.
- Clean-CDS allocation delta: -0.1549.
- uORF CPM log2 ratio: +2.0325; clean-CDS CPM log2 ratio: -0.0565.
- Replicate allocation ranges are completely separated by 0.1280.
- uORF phase-0 fraction: 0.873 versus 0.817.
- Coordinate-shift minimum absolute delta: 0.1465.
- The worst single-position deletion leaves an allocation delta of -0.0819;
  removal of the three highest pooled positions (35, 38, and 23) leaves
  -0.0444.

CN34-LM1a versus CN34-parental:

- 2 versus 2 biological replicates.
- Union-uORF counts: 349 versus 347.
- Unique clean-CDS counts: 1338 versus 2775.
- Clean-CDS allocation delta: -0.0957.
- uORF CPM log2 ratio: +0.9525; clean-CDS CPM log2 ratio: -0.0827.
- Replicate allocation ranges are completely separated by 0.0704.
- uORF phase-0 fraction: 0.842 versus 0.706.
- Coordinate-shift minimum absolute delta: 0.0800.
- This arm is more dependent on the position-35 pause: the worst
  single-position deletion reduces the allocation delta to -0.0224, and
  deleting the pooled top three positions reduces it to -0.0028.

Both metastatic comparisons retain a conservative shape warning because uORF1
contains a dominant, reproducible phase-0 pause. MDA-LM2 is formally classified
`C_shape_flag_leaveout_robust`; CN34 remains
`D_spike_or_narrow_union_signal` and is supportive but not independently
claim-ready.

### Later independent MDA-LM2 dataset does not replicate the state effect

A later Goodarzi MDA-parental/MDA-LM2 Ribo-seq experiment, GSE186639 /
PRJNA774895, was present in the atlas but absent from the initial focused
packet. It is an important negative result:

- 2 versus 2 biological replicates;
- union-uORF counts: 243 MDA-LM2 versus 342 parental;
- clean-CDS counts: 2123 versus 2071;
- clean-CDS allocation delta: +0.0390, opposite to the original dataset;
- uORF CPM log2 ratio: -0.7798;
- coordinate-shift minimum absolute delta: 0.0133;
- frame-0 usage: 80.57% versus 63.29%;
- top read length: 21 nt in both MDA-LM2 libraries versus 26 nt in both
  parental libraries;
- classification: `D_frame_imbalance_risk`.

The within-condition replicates are consistent, but condition is inseparable
from footprint length/frame composition. This dataset cannot validate either
direction without raw length-specific offsets. More importantly, it means the
original MDA-LM2 and CN34 results must not be presented as independent proof of
a general metastatic-state switch. The evidence packet now includes this
nonreplication explicitly.

### Sequence-level interpretation

EIF4G2 uORF1 is 17 codons including its stop:

`ATG GAG GTG GCA GCG GGT ACC GAG TGG CGG CTG CAG CAG CGA CTC CTC TGA`

It contains two Glu codons, both GAG: 2/17 codons (11.8%), compared with 72/908
codons (7.9%) in the annotated main CDS. The source paper explicitly reports
that tRNA-Glu(UUC)-family tRNAs can wobble-decode both GAA and GAG and that GAR
content predicts increased ribosome occupancy after tRNA-Glu(UUC)
overexpression. The EIF4G2 result is nevertheless not a trivial replication of
that gene-level observation: it is a redistribution between a short uORF and
the same transcript's non-AUG main CDS, and the source article does not mention
EIF4G2.

The working mechanistic model is that increased tRNA-Glu availability changes
ribosome flux through the Glu-containing EIF4G2 uORF1 and amplifies occupancy
at a pre-existing phase-0 pause near codon 5, altering the fraction of ribosomes
that reach or initiate at the downstream GUG main start. This is a hypothesis,
not yet proof of causality at the two Glu codons. The null transcriptome-wide
codon-dose test makes an indirect, locus-specific response at least as plausible
as direct decoding at the two GAG codons.

The perturbation also does not distinguish mature-tRNA decoding from products
derived from the overexpressed tRNA. tRNA-Glu-derived fragments can have
post-transcriptional regulatory activity, including reported CDK1 regulation in
other biological systems, although those fragment studies do not describe this
human leader-uORF effect. Small-RNA/Northern profiling and a construct design
that separates mature tRNA abundance, anticodon decoding, and cleavage-derived
fragments are therefore necessary before naming the molecular route.

The original Cell source is
<https://pmc.ncbi.nlm.nih.gov/articles/PMC4915377/>. It establishes the tRNA
perturbations, metastatic phenotypes, codon mechanism, and the two metastatic
derivative models, but does not report this EIF4G2 uORF switch. A later breast
cancer study independently reports that EIF4G2/DAP5 expression and activity
are associated with and required for breast cancer EMT/metastasis:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC10895648/>. That makes the intersection
biologically consequential, while also requiring a careful novelty statement.

Two newer papers further constrain novelty. `ribofootPrinter` used EIF4G2 as a
worked example of a densely translated uORF-containing leader, but did not test
tRNA-specific regulation:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC12458181/>. More decisively, the 2025
mitotic-arrest study detected an EIF4G2 uORF-derived HLA peptide and validated
the leader/start-codon dependence with reporters:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC12391379/>. Therefore neither basal
uORF translation, stress inducibility in general, nor peptide production is new.
The searched literature still does not report the tRNA-Glu-specific branch
switch or compare it with tRNA-Arg.

A 2026 pancreatic cancer study reports the opposite organismal direction from
breast cancer—EIF4G2 loss accelerates pancreatic tumor progression and
metastasis—showing that DAP5 is a context-dependent cell-state regulator rather
than a universal metastasis driver:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC13276777/>. This further argues for a
precise breast/tRNA/uORF claim instead of a pan-cancer narrative.

### Focused experimental falsification

The shortest decisive follow-up would use the EIF4G2 leader reporter with its
native GUG main start and compare:

1. wild-type uORF1;
2. uORF1 AUG-to-nonstart mutation;
3. synonymous GAG-to-GAA and/or non-Glu substitutions that preserve peptide
   properties as far as possible;
4. mutation of the codon-5 pause site;
5. tRNA-Glu(UUC), tRNA-Arg(CCG), and matched vector backgrounds.
6. mature-tRNA and tRNA-fragment quantification, ideally with
   decoding-competent/cleavage-altered and anticodon controls.

Reporter output plus targeted ribosome profiling or toeprinting would
distinguish altered uORF initiation, altered elongation/pause occupancy, and
changed downstream GUG initiation. EIF4G2 protein measurement in the existing
stable lines would establish whether the redistribution changes functional
DAP5 output.

## Strict validation design

The focused validator is
`scripts/dominant_rdg_union_replicate_validation.R`, with reusable helpers in
`R/union-validation.R` and tests in
`tests/testthat/test-union-validation.R`.

It deliberately applies more stringent gates than the original grouped model:

- exact study/cell/context case-control recovery;
- SRA runs collapsed to biological units by explicit replicate labels when
  available, then Experiment, BioSample, sample name, or Run as conservative
  fallbacks;
- union of all annotated leader/overlapping uORFs, without double-counting;
- clean CDS made disjoint from the uORF union;
- minimum counts in both arms;
- broad positional support rather than a one-position spike;
- phase-0 support on unambiguous uORF positions;
- frame-QC balance between groups;
- direction stability after shifting region coordinates by -2 to +2 nt;
- complete separation of biological-replicate allocation ranges for the
  strongest `A` class.

## Exact EIF4G2 transcript geometry

Transcript: ENST00000339995.

- Leader length: 308 nt.
- uORF1: transcript positions 23--73, 51 nt, ATG, moderate Kozak context.
- uORF2: transcript positions 128--220, 93 nt, TTG, weak Kozak context.
- Main clean CDS: positions 309--3032, 2724 nt.
- The annotated main start is the known non-AUG GUG start of EIF4G2.
- The uORFs do not overlap each other or the main CDS.
- uORF-union length: 144 nt; clean-CDS length: 2724 nt.

Across the compendium, EIF4G2 is already unusually leader-heavy: its global
primary clean-CDS usage is approximately 0.138, versus a count-only expected
value of approximately 0.453. Therefore, the candidate effect is a further
stress-dependent redistribution within an already uORF-rich transcript, not
the appearance of a few isolated leader reads.

## UV result: GSE141459 / PRJNA593553

Important metadata correction: the compendium labels this context as HEK293 /
kidney, but the raw GEO record describes HeLa monosome Ribo-seq. The experiment
used a broad 15--34-nt footprint selection and is explicitly designed to
capture UV-induced short footprints.

Biological-unit validation:

- 2 UV and 2 control biological units.
- Union-uORF counts: 456 UV versus 239 control.
- Unique clean-CDS counts: 624 UV versus 1819 control.
- Clean-CDS allocation delta: -0.3061.
- uORF CPM log2 ratio: +1.9709, about a 3.92-fold increase.
- clean-CDS CPM log2 ratio: -0.341.
- uORF phase-0 fraction: 0.500 UV versus 0.485 control.
- Biological-unit allocation range separation: 0.2513.
- Coordinate-shift minimum absolute delta across -2 to +2 nt: 0.2796.
- Classification: `A_union_validated_range_separated`.

Feature decomposition:

- uORF1 counts: 377 UV versus 187 control.
- uORF1 raw-density log2 ratio: +1.008.
- uORF1 phase-0 fraction: 0.525 UV versus 0.481 control.
- uORF1 start-9-nt fraction: 0.385 UV versus 0.337 control.
- uORF1 stop-6-nt fraction: 0.138 UV versus 0.257 control.
- uORF2 counts: 79 UV versus 52 control, with weaker phase evidence
  (0.380 versus 0.500); uORF2 is not part of the robust claim.

The original study reported UV-triggered short footprints, ribosome
collisions, ZAK signaling, and the integrated stress response. It did not, as
far as the current literature search found, report an EIF4G2 leader-to-CDS
switch. Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7384957/> and GEO
record <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141459>.

### Raw read-length-resolved UV falsification result

All four SRA FASTQs were independently reprocessed with the documented random
base/adaptor structure retained until trimming, uniquely aligned to GRCh38,
and checked with `samtools quickcheck`. P-site offsets were inferred separately
for every run and read length against 2,000-plus non-overlapping CDS-start
calibration windows. A length could enter the decisive strata only if it was
phase-0 enriched in all four libraries. The shared well-phased lengths are 17,
19--22, 24--31, and 33 nt. One physically impossible detector result
(16 nt with a -18-nt start offset in control replicate 2) was discarded before
target mapping; that length is absent from every shared evidence category.

The exact raw result is stronger than the original pooled track:

- All shared well-phased lengths: uORF 1035 UV versus 389 control; clean CDS
  1081 versus 2279; allocation 0.5109 versus 0.8542; delta -0.3433.
- Canonical 28--32-nt shared well-phased lengths: uORF 438 versus 162; clean
  CDS 576 versus 1621; allocation 0.5680 versus 0.9091; delta -0.3411.
- Short 20--23-nt shared well-phased lengths: uORF 360 versus 143; clean CDS
  274 versus 344; allocation 0.4322 versus 0.7064; delta -0.2742.
- All-shared per-run allocation is 0.4954 and 0.5243 in UV versus 0.8088 and
  0.8699 in controls. Canonical-only and short-only ranges are also completely
  separated.
- uORF phase 0 is dominant in both arms of every passing stratum. For the
  all-shared pool it is 0.4850 UV versus 0.4807 control; for canonical reads it
  is 0.5685 versus 0.5062.
- Across -2 to +2 nt coordinate perturbations, the minimum absolute deltas are
  0.3133 all-shared, 0.3033 canonical, and 0.1906 short.
- Deleting the single position that most weakens the signal leaves deltas of
  -0.3299, -0.3017, and -0.2615, respectively.
- Greedily deleting the three positions that most weaken each signal still
  leaves -0.3115 all-shared, -0.2780 canonical, and -0.2380 short.

Nearly every exact read length with usable target counts shifts in the same
direction; the effect is not created by mixing a UV-enriched short class with a
control-enriched canonical class. This closes the principal raw-data artifact
gate for the UV recurrence. The raw outputs are in
`results/dominant_rdg_eif4g2_uv_length_resolved/`.

## ABCE1-loss recurrence: GSE144163 / PRJNA602917

This is an independent HCT116 experiment with two ABCE1-targeting sgRNAs and
two non-targeting controls, treated as two true biological units per arm.

Biological-unit validation:

- 2 ABCE1-KO and 2 control biological units.
- Union-uORF counts: 750 KO versus 191 control.
- Unique clean-CDS counts: 1335 KO versus 981 control.
- Clean-CDS allocation delta: -0.1967.
- uORF CPM log2 ratio: +1.7372, about a 3.33-fold increase.
- clean-CDS CPM log2 ratio: +0.1845.
- uORF phase-0 fraction: 0.501 KO versus 0.576 control.
- Biological-unit allocation range separation: 0.1328.
- Coordinate-shift minimum absolute delta across -2 to +2 nt: 0.1663.
- Classification: `A_union_validated_range_separated`.

Individual units:

- sgABCE1 SRR10959087: uORF 272, clean CDS 584, clean-CDS allocation 0.6822.
- sgABCE1 SRR10959088: uORF 478, clean CDS 751, allocation 0.6111.
- Control SRR10959089: uORF 123, clean CDS 542, allocation 0.8150.
- Control SRR10959090: uORF 68, clean CDS 439, allocation 0.8659.

Feature decomposition:

- uORF1 counts: 598 KO versus 117 control.
- uORF1 raw-density log2 ratio: +2.344.
- uORF1 phase-0 fraction: 0.599 KO versus 0.752 control.
- uORF1 start-9-nt fraction: 0.378 KO versus 0.368 control.
- uORF1 stop-6-nt fraction: 0.324 KO versus 0.222 control.
- uORF2 counts: 152 KO versus 74 control, but phase evidence is poor
  (0.118 versus 0.297); uORF2 is not part of the robust claim.

The increased uORF1 stop-end occupancy under direct ABCE1 loss is consistent
with a termination/recycling defect. UV does not share that exact within-uORF
topology, which is useful: the recurrence is not merely the same narrow peak
reproduced by both studies.

The original ABCE1 study focused on ribosome recycling, 3' UTR ribosomes, NMD,
and iron/ROS biology, and did not identify EIF4G2 in the current search. GEO:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144163>.

A 2026 Nucleic Acids Research reanalysis reported that ABCE1 depletion can
increase stop occupancy and can unexpectedly preserve/increase main translation
for transcripts with very short start-stop elements. This is mechanistically
adjacent and reduces the breadth of any novelty claim, but the searched article
does not establish the exact 51-nt EIF4G2 uORF1 recurrence across UV and ABCE1
loss: <https://academic.oup.com/nar/article/54/13/gkag627/8725953>.

## Why this is stronger than earlier leads

The original 95-row review queue collapsed to 60 candidate directions and had
no claim-ready row. Several initially exciting leads failed on inspection:

- PPP1R15A after EV71 infection survived exact union validation, but the
  GADD34/PPP1R15A-uORF stress mechanism in EV71 has already been published, so
  novelty is inadequate.
- ATF4/TRIB3 in PRJNA957283 reflected a mislabeled context: the experiment was
  doxycycline-induced SAMD9 R1293W ribotoxic stress, not infection, and the
  mechanism is already central to the source paper.
- HSPA5, LDHA, and EIF4G1 UV redistributions failed the strict uORF phase gate.
- HBS1L survived the UV gate but is close to the count threshold and is much
  weaker than EIF4G2.
- PRJNA401882 lanes initially inflated replicate counts; collapsing by SRA
  Experiment revealed only two biological replicates per arm.

EIF4G2 remains strong because it has hundreds of uORF counts in both
experiments, two biological units in every arm, range-separated replicate
effects, broad positional support, uORF phase coherence, coordinate-shift
stability, a defined single ATG uORF1 carrying most of the signal, and
independent recurrence under direct loss of a ribosome-recycling factor.

## Literature boundary

Known biology that must not be claimed as new:

- EIF4G2/DAP5 translation uses a non-AUG GUG main start and is regulated by a
  translated upstream ORF.
- EIF4G2 translation can resist severe elongation stress, and ribosome queuing
  can enable its non-AUG initiation.
- UV induces short ribosome footprints/collisions and stress signaling.
- ABCE1 loss perturbs ribosome recycling and stop-region occupancy.
- EIF4G2 uORF translation produces an HLA-presented peptide and is inducible
  by mitotic arrest in a leader/start-codon-dependent reporter.

Candidate new intersection after the raw UV and current literature audits:

> tRNA-Glu(UUC), but not the matched metastasis-promoting tRNA-Arg(CCG),
> produces a large, phase-coherent, replicate-separated redistribution into
> the exact immunogenic EIF4G2 ATG uORF1 while the dominant EIF4G2 mRNA isoform
> remains stable. UV reproduces an even larger switch independently in both
> short and canonical well-phased footprints, and ABCE1 loss provides a second
> recycling perturbation. Neither recurrence establishes how tRNA-Glu causes
> the switch.

Repeated exact searches through 2026-08-18 found no primary report connecting
tRNA-Glu(UUC) abundance to translation of the EIF4G2 transcript's immunogenic
ATG uORF. The adjacent components are known separately: Goodarzi et al. showed
tRNA-Glu/Arg-driven metastasis and codon-selective expression; DAP5/eIF4G2 is a
known regulator of other uORF-rich transcripts; and the 2025 immunogenic-uORF
paper established this EIF4G2 peptide and its mitotic inducibility. The
defensible novelty is their intersection, not any component alone.

Useful EIF4G2 reviews:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC9945437/> and
<https://pmc.ncbi.nlm.nih.gov/articles/PMC9977666/>.

## Work in progress and recovery information

### 2026-08-19 manuscript and raw tRNA checkpoint

A Nature-style working manuscript is now maintained in `manuscript/`. Its
Markdown sources separate the article, Methods, main-figure plan and legends,
Extended Data, cover letter, and an explicit claims/evidence/gaps ledger. Run
`scripts/build_eif4g2_manuscript.sh` to audit the central numbers against the
saved evidence tables and generate the Google-Docs-ready file
`manuscript/generated/eif4g2_working_manuscript.docx`. The generated Word file
is intentionally ignored; the Markdown files are canonical. The draft is
explicit that causal reporter, DAP5-protein and HLA-output experiments remain
missing.

The final unresolved computational caveat was tested directly by raw,
length-preserving reconstruction of all tRNA-Glu and matched tRNA-Arg RPF
libraries and passed the frozen gate described below. The exact 16-run manifest is encoded in
`scripts/download_goodarzi_trna_rpf_fastq.sh` and
`scripts/dominant_rdg_eif4g2_trna_length_resolved.R`. Runtime FASTQs and BAMs
are outside the repository at:

- `/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/Goodarzi_trna_fastq`
- `/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/Goodarzi_trna_length_resolved`

The downloader uses resumable ENA transfers and verifies the exact published
MD5 plus gzip integrity. Preprocessing uses the ARTseq/Illumina adapter
`AGATCGGAAGAGCACACGTCT`, retains 15--40-nt inserts, and uniquely aligns to the
same GRCh38 STAR index used for the UV audit. At the first pilot checkpoint,
SRR3131003 retained 8,769,318 of 9,593,373 reads after trimming (91.4%); the
major biological insert lengths were 28--34 nt. STAR uniquely aligned
1,452,176 reads (16.6%) while 75.1% were multimapping, as expected for these
short footprints. The indexed pilot BAM passed `samtools quickcheck`.
Independent raw P-site calibration of this control lane recovered 185 EIF4G2
union-uORF and 907 clean-CDS reads, for clean-CDS allocation 0.8306 and uORF
phase-0 fraction 0.778. The prior compendium track for the same run contained
176 and 815 reads, allocation 0.8224 and phase-0 fraction 0.773. Thus the new
pipeline closely reconstructs the target geometry without copying the old
P-site calls; this single control lane does not test the case--control effect.

The next preregistered-style checkpoint added both condition-matched lanes 1
and 2 from biological replicate 1. Across all independently calibrated lengths,
the two raw tRNA-Glu lanes had clean-CDS allocations 0.705 and 0.698, whereas
the matched controls were 0.831 and 0.848. Aggregated only as a technical-lane
audit, this is 0.701 versus 0.839 (delta -0.1376); the old processed-track
subset was 0.681 versus 0.831 (delta -0.1503). Raw uORF phase-0 fractions were
0.784 and 0.791 in the case lanes and 0.778 and 0.730 in controls. This early
checkpoint reproduces the direction, magnitude, separation and phase within a
fully length-design-matched subset, but it is not the final claim: replicate
2, lane 3, strict matched-length intersections and tRNA-Arg remain pending.

After case replicate 2 became available, lane 1 supplied a stricter complete
two-versus-two biological-replicate checkpoint. Read lengths 21, 23, 24,
26--31 nt passed the independent CDS-start phase gate in all four lane-1 runs.
Using only that intersection, case allocations were 0.691 and 0.731 and control
allocations were 0.826 and 0.840. Aggregated counts were 188 uORF/468 CDS in
case and 376/1,892 in control, giving 0.713 versus 0.834 (delta -0.1208), with
non-overlapping biological-unit ranges and uORF phase-0 fractions 0.766 and
0.750. This is the first raw checkpoint that already satisfies independent
length calibration and true biological replication within one matched lane;
full all-lane coordinate/deletion gates and tRNA-Arg specificity are still
pending.

Lane 2 independently passed the same complete checkpoint. Its all-four-run
phase intersection was 21, 24 and 26--31 nt. Case biological-unit allocations
were 0.660 and 0.703 versus controls 0.842 and 0.862. Aggregated counts were
215 uORF/466 CDS in case and 321/1,864 in control, giving 0.684 versus 0.853
(delta -0.1688), again with disjoint ranges and phase-0 fractions 0.791 and
0.751. Thus both independently processed 28-nt-dominant lane groups reproduce
the effect in both biological replicates using their own all-run calibrated
length intersections. The two lane-specific deltas bracket the original full
effect (-0.1208 and -0.1688 versus -0.1371); lane 3, joint robustness and
tRNA-Arg still remain decisive pending checks.

On 2026-08-19, before the remaining lane-3 case profiles and any raw
tRNA-Arg profiles were inspected, the final validator was tightened in two
ways. First, phase dominance, minimum counts, allocation direction and
biological-unit range separation must now pass for the exact 51-nt EIF4G2
uORF1 itself; the union of uORF1 and uORF2 is no longer sufficient. Second,
the all-length, 20--23-nt and 28--32-nt lane-matched categories each receive
their own lane-coverage gate, so a stratum cannot borrow support from a
different read-length class in a missing lane. At this checkpoint nine of the
twelve tRNA-Glu runs have complete per-run profiles; the lane-3 control from
replicate 2 is being rebuilt, and all remaining downloads are integrity
checked before promotion. These stricter gates are recorded now to avoid
post-result relaxation.

The next completed profile, lane-3 tRNA-Glu case replicate 1 (SRR3131023),
contained 329 exact uORF1, 9 uORF2 and 909 clean-CDS P sites after independent
raw calibration. Exact uORF1 phase 0 was 0.781. Its clean-CDS-versus-uORF1
allocation was 0.734, below both independently rebuilt lane-3 controls (0.835
and 0.867); the corresponding union allocation was 0.729. Thus every lane with
complete case data still reproduces the direction. The second lane-3 case,
the combined robustness gates and raw tRNA-Arg profiles remain pending and no
final gate is assigned at this checkpoint.

All 12 raw tRNA-Glu/control profiles are now complete. Across all independently
calibrated P sites, case contained 1,188 exact uORF1, 37 uORF2 and 3,267
clean-CDS reads, whereas control contained 2,057, 70 and 11,963. Clean-CDS-
versus-uORF1 allocation was 0.7333 in case and 0.8533 in control (delta
-0.1199); union allocation was 0.7273 versus 0.8490 (delta -0.1217). Every
case lane remained below every control lane, including the final lane-3 case
replicate (uORF1 allocation 0.741), and exact uORF1 was phase-0 dominant in
all 12 runs (per-run fractions 0.728--0.847). The raw biological-unit ranges
are also separated. These all-calibrated values are descriptive until the
predeclared lane-shared/read-length, coordinate, position-deletion and raw
tRNA-Arg specificity gates finish.

The decisive validator calibrates P sites independently for each run and read
length against the same 4,670 external CDS-start windows. A read length may
enter a study-wide tRNA-Glu or tRNA-Arg category only if it passes the phase
gate in every run of that perturbation. Because the tRNA-Glu material was
deliberately split into three lanes with different dominant lengths, a second
predeclared category admits a length only when it passes in all four
case/control runs from the same lane, requires every lane to contribute at
least one such length, and only then aggregates lanes into the two biological
units. This keeps length composition condition-matched without demanding that
different technical lanes share an identical footprint distribution. The
validator also reports exact-length and individual-lane contrasts, requires
uORF phase 0 in both arms, and repeats the coordinate and worst-case
position-deletion challenges. Its source parses successfully but should not be
interpreted until every BAM and output table is complete.

Safe recovery commands are:

```bash
GOODARZI_TRNA_DOWNLOAD_TMPDIR=/dev/shm \
  scripts/download_goodarzi_trna_rpf_fastq.sh
GOODARZI_TRNA_THREADS=8 scripts/process_goodarzi_trna_length_resolved.sh
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt \
  env -u LC_ALL Rscript --vanilla \
  scripts/dominant_rdg_eif4g2_trna_length_resolved.R
```

Both transfer and preprocessing scripts are resumable: validated FASTQs,
trimmed FASTQs, indexed BAMs and per-run P-site checkpoints are reused. No
claim should be upgraded merely because the pooled raw result agrees; the
study-intersection length categories and matched tRNA-Arg negative control are
the predeclared decision surfaces.

Integrity incident retained for provenance: the first SRR3131015 transfer had
the exact ENA byte count (514,107,844 bytes) but MD5
`73f7d9edc44e47c7a92a307f5026c5f6`, which disagreed with ENA's primary record
`d9d283f1ebae4951c28600089e275d22`. It was not processed or silently accepted;
the file was preserved with suffix `.md5_failed_73f7d9edc44e47c7`. A second
direct-to-S-drive transfer produced a different bad MD5
(`c7b211ce3edcab5dee58a812967054a2`) and failed gzip CRC and length checks; it
was preserved separately. Downloading the same official ENA object to the
memory-backed `/dev/shm` filesystem produced the declared MD5 and passed gzip,
and its subsequent copy to the S drive passed both checks again. This rules out
a persistently corrupt ENA object and is consistent with failure along the
direct concurrent-write path, but it does not identify the hardware or network
cause. The downloader now writes to `.incomplete`, validates MD5 and gzip,
and promotes a file to its final FASTQ name only after both checks pass; it also
supports `GOODARZI_TRNA_DOWNLOAD_TMPDIR=/dev/shm` so future transfers can be
verified before copying to the runtime drive.

### Completed raw UV checkpoint

The four GSE141459 FASTQ files were downloaded outside the repository to:

`/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/GSE141459_fastq`

Runs: SRR11569080, SRR11569081, SRR11569082, and SRR11569083.

Their verified compressed byte sizes are 1,077,274,739; 1,407,487,328;
1,898,418,828; and 1,514,124,242, respectively. Downloads used `wget -c`, so
an interrupted transfer can be resumed safely.

The library structure inferred from GEO and raw records is:

- remove the first four random 5' nucleotides;
- locate the constant 3' adapter `CACTCGGGCACCAAGGA`;
- remove the six random nucleotides immediately preceding that adapter;
- retain 15--34-nt biological inserts.

Completed decisive checks:

1. Trimmed and uniquely aligned every raw run with read length retained.
2. Inferred and independently validated P-site offsets per run and length.
3. Quantified EIF4G2 uORF1, uORF2, and clean CDS separately in short
   (especially 20--23 nt) and canonical footprint-length strata.
4. Required and observed the effect within both short and canonical
   well-phased strata in both biological replicates.
5. Confirmed exact-position, coordinate-shift, and worst-case one/three-position
   deletion robustness using unique alignments only.
6. Searched additional EIF4G2 contrasts and recycling perturbations, finding
   ABCE1 recurrence and a later metastatic-state non-replication.

Still required experimentally are a tRNA-Glu/EIF4G2 leader reporter, uORF-start
and GAG/GAA codon mutants, mature-tRNA versus tRNA-fragment measurements, DAP5
protein output, and HLA-peptide quantification.

Length-preserving preprocessing is implemented in
`scripts/process_gse141459_length_resolved.sh`. It reproduces the source
library structure with Skewer 0.2.2, removes the four 5' random nucleotides,
retains 15--34-nt biological inserts, aligns uniquely to GRCh38 with STAR, and
indexes each coordinate-sorted BAM. Runtime outputs are written to:

`/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/GSE141459_length_resolved`

The process can be safely resumed, or restricted to specified runs, for
example:

```bash
GSE141459_RUNS="SRR11569080 SRR11569081" \
  scripts/process_gse141459_length_resolved.sh
```

The raw evidence gate is implemented in
`scripts/dominant_rdg_eif4g2_uv_length_resolved.R`, with testable helpers in
`R/length-resolved-validation.R` and unit tests in
`tests/testthat/test-length-resolved-validation.R`. It uses ORFik's
CDS-start-based offset detector separately for each run and read length, then
requires phase-0 enrichment in all four libraries before a length can enter
the shared well-phased evidence strata. It quantifies exact lengths,
UV-associated 20--23-nt footprints, canonical 28--32-nt footprints, the full
calibrated pool, and coordinate-shift sensitivity. The final gate also requires
the EIF4G2 uORF union itself—not only the genome-wide calibration set—to be
phase-0 dominant in both arms of a passing length category. Run it only after
all four BAMs exist:

```bash
RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt \
  env -u LC_ALL Rscript --vanilla \
  scripts/dominant_rdg_eif4g2_uv_length_resolved.R
```

As of the latest 2026-08-18 raw-processing checkpoint, all four FASTQs and all
four indexed BAMs are complete, the preprocessing script passes `bash -n`, and
`samtools quickcheck` reports no BAM errors. A rerun can be restricted to one
interrupted run with:

```bash
GSE141459_RUNS="SRR11569082" GSE141459_THREADS=8 \
  scripts/process_gse141459_length_resolved.sh
```

The length-resolution helper tests, including impossible-offset and positional
leaveout regressions, pass. The union-validation helper tests pass 14/14. The
complete package suite passes 66/66 after the raw UV validator and leaveout
helper changes. The broad pipeline invariant suite also passes, including
script parsing, Shiny startup, registry/output assumptions, and key model
tables.

Current focused outputs are under ignored local result directories:

- `results/dominant_rdg_eif4g2_recycling_stress_evidence/`
- `results/dominant_rdg_eif4g2_trna_glu_lane_collapsed/`
- `results/dominant_rdg_eif4g2_trna_arg_negative_control/`
- `results/dominant_rdg_eif4g2_mda_lm2_validation/`
- `results/dominant_rdg_eif4g2_cn34_lm1a_validation/`
- `results/dominant_rdg_eif4g2_evidence_packet/`
- `results/dominant_rdg_eif4g2_uv_length_resolved/`
- `results/dominant_rdg_union_replicate_validation_eif4g2_recurrence/`
- `results/dominant_rdg_union_replicate_validation_eif4g2_top30/`

The combined evidence packet is generated by
`scripts/dominant_rdg_eif4g2_evidence_packet.R`. It contains one contrast-level
CSV, biological-unit and uORF1 tables, the complete position-deletion audit,
raw UV length/run/position-leaveout tables, and PNG/PDF overview figures. The
raw falsification figure is
`results/dominant_rdg_eif4g2_evidence_packet/eif4g2_uv_length_falsification.png`.

Do not describe this as a demonstrated tRNA mechanism. The computational
artifact gates and current novelty search now support an unusually strong,
experiment-ready candidate; causal reporter/protein/immunopeptidomic work is
still required for a publication-level biological claim.

### 2026-08-19 rapid laboratory decision plan

UV is now explicitly framed as an orthogonal confirmation of EIF4G2 branch
responsiveness, not as an equal primary story and not as evidence that UV and
tRNA-Glu share a mechanism. The cheapest causal decision is split into two
small screens: MDA-MB-231 vector/tRNA-Glu/tRNA-Arg crossed with native
EIF4G2-leader WT/uORF1-start mutant; and HeLa untreated/UV crossed with the same
reporter pair and a verified UV stress-engagement marker. The first screen
decides the tRNA paper. The second validates the locus and assay; it cannot
rescue tRNA failure. New Ribo-seq, endogenous editing and HLA work are deferred
until the tRNA-by-start interaction passes with verified mature tRNA and stable
reporter RNA. Full go, stop, redesign and spending criteria are saved in
`manuscript/eif4g2_lab_go_no_go_plan.md` and incorporated into the generated
Google Docs manuscript bundle.

### 2026-08-19 completed raw tRNA gate

All 16 raw Goodarzi tRNA libraries completed the frozen validator. No final
criterion was relaxed after the last lane-3 case or either tRNA-Arg target
profile became available.

The decisive tRNA-Glu lane-matched pool used lengths that passed independent
CDS-start calibration in all four case/control runs within each technical
lane, with every lane required to contribute. It contained 1,028 exact uORF1,
28 uORF2 and 2,568 clean-CDS P sites in case, versus 1,870, 65 and 10,488 in
control. Union clean-CDS allocation was 0.70861 versus 0.84424 (delta
-0.13563); exact-uORF1 allocation was 0.71413 versus 0.84868 (delta -0.13455).
The two case biological units were 0.70291 and 0.71303 for the union and
0.70783 and 0.71901 for exact uORF1. The controls were 0.82897 and 0.85564 for
the union and 0.83431 and 0.85938 for exact uORF1. Thus both metrics have
complete biological-unit range separation. Exact uORF1 was phase-0 dominant
in every one of the 12 tRNA-Glu/control raw runs.

The three independently matched lanes gave union deltas -0.12080, -0.16880
and -0.13019. The eight study-wide tRNA-Glu shared lengths were 21, 24 and
26--31 nt and produced a -0.14146 shift. Canonical 28--32-nt footprints gave
-0.14453. The short 20--23-nt pool gave -0.05872 but failed its exact-uORF1
phase and measurement gates; it is explicitly non-evidence. Coordinate shifts
from -2 to +2 nt retained at least 0.11605 absolute effect, deletion of any
single position retained at least 0.10027, and greedy deletion of positions
35, 32 and 23 retained -0.06030.

The identically processed raw tRNA-Arg lane-matched pool gave +0.01182 rather
than the Glu direction. The Glu-minus-Arg contrast was -0.14745 and the frozen
specificity gate passed. tRNA-Glu also passed the study-shared and canonical
specificity gates; tRNA-Arg did not pass a positive evidence gate. This closes
the raw read-length/P-site/lane/peak artifact question for the public source
experiment. It does not establish causal tRNA action, distinguish mature tRNA
from fragments, or measure DAP5 or peptide--HLA output.

Final tables and the four-panel falsification plot are under
`results/dominant_rdg_eif4g2_trna_length_resolved/`. The machine-readable plot
input is `eif4g2_trna_length_figure_source_data.csv`. The plot deliberately
shows non-passing exact-length points; aggregate claims use only the frozen
shared-length gates.
