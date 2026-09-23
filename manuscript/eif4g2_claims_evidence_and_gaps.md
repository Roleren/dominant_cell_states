# Claims, evidence and submission gates

Updated: 2026-08-19

This document is the paper's anti-overclaiming ledger. A sentence can enter the
title, summary or cover letter only at the level supported below. “Pass” means
the saved analysis supports the statement; it does not mean the full paper is
submission-ready.

## Claim hierarchy

| Claim | Current status | Strongest evidence | What remains unresolved | Allowed wording now |
|---|---|---|---|---|
| tRNA-Glu(UUC) overexpression redistributes EIF4G2 occupancy from clean CDS towards the uORF union | Strong pass as a public-data observation | Raw lane-matched delta −0.136; `n = 2` versus `n = 2`; separated unit ranges; study-shared −0.141; canonical −0.145 | Independent newly generated perturbation | “tRNA-Glu overexpression was associated with a 0.136 clean-CDS allocation decrease in independently calibrated raw footprints” |
| The tRNA response is selective relative to tRNA-Arg(CCG) | Strong pass within the matched public experiments | Identically processed raw tRNA-Arg delta +0.012; Glu-minus-Arg −0.147; specificity gate passes; feature-wide EIF4G2 rank 1 | The two stable-line studies are separate experiments; matched vector/construct work in one factorial experiment | “The matched tRNA-Arg perturbation did not reproduce the switch” |
| The exact ATG uORF1 is the main EIF4G2 leader feature involved | Strong pass as occupancy assignment | Raw exact-uORF1 delta −0.135; 1,028 versus 1,870 case/control P sites; phase-0 dominant in all 12 Glu/control runs | Initiation-site mapping in a new experiment; uORF1 start mutation | “Occupancy maps predominantly to the annotated ATG uORF1” |
| The tRNA effect is not caused solely by lane/read-length/P-site or one peak | Strong pass against tested artifacts | Raw matched lanes −0.121, −0.169 and −0.130; ±2 nt minimum 0.116; worst single deletion 0.100; greedy-three delta −0.060 | A new experiment is still required to exclude source-study-specific biology | “The effect survived raw lane, read-length, frame, coordinate and position-deletion challenges” |
| EIF4G2 RNA abundance or major isoform replacement explains the shift | Evidence against | NM_001418 total-RNA log2 change −0.038; isoform fractions change only 0.025–0.028 | Long-read or isoform-targeted RNA measurement | “A commensurate RNA or dominant-isoform change was not detected” |
| The two uORF1 GAG codons directly sense tRNA-Glu | Not established | Sequence contains two GAG codons | Genome-wide codon-dose test is null; requires synonymous and non-Glu reporter mutants | “The GAG codons are a reporter-level hypothesis” |
| Mature tRNA decoding, rather than a tRNA-derived fragment, causes the shift | Not established | Original perturbation increases tRNA-Glu construct abundance | Mature and fragment species not separated | “The perturbation does not distinguish mature-tRNA and fragment-mediated routes” |
| UV reproduces the branch shift independently of short-footprint/P-site mixing | Strong pass | Raw all-shared −0.343; canonical −0.341; short −0.274; all with range, coordinate and deletion robustness | Different cell type and perturbation; no causal link to tRNA-Glu | “UV independently reproduced the occupancy redistribution in short and canonical footprints” |
| ABCE1 loss reproduces the branch shift | Pass as recurrence | Delta −0.197; separated ranges; exact-union class A | No raw length reprocessing; different downstream CDS response | “ABCE1 loss produced a second recycling-related recurrence” |
| tRNA-Glu, UV and ABCE1 act through the same molecular mechanism | Not established | Directionally concordant allocation | Heterogeneous cells/perturbations; UV and ABCE1 topology and CDS outputs differ | “The perturbations converge at the same branch” |
| tRNA-Glu increases DAP5 protein | Not measured | Clean-CDS CPM decreases in the public tRNA contrast; no protein data | Quantitative endogenous DAP5 protein | Do not claim |
| tRNA-Glu increases the EIF4G2-derived HLA peptide | Not measured | The exact uORF is independently known to yield an HLA ligand under mitotic arrest | Targeted immunopeptidomics or peptide–HLA assay under tRNA perturbation | “The affected uORF has an independently established immunogenic output” |
| This is a general metastatic-state switch | Rejected | Original MDA-LM2 and CN34 support | Later MDA-LM2 is weak/opposite with length/frame imbalance | “A later dataset did not support a general metastatic-state effect” |
| EIF4G2 is the only relevant tRNA-Glu leader target | Not supported | EIF4G2 is the cleanest strict specificity outlier | CDK1 has a larger union shift but ambiguous overlapping-uORF architecture | “EIF4G2 is the strongest resolved non-overlapping leader-uORF case” |

## Submission blockers

A Nature submission should not be sent until the first five gates pass. Gates
6–8 are publication-quality requirements rather than mechanistic experiments.

### Rapid investment triage — before the expensive gates

Run two small, separately interpretable reporter screens before new Ribo-seq,
endogenous editing or HLA work: MDA-MB-231 vector/tRNA-Glu/tRNA-Arg crossed
with native EIF4G2-leader WT/uORF1-start mutant, and HeLa untreated/UV crossed
with the same reporter pair. The tRNA screen decides whether the central story
is causal; the UV screen validates the branch and assay in the source UV cell
context. UV success cannot rescue tRNA failure. Exact go/stop criteria and the
spending ladder are in `eif4g2_lab_go_no_go_plan.md`.

### Gate 1 — Independent factorial perturbation

Generate new MDA-MB-231 perturbations in one experiment with vector,
tRNA-Glu(UUC) and tRNA-Arg(CCG), using at least three independently cultured
biological replicates. Verify perturbation abundance and cell state. Avoid
treating stable-line clones as biological replication; if stable lines are
used, include multiple independently derived clones and experiments.

**Pass:** endogenous or reporter allocation changes selectively under
tRNA-Glu, with an estimated perturbation interaction consistent with the public
data and confidence intervals excluding a biologically trivial effect defined
before data collection.

**Kill:** tRNA-Glu does not separate from tRNA-Arg/vector in a technically
adequate experiment, or the result follows clone identity rather than tRNA.

### Gate 2 — uORF1-start dependence

Cross the tRNA perturbations with an EIF4G2 leader reporter retaining the
native downstream GUG and an otherwise matched uORF1 ATG-to-nonstart mutant.
Measure uORF translation and downstream output separately. Analyse the
construct × tRNA interaction.

**Pass:** the tRNA-Glu-specific branch response requires the intact uORF1
start, while baseline reporter integrity and RNA abundance remain acceptable.

**Kill:** the effect persists unchanged after uORF1 initiation is abolished,
indicating that the counted uORF signal is not the functional mediator.

### Gate 3 — Direct-codon versus indirect-route discrimination

Test GAG→GAA synonymous substitutions, non-Glu substitutions designed to
minimize peptide/structure disruption, and mutation of the codon-5 pause
neighbourhood. Quantify reporter RNA and, where possible, RNA structure.

**Pass for direct decoding:** the tRNA-Glu interaction changes predictably with
codon identity while tRNA-Arg and RNA abundance controls remain stable.

**Pass for an indirect model:** uORF-start dependence is real but codon
substitution does not explain it; the paper must then abandon a GAG-sensor
claim and identify the upstream route.

**Kill:** all effects track nonspecific reporter expression, RNA abundance or
construct destabilization.

### Gate 4 — Mature tRNA versus tRNA-derived fragments

Measure mature tRNA-Glu and derived fragments by Northern blot and/or validated
small-RNA methods. Compare wild-type, decoding/anticodon-altered and
cleavage-altered constructs while controlling total expression and charging
where feasible.

**Pass:** one molecular species and construct property tracks the branch
interaction reproducibly.

**Kill for mature-decoding claim:** cleavage products, rather than mature or
charged tRNA, explain the phenotype. This would not kill the branch but would
change the mechanism and title.

### Gate 5 — Endogenous functional output

Measure endogenous EIF4G2 uORF/CDS occupancy and DAP5 protein in the new
perturbations. Then quantify the exact EIF4G2-uORF-derived HLA ligand using
targeted immunopeptidomics or an appropriately validated peptide–HLA assay.

**Pass for translation claim:** endogenous branch redistribution and a
corresponding DAP5 protein consequence are observed.

**Pass for immunogenic claim:** the exact endogenous peptide–HLA product changes
under tRNA-Glu relative to tRNA-Arg/vector and depends on the uORF where the
system permits genetic testing.

**Kill for immune claim:** ribosome occupancy changes without a reproducible
change in peptide–HLA display. The paper can retain a translation branch but
must remove “exposes,” “enhances presentation” and analogous immune language.

### Gate 6 — Raw tRNA study reconstruction

**Status: passed on 2026-08-19.** All 16 public FASTQs were integrity checked,
trimmed with biological insert length retained, uniquely aligned and calibrated
independently by run and length. The lane-matched tRNA-Glu union effect was
−0.136 and the exact-uORF1 effect was −0.135, with disjoint biological-unit
ranges. The study-shared and canonical effects were −0.141 and −0.145. The
identically processed tRNA-Arg effect was +0.012, giving a −0.147 specificity
contrast. Coordinate, single-position and greedy-three-position challenges
passed. The sparse short-read stratum failed its own measurement and target-
phase gates and was not counted as support.

Reprocess the public tRNA-Glu and tRNA-Arg FASTQs with biological insert length
retained, unique alignment and per-run/per-length P-site offsets. Repeat the UV
style exact-length, coordinate and deletion analysis.

**Pass:** the tRNA-Glu-specific effect recurs within adequately counted,
condition-matched length classes and remains absent in tRNA-Arg.

**Kill or downgrade:** the effect is confined to a condition-linked length with
an incompatible P-site or disappears under independently calibrated offsets.

### Gate 7 — Frozen reproducibility package

Release a tagged repository, dependency lock/container, run manifest,
checksums, source-data tables, complete output data dictionary and DOI-minted
archive. Every main and Extended Data panel should rebuild from one documented
command without private manual edits.

### Gate 8 — Reporting and authorship completion

Complete the Nature life-sciences reporting summary, exact author list and
contributions, affiliations, funding, competing interests, ethics statement if
new work requires it, data/code DOIs and the final generative-AI disclosure.

## Experimental priority order

1. Run the small MDA-MB-231 tRNA × WT/start-mutant screen. If the verified
   tRNA-Glu-by-uORF-start interaction fails, stop the mechanism paper.
2. In parallel, run the HeLa UV × WT/start-mutant screen as an orthogonal
   branch/assay validation. Do not make it a condition for the main tRNA claim.
3. Only after the tRNA screen passes, expand Gates 1 and 2 to a powered
   independent experiment and add an uORF-specific readout.
4. In the same reporter backbone, run Gate 3 while building the mature/fragment
   measurements for Gate 4.
5. Confirm the endogenous branch and DAP5 protein before interpreting reporter
   magnitude; proceed to targeted HLA work only after that link is real.

## Language blacklist until the corresponding gate passes

- “tRNA-Glu directly decodes the EIF4G2 uORF”
- “the GAG codons sense tRNA-Glu abundance”
- “tRNA-Glu increases DAP5 expression”
- “tRNA-Glu exposes/increases an EIF4G2 neoantigen”
- “a universal collision/recycling sensor”
- “a metastasis-specific switch”
- “causal” when referring only to reanalysis of the public overexpression data
- “replicated six times” or any wording that treats sequencing lanes as `n`

## Pre-submission go/no-go decision

The central paper survives if independent tRNA-Glu perturbation requires uORF1
and changes at least one endogenous functional output. Direct GAG decoding is
optional; a fragment-mediated or indirect upstream route could be equally
novel if experimentally resolved. The Nature-scale immune framing survives only
if the exact endogenous peptide–HLA output changes. Without that endpoint, the
work can still become a strong translation-mechanism paper, but the title,
summary and cover letter should be narrowed accordingly.
