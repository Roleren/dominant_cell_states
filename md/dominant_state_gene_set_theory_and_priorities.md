# Which dominant states to track, and why: theory, compendium reality, and priority order

Updated: 2026-09-24

## Purpose

`dominant_state_module_definitions.R` currently scores nine states (eight
biological programs plus baseline) and, as of this session, a tenth
(post-viral fatigue). This document answers a question the code does not:
**why these states and not others, and in what order should new ones be
added?** The answer has to satisfy three separate constraints at once, and a
good candidate state should score well on all three, not just one:

1. **Theory**: the state must plausibly be regulated at the level of
   translation -- ribosome allocation between a clean CDS and a uORF/ouORF/
   NTE/iORF branch -- not purely at the level of transcription. This atlas
   cannot say anything interesting about a state that only changes mRNA
   abundance; that is what every other RNA-seq compendium already does.
2. **Nature-worthiness**: the candidate mechanism should be specific enough
   to produce a clean, falsifiable, single-gene-or-small-module story (like
   EIF4G2 uORF1), not just another "this pathway is stress-responsive"
   restatement. Independent literature that could orthogonally validate a
   finding (as Kowar et al. 2025 did for EIF4G2) is a strong bonus, not a
   requirement.
3. **Feasibility**: the compendium must actually contain enough samples in
   the relevant perturbation/condition to detect the effect. A theoretically
   perfect candidate with six samples in the whole compendium is not a
   near-term project.

## Theory: what is actually known to be translationally regulated

Four mechanisms have decades of biochemical characterization behind them as
*bona fide* translational (not transcriptional) control systems, meaning a
uORF/CDS allocation shift in these pathways has strong prior mechanistic
support rather than needing to be discovered from scratch:

**1. The integrated stress response (ISR).** Four kinases (PERK, GCN2, PKR,
HRI) converge on phosphorylating eIF2-alpha in response to ER stress, amino
acid starvation, dsRNA/viral infection, and heme/oxidative stress
respectively. Phosphorylated eIF2-alpha reduces global cap-dependent
initiation while paradoxically *increasing* translation of uORF-containing
transcripts like ATF4, ATF5, and PPP1R15A/GADD34, via a well-characterized
uORF-reinitiation mechanism. This is the single best-characterized
translational switch in mammalian biology and is exactly what this atlas's
clean-CDS/uORF allocation metric is built to detect. It is also already the
project's strongest validated category (EIF4G2, and the viral ISR anchors
DDIT3/PPP1R15A/ATF4).

**2. The mTORC1/TOP-mRNA program.** mRNAs with a 5' terminal oligopyrimidine
(TOP) motif -- overwhelmingly ribosomal proteins and translation factors --
are selectively translated in an mTORC1-activity-dependent manner via
LARP1 and 4E-BP/eIF4E. This is the mechanistic basis of the existing
"mTOR / Translation capacity (TOP program)" module, and it is a second
completely independent, well-established translational control axis.

**3. Ribosome rescue and recycling (RQC/no-go decay).** When a ribosome
stalls or fails to terminate (truncated mRNA, strong secondary structure,
collision), Pelota (PELO)/HBS1L and ABCE1 rescue and recycle it. Loss of this
machinery is mechanistically expected to change the fate of transcripts
carrying uORFs and overlapping ORFs, because failure to recycle after a uORF
directly determines whether the ribosome reinitiates at the downstream CDS.
This is not a hypothesis for this project specifically: it is the literature
mechanism, and this project has already found empirical support for it --
ABCE1 loss is one of the two orthogonal recurrences behind the EIF4G2 uORF1
finding (alongside UV stress). That result was found essentially by accident
while investigating EIF4G2; it justifies treating ribosome
rescue/recycling as its own dominant-state module rather than a one-off
confirmation.

**4. m6A-dependent translational control.** METTL3/METTL14 (writers),
FTO/ALKBH5 (erasers), and YTHDF1-3 (readers) modulate translation efficiency
and can influence start-codon selection and uORF bypass. This mechanism is
newer and less exhaustively characterized than ISR or TOP-mRNA control, which
raises novelty upside but also raises the bar for what would count as a
clean, defensible finding -- a m6A-machinery perturbation that
redistributes ribosome occupancy from a clean CDS into a uORF branch would be
a genuinely new observation, not a confirmation of known biology.

Two more mechanisms are worth tracking as diagnostic, secondary candidates,
not primary modules yet: **IRES-mediated cap-independent translation**
(MYC, VEGFA, and other genes with documented IRES elements continue
translating when global cap-dependent initiation is suppressed -- during
mitosis, apoptosis, or severe stress -- which could show up as an
anti-correlated clean-CDS signal relative to the rest of the transcriptome
during exactly the conditions this atlas already screens), and **stress
granule dynamics** (G3BP1/G3BP2-nucleated mRNP condensation sequesters
untranslated mRNAs; whether this changes measured ribosome occupancy on
specific uORF-containing transcripts differently than on clean-CDS-only
transcripts is untested here and would be a new observation either way).

What this theory section deliberately excludes: pure transcriptional
programs (EMT, senescence/SASP as currently scored) are useful *covariates*
and *context* for interpreting translational findings, but a clean-CDS
allocation shift in an EMT signature gene is not, by itself, evidence of
translational control -- it could simply track total mRNA abundance change.
These modules should stay in the atlas as context/confounder controls, not be
mistaken for primary discovery targets.

## Compendium reality: what ~3,857 human samples actually support

Cross-referencing `all_samples-Homo_sapiens`'s 3,857 runs against
`metadata_done_samples_extended_qc.csv`'s `CONDITION`, `INHIBITOR`, and `GENE`
columns gives a direct, empirical feasibility check for each mechanism above
(all counts are runs within this exact compendium subset, not the full
13,846-row cross-organism sheet):

| Mechanism / perturbation | Approx. samples | Feasibility |
|---|---|---|
| Viral infection (`CONDITION` contains infected/infect/infe) | ~252 | Strong -- already the project's proven category |
| Gene knockdown/knockout, any target (`CONDITION` = KD/KO/DKO) | ~594 | Strong as a general resource; power depends on the specific target gene below |
| Ribosome rescue/recycling (`GENE` = PELO, ABCE1, eIF2A, eIF2D, combined) | ~62 | Moderate, and already mechanistically connected to a real finding (EIF4G2 x ABCE1) |
| m6A machinery (`GENE` = METTL3/14/16/1/5/8, FTO, ALKBH3/5, YTHDF1, combined) | ~44 | Moderate; fragmented across many small perturbation sets rather than one large one |
| Stress granules (`GENE` = G3BP1/G3BP2 and variants, combined) | ~39 | Moderate |
| Amino acid/serum starvation (`CONDITION` contains starv/deprived) | ~62 | Moderate; the classical GCN2 arm of ISR, currently not separated from the general ISR module |
| Hypoxia (`CONDITION` = hypoxia) | 6 | Weak -- existing Hypoxia/HIF1A module cannot be freshly validated from this subset alone; treat as a context covariate, not a discovery target, until more hypoxia runs are added |
| Heat shock (`CONDITION` = HS) | 7 | Weak, same caveat as hypoxia |
| ER stress drugs specifically (tunicamycin/thapsigargin `CONDITION`) | 0 | Absent as a clean perturbation; HSPA5/XBP1 can still be scored as passive markers inside other conditions, but a dedicated UPR discovery push is not supported by this compendium as-is |
| tRNA pool / charging perturbations (`GENE` contains tRNA/synthetase) | 0 | Absent -- confirms the EIF4G2 tRNA-Glu mechanism is not compendium-native; it required going to one specific external study (Goodarzi et al.) and will not generalize to other genes without similarly sourcing external data |

Tissue/cell-type coverage is broad enough for the atlas's existing
context-interaction and transportability layers to be meaningful (kidney,
blood, lung, brain, cervix, breast, bone, liver, embryonic, colon, muscle all
have 60+ samples; cell type is dominated by epithelial at 57%, so any
tissue-specific claim should be checked against that skew before being
trusted as tissue-general).

## Prioritized candidate list

Ranked by the product of theoretical plausibility, nature-worthiness upside,
and compendium feasibility -- not by any one criterion alone.

1. **Ribosome rescue/recycling (PELO/HBS1L/ABCE1) as its own dominant-state
   module.** This is the highest-priority *new* addition. Theory is strong
   (direct mechanistic link between recycling failure and uORF reinitiation
   fate), feasibility is real (~62 samples, not huge but genuine, and
   concentrated in a small number of clean genetic perturbations rather than
   diffuse conditions), and it is not starting from zero -- it already has
   one confirmed hit (EIF4G2 x ABCE1) that this module would be built to
   generalize beyond a single gene. Should be added to
   `dominant_state_module_definitions.R` as e.g. `"Ribosome rescue / RQC"`,
   scored from PELO/ABCE1/HBS1L knockdown-or-loss samples versus matched
   controls, and cross-referenced against the existing atlas cards to see if
   other genes besides EIF4G2 show the same uORF-ward redistribution under
   the same perturbation class.
2. **Continue post-viral fatigue / ISR as the flagship discovery track.**
   Already covered in depth in `dominant_cell_state_iteration_review_2026-09-24.md`;
   restated here only to note that it sits at the top of this list too, on
   pure theory-times-feasibility grounds, independent of how much work has
   already gone into it.
3. **m6A-machinery-dependent translational routing.** Higher novelty upside
   than ISR/RQC (less exhaustively characterized as a uORF-routing
   mechanism specifically), moderate feasibility (~44 samples spread across
   many small perturbations, which will limit per-target statistical power
   but is enough for an exploratory screen). Worth a dedicated
   `dominant_rdg_*` audit script analogous to the EIF4G2 investigation,
   explicitly looking for genes whose clean-CDS/uORF allocation shifts under
   METTL3/METTL14 loss or FTO/ALKBH5 perturbation.
4. **Stress granule dynamics (G3BP1/G3BP2).** Similar profile to m6A:
   moderate feasibility (~39 samples), real but less-established
   translational mechanism specifically for uORF routing (most G3BP
   literature is about mRNA sequestration/condensation, not uORF
   reinitiation per se, so a positive result here would be more novel but
   also needs a more careful mechanistic argument before claiming it is a
   "translation" effect rather than an "mRNA availability" effect).
5. **Amino-acid/serum starvation as an ISR sub-axis, separated from the
   general ISR module.** ~62 samples is enough to ask whether the GCN2 arm
   of ISR (starvation) produces the same uORF-routing signature as the
   PERK/PKR arms (ER stress/viral) already captured in the general ISR
   module, or whether it is mechanistically distinguishable. This is a
   refinement of an existing module rather than a new one, and is cheap to
   test now that the fatigue module's gene set is restored (see below).
6. **IRES-dependent survival translation as a diagnostic overlay, not a
   standalone module yet.** No dedicated perturbation exists in the
   compendium to isolate this cleanly; the practical path is to check
   whether known IRES-containing genes (MYC, VEGFA, XIAP) already show
   atypical clean-CDS behavior in the stress/mitotic conditions the atlas
   already screens, as a retrospective look rather than a new data
   acquisition target.
7. **Do not invest further in hypoxia- or heat-shock-specific discovery
   from this compendium alone.** Both mechanisms are theoretically sound
   translational-control candidates (HIF-dependent internal ribosome entry,
   heat-shock-dependent uORF usage in HSPA-family transcripts are both
   documented), but 6 and 7 samples respectively cannot support a
   compendium-scale claim. Keep both as context/covariate modules only
   unless more runs are added.
8. **Do not expect the tRNA-pool/codon-optimality mechanism (the EIF4G2
   story's own category) to generalize compendium-wide.** Confirmed absent
   from this compendium's own perturbation set. Any follow-on tRNA-pool
   story will need to source its own external dataset the way the Goodarzi
   study did for EIF4G2, not rely on `all_samples-Homo_sapiens` alone.

## A concrete process gap this inspection surfaced

While cross-checking module #2 above, `results/curated_inputs/` did not exist
on this server at all -- meaning `postviral_fatigue_dominant_state_genes.csv`
(the actual gene list behind the post-viral fatigue module) was missing, and
`load_dominant_state_modules()` was silently returning the module with zero
genes (matching the `"No post-viral fatigue genes found for scoring tiers:
conservative"` warning seen in pipeline logs this session). This file has
been reconstructed from the conservative-tier gene list documented in
`postviral_fatigue_translation_note.md` and placed back at
`results/curated_inputs/postviral_fatigue_dominant_state_genes.csv`
(26 up-genes, 17 down-genes, verified to load correctly). The same is true in
principle for `expanded_oxphos_marker_genes.csv` and
`curated_additional_gene_priorities.csv`, though those two degrade
gracefully to a smaller hardcoded default rather than silently zeroing out an
entire module, so they are lower urgency. The general lesson: anything under
`results/curated_inputs/` is actually *curated source content*, not a
regenerable pipeline output, despite living inside the gitignored `results/`
tree -- it should be tracked in git (or moved outside `results/`) so it
cannot silently vanish on a fresh checkout the way it did here. See
`dominant_cell_state_iteration_review_2026-09-24.md` for the broader pattern
of curated-content-living-in-a-gitignored-directory issues found this
session.
