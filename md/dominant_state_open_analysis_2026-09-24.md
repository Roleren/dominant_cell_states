# Open analysis of the dominant-state results: what the numbers actually show

Updated: 2026-09-24

## Purpose

This is a free exploration of the just-completed full-server atlas run, not
tied to a specific story. It asks: looking at the real per-sample score
matrix, the score audit, the uORF structure-rule outputs, the transport
robustness table, the pooling-bias audit, and the auxiliary-state screen
with fresh eyes, what shows up? Seven findings below are substantial enough
to change how parts of the atlas should be read; the rest are smaller,
concrete observations and one piece of analysis that didn't exist and has
been added (the code is in this note; it should be promoted to a proper
script if it holds up on review).

## Finding 1: the raw "dominant state" winner label is not comparable across states, and misrepresents the compendium

Across all 3,857 samples, the raw per-state scores in
`human_dominant_cell_states_clean_cds.csv` have wildly different means:
Interferon averages **-6.00**, ISR **-4.29**, Proliferation **-3.75**, while
OXPHOS is the only state with a **positive** mean (**+0.98**). Since
`dominant_state` is computed as an argmax across these raw scores
(`dominant_cell_states_fst_clean_cds.R:879-885`, `.SD[which.max(score)]`),
this guarantees OXPHOS wins almost by construction: **OXPHOS is the "winner"
for 2,447 of 3,857 samples (63.4%)**, while ISR -- the atlas's own flagship,
best-validated state -- wins for **exactly one sample**, and Interferon and
Post-viral fatigue win for **zero**.

The mechanism is in the scoring code itself:
`gene_z <- sweep(expr_norm, 2, baseline, "-")` (`dominant_cell_states_fst_clean_cds.R:857`)
subtracts a **per-sample** housekeeping-gene reference level from every
gene, then each state's score is the mean of its own marker genes' values in
that same units. This makes a state's score reflect "expression level
relative to housekeeping genes" -- and OXPHOS's marker genes (mitochondrial
complex subunits) are, in most cell types, expressed at levels comparable to
or above average housekeeping genes as a matter of routine cell biology, not
because a sample is in a distinctive "OXPHOS state." Interferon- and
ISR-target genes are, in the overwhelming majority of samples that are not
actively infected or stressed, expressed well below housekeeping level, so
they structurally lose the argmax contest even in the specific samples
where they are, in fact, the most *activated relative to their own normal
range* -- which is the question a "dominant state" label should answer.

**This is not a new problem this note is introducing -- it was already
independently found and partially addressed.** `dominant_state_score_audit.R`
computes a percentile-rank calibration (`percentile_z`, by `dominant_state`)
and reports both a `Raw signed score` winner count (matching the numbers
above exactly) and a `Within-state percentile` winner count. That
calibrated count is far more balanced: OXPHOS drops to 636 (16.5%), ISR
rises to 293, Interferon to 353, and Post-viral fatigue to 139 -- all in
the same 8-16% range. An independent z-score-based recalibration, computed
fresh for this note directly from the raw score matrix, agrees closely
(OXPHOS 508, ISR 311, Interferon 329, Fatigue 123, and every state landing
in a 300-550-sample range) despite using a different calibration method.
**Two independent methods agreeing this closely is good evidence the
correction is real, not an artifact of either method's specific choices.**
Overall, **70.8% of all 3,857 samples get reassigned to a different
"dominant state" once calibration is applied.**

**Correction (verified 2026-09-24, after an initial fix attempt that turned
out to target the wrong file):** the first pass at scoping this grepped for
per-sample filtering by the raw `dominant_state` label and found (and
patched) `case_dominant_call_fraction`/`control_dominant_call_fraction` in
`dominant_cell_state_control_aware_enrichment.R`. Running the full pipeline
end-to-end to verify that patch showed the patch had **no effect** on the
actual output. Tracing why revealed the real story is bigger than first
scoped: `summarise_direct()` in that same script never propagates
`case_dominant_call_fraction`/`control_dominant_call_fraction` out of the
raw per-stratum `contrasts` table into `condition_summary`/`gene_summary` --
so patching them only affected an internal column that gets written to
`dominant_state_control_aware_condition_contrasts.csv` and then dropped, not
anything downstream. The `dominant_call_fraction_in`/`_out`,
`dominant_calls_in`/`_out`, and `fisher_odds_ratio`/`fisher_greater_p`
columns actually visible in the final review/summary tables come from a
**different, currently-active pipeline script**:
`dominant_cell_state_adjusted_enrichment.R`
(`enrich_one_field_adjusted()`, `call_in <- sum(dt$dominant_state[in_level]
== state, ...)`), which runs one step earlier in the pipeline and feeds
`dominant_state_clean_cds_enrichment_all_main.csv` /
`_cell_line.csv` / `_condition.csv` / `_gene.csv` / `_tissue.csv` /
`_interactions.csv` and their heatmaps -- the whole "adjusted enrichment"
analysis layer, not one diagnostic pair in one script. That script's Fisher
test (`safe_fisher_greater(call_in, no_call_in, call_out, no_call_out)`) is
an actual significance test on raw-label calls, not just a descriptive
fraction, so this was silently distorting more than headline percentages:
for low-baseline-rate states (ISR, Interferon, Fatigue), `call_in`/`call_out`
were both ~0 almost everywhere, making that Fisher test underpowered and
its enrichment/depletion calls across CELL_LINE/TISSUE/GENE/CONDITION
largely uninformative for exactly those states, while OXPHOS's
near-universal "win" produced enrichment signal driven mostly by which
levels happen to have slightly more or fewer samples, not by real biology.
**This has now been fixed at the actual source** (both
`dominant_cell_state_adjusted_enrichment.R` and, redundantly but harmlessly,
`dominant_cell_state_control_aware_enrichment.R` now merge in
`human_dominant_state_activity_metadata.fst`'s calibrated
`dominant_state_top` and use it for all call-based columns), and verified
directly: before the fix, `dominant_call_fraction_in` averaged 0 for ISR,
Interferon, and Fatigue and 0.70 for OXPHOS across
`dominant_state_clean_cds_enrichment_all_main.csv`; after the fix, every
state falls in a realistic 0.06-0.18 range (OXPHOS 0.18, mTOR 0.14, ...,
Fatigue 0.057, Interferon 0.066, ISR 0.084) with no state at exactly zero.

**Two more instances surfaced from a systematic follow-up grep** (searching
every production script for a raw `dominant_state == <something>`-style
comparison, not just the two already found), and both are bigger than the
enrichment tables:

1. **`dominant_ribosome_decision_graphs.R`, the RDG atlas builder itself.**
   `build_branch_terms()` always constructs exactly 9 fixed
   `"sample_dominant_call"` branches, one per state (`branch_source =
   "sample_dominant_call"`, `field = "dominant_state_call"`) -- these are
   not conditional on any enrichment significance test, they are always
   present. `branch_runs()` then matched each branch's membership via
   `context$dominant_state == branch$dominant_state`, i.e. the same raw
   argmax. In practice this meant the OXPHOS `sample_dominant_call` branch's
   "in" group was really just **63% of the entire compendium** (not an
   OXPHOS-specific subset), while the ISR/Interferon/Post-viral-fatigue
   `sample_dominant_call` branches had ~0-1 "in" samples each -- almost
   certainly falling below `min_branch_samples` and silently vanishing from
   the atlas rather than appearing as obviously wrong. Any per-feature
   branch-level score attributed to "OXPHOS, sample-dominant-call" in the
   atlas was therefore mostly measuring the population average dressed up
   as state-specific evidence. Fixed the same way, in `build_sample_context()`
   and `branch_runs()`.
2. **`postviral_fatigue_pipeline_summary.R`, the fatigue module's own
   dedicated report.** `winner_counts <- state_dt[, .N, by =
   dominant_state]` fed directly into a headline metric,
   `fatigue_winner_samples`, in that report's own `summary_metrics` table --
   meaning the fatigue module's own summary document was reporting (or very
   close to reporting) the same near-zero raw winner count this note
   already flagged as misleading in Finding 1, inside the one document
   whose entire purpose is to argue the fatigue state is real and worth
   keeping. Fixed the same way.

Both fixes follow the identical pattern used for `adjusted_enrichment.R`:
merge in `human_dominant_state_activity_metadata.fst`'s calibrated
`dominant_state_top` and use it wherever raw-label membership/counting is
needed. The corresponding pipeline dependency declarations
(`pipeline_step_dependencies()` in `run_dominant_cell_state_pipeline.R`) were
updated to include `outputs$state_activity_metadata` for `rdg_graphs` and
`adjusted_enrichment` (matching the earlier fix already applied to
`control_aware`), so cache invalidation stays correct if the calibrated file
ever changes.

**Verified end-to-end**, via a full `--force-step rdg_graphs` pipeline run:
`rdg_graphs` and `fatigue_summary` both completed without error, and the
resulting numbers are exactly the expected calibrated distribution.
`rdg_branch_feature_scores.csv`'s `sample_dominant_call` branches now show
real, non-degenerate "in"-group sample counts for every state -- OXPHOS 636,
mTOR 568, Proliferation 527, Senescence 487, EMT 453, Hypoxia 401,
Interferon 353, ISR 293, Post-viral fatigue 139 -- instead of OXPHOS's
former ~2,447 (63% of the whole compendium) and three states at ~0. The
fatigue report's own `fatigue_winner_samples` metric now reads **139**
(matching the calibrated count exactly) instead of the near-zero raw
number.

Four locations fixed in total: `dominant_cell_state_control_aware_enrichment.R`,
`dominant_cell_state_adjusted_enrichment.R`, `dominant_ribosome_decision_graphs.R`,
`postviral_fatigue_pipeline_summary.R`. There is also already a fully
correct, calibrated per-sample winner sitting in
`human_dominant_state_activity_metadata.fst`
(`dominant_state_activity_metadata_fst.R`) that all four fixes now draw
from. The durable lesson from this whole exercise, worth stating plainly:
**a grep-based scoping pass that finds the right symptom does not
guarantee it finds every place the symptom's *output* actually surfaces**,
and neither does a first follow-up grep -- verifying a fix means running
the real pipeline and checking real numbers change in the expected
direction, twice over in this case, not trusting that grep found everything
the first time. Before trusting any future call-based summary column in
this atlas, check whether its source script computes a raw
`dominant_state == state`-style comparison anywhere in its own code, not
just whether the specific column name looks familiar from this note.

## Finding 2: the newly-restored Post-viral fatigue state is not yet demonstrated to be independent of Interferon/ISR -- but it isn't purely definitional either

This matters because the fatigue module's gene list was reconstructed this
session (it did not exist on this server before), so its validity had not
yet been checked against real data. Its correlation with other states
across all samples: Interferon **r=0.80**, ISR **r=0.57**, Hypoxia
**r=0.43**, OXPHOS **r=-0.36**. An r=0.80 correlation with an existing
state is a real concern for whether "Post-viral fatigue" earns its place
as a distinct dominant state rather than a relabeling of Interferon.

Checking the gene lists explains most, but not all, of this: **7 of
Interferon's 8 up-genes** (CXCL10, IFIT1, IFIT3, ISG15, MX1, OAS1, STAT1)
and **9 of ISR's 10 up-genes** (ASNS, ATF3, ATF4, DDIT3, DDIT4, HERPUD1,
HSPA5, PPP1R15A, XBP1) are directly reused verbatim in the fatigue module's
26-gene up-list. Only 10 of the 26 up-genes are unique to fatigue (CCL2,
GDF15, IFI27, IFITM1, IL6, IRF7, NFKBIA, SERPINE1, TNFAIP3, TRIB3) and none
of the 17 down-genes overlap Interferon or ISR (which have no down-arm at
all). A high correlation driven mostly by literal shared marker genes would
make "fatigue is a distinct state" a weak claim on its own.

To check whether anything survives beyond that overlap, this note computed
(fresh, not previously in the pipeline) a **residual fatigue score using
only the 10 non-shared up-genes and the 17 down-genes**, via the same
log2-then-z-score-then-mean approach the pipeline itself uses, from
`human_dominant_cell_states_clean_cds_gene_expression.csv`. That residual
score still correlates with the full official fatigue score at r=0.94 (so
it is a reasonable proxy), and its own correlations with the other states
are: Interferon **r=0.69**, ISR **r=0.47**, OXPHOS **r=-0.39**, Hypoxia
**r=0.45**. The correlation drops modestly (0.80 to 0.69) but stays large
even using genes that share *zero* names with Interferon or ISR. **This
means the correlation is not purely definitional** -- CCL2, IL6, NFKBIA,
TNFAIP3 (classic NF-kB/inflammatory response genes) and IRF7/IFI27/IFITM1
(interferon-adjacent but not literally Interferon-module genes) genuinely
do co-vary with interferon/ISR activation in this compendium, which is
exactly what the post-viral-fatigue literature hypothesis predicted
(chronic immune/stress co-activation, not either alone).

The honest characterization, then: **Post-viral fatigue, as currently
scored, functions less like an independent tenth axis and more like a
specific, NF-kB/inflammatory-flavored refinement of the Interferon/ISR
axis, whose most distinguishing feature relative to plain Interferon is a
somewhat stronger anti-correlation with OXPHOS** (-0.39 vs Interferon's own
-0.21). That is a real, falsifiable, more precise claim than "fatigue is a
new dominant state," and it is the claim this atlas's actual data supports
right now. Whether that refinement is worth a dedicated state or should be
folded into Interferon as a sub-score is a genuine open question -- not
resolved here, but now precisely statable instead of assumed.

## A consensus-candidate observation worth a closer look: LDHA under ISR, not hypoxia

The top-ranked row in `dominant_rdg_consensus_candidates.csv` by consensus
score is not one of the usual viral-ISR anchors -- it is **LDHA, under the
`isr_er_translation_stress` design family**, with 7 supporting studies at
`B_replicated_branch_consensus` tier, and a clean-CDS change of **-13.9%**
(a *decrease*). LDHA (lactate dehydrogenase A) is a canonical HIF1A/hypoxia
transcriptional target, not a typical ISR gene, and its atlas row here is
specifically filed under ISR/ER/translation-stress conditions, not hypoxia,
with translation going *down* rather than the up-regulation that would be
expected from HIF-driven transcriptional induction. TFAM shows a related,
smaller pattern -- present under both `isr_er_translation_stress` (+3.1%)
and `nutrient_starvation` (+1.9%), both modest but positive, both
well-replicated (7 and 12 studies).

The most parsimonious reading is not that LDHA is being transcriptionally
repressed -- HIF-driven LDHA induction under hypoxia is extremely
well-established -- but that this atlas may be catching an **uncoupling
between transcriptional induction and translational output specifically
under ISR/translation-stress conditions**: LDHA mRNA may still be induced
at the transcript level while its *ribosome allocation* is redirected away
from the clean CDS when the ISR is simultaneously active, i.e. a
translational throttle layered on top of transcriptional induction. That
is exactly the kind of buffering/compensation signature this atlas's
`allocation_output_disagreements` layer is built to catch (see
`dominant_cell_state_iteration_review_2026-09-24.md` and
`uorf_to_cds_effects_and_distance_patterns.md`), but LDHA does not
currently appear in that specific table -- it may be worth adding LDHA to a
future active-review batch specifically to check whether its clean-CDS
output and total mRNA-level output move in opposite directions under these
conditions, since if so, this would be a clean, mechanistically specific,
and so far unflagged example of transcription/translation uncoupling in the
atlas's own top-ranked candidate.

## Finding 3: six genes have far less clean-CDS usage than their uORF count predicts -- DDIT3 (the atlas's own #1 flagship gene) most of all

`dominant_uorf_structure_rules` fits a simple count-only model
(`expected_clean_cds_usage_from_uorf_count`, from
`uorf_structure_multivariable_coefficients.csv`'s `count_only` term:
each additional uORF predicts roughly a fixed drop in clean-CDS usage) and
flags genes whose *actual* clean-CDS usage falls well below what their
uORF count alone would predict, on the logic that a residual this large is
more likely to mean "this gene has an uORF-like repressive element that
isn't in the annotation yet" than "this gene is simply an outlier." Six
genes clear the high/medium bar in
`uorf_structure_hidden_uorf_candidates.csv`, ranked by the size of the gap:

| gene | uORFs (all leader) | actual clean-CDS usage | count-model expected | gap |
|---|---|---|---|---|
| **DDIT3** | 1 | 15.3% | 97.9% | **-82.6 pts** |
| ZEB2 | 3 | 16.9% | 83.9% | -67.0 pts |
| TNFAIP3 | 2 | 27.4% | 93.9% | -66.5 pts |
| RRM1 | 4 | 23.6% | 63.8% | -40.2 pts |
| VEGFA | 4 | 28.0% | 63.8% | -35.8 pts |
| STAT1 (medium) | 5 | 19.1% | 37.3% | -18.1 pts |

None of these six have a predicted overlapping uoORF -- the model's
"expected" number assumes their uORFs behave like the typical
leader-only, non-overlapping case (see `uorf_structure_rule_interpretation.csv`
row 4: "Low clean-CDS usage in genes without predicted uoORFs should be
treated as an annotation-quality signal... A low clean-CDS residual is not
proof of a missed uoORF; it could also reflect RNA structure, poor uORF
prediction, alternate transcript choice, or other translated leader
features"). That caution matters most for the top hit specifically:

**DDIT3 (CHOP)** is the single largest residual in the entire table, and it
is also this atlas's own #1 consensus-ranked flagship gene (see
`dominant_rdg_consensus_candidates.csv` discussion above and
`uorf_architecture_score_genetics_case.md`). But a literature check
([eLife 2015, "Translation of 5' leaders is pervasive in genes resistant
to eIF2 repression"](https://elifesciences.org/articles/03971);
[PubMed 26817837, "Ribosome Elongation Stall Directs Gene-specific
Translation in the Integrated Stress Response"](https://pubmed.ncbi.nlm.nih.gov/26817837))
shows CHOP's single leader uORF is already known to be mechanistically
unusual: translational control here does not work by the standard
"short uORF -> ribosome reinitiates downstream" logic the count-only model
implicitly assumes, but by **elongating-ribosome stalling within the uORF
itself**, which blocks reinitiation at the main CDS under low-ISR
conditions and is relieved (via eIF2alpha phosphorylation increasing leaky
scanning) under high-ISR conditions. A single-uORF gene built on an
elongation-stall mechanism would be expected to show *much* lower baseline
clean-CDS usage than a generic "one short uORF, simple reinitiation" gene
with the same uORF count -- exactly the direction and rough magnitude of
DDIT3's residual here. So DDIT3's flag is likely **not** a missed second
uORF so much as **the count-only model's core assumption (all uORFs behave
like generic short reinitiation-permissive elements) breaking down for a
mechanistically atypical, well-characterized case** -- a genuine finding
about the model's limits, not necessarily an annotation gap to fix by
adding a uORF.

The other five have no comparably specific literature mechanism checked
yet and remain plausible annotation-gap candidates in the ordinary sense
(unannotated overlapping/near-CDS uORF, alternate 5'UTR isoform, or
structured-RNA reinitiation block) -- worth prioritizing for manual
browser review in roughly this order, since **TNFAIP3** is one of the 26
genes just added via the Post-viral fatigue module (so cleaning up its
structure would directly firm up a module this session already flagged as
needing independence evidence, see Finding 2), and **RRM1**/**VEGFA** are
both consensus-tier RDG genes with active therapeutic-perturbation
predictions downstream, meaning a corrected clean-CDS baseline here would
propagate to already-published-looking numbers elsewhere in the atlas.

## Finding 4: in the transport-model-agreement table, DDIT3 and PPP1R15A are more robust than ATF4 itself, and only under viral infection

`dominant_rdg_model_agreement_transport_gene_design.csv` runs a
leave-one-axis-out robustness check per gene x design-family row
(`transport_stability_score`, with the supporting/fragile axis broken out
by name: cell line, tissue, footprint length, protocol, author,
publication). Restricting to DDIT3, PPP1R15A, and ATF4 -- the three core
ISR effectors this atlas already treats as its central story -- surfaces a
pattern worth naming: the two **downstream** effectors have the highest,
most broadly cross-axis-robust consensus scores of the three, and only in
one design family, while ATF4 itself -- nominally the pathway's master
regulator -- never reaches that top tier anywhere.

| gene | design family | consensus class | stability | fragile on |
|---|---|---|---|---|
| PPP1R15A | viral_infection | support_robust | **0.935** | (none) |
| DDIT3 | viral_infection | support_robust | **0.910** | (none) |
| ATF4 | viral_infection | direction_robust | 0.780 | footprint length, tissue |
| ATF4 | hypoxia_oxidative_stress | direction_robust | 0.760 | footprint length |
| ATF4 | isr_er_translation_stress | direction_robust | 0.590 | protocol |
| ATF4 | nutrient_starvation | **fragile** | 0.180 | cell line, footprint length, tissue |

DDIT3 and PPP1R15A's `viral_infection` rows are robust against every axis
tested (cell line, footprint length, *and* tissue context all listed as
"support_robust", none fragile) -- about as strong a consensus signal as
this table produces anywhere. ATF4, by contrast, tops out at
"direction_robust" (the LOO check agrees on the *sign* of the effect but
not reliably on its *magnitude/support*) in every design family it appears
in, and is outright fragile under nutrient starvation -- arguably the
textbook condition for ATF4 induction via uORF2 reinitiation. All three
genes are `non_consensus_review` (stability 0) in every other design
family shown (tumor/cancer context, hypoxia for DDIT3/PPP1R15A, DNA
damage/senescence, etc.), meaning this robust signal is specific to viral
infection, not a general ISR-wide property of these genes.

Read cautiously, this does not mean ATF4's biology is less real or less
central -- ATF4 protein induction is one of the best-characterized events
in the ISR literature, and a single gene's uORF-reinitiation control is
mechanistically more von Neumann-sensitive to exact conditions (start-codon
context, uORF spacing, eIF2-GTP levels) than a downstream transcriptional
target's induction, which could plausibly make ATF4's *ribosome-profiling*
signal noisier across studies even when the underlying pathway activation
is completely consistent. The more defensible reading is narrower: **in
this specific compendium, at the ribo-seq clean-CDS-usage level, DDIT3 and
PPP1R15A currently make a more reproducible case for a robust
viral-infection-specific transport effect than ATF4 does** -- which is
useful to know before treating "ATF4 up" as the most citable single-gene
claim from this atlas's viral-infection story; DDIT3 or PPP1R15A may
actually be the safer flagship gene for that specific design family,
with ATF4 kept as the mechanistic explanation rather than the most robust
individual data point.

## Finding 5: two independent methods agree that most explainable dominant-state variance is cell-line/batch identity, not treatment -- and this varies sharply by state

This combines an unsupervised result (PCA + clustering directly on the 9
state scores, `dominant_state_latent_program_*`) with a supervised one
(regularized regression of each state score against all available
metadata, `dominant_state_glmnet_metadata_*`) that were computed completely
independently of each other and of this note, and that agree closely.

**The unsupervised side first, because it's also good news for the module
design.** PCA on the 3,857 x 9 state-score matrix finds clean, interpretable
top axes: **PC1 (38.0% of variance)** is a coherent "stress/damage vs.
growth-biosynthesis" contrast -- Post-viral fatigue, Senescence, Interferon,
and Hypoxia all load positively, while mTOR, Proliferation, and OXPHOS load
negatively, exactly the grouping a biologist would draw by hand. **PC2
(20.1%)** contrasts Proliferation/Hypoxia/mTOR/ISR against
Senescence/OXPHOS/EMT, a plausible cell-cycle-active vs.
arrested/differentiated axis. Together the top two PCs already capture 58%
of all variance in just 9 numbers -- the modules are not scoring 9
independent, arbitrary things, they cohere into a small number of
real biological axes.

But when samples are then clustered in this same latent space
(`dominant_state_latent_program_clusters.csv`, 10 clusters), the clusters
that fall out are almost entirely explained by **which cell line or study a
sample came from**, not by treatment. Seven of the ten clusters have their
single strongest metadata association be a specific `CELL_LINE`, `study`,
or `AUTHOR` term, several with enormous odds ratios: `CELL_LINE=LCL`
OR=2050 (q=5.6e-244) for the 241-sample cluster LC1; `study=PRJNA707431`
OR=Inf (q=1.1e-82) for the 614-sample LC6 (all NSCLC/A549/Luan-authored);
`CELL_LINE=HBEC` OR=Inf (q=2.2e-47) for LC10. In plain terms: a large
fraction of this atlas's unsupervised structure is "which lab/cell-line
generated this sample," not "what state is this sample in."

**The supervised side quantifies exactly how much.**
`dominant_state_glmnet_metadata_model.R` fits a separate regularized
(elastic-net) regression per state, predicting the state score from every
available metadata field, field interaction, and its own stability
selection across resamples, and classifies each selected term into one of
four roles: `perturbation` (a condition/treatment effect that generalizes
across identity, e.g. "infected" as a main effect), `biology` (pure
identity, e.g. a specific cell line or tissue), `biology_perturbation` (a
condition effect specific to one identity context, e.g. "tumor in kidney"),
or `nuisance_or_design` (batch/protocol/author/technical terms). Across all
9 states combined (2,839 selected terms total), weighting by effect
magnitude (`coefficient_abs`): **`nuisance_or_design` accounts for 49.0% of
total explanatory magnitude, `biology` (identity) for 26.2%,
`biology_perturbation` for 10.4%, and pure, generalizable `perturbation`
terms for only 14.4%.** By simple term count the split is similar (44.7% /
21.0% / 18.0% / 16.3%). Put plainly: **across this whole atlas, only
roughly one-seventh of what predicts a dominant-state score is a
treatment effect that would be expected to replicate in a new cell line**;
the rest is technical/batch structure or context-specific identity.

**This varies sharply and informatively by state**, which is the most
actionable part:

| state | nuisance_or_design | perturbation |
|---|---|---|
| Proliferation / Cell cycle | **74.2%** | 6.3% |
| Interferon / antiviral | **66.5%** | 8.0% |
| ISR / ER stress (ATF4 axis) | 56.6% | 14.1% |
| mTOR / Translation capacity | 48.1% | 15.0% |
| Hypoxia / HIF1A program | 41.3% | 18.8% |
| Post-viral fatigue / ribosome stress | 41.6% | 20.1% |
| OXPHOS (mitochondrial) | 40.6% | 17.6% |
| EMT / Mesenchymal shift | 34.1% | 16.0% |
| Senescence / SASP | 31.0% | **21.3%** |

Two results stand out as genuinely surprising. **Proliferation** is the
single most batch-confounded state in the atlas (74.2% nuisance,
6.3% perturbation) -- its score is barely predictable from what was done to
a sample beyond which batch/author produced it, which is a real caution
against citing "proliferation increased under X" claims from this atlas
without checking this table first. **Interferon** is nearly as bad (66.5%
nuisance, 8.0% perturbation) -- surprising for what should be the most
directly infection-driven state in the panel; the perturbation terms that
*do* survive are biologically sensible (`FRACTION=SARS2`,
`CONDITION=cytokines`), so the signal that exists is real, there is just
comparatively little of it relative to batch structure. By contrast,
**Senescence** (21.3%) and **Post-viral fatigue** (20.1%) have the highest
treatment-attributable fractions in the panel -- an independent, and
reassuring, piece of evidence for the fatigue module's validity alongside
the residual-correlation check in Finding 2, since it means fatigue's score
is comparatively *more* driven by real condition effects than most of the
established states it was checked against, not less.

**Practical takeaway:** before citing any "state X is elevated/depressed
under condition Y" claim from this atlas, check `feature_role` for that
state/term combination in
`dominant_state_glmnet_metadata_selected_term_scores.csv`. For Proliferation
and Interferon specifically, treat such claims as needing direct replication
across multiple cell lines/studies before trusting them at face value; for
Senescence and Post-viral fatigue, the atlas's own regression evidence is
comparatively more supportive of a real, generalizable treatment effect.

## Finding 6: the uORF "on/off switching" catalog looks like broad multi-gene evidence but is actually 43% one single osmotic-stress study

`dominant_ouorf_cds_coupling_analysis.R` builds a catalog of genes whose
uORF translation appears to switch qualitatively on/off between case and
control (as opposed to a graded quantitative shift) --
`dominant_ouorf_cds_coupling_onoff_summary.csv` has 28 distinct
gene/uORF-feature rows that clear its `ouorf_onoff_review` evidence bar,
spanning genes as varied as IGFBP5, TRIB3, BNIP3, CDK1, DDIT4, MYC, RPLP0,
TFAM, LMNB1, SNAI2, TBP, MKI67, SLC2A1, ATF3, IFIT1, and EGLN3. Read at
face value, a list spanning that many unrelated genes and pathways looks
like broad evidence for a general phenomenon (uORFs frequently switch
on/off rather than shift gradually).

Checking `strongest_onoff_context` (which study/cell-line/condition
produced each row's strongest evidence) tells a narrower story: **12 of the
28 rows (42.9%) -- including 12 of the top 20 by evidence score (60%) --
all trace to the exact same single experiment: `PRJNA976655-homo_sapiens`,
HEK293 kidney cells, `NaCl_wash` (hyperosmotic stress) vs. WT.** IGFBP5,
TRIB3, BNIP3, CDK1, DDIT4, MYC (x2 features), TBP, MKI67, and SLC2A1 all
show their single strongest on/off evidence in this one study. The second
most-represented study, `PRJNA472636` (MCF10A breast, nutrient starvation),
contributes another 5 rows. So roughly two-thirds of this catalog's
apparent gene diversity collapses to two single-study, single-perturbation
sources.

This does not mean the phenomenon is fake -- it is entirely plausible that
acute hyperosmotic stress genuinely triggers unusually clean, fast,
binary-looking uORF engagement switches across many genes at once in this
one well-powered experiment (a real, single-condition biological effect,
consistent with how sharply and rapidly osmotic stress is known to remodel
translation). But it does mean the catalog should not be read as "N
independent genes across the compendium show uORF on/off switching" --
it should be read as "one deeply-sampled osmotic-stress experiment
supplies most of the evidence, plus a smaller nutrient-starvation
contribution, with only a handful of rows independently supported
elsewhere." This is the same underlying pattern already flagged for
`tumor_cancer_context` in the pooling-bias audit (large apparent row counts
masking a narrow, non-independent evidential base) and for the lung/viral
"off-frame" cross-reference this session initially over-read before
checking -- a third, independently-discovered instance of the same
meta-lesson: **row count and gene-name diversity in this atlas's tables are
not a reliable proxy for the number of independent experiments actually
supporting a pattern; always check `strongest_onoff_context` (or the
equivalent study/context column) before citing breadth.** Concretely, any
future claim built from this table naming more than one or two of
IGFBP5/TRIB3/BNIP3/CDK1/DDIT4/MYC/TBP/MKI67/SLC2A1 as independent evidence
for a general on/off-switching phenomenon should first confirm they are not
all restatements of the same PRJNA976655 NaCl-wash experiment.

## Finding 7: the atlas's own held-out cross-validation shows its flagship hierarchical model predicts the *direction* of a gene's change worse than a coin flip, for most design families

`dominant_rdg_mechanistic_motif_model.R` re-evaluates the pre-existing
`hierarchical` effect model's already-computed **leave-one-study-out**
held-out predictions (`dt$heldout_log2fc` vs `dt$predicted_log2fc`,
`n_studies` = number of held-out studies), and separately fits two new
models of its own (`motif_model`, `motif_design_model`, incorporating
uORF/sequence-motif features) evaluated by a *different*, stricter scheme:
**leave-one-gene-out** (`block_col = "gene_symbol"` in
`fit_blocked_direct_model()`), i.e. predicting a gene the model has never
seen at all, from generic motif/context features alone. Both schemes report
a `direction_accuracy` = the fraction of held-out predictions whose *sign*
(up vs. down) matches the observed sign of the actual log2FC -- a metric
where 50% is what random guessing would achieve.

**Aggregated across all 11,808 held-out predictions, the flagship
`hierarchical` model gets the direction right only 37.3% of the time --
meaningfully worse than a coin flip** -- and even the best-performing
variant, `motif_model`, reaches only 46.1%, still below chance. Broken down
by design family, the picture is not uniformly bad, which is the important
part:

| design family | n | hierarchical | motif_model | motif_design_model |
|---|---|---|---|---|
| acute_interferon_cytokine | 328 | 25.0% | 26.2% | 25.9% |
| other_non_control | 1,886 | 29.9% | 37.5% | 37.4% |
| nutrient_starvation | 984 | 36.0% | 39.6% | 39.5% |
| isr_er_translation_stress | 574 | 35.9% | 40.9% | 39.4% |
| drug_or_stimulus | 1,230 | 37.2% | 43.6% | 42.9% |
| tumor_cancer_context | 4,264 | 36.3% | 47.6% | 47.4% |
| loss_of_function | 820 | 41.0% | **54.4%** | **54.4%** |
| fraction_or_protocol | 164 | 42.7% | **62.8%** | **61.0%** |
| genetic_variant_or_mutant | 164 | 29.9% | **55.5%** | **57.3%** |
| hypoxia_oxidative_stress | 328 | **52.4%** | **65.5%** | **58.2%** |
| viral_infection | 1,066 | **52.8%** | **56.5%** | **53.9%** |

Two of the largest, most-cited design families in this atlas's own
headline story --**tumor_cancer_context** (the single largest family by
data volume, already flagged in the "smaller findings" section below as
the highest pooling-hazard family too) and **acute_interferon_cytokine** --
are among the worst performers here, with `acute_interferon_cytokine`
essentially failing outright (25-26% across all three models, i.e. *wrong*
direction about three times out of four). By contrast, **viral_infection**
and **hypoxia_oxidative_stress** are the two families where every model
variant clears 50%, and the motif-aware models clear it comfortably
(56.5-65.5%) -- consistent with Finding 4's independent observation that
viral-infection evidence in this atlas tends to be unusually robust, and a
new, convergent piece of evidence for the same conclusion from a completely
different analysis (predictive cross-validation rather than LOO-axis
robustness).

This does not mean the atlas's descriptive findings (e.g. "gene X's
clean-CDS usage differs between case and control in study Y") are wrong --
those are direct, in-sample observations, unaffected by this. What it means
is that the *predictive/generalization* claims built on top of these models
-- most importantly `dominant_rdg_therapeutic_perturbation_predictions.csv`
and any other output that presents a *predicted* effect for a
gene/condition combination not directly observed -- should be trusted in
rough proportion to how well the relevant design family scores here, not
uniformly. A predicted effect for a **viral_infection** or
**hypoxia_oxidative_stress** gene/condition pair has real, cross-validated
support behind its direction; a predicted effect for
**acute_interferon_cytokine** or **tumor_cancer_context** currently does
not, and should be treated as a hypothesis-generating guess rather than a
model-validated prediction until the underlying models improve for those
families specifically. Two caveats worth carrying forward: (1) the
`hierarchical` and `motif_model`/`motif_design_model` numbers are not
strictly apples-to-apples, since they use different held-out schemes
(leave-one-study-out vs. leave-one-gene-out) answering subtly different
generalization questions; and (2) `sign()` on a log2FC near zero is
inherently noisy, so some of the below-chance accuracy in small-effect
design families may reflect genuinely small, hard-to-call effects rather
than a model that is actively anti-predictive.

## Smaller findings from the auxiliary-state screen

`dominant_state_auxiliary_module_screen.R` already independently proposed
several candidate states beyond the current ten, worth surfacing because
some connect directly to other threads in this project:

- **"Ribosome quality control / collision"** (ZNF598, GIGYF2, EIF2A, PELO,
  HBS1L, ASCC3, NEMF, LTN1, GCN1, EIF2AK4) is already registered as a
  candidate, currently blocked purely by `needs_fst_marker_expansion` (0 of
  10 marker genes have FST coverage yet). This is essentially the same
  module independently recommended in
  `dominant_state_gene_set_theory_and_priorities.md` on theoretical and
  compendium-feasibility grounds (ribosome rescue/RQC, motivated there by
  the EIF4G2 x ABCE1 recurrence) -- two independent routes into the same
  recommendation is a good sign it is worth prioritizing once FST pages are
  expanded.
- **"DecodeME GWAS-prior expression layer"** (BTN2A2, OLFM4, RABGAP1L,
  **ZNFX1**, FBXL4, CA10, ARFGEF2, CSE1L) is blocked for the identical
  reason. ZNFX1 is one of the four genes flagged earlier this session as
  unresolvable in `manual_translons` (not in the current clean-CDS marker
  panel) -- the same root cause (marker-panel/FST coverage) blocks two
  unrelated analyses at once, which sharpens the case for prioritizing that
  expansion over other options.
- **"Host mRNA export / initiation shutoff"** (NUP98, RAE1, NXF1, ALYREF,
  EIF4G1, **EIF4G2**, EIF4E, EIF4A1, EIF4A2, PABPN1) is a registered
  pipeline candidate, correlated with Post-viral fatigue (r=0.76). EIF4G2
  appearing here is a coincidental but genuine cross-reference back to this
  project's own EIF4G2 story -- worth a note if this module is ever
  developed further, since it means the atlas's own compendium-wide screen
  independently flagged EIF4G2's functional neighborhood as biologically
  interesting, from a completely different angle (nuclear export/shutoff
  machinery) than the uORF-specific finding.
- **Tumor/cancer context is simultaneously the atlas's largest and its
  most confounded design family.** `pooling_bias_report.md`'s
  "Top Pooling-Hazard Design Families" table ranks
  `tumor_cancer_context` #1 by `pooling_hazard_score` (0.746, versus 0.615
  for the next-highest, `viral_infection`) *and* it independently has by
  far the most supporting data of any family (9,430
  `branch_context_rows`, roughly 6x `viral_infection`'s 1,148 and 20x most
  others). That combination is worth flagging on its own: it is the
  opposite of the usual "more data means more confidence" intuition --
  here, more data has accumulated precisely because tumor/cancer studies
  pool together many different tissues (kidney, bone, breast, brain,
  colon per `top_tissues`) and case conditions (KD/KO/treated/OE/DKO --
  a heterogeneous grab-bag rather than one consistent perturbation) under
  a single design-family label. Any consensus claim in this atlas that
  leans on `tumor_cancer_context` support should be treated as recurrence
  evidence across a loose, biologically heterogeneous bucket, not as a
  well-replicated single-condition result, regardless of how large its
  study count looks -- exactly the distinction the report's own "Claim
  Rules" section already recommends, but worth restating here since this
  family's sheer row-count could otherwise read as its strongest, not its
  weakest, evidence.
- **"Antigen presentation / immune composition"** (TAP1, TAP2, PSMB8,
  PSMB9, CIITA, HLA-A/B/C) shows the single largest effect size in the
  entire screen (Cohen's d = 4.49, q = 3.4e-204), but its driving metadata
  axis is `CELL_LINE = LCL` (lymphoblastoid cell lines), not a treatment or
  condition. LCLs are B-cell-derived and constitutively express antigen
  presentation machinery as a matter of cell identity, not induced state.
  This candidate should not be promoted as a dynamic "dominant state"
  without first confirming the effect survives *within* non-LCL contexts
  too -- as scored now, it is at real risk of being a cell-type-composition
  confound dressed as a biological state, and should be flagged as such if
  it moves further through the pipeline.

## What was added versus what already existed

**Code fixed and end-to-end verified this session** (Finding 1): the raw,
miscalibrated `dominant_state` label was propagating uncorrected into four
active production scripts --
`dominant_cell_state_control_aware_enrichment.R`,
`dominant_cell_state_adjusted_enrichment.R` (the actual source of the
enrichment-table call-fraction/Fisher-test columns; the first fix attempt
targeted the wrong file and was inert until this was found),
`dominant_ribosome_decision_graphs.R` (the RDG atlas's always-present
`sample_dominant_call` branches), and `postviral_fatigue_pipeline_summary.R`
(the fatigue module's own headline `fatigue_winner_samples` metric). All
four now merge in `human_dominant_state_activity_metadata.fst`'s calibrated
winner instead, the pipeline's cache-dependency graph was updated to match,
and a full `--force-step` pipeline rebuild confirmed the corrected numbers
propagate through (verified concrete before/after counts are in Finding 1).

**New analysis, computed fresh for this note, not yet standalone scripts**:
the recalibrated winner-count comparison and residual-fatigue-score check
(Finding 2), the hidden-uORF-candidate literature cross-check (Finding 3),
the transport-robustness comparison across DDIT3/PPP1R15A/ATF4 (Finding 4),
the PCA/latent-cluster vs. glmnet-metadata cross-validation of the
cell-line-vs-treatment confound (Finding 5), the single-study concentration
check on the uORF on/off-switching catalog (Finding 6), and the
held-out-cross-validation direction-accuracy read-out of the atlas's own
predictive models (Finding 7). Each is a one-off validity audit built from
already-computed atlas outputs rather than new infrastructure -- worth
keeping as documented findings rather than recurring pipeline steps unless
the underlying module/model definitions change again and need re-checking.
If any of Findings 5-7 turn out to be load-bearing for how results get
cited going forward, the natural next step would be promoting the
`feature_role` breakdown (Finding 5) and the direction-accuracy-by-design-family
table (Finding 7) into a small, explicit "confidence" column joined onto
the relevant output CSVs, rather than leaving them as prose caveats only
discoverable by reading this note.
