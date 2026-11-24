# Gene-set revision plan: what the correlation structure says should change

Updated: 2026-09-25

## Purpose and a methodological rule to apply throughout

This is a concrete, actionable plan -- not just observations -- for which
genes in `dominant_state_module_definitions.R`'s 9 modules should be
moved, removed, added, re-signed, or left alone, built from a per-gene
"leave-one-out cohesion" audit (each up-gene's mean correlation with the
*rest* of its own module's up-genes; each down-gene's mean correlation with
its module's up-gene set, which should be negative) across the full
3,857-sample compendium. It also proposes concrete master-metadata-sheet
additions. It builds directly on
`dominant_state_pca_umap_correlation_patterns.md`'s two-supercluster
finding and digs one level deeper, to individual genes.

**Rule to apply before changing anything**: this atlas measures **ribosome
occupancy / clean-CDS relative usage**, not mRNA abundance. A gene's
GO/literature-assigned direction is almost always established at the
**transcript level**. When an empirical correlation contradicts a gene's
assigned direction, the first hypothesis to rule out is **transcription/
translation decoupling** -- e.g. a gene whose own translation is
uORF-gated (like `ATF4`) can show a real, mechanistically-correct
stress response as a *change in the fraction of its mRNA being
translated* without its raw ribosome-footprint count (this atlas's
`clean_cds` metric, expressed relative to a per-sample baseline) moving
the same way a bulk-RNA-seq study would show. Only when translation-level
reasoning is checked and doesn't explain the mismatch -- or when a cleaner
explanation is a different confound (cell-line identity, mixed timescales,
a normalization-reference problem) -- should a gene actually be moved or
removed. Every recommendation below states which of these applies.

## Method

For each of the 9 modules' up-gene list, and each gene `g` in it: compute
`mean(cor(g, other up-genes in the same module))` using the 129x129
Spearman correlation matrix from
`dominant_state_pca_umap_correlation_patterns.md`. For each down-gene,
compute its mean correlation with the module's up-gene set (expected:
negative). Full numbers in the companion table
(`module_gene_audit.csv` / `module_cohesion_summary.csv`, scratch outputs,
not yet promoted to `scripts/`).

## Module-by-module findings and recommendations

### Baseline (control) -- the normalization reference set itself has two weak genes

Ranked by cohesion with the other 15 housekeeping genes: `RPS27A` (0.39),
`RPS18` (0.38), `RPL27` (0.36), `RPL30` (0.35), `EEF1A1` (0.33),
`RPL13A` (0.32), `PPIA` (0.31), `PSMB4` (0.28), `TPT1` (0.27),
`GAPDH` (0.27), `HPRT1` (0.22), `YWHAZ` (0.14), `TBP` (0.13),
`ACTB` (0.10), **`UBC` (-0.014), `B2M` (-0.040)**.

`B2M` and `UBC` are the two clear outliers, and both have a specific,
known reason to be unstable reference genes, not just noise: **`B2M`
(beta-2-microglobulin) is a core MHC class I component and a
well-documented interferon-stimulated gene**, and **`UBC` (ubiquitin C)
is itself stress-inducible** (proteotoxic/heat-shock stress increases free
ubiquitin demand). Both of these are real, literature-consistent reasons
these genes move with condition rather than staying flat -- this is not a
transcription-vs-translation subtlety, it is a straightforward case of
"these were never good housekeeping genes to begin with," and their
current low cohesion in *this atlas's own translation data* directly
confirms that generic literature-housekeeping status does not guarantee
translational stability either.

**This matters more than a typical single-gene flag**, because `Baseline`
is not just another state -- it is the **per-sample normalization
reference** (`gene_z <- sweep(expr_norm, 2, baseline, "-")`, see Finding 1
of the open-analysis note). If `B2M` and `UBC` are not stable, they
introduce sample-specific noise directly into the reference level every
other state's score is computed against.

**Recommendation**: remove `B2M` and `UBC` from the `Baseline` up-gene list
(drop to 14 genes). `ACTB`/`TBP`/`YWHAZ` (0.10-0.14) are weaker than the
core ribosomal/translation-factor genes but not clearly wrong-signed --
keep them, but do not add more genes from this same weak tier without
checking cohesion first. **Priority: high** -- this is the one change in
this note that could shift every other module's scores, not just its own.

### OXPHOS (mitochondrial) -- the up-genes are solid, but its own glycolysis "down-genes" are backwards

Up-genes are cohesive and clean: `COX5A` (0.43), `ATP5MC1` (0.38),
`SDHB` (0.37), `PHB2` (0.34), `NDUFA1` (0.30), `ATP5F1A` (0.28),
`UQCRC1` (0.27), `NDUFS1` (0.22) -- no changes needed here.

The down-genes are the problem: **`LDHA` (+0.183) and `PDK1` (+0.134)
correlate *positively*, not negatively, with the OXPHOS up-gene set**
(`SLC2A1` at -0.036 is weak but at least correctly signed). Both `LDHA`
and `PDK1` are simultaneously assigned as *up*-genes for `Hypoxia`, where
they cohere well and correctly (0.175 and 0.171, mid-pack among Hypoxia's
own up-genes) -- so the issue is specifically their reuse as an
OXPHOS-down marker, not their quality as hypoxia markers.

Applying the translation-first rule: the canonical Warburg-effect logic
(glycolysis genes transcriptionally up, OXPHOS genes transcriptionally
down, under HIF1A) is well-established **at the mRNA level**. This
atlas's own `dominant_rdg_consensus_candidates.csv` already independently
flagged `LDHA` (this atlas's #1 consensus-ranked gene) as a case where
clean-CDS usage moves in the opposite direction from what HIF-driven
transcriptional induction would predict, hypothesizing a **translational
throttle on top of real transcriptional induction** -- i.e., `LDHA` mRNA
may genuinely go up under hypoxia (literature is not wrong at the
transcript level) while its *ribosome allocation specifically to the
clean CDS* does not track that increase 1:1, and may instead track
"how much cap-dependent/uORF-gated translation is generally permitted in
this sample" -- a variable that also happens to correlate with how well
OXPHOS genes translate. That would produce exactly this pattern: `LDHA`/
`PDK1` positively correlated with OXPHOS at the translation level despite
being real, opposite-direction hypoxia/Warburg genes at the transcript
level.

**Recommendation**: do not treat this as "the gene set is wrong" and
delete `LDHA`/`PDK1` outright. Instead: (1) remove them specifically from
`OXPHOS`'s *down*-gene role only (keep them as `Hypoxia` up-genes, where
they are empirically correct); (2) add a note to
`dominant_rdg_consensus_candidates.csv`'s LDHA row and to the module
definitions file cross-referencing this translation/transcription
distinction, so a future reader does not "fix" `LDHA` back to a
transcript-level assumption without checking translation data first; (3)
the concrete follow-up that would settle this definitively: compare
`LDHA`/`PDK1` clean-CDS usage against matched bulk-RNA-seq (or total,
non-clean-CDS-normalized ribo-seq footprint count) for the same
samples/studies where available, to directly test the "transcript up,
translation-allocation flat/down" hypothesis rather than inferring it.
**Priority: medium** -- affects one gene's role in one module, not the
normalization backbone.

### Senescence / SASP -- CDKN2A does not covary with the rest of the module, plausibly a timescale mismatch, not a translation artifact

Up-genes: `SERPINE1` (0.42), `IGFBP7` (0.40), `TGFBI` (0.40),
`TIMP1` (0.39), `CCND1` (0.37), `IL6` (0.33), `CDKN1A` (0.31),
`CXCL8` (0.28) form a cohesive core. `IGFBP5` (0.075) is weak (see also
Finding 3 of the open-analysis note, which flagged `IGFBP5` as a "hidden
uORF candidate" with unexplained low clean-CDS usage -- its weak module
cohesion here is a second, independent signal that something about
`IGFBP5` is idiosyncratic, not just noise). **`CDKN2A` (p16) is the clear
outlier at -0.124**, and the down-gene `LMNB1` is correctly, strongly
negatively signed (-0.241, a real, validated senescence marker).

Here the translation-vs-transcription lens does *not* obviously apply
(nothing about p16 is known to be specifically uORF-gated in a way that
would explain this), so the more likely explanation is a **timescale/
mechanism mismatch already recognized in the literature**: `CDKN2A`
accumulates slowly via chromatin derepression over many cell divisions
(replicative or oncogene-induced senescence, a process of days-to-weeks),
while `IL6`/`CXCL8`/`SERPINE1`/`TIMP1`/`IGFBP7`/`CCND1`/`CDKN1A` (the
acute SASP/cell-cycle-arrest-onset genes) can be induced within hours by
many acute stresses that are not necessarily committed, chromatin-level
senescence. This compendium is dominated by short-timescale experimental
perturbations (drug treatment, knockdown, acute stress -- see the
`design_family` breakdown in the open-analysis note), so it would see the
fast SASP program far more often than the slow chromatin-commitment
program, which would produce exactly this low-CDKN2A-cohesion pattern
even if both genes are "real" senescence markers on different clocks.

**Recommendation**: do not remove `CDKN2A` (it is too well-established a
marker to drop from the panel), but (1) stop treating it as
interchangeable with the acute SASP genes in any score that assumes they
move together; consider a `senescence_acute_sasp_score` (the 8-gene
cohesive core) vs. a separate `senescence_commitment_score` (`CDKN2A`,
plus `LMNB1` loss) as two sub-scores rather than one blended module,
mirroring the split already worth doing for ISR (below); (2) this is
directly testable with the metadata addition proposed below (an
induction-method/timescale field): if the hypothesis is right, `CDKN2A`
should show good cohesion with the SASP genes specifically within
replicative/long-timescale senescence experiments, and poor cohesion
within short-timescale acute-stress experiments. **Priority: medium.**

### ISR / ER stress (ATF4 axis) -- a real internal structure, not a binary split, but still worth two sub-scores

All 10 up-genes show *positive* cohesion with each other (no sign
problems), but with a clear internal gradient: `DDIT3` (0.31),
`HERPUD1` (0.28), `ATF4` (0.26), `ATF3` (0.25), `DDIT4` (0.24),
`PPP1R15A` (0.23), `XBP1` (0.21), `CEBPB` (0.18), `ASNS` (0.15),
`HSPA5` (0.09, the weakest link). Read in isolation, this module is
reasonably cohesive. The sharper finding from
`dominant_state_pca_umap_correlation_patterns.md`'s whole-panel
clustering is that **when compared against all 129 genes, not just its
own 10, `ATF4`/`ATF3`/`ASNS` correlate *even more strongly* with
proliferation/translation-machinery genes than they do with the rest of
the ISR set**, pulling them into the opposite supercluster during global
clustering, while `DDIT3`/`DDIT4`/`TRIB3`/`PPP1R15A`/`HSPA5`/`XBP1` stay
firmly in the stress/damage supercluster.

Applying the translation-first rule, this is a strong candidate for a
genuine, literature-grounded, translation-specific explanation: **`ATF4`
is the textbook example of a gene whose own stress response is a
*translation-efficiency switch* driven by its 5' leader uORFs**, not
primarily a transcriptional one -- under stress, a larger *fraction* of
existing `ATF4` mRNA gets translated, without necessarily a proportional
change in total mRNA level. `ATF3`/`ASNS` are direct ATF4 transcriptional
targets whose induction is comparatively more transcription-driven. If
this atlas's clean-CDS metric (ribosome allocation relative to a
per-sample baseline) responds more sensitively to *transcriptional*
induction (more total mRNA to allocate ribosomes to) than to a pure
translation-efficiency switch on an unchanged mRNA pool, `ATF4` could show
a smaller apparent shift than heavily transcriptionally-induced genes like
`DDIT3`/`PPP1R15A` -- exactly the pattern observed. This is a genuinely
useful, checkable hypothesis, not just an excuse: it predicts `ATF4`
should show a *bigger* signal in the atlas's uORF-specific outputs
(`dominant_uorf_regulation_*`, `uorf_structure_rules`) than in the plain
clean-CDS score, since those outputs are built to detect exactly this kind
of translational-efficiency change.

**Recommendation**: (1) keep all 10 genes in the module (none are
wrong-signed), but split reporting into two sub-scores --
`isr_atf4_direct` (`ATF4, ATF3, ASNS`) and `isr_terminal_upr`
(`DDIT3, PPP1R15A, DDIT4, HERPUD1, HSPA5, XBP1, CEBPB`) -- and check
whether they diverge under specific conditions (prolonged vs. acute ER
stress in the same study is the cleanest test); (2) directly test the
uORF-translation-efficiency hypothesis for `ATF4` specifically by cross-
referencing its clean-CDS score trajectory against its own uORF-usage
metrics in `dominant_uorf_regulation_edges.csv` for the same samples; (3)
this split also plausibly explains part of Finding 7's weak direction-
accuracy for `isr_er_translation_stress` (a model trying to predict one
direction for a module secretly containing two sub-programs would
struggle) -- worth re-running that model split by sub-score once (1) is
implemented. **Priority: medium-high** -- directly actionable and
testable, and connects two other open threads (Finding 3's uORF work,
Finding 7's model weakness).

### EMT / Mesenchymal shift -- solid module, one weak up-gene, down-genes correctly signed but two are weak

Mesenchymal up-genes are cohesive: `MMP2` (0.48), `POSTN` (0.39),
`FN1` (0.38), `CDH2` (0.38), `VIM` (0.38), `THY1` (0.34),
`SNAI2` (0.30), `ZEB1` (0.30), `ZEB2` (0.27) -- only **`SNAI1` (0.13)**
is notably weaker. `SNAI1` (Snail) and `SNAI2` (Slug) are both classic EMT
transcription factors but have different kinetics in the literature
(`SNAI1` a faster, more immediate-early TGF-beta response; `SNAI2` more
stable/sustained) -- plausibly a real, mechanistic kinetics difference
rather than an error, similar in spirit to the Senescence timescale issue.
Epithelial down-genes are all correctly negatively signed (`CDH1` -0.180,
`EPCAM` -0.168, `KRT8` -0.092, `KRT18` -0.019, `KRT19` -0.014), though the
two keratins are weak enough to add little information beyond `CDH1`/
`EPCAM`.

**Recommendation**: no removals. Consider de-weighting `SNAI1` relative to
`SNAI2` if the module score is ever changed from a simple mean to a
weighted one; this is low priority. The more important EMT issue is
**composition, not gene selection** -- see
`dominant_state_pca_umap_correlation_patterns.md`'s finding that the
mesenchymal pole correlates with stromal/vascular/macrophage genes
(`VWF`, `C1QA/B/C`) at the whole-panel level, which is a sample-level
metadata problem (addressed below), not a gene-list problem.

### Proliferation, mTOR, Interferon, Hypoxia -- no changes needed

`Proliferation` (mean up-up cohesion 0.650) and `Interferon` (0.542) are
the two most internally cohesive modules in the entire panel -- leave
them exactly as they are. `mTOR` (0.453) is also solid. `Hypoxia`'s
up-genes (0.185 mean, `CA9` 0.25 down to `NDRG1` 0.148) are all positively
and fairly evenly cohesive with no outliers worth flagging -- weaker
cohesion than Proliferation/Interferon simply reflects a more diverse set
of downstream HIF1A target pathways (angiogenesis, glycolysis, pH
regulation), not a data quality problem.

### Post-viral fatigue / ribosome stress -- inherits ISR's structure, not a new issue

Fatigue shares most of its up-genes with ISR/Interferon (already
quantified in the open-analysis note's Finding 2), so it inherits the same
`ATF4`/`ATF3`/`ASNS`-vs-terminal-UPR structure discussed above rather than
presenting a new problem. Its down-genes (`ATP5F1A, ATP5MC1, COX5A, CS,
EEF2, EIF4B, NDUFA9, NDUFS1, PABPC1, PPARGC1A, RPL32, RPLP0, RPS3, RPS6,
SDHB, TFAM, UQCRC1`) are mostly OXPHOS/translation-machinery genes
reused with the opposite sign, which is coherent with fatigue being framed
as a hypometabolic/hypo-translational state -- no changes recommended
here beyond what the ISR split above already proposes.

## Do we need more genes? Corrected after implementation: no, not from this set

**Update after implementation (see
`dominant_state_gene_justifications/README.md`): this section's original
recommendation was proposed, then reverted the same session.** 14 genes
are already present in the atlas's 129-gene FST/clean-CDS panel but are
not assigned to *any* of the 9 scored modules. The original guess was that
these were orphaned candidates from a still-blocked "antigen presentation"
idea. Checking their actual source instead
(`results/curated_inputs/postviral_fatigue_dominant_state_genes.csv`)
found they are in fact **all 14 of that file's own `Tier = exploratory`
genes** (vs. 43 `Tier = conservative` genes that are actually scored),
and every one of them -- including the four this note originally proposed
promoting -- already carries a curator-written `Rationale` ending in
"tissue-composition risk". The atlas's own prior curation had already
answered this question with better information than the correlation
numbers below can provide on their own; the numbers are kept below as a
record of what motivated the (reverted) proposal, not as a live
recommendation. **The corrected answer to "do we need more genes" is: not
from this already-curated exploratory set -- if more genes are wanted,
they should come from outside it, or the exploratory tier itself should
be revisited on its own biological merits, not promoted piecemeal based on
correlation strength.**

| gene | best-correlated module | r | second-best |
|---|---|---|---|
| ICAM1 | Interferon | 0.45 | Post-viral fatigue (0.30) |
| CFB | Interferon | 0.40 | Senescence (0.32) |
| C3 | Interferon | 0.37 | Senescence (0.27) |
| CD74 | Interferon | 0.35 | Post-viral fatigue (0.18) |
| VCAM1 | Interferon | 0.34 | EMT (0.21) |
| MS4A1 | Interferon | 0.31 | Post-viral fatigue (0.14) |
| HLA-DRA | Interferon | 0.27 | Post-viral fatigue (0.12) |
| CD79A | Interferon | 0.24 | Post-viral fatigue (0.13) |
| F3 | Senescence | 0.41 | EMT (0.27) |
| C1QB | Interferon | 0.15 | Post-viral fatigue (0.07) |
| C1QC | Interferon | 0.13 | EMT (0.07) |
| HLA-DRB1 | Proliferation (weak/noise) | 0.15 | Interferon (0.13) |
| C1QA | Interferon | 0.10 | EMT (0.04) |
| VWF | Senescence (weak) | 0.11 | Interferon (0.11) |

**Originally proposed** (then reverted): `ICAM1`, `CFB`, `C3`, `CD74`,
`VCAM1` correlate with `Interferon` at 0.34-0.45, and `F3` with
`Senescence` at 0.41 -- comparable to or better than some of those
modules' existing weaker members. The original reasoning tried to split
this 14-gene set into a "safer" subset (`ICAM1`, `CD74`, `VCAM1`, `F3`)
and a composition-confound subset (the complement/MHC-II/B-cell genes),
but the curator rationale found after implementation tags **all 14**,
including the "safer" four, as tissue-composition risk -- `ICAM1`/`VCAM1`
specifically as endothelial adhesion markers, `CD74` as an EBV/B-cell
context marker, `F3` as a coagulation marker. All 14 should be treated as
composition markers, not folded into any dynamic state, unless
re-evaluated with per-cell-line/tissue-composition-controlled evidence
rather than whole-compendium correlation alone. See
`dominant_state_gene_justifications/03_interferon_add_ICAM1_CD74_VCAM1.md`
and `04_senescence_add_F3.md` for the full record.

## Metadata master-sheet additions this analysis motivates

Three separate lines of evidence this session (Finding 5's cell-line-
dominance result, the UMAP islands, and the CDKN2A/SASP timescale
hypothesis above) all point to the same gap: **the master metadata sheet
does not currently carry the fields needed to correct for, or test,
these confounds directly.** Concrete additions:

1. **A persisted `design_family`/`perturbation_class` column on the master
   sheet itself**, computed once by a canonical function and stored, not
   re-derived independently inside `dominant_cell_state_adjusted_
   enrichment.R`, `control_aware_enrichment.R`, `dominant_condition_
   families.R`, etc. separately. Right now the same condition string can
   in principle be classified slightly differently by different scripts'
   own copies of this logic (a real consistency risk, not yet confirmed to
   have caused an actual discrepancy, but worth closing given how many
   scripts each carry their own version).
2. **A senescence-induction-method / experiment-timescale field**
   (`replicative`, `oncogene_induced`, `drug_induced_acute`, `stress_acute`,
   `not_senescence_related`, with a rough duration bucket if extractable
   from the condition string) -- this is the direct, falsifiable test for
   the CDKN2A-vs-SASP timescale hypothesis above: split samples by this
   field and check whether CDKN2A's cohesion with the SASP core genes
   improves specifically within the `replicative`/`oncogene_induced`,
   long-timescale subset.
3. **A per-cell-line empirical state-calling prior**, precomputed once
   (e.g. from `human_dominant_state_activity_metadata.fst` itself, as
   already sketched informally in
   `dominant_state_pca_umap_correlation_patterns.md`'s A549/HEK293/K562
   table) and stored as a lookup the enrichment/control-aware scripts can
   join against, so a future analysis can *report* "this effect exceeds
   what this cell line normally shows for this state" rather than needing
   to separately re-derive that baseline each time. This directly
   operationalizes the caution already raised for A549 (73.1% baseline
   Senescence-calling) and LCL (antigen-presentation confound).
4. **A gene-level sidecar table** (not a sample metadata column, but a
   parallel small reference table keyed by `gene_symbol`) recording, per
   marker gene: known uORF/translational-control status (yes/no/unknown),
   and which regulatory level (transcriptional, translational, both)
   its literature-established direction was established at. This is the
   single most useful artifact this plan's own reasoning process needed
   and had to reconstruct informally each time (for `ATF4`, `LDHA`,
   `PDK1`) -- formalizing it would make future gene-direction disputes
   fast to check rather than requiring a fresh literature dive each time.

## Summary action list

| action | genes/fields | priority | status |
|---|---|---|---|
| Remove from `Baseline` up-genes | `B2M`, `UBC` | **High** (affects normalization) | **Implemented** |
| Remove from `OXPHOS` down-genes only (keep as `Hypoxia` up-genes) | `LDHA`, `PDK1` | Medium | **Implemented** |
| ~~Add as new `Interferon` up-genes~~ | `ICAM1`, `CD74`, `VCAM1` | -- | Reverted (curator-tagged composition risk) |
| ~~Add as new `Senescence` up-gene~~ | `F3` | -- | Reverted (curator-tagged composition risk) |
| Do NOT add to any dynamic state | all 14 `Tier = exploratory` genes in `postviral_fatigue_dominant_state_genes.csv` | -- | Confirmed as designed |
| Split into two sub-scores | ISR -> `isr_atf4_direct` / `isr_terminal_upr` | Medium-high | **Implemented** (standalone script, see justifications/05) |
| Split into two sub-scores | Senescence -> acute SASP core / commitment (`CDKN2A`+`LMNB1`) | Medium | **Implemented** (standalone script, see justifications/05) |
| New master-sheet column | canonical persisted `design_family` | Medium | Not yet implemented |
| New master-sheet column | senescence induction-method/timescale | Medium (testable, falsifiable) | Not yet implemented |
| New reference artifact | per-cell-line state-calling prior lookup | Medium | Not yet implemented |
| New reference artifact | gene-level transcription-vs-translation sidecar table | Medium | Not yet implemented |

Two of the ten actions above are implemented and verified (see
`dominant_state_gene_justifications/`); two were proposed and then
reverted the same session after a follow-up check found stronger contrary
evidence than the correlation analysis that motivated them; the rest
remain plans pending review.
