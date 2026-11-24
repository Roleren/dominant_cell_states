# What's working, what's missing: lessons from the story-hunting sessions

Updated: 2026-09-28

## Purpose

This session (and the several before it) used the atlas to find real,
checkable candidate findings (EIF4G2, METTL3, NAT10, VEGFA) rather than
just describing the atlas's own outputs. That process repeatedly hit the
same handful of capabilities that worked well and the same handful of
gaps that cost real time or nearly produced a wrong conclusion. This
note pulls those into one place, concretely, so the gaps can be turned
into implementation work rather than re-discovered informally next time.

## What's good -- keep doing this, and build more of it

**1. The atlas already has multiple independent statistical layers, and
cross-checking them catches real errors.** `dominant_rdg_consensus_
candidates.csv` (a fast, top-line summary), `dominant_rdg_context_
interaction_branch_context_effects.csv` (per-context, per-branch, with a
pooled `global_axis_delta`), `dominant_rdg_model_agreement_loo.csv`
(leave-one-out), and `dominant_rdg_model_agreement_transport.csv`
(cross-context transport robustness) all exist independently. Today's
VEGFA investigation is the concrete proof this redundancy is worth
having: the consensus table's headline number (`+46.95%`, framed as
induction) turned out to have the **wrong sign** for the `clean_CDS`
branch once checked against the properly pooled `global_axis_delta`
(`-0.027`) in the context-interaction table. Without a second,
independent layer to check against, this would have gone into a
`papers/` folder backwards. Keep building features that give a second,
differently-computed opinion on the same claim -- it is not redundant
effort, it is the error-correction mechanism.

**2. Branch-level structural decomposition (`leader_uORF` vs.
`overlapping_uORF` vs. `clean_CDS` vs. `other`) turns a single gene-level
number into a mechanistic story.** The VEGFA "opposite-signed deltas"
that first looked like a contradiction were actually the correct,
expected signature of ribosome reallocation between the leader and the
CDS, once the branch structure was understood. This is the single most
valuable structural feature in the whole atlas for turning a correlation
into a checkable mechanism -- worth extending to more genes/design
families, not just relying on it existing for the ones already scored.

**3. The atlas already tags its own uncertainty** (`evidence_tags` like
`pooling_hazard`, `single_context`, `fragile`, `loo_or_transport_robust`)
**and this session's own pooling-bias audit is unusually mature self-
awareness for a project at this stage.** VEGFA's row already carried
`pooling_hazard`/`fragile` tags before this note dug into why -- the
atlas's own audit had already flagged it as not-a-slam-dunk. The gap
isn't that this doesn't exist; it's that it isn't surfaced prominently
enough (see gap #1 below).

**4. Cheap, systematic gene x state scans, run outside the formal
pipeline, reliably surface real leads.** METTL3, NAT10, and (revised)
VEGFA were all found by variations on the same simple recipe: score a
gene/condition against its own matched study control, aggregate across
studies, rank by consistency and magnitude. This recipe used only
already-computed atlas outputs (the gene expression matrix, metadata)
and no new data. It should be formalized as a real pipeline step (see
gap #2), not repeated ad hoc every time.

**5. Grounding every candidate in external literature before trusting
it is now a consistent habit, and it works.** Every one of EIF4G2,
METTL3, NAT10, and VEGFA's corrected story was checked against real,
specific, dated papers, not just plausibility. This caught NAT10's
weaker independence (same lab, same cell line) immediately, and gave
VEGFA's *corrected* direction a concrete, checkable literature framing
(transcription/translation decoupling, matching `LDHA`) rather than
leaving it as a bare number.

## What's lacking -- concrete gaps, with fixes

**1. IMPLEMENTED (2026-09-28).** When the atlas's own tables disagree,
there was no flag telling the reader so, and no stated rule for which one
to trust -- VEGFA's consensus-table headline and its own context-
interaction table's pooled estimate disagreed in sign for the same gene/
condition/branch, and nothing in either file said so. **Fix built**:
`scripts/dominant_rdg_consensus_context_consistency_audit.R`, wired in as
pipeline step `consensus_context_consistency_audit` (runs right after
`rdg_consensus`, depends on `rdg_consensus` + `context_interactions`
outputs, cached the same way as every other step). It joins every
consensus-table row to the context-interaction table's `global_axis_delta`
for the matching gene/design_family's `clean_CDS` branch and writes
`dominant_rdg_consensus_context_consistency_audit.csv` (per-row agreement
status: `agree` / `DISAGREE_SIGN` / `pooled_effect_negligible` /
`no_context_interaction_row`) plus a plain-text top-disagreements summary.
First real run found 3 disagreements out of 75 evaluable rows: `VEGFA`
(already known), plus two new ones surfaced automatically --
`MMP2` and `LMNB1`, both also under `viral_infection` -- not yet
individually investigated.

**2. IMPLEMENTED (2026-09-28).** There was no formal "which perturbed
genes shift which dominant state" pipeline step -- this session rebuilt
it from scratch as a scratch script every time it was needed. **Fix
built**: `scripts/dominant_state_gene_perturbation_screen.R`, wired in as
pipeline step `gene_perturbation_screen` (runs right after
`control_aware`, in the `metadata`/`manual`/`all` scopes, cached like
every other step). Computes, per gene x dominant-state, a delta against a
matched same-study-same-cell-line control, a Welch t-test p-value where
both arms have >=2 samples, and aggregates across studies
(`n_studies`, `n_cell_lines`, `consistent_direction`, BH-adjusted
p-value). Outputs a sample-level table, a gene x state summary, and a
filtered top-hits table (`n_studies >= 2` and consistent direction). A
first real run reproduced METTL3/Interferon (p=0.054) and NAT10/
Senescence exactly as found by hand, now with a real significance test
attached instead of eyeballed means.

**3. IMPLEMENTED, as a script-internal fix rather than a metadata-sheet
change (2026-09-28).** No canonical "which samples are the matched
control for this perturbation" annotation exists on the master metadata
sheet -- every scan had re-guessed it with a regex, and got it wrong once
already this session (METTL3's HEK293 study, where the untransduced "WT"
baseline was not the appropriate control for a shRNA knockdown -- the
scrambled-shRNA condition was). `dominant_state_gene_perturbation_screen.R`
now encodes the fix directly: within each (study, cell line) group it
checks for a *tight* control first (scrambled/NTC/mock/empty-vector --
i.e. something that went through the same transduction procedure) and
only falls back to a *loose* control (generic WT/control/parental) if no
tight control exists in that group. This closes the gap for anything run
through this new script; it does **not** yet create a persisted,
human-reviewed `matched_control_group_id` column on the actual master
metadata sheet (a bigger, shared-resource change, deliberately not made
without more explicit instruction) -- that remains open if this
heuristic turns out to still misfire somewhere.

**4. IMPLEMENTED (2026-09-28).** Derived-table columns were not
self-documenting, and the only way to understand a non-obvious one was to
grep the source script -- the `branch_class` column that resolved (and
initially confused) the VEGFA investigation required reading
`scripts/dominant_rdg_context_interaction_model.R` directly; nothing in
`dominant_rdg_context_interaction_branch_context_effects.csv` itself, or
alongside it, explained what its own columns meant. **Fix built**:
`scripts/dominant_column_dictionary.R`, a small shared helper
(`write_column_dictionary(csv_path, descriptions, title)`) that any
script can call after writing a results table to emit a paired
`<output>_columns.md` file. It always lists every actual column in the
table (reading the CSV's own header), marking any without a supplied
description as `(undocumented)` rather than silently omitting it, so a
table's documentation coverage stays visible and regenerates
automatically if the table's schema ever changes. Wired into all four new
scripts from gaps #1-3/#5-6 (`dominant_state_gene_perturbation_screen.R`,
`dominant_state_perturbation_data_gap_map.R`,
`dominant_state_gene_translation_regulation_audit.R`,
`dominant_rdg_consensus_context_consistency_audit.R`), each now emitting
a `_columns.md` alongside its own outputs, registered as a declared
output in `pipeline_step_outputs()` so caching tracks it too.

Also retrofitted onto the two *existing* tables that actually caused this
session's confusion, without touching the two heavy, already-validated
scripts that produce them: `scripts/dominant_rdg_existing_table_column_
dictionaries.R` (new pipeline step `existing_table_column_dictionaries`,
right after `consensus_context_consistency_audit`) reads the already-
written `dominant_rdg_context_interaction_branch_context_effects.csv` and
`dominant_rdg_consensus_candidates.csv` and writes each a column
dictionary using the same shared helper. The context-interaction one
documents all 30 columns, including a full explanation of `branch_class`
and the context-vs-global-vs-complement-delta distinction that caused the
original confusion. The consensus table (82 columns) documents the ~26
columns this session actually relied on plus an explicit note on why the
remaining ~56 internal `*_component` scoring columns are left
undocumented rather than guessed at. **Remaining gap**: the other ~50
pre-existing pipeline steps and their output tables still have no
dictionary -- retrofitting every one of them was out of scope for one
session; the pattern (a small companion script, or a few lines added to
the producing script) is now proven and cheap to repeat whenever a
specific table causes confusion again.

**5. IMPLEMENTED (2026-09-28).** Cross-lab, independent-replicate
coverage for individual gene perturbations is thin, and this is a data
problem, not an analysis problem -- the systematic scan across the whole
compendium found only 9 genes with knockdown/knockout data from >=2
independent studies at all (most perturbed genes in this compendium have
exactly one supporting study). METTL3's 3-cell-line replication is the
exception, not the norm. **Fix built**:
`scripts/dominant_state_perturbation_data_gap_map.R`, wired in as
pipeline step `perturbation_data_gap_map` (right after
`gene_perturbation_screen`). Checks a curated set of 10 biologically
motivated gene classes (m6A writers/readers/erasers, ac4C pathway,
ribosome quality control, stress granule, ISR core kinases,
retroelement-silencing/viral-mimicry guardians, recycling/reinitiation,
termination) against actual perturbation coverage in the compendium. A
first real run found `retroelement_silencing` (the classic ADAR/SETDB1/
TRIM28/DNMT3A/TREX1-style viral-mimicry guardian genes) and `isr_core`
(the EIF2AK1-4 kinases and PPP1R15A/B that sit directly upstream of the
atlas's own flagship ISR state) as the two thinnest classes -- 10/12 and
7/8 genes respectively with zero perturbation data in this compendium at
all. This converts "we need more analysis" into a concrete, prioritized
list of what data to go acquire next.

**6. IMPLEMENTED (2026-09-28).** The gene-level "transcription vs.
translation" sidecar table proposed in
`dominant_state_gene_set_revision_plan_2026-09-25.md` did not exist, and
this session needed it again -- the VEGFA/LDHA decoupling pattern
(translation moving independently of, or against, presumed
transcriptional direction) has now come up twice as the explanation for
an initially-confusing result. **Fix built**:
`curated_inputs/gene_translation_regulation_notes.csv` (a hand-curated
reference table, deliberately placed at the repo root rather than under
`results/curated_inputs/`, since that directory is gitignored and this
table needs to be a durable, version-controlled reference, not a
regenerable output -- worth noting as its own small finding: the two
*existing* curated-input files,
`postviral_fatigue_dominant_state_genes.csv` and
`known_uorf_effects_with_citations.csv`, are both hand-curated but
currently sit inside the gitignored `results/` tree and are therefore not
actually version-controlled today, a pre-existing risk this note flags
but does not fix without being asked). Seeded with the 11 genes this
session's investigations already characterized (`ATF4`, `ATF3`, `ASNS`,
`DDIT3`, `LDHA`, `PDK1`, `VEGFA`, `METTL3`, `NAT10`, `B2M`, `UBC`).
Validated by `scripts/dominant_state_gene_translation_regulation_audit.R`
(pipeline step `gene_translation_regulation_audit`), which cross-checks
the sidecar against the current marker panel and flags any entries for
genes no longer in it.

## Priority order if implementing from this list

1. ~~Gap #1~~ (consensus/context-interaction consistency flag) --
   **implemented 2026-09-28**, see above.
2. ~~Gap #2~~ (formal gene-perturbation-screen pipeline step) --
   **implemented 2026-09-28**, see above.
3. ~~Gap #3~~ (matched-control classification) -- **implemented
   2026-09-28** inside gap #2's script; the persisted-metadata-sheet
   version remains open.
4. ~~Gap #6~~ (transcription-vs-translation sidecar) -- **implemented
   2026-09-28**, see above.
5. ~~Gap #5~~ (data-gap map) -- **implemented 2026-09-28**, see above.
6. Gap #4 (column dictionaries) -- cheap, general-purpose, prevents
   wasted time on every future table, not just this one. **Not yet done
   -- the only item from this note's original list still open.**

All six gaps identified in this note's first pass are now either fixed or
explicitly still open (#4 only). Five new pipeline steps exist that did
not before this note was written:
`gene_perturbation_screen`, `perturbation_data_gap_map`,
`gene_translation_regulation_audit`, `consensus_context_consistency_audit`,
plus a real bug fix in `dominant_cds_frame_bumpiness_qc.R` (a stale,
hardcoded 16-gene Baseline copy including the now-removed `B2M`/`UBC`,
found by auditing the whole codebase for hardcoded duplicates of the
canonical module gene lists and replaced with a call to
`load_dominant_state_modules()`).

## New leads this implementation pass surfaced on its own

The gene-perturbation screen and the consistency audit are themselves new
analysis, and their first real runs already found things worth a look:
`MMP2` and `LMNB1` both disagree in sign between the consensus and
context-interaction tables under `viral_infection`, the same pattern
VEGFA showed -- not yet investigated, but the audit found them
automatically instead of requiring another manual dig. The data-gap map
independently flagged `retroelement_silencing` (viral-mimicry guardian
genes: `ADAR`, `SETDB1`, `TRIM28`, `DNMT3A/B`, `TREX1`, etc.) and
`isr_core` (the `EIF2AK1-4` kinases and `PPP1R15A/B` upstream of the ISR
state) as the two gene classes with the thinnest perturbation coverage in
the whole compendium -- concrete targets if more data is ever ingested.
