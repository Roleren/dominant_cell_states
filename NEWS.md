# Dominant Cell States Analysis NEWS

This file records implementation history, motivation for pipeline changes, and iteration-level notes. The scientific interpretation now lives in `md/dominant_cell_state_scientific_note.md`.

## 2026-07-10: R package scaffold and ignored results tree

Motivation: after splitting `dominant_cell_states` into a standalone checkout,
the repo needed normal R package structure and a clean separation between
source-like files and multi-GB local outputs.

Changes:

- Moved result-bound curated/app CSVs into `results/curated_inputs/`, moved
  manual-translon candidate CSVs into `results/manual_translons/`, and kept
  only `manual_translons/translons_to_be_added_to_manual.txt` as source.
- Moved the Shiny app payload to `inst/shiny/rdg_atlas_browser/` and added the
  exported, documented `run_dcs_app()` package launcher.
- Added `DESCRIPTION`, `NAMESPACE`, `R/`, and `tests/testthat/`.
- Added package helpers for project-root discovery, `results/` paths, local
  RiboCrypt discovery through `RIBOCRYPT_REPO`, pipeline loading, and cleanup
  of stale legacy result symlinks.
- Moved generated outputs, FST page manifests, static assets, readlength inputs,
  and result-sanity outputs under ignored `results/`.
- Kept app-derived marker/prior CSVs, known-uORF citations, and browser
  validation notes under ignored `results/curated_inputs/`.
- Updated old script startup code so a new analysis-package `DESCRIPTION` is
  not mistaken for the RiboCrypt package root.

Results:

- `git status --short` now surfaces source-like files while generated artifacts
  are ignored under `results/`.
- The pipeline writes generated outputs to `results/` and does not recreate
  root-level output symlinks.
- FST pages, ORFik annotations, BAMs, reference indexes, and other Bio_data
  artifacts remain external runtime inputs and should not be committed.

## 2026-07-10: Standalone repo handoff documentation

Motivation: the analysis is being split out of the RiboCrypt source tree into
`/home/roler/Desktop/forks/dominant_cell_states`. Future Codex sessions started
from that location need local context, startup commands, and Git hygiene without
depending on chat history.

Changes:

- Added `README.md` with the standalone path, migration procedure, key files,
  RiboCrypt dependency, current merged-translon review target, and typical run
  commands.
- Added `AGENTS.md` with local agent instructions for loading RiboCrypt,
  mounting/using S-drive data, output/Git policy, current algorithm state, and
  verification commands.
- Added a standalone `.gitignore` that excludes generated output directories,
  caches, static artifacts, and large binary/data products.
- Updated `scripts/run_dominant_cell_state_pipeline.R` so the RiboCrypt package
  root and `dominant_cell_states` analysis directory are discovered separately;
  standalone runs can use `RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt`.
- Updated the scientific note with the standalone repository handoff, current
  2026-07-10 result summary, standalone run pattern, and merged-translon gate
  counts.
- Updated the browser tutorial with standalone startup instructions and the
  current merged-replicate translon-shift candidate table as the first review
  entry point.

Results:

- A future agent starting from `/home/roler/Desktop/forks/dominant_cell_states`
  should read `AGENTS.md`, `README.md`, the scientific note, tutorial, and this
  NEWS file before editing.
- The docs now explicitly state that generated outputs are local artifacts and
  should not be committed blindly to GitHub.
- The docs preserve the key current algorithmic state: 95 raw review rows
  collapse to 60 candidate-direction rows, with no claim-ready rows and many
  viral leads routed to protocol/start-site/pooling audit.
- Standalone-style dry run with `RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt`
  resolved separate RiboCrypt and analysis paths successfully.
- Full pipeline invariant suite passed after the root-discovery change.

## 2026-07-09: Merged-replicate translon-shift claim gate

Motivation: the main discovery target is now condition-level merged replicate cases where clean-CDS allocation changes relative to other translons, especially leader and overlapping uORFs. The pipeline needed a dedicated gate that separates strong table-level clean-CDS/uORF opposition from actual claim readiness.

Changes:

- Added `scripts/dominant_rdg_merged_replicate_translon_shift.R`.
- Added the `merged_replicate_translon_shift` pipeline step after `rdg_evidence_gaps` and before `rdg_review_batches` in the `rdg`, `manual`, and `all` scopes.
- Added declared outputs and cache dependencies for joint branch allocation, atlas cards, and evidence gaps.
- Updated the pipeline invariant test to require the new step and conservative claim/readiness gates.
- Updated the scientific note with the merged-replicate pooling interpretation and next-model instructions.
- Added a collapsed candidate queue and pooling-scope summary so review starts from gene/transcript/design/shift-direction candidates instead of duplicated context rows.
- Merged allocation-versus-output audit context into the candidate queue when available.

Results:

- Wrote outputs under `dominant_rdg_merged_replicate_translon_shift/`.
- Current run: 22,896 joint rows, 15,840 exact same-condition merged rows, 8,938 clean-CDS/uORF opposition rows, and 95 inspectable review-queue rows.
- The 95 review rows collapse to 60 candidate-direction rows, including 8 multi-study exact-recurrence candidates, 1 multi-condition exact-recurrence candidate, and 18 candidates with direction conflicts elsewhere in the same gene/transcript/design.
- The top collapsed recurrent viral/stress candidates include `PPP1R15A | viral_infection`, `ATF3 | nutrient_starvation`, `TFE3 | viral_infection`, `DDIT3 | viral_infection`, and `ATF4 | viral_infection`.
- The queue retains viral ISR/antiviral anchors such as `PPP1R15A`, `DDIT3`, and `ATF4`, but marks the current viral rows as protocol/pooling audit leads rather than claim-ready rows.
- No rows are currently promoted to `claim_candidate_exact_merged_review_ready`, which is expected until browser-label, protocol/pooling, static coverage, and robustness flags are resolved.

## 2026-06-24: Validation-label bridge for motif and review calibration

Motivation: the feedback-ingest step created a stable readback path, but the current atlas also has older curated browser-validation notes. The next useful layer is a single label product that combines curated browser labels with future active-review feedback and then shows which label classes are still missing for motif/context calibration.

Changes:

- Added `scripts/dominant_rdg_validation_label_bridge.R`.
- Added the `rdg_validation_labels` pipeline step after `rdg_review_feedback` in the `rdg`, `manual`, and `all` scopes.
- Added declared outputs and cache dependencies for curated browser notes, motif-prior calibration, review feedback, and the linked review template.
- Updated the pipeline invariant test to require the new step.
- Updated the scientific note and next-model instructions to use the label-gap queue after each review pass.

Results:

- Unified the existing curated browser labels into `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_set.csv`.
- Current seed label set: 2 positive labels over 2 gene-designs, `DDIT3 | viral_infection` and `ATF3 | viral_infection`.
- No active-review feedback labels exist yet, so the balance table correctly reports missing weak positives, negatives, confounded rows, measurement-problem rows, and negative-control labels.
- Wrote a 970-row motif overlay and a 129-row label-gap queue; 126 active-review rows remain unlabeled after matching curated browser notes.
- The top gap rows now prioritize ATF4 and PPP1R15A viral-infection counterfactual/browser checks, then IFIH1/VEGFA and early negative-control rows.
- Pipeline invariant tests passed.
- Result sanity now reads 468 CSV files with zero read errors and no warnings/failures.

## 2026-06-24: Review-feedback ingest for linked active-review batches

Motivation: the linked active-review sheet made browser inspection practical, but the pipeline still had no stable readback step that turns completed manual review rows into label tables for calibration, browser-note curation, or later model fitting.

Changes:

- Added `scripts/dominant_rdg_review_feedback_ingest.R`.
- Added the `rdg_review_feedback` pipeline step after `rdg_review_links` in the `rdg`, `manual`, and `all` scopes.
- Added declared outputs and cache dependencies for the feedback step.
- Updated the pipeline invariant test to require the new step.
- Updated the result-sanity audit so empty feedback label/training/browser-note candidate files are expected until the linked review template contains manual labels.
- Updated the scientific note and next-model instructions with the feedback-readback workflow.

Results:

- The real linked template currently has 129 review rows and zero completed manual-feedback rows.
- Wrote stable feedback outputs under `dominant_rdg_review_feedback/`: header-only label, training, and browser-note candidate files; a 15-row summary; a six-row batch summary; and a Markdown report.
- A synthetic filled-template check produced 3 normalized labels and 3 training rows for positive, negative-control-pass, and absent-browser-signal cases.
- Pipeline invariant tests passed.
- Result sanity now reads 462 CSV files with zero read errors, no warnings/failures, and the three empty feedback CSVs annotated as expected.

## 2026-06-23: Direct Observatory links for review batches

Motivation: the active-review batches were balanced, but still required manual reconstruction of case/control subsets in RiboCrypt. The next practical step was to attach direct Observatory URLs and compact review cards to every selected row.

Changes:

- Added `scripts/dominant_rdg_review_batch_links.R`.
- Added the `rdg_review_links` pipeline step after `rdg_review_batches` in the `rdg`, `manual`, and `all` scopes.
- Added declared outputs and cache dependencies for the review-link step, including the metadata CSV.
- Updated the pipeline invariant test to require the new step.
- Updated the scientific note and next-model instructions to use the linked review template.

Results:

- Wrote direct Observatory URLs for all 129 active-review rows.
- With `/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv` mounted, 112 rows used exact context case/control subsets and 17 used study+condition fallback.
- Zero rows fell back to gene-only links.
- The preferred manual file is now `dominant_rdg_review_links/dominant_rdg_review_batch_template_with_links.csv`.
- A compact clickable HTML card page is available at `dominant_rdg_review_links/dominant_rdg_review_cards.html`.

## 2026-06-23: Active-review batches from evidence gaps

Motivation: the evidence-gap queue made the missing evidence explicit, but 1,680 rows are still too many for manual browser review. The next useful product is a balanced active-learning review plan with positives, counterfactual candidates, context audits, measurement-risk rows, replication follow-ups, and negative controls.

Changes:

- Added `scripts/dominant_rdg_active_review_batches.R`.
- Added the `rdg_review_batches` pipeline step after `rdg_evidence_gaps` in the `rdg`, `manual`, and `all` scopes.
- Added declared outputs and cache dependencies for the review-batch step.
- Updated the pipeline invariant test to require the new step.
- Updated the scientific note with the active-review-batch concept and next-review instructions.

Results:

- Wrote 129 active-review rows across 61 genes in 6 batches.
- Batch composition: 22 candidate-positive browser checks, 6 counterfactual-readiness checks, 28 context/confounder audits, 15 measurement/translon reviews, 28 replication/transport follow-ups, and 30 negative controls.
- The counterfactual overlay intentionally reuses ATF4, PPP1R15A, DDIT3, ATF3, VEGFA, and IFIH1 viral rows because predictive-readiness review is a different question from first-pass browser validation.
- Main outputs live in `dominant_rdg_review_batches/`; the preferred manual label sheet is `dominant_rdg_review_batch_template.csv`.

## 2026-06-19: RDG evidence-gap acquisition queue

Motivation: the atlas had many effect estimators, consensus rows, and validation dossiers, but it did not clearly answer what evidence should be acquired next. Another broad screen would mostly add volume; the more useful layer is a triage policy that separates browser validation, translon/measurement review, context replication, metadata-confounder audit, perturbation follow-up, and negative-control review.

Changes:

- Added `scripts/dominant_rdg_evidence_gap_prioritization.R`.
- Added the `rdg_evidence_gaps` pipeline step after `rdg_validation_dossiers` in the `rdg`, `manual`, and `all` scopes.
- Added declared outputs and cache dependencies for the evidence-gap step.
- Updated the pipeline invariant test to require the new step.
- Updated the scientific note with acquisition-priority and counterfactual-readiness interpretation.

Results:

- Ranked 1,680 atlas rows across 120 genes and 14 design families.
- Produced 71 high-value acquisition rows and 4 near-counterfactual candidate rows.
- Assigned 100 rows to browser validation, 7 to translon/measurement review, 668 to context/metadata audit, 904 to replication/robustness follow-up, and 1 to perturbation follow-up.
- Wrote a 250-row negative-control table for falsification review.
- Main outputs live in `dominant_rdg_evidence_gaps/`.

## 2026-06-18: Canonical FST page-list product

Motivation: the analysis directory had accumulated separate core, post-viral, expanded-OXPHOS, and literature-gap FST page lists. Those were useful while building the atlas, but they made the current product harder to use when preparing refreshed FST pages.

Changes:

- Replaced the three old FST-readiness/page-planning steps with one pipeline step: `fst_pages`.
- Made `scripts/dominant_cell_state_required_fst_pages.R` the single page-list generator for the full atlas target set.
- Renamed the curated auxiliary-prior input from `literature_gap_fst_page_priorities.csv` to `curated_additional_gene_priorities.csv`.
- Updated `scripts/dominant_cell_state_pack_fst_pages.R` so the default input is `required_fst_pages.csv`.
- Removed obsolete subset page-list scripts and generated files for post-viral-only, expanded-OXPHOS-only, literature-gap-only, and old dominant-cell-state page batches.

Results:

- Canonical product files are now `required_fst_pages.csv`, `required_fst_pages.txt`, `required_fst_page_basenames.txt`, `required_fst_pages_by_gene.csv`, `required_fst_selected_transcripts.csv`, and `required_fst_pages_summary.txt`.
- The regenerated page list contains 290 target genes, 290 selected transcripts, and 252 unique FST pages: 129 forward and 123 reverse.
- The selected transcript audit reports 283 canonical selections and 44 selected transcripts with at least one CDS-overlapping uORF.

## 2026-06-16: RDG table-only pipeline mode

Motivation: the expanded RDG graph step was spending most of its runtime regenerating overlay figures even when the analysis only needed refreshed tables, manifests, app inputs, or downstream model inputs. This made small iterations slow and encouraged unnecessary full-plot reruns.

Changes:

- Added `rdg_plot_mode` to `run_dominant_cell_state_pipeline()`.
- Added CLI support for `--skip-rdg-plots` and `--rdg-plot-mode all|tables`.
- In table-only mode, `dominant_ribosome_decision_graphs.R` regenerates RDG CSV/table outputs and relative-usage matrices but skips per-gene and per-state overlay plotting.
- Existing PNG/PDF plot manifests are preserved by scanning already written RDG output folders, so the app and downstream tables can still link to existing figures.
- Added `.runtime/rdg_plot_mode.txt` as an explicit cache dependency for the `rdg_graphs` step. This prevents full-plot and table-only modes from sharing an invalid cached fingerprint.
- Extended the internal invariant test to check CLI parsing and the RDG plot-mode cache dependency.

Results:

- A forced manual-scope `rdg_graphs` table-only run completed in about 204 seconds.
- The regenerated figure manifests referenced 288 existing general RDG figure files and 2,590 existing state-specific figure files.
- The relative-usage matrix manifest referenced 9 existing PNG and 9 existing PDF files.
- `dominant_cell_states/tests/integration/test_pipeline_invariants.R` passed after the change.

Interpretation:

- Table/model iterations can now avoid the expensive all-overlay redraw unless figures need to be refreshed.
- Use full plotting when manual translons, plotting code, colors, layout, or figure content need visual regeneration.
- Use table-only mode for scoring/model/app-data iterations where existing figures remain acceptable.

## 2026-06-16: Manual PELO/HK2 translon refresh and full manual-scope rerun

Motivation: the manual translon list was updated again, including new PELO and HK2 features that are not present in the current ORFik/Transcode union. The goal was to verify that the manual list is automatically incorporated and to rerun the downstream atlas layers with the refreshed RDG structures.

Changes:

- Validated `manual_translons/translons_to_be_added_to_manual.txt` through the pipeline manual-translon step.
- Selected 27 manual translons across 18 genes.
- Added PELO manual features on `ENST00000274311`: two leader uORFs and one overlapping uORF.
- Added HK2 manual features on `ENST00000290573`: leader/overlapping candidates plus one manual NTE-like feature.
- Reran `run_dominant_cell_state_pipeline(scope = "manual", force = FALSE, cache_version = "literature_gap_v1")`.
- Updated the scientific note with the refreshed atlas scale, PELO/HK2 interpretation, and current graph/motif model conclusions.

Results:

- Pipeline status: 47 successful steps and 2 expected cache-skipped FST-readiness checks; no failed registered steps.
- Manual translons: 27 selected; PELO has 3 selected manual features and HK2 has 4.
- The clean-CDS marker panel still measured 290 marker genes.
- Atlas scale: 1,680 gene-design rows over 120 genes and 14 design families.
- Consensus review queue: 1,285 rows.
- Validation dossiers: 184 primary review packets across 72 genes, with 380 matched branch contexts.
- Top next-model candidates are now IFIT1, DDIT4, EGLN3, MKI67, ATF3, CDK1, PELO, LMNB1, RPLP0, MMP2, HK2, and HSPA5.
- PELO is now rank 7 in the next-model table. Its viral row is context-specific review, not robust viral consensus.
- HK2 is rank 11. Its strongest current row is ISR/ER translation stress replicated robust consensus; its viral row is weak.

Interpretation:

- Manual translon recovery can materially change candidate ranking and branch interpretation.
- PELO is a ribosome-rescue/RQC review hypothesis, not yet a robust viral discovery.
- HK2 is more likely a glycolytic stress/hypoxia/ISR review hypothesis than a viral candidate in the current atlas.
- The strongest viral consensus biology remains DDIT3, PPP1R15A, ATF4, IFIH1, VEGFA, and the newer BBC3 branch-consensus row, with DDIT3 as the browser-supported anchor.
- The newest graph/motif tests strengthen the negative lesson: raw geometry and fixed hand-built motifs are useful annotations, but not yet predictive enough for counterfactual RDG forecasting.

## 2026-06-16: Literature-gap FST expansion and full downstream rerun

Motivation: the literature review identified missing marker genes for viral/post-viral fatigue, ISR kinase/initiation control, RNA sensing, host translation shutoff, ribosome quality control, antigen/complement context, DNA-damage/p53 context, and DecodeME-prior expression layers. The user downloaded the missing FST pages, so the pipeline could move from "page planning" to a real expanded run.

Changes:

- Added `literature_gap_fst_ready` as a first-class pipeline step.
- Expanded dominant-state target gene resolution to include literature-gap genes and aliases such as `CHOP -> DDIT3` and `DDX58 -> RIGI`.
- Added auxiliary modules for ISR kinase/initiation control, RNA sensing/RNase L/stress granules, host mRNA export/initiation shutoff, ribosome quality control/collision, and DecodeME-prior genes.
- Hardened translon category mapping and uORF feature measurement against invalid transcript-coordinate mappings.
- Fixed `dominant_uorf_regulation_inference.R` so mapped feature positions are restricted to finite 1-based coverage-matrix rows before `colSums()`.
- Reran `run_dominant_cell_state_pipeline(scope = "manual", force = FALSE, cache_version = "literature_gap_v1")` after the fix. Cached early steps were reused; downstream RDG/model/atlas layers reran.

Results:

- Literature-gap FST readiness: 118 required page rows, all local files present.
- Clean-CDS marker panel: 290 measured marker genes.
- Pipeline status: 13 cached steps and 36 successful downstream steps, no failed step.
- RDG plotting was the runtime bottleneck, about 22 minutes for the expanded state-specific figure set.
- Consensus atlas: 1,680 atlas rows, 120 genes, 130 model-agreement consensus rows, six viral consensus rows, and two viral transport-support-robust rows.
- Viral-infection review leads: DDIT3, ATF3, PPP1R15A, IFIH1, MYC, ATF4, MRPS12, VEGFA, TFE3, GABARAPL1, PPP1R15B, RRM1, G3BP1, DDIT4, and MAVS.
- Top consensus review queue: DDIT3 viral infection, PPP1R15A viral infection, LDHA ISR/stress, ATF4 viral infection, TFAM ISR/stress, EIF4G1 ISR/stress, IFIH1 viral infection, DDIT4 ISR/stress, LMNB1 ISR/stress, and CDK1 ISR/stress.
- Therapeutic hypothesis layer: ISR/eIF2B-ATF4 attenuation is the top preclinical perturbation hypothesis; JAK/STAT IFN attenuation is the top clinical-adjacent biomarker/stratification hypothesis.

Interpretation:

- The new FST pages materially expanded the atlas and strengthened the viral ISR/antiviral discovery test.
- DDIT3 is still the strongest browser-supported anchor. PPP1R15A, ATF4, IFIH1, ATF3, and VEGFA are the next viral candidates for browser validation.
- This run supports candidate discovery and preclinical hypothesis generation, not a medical-trial claim.
- Next engineering pass should add a pipeline flag to regenerate all RDG plots only when needed, because table/model reruns are much faster than full figure regeneration.

## 2026-06-12: RDG validation dossiers

Motivation: the consensus atlas solved prioritization, but browser review still required manually joining ranked candidates, branch contexts, warnings, and expected signal direction across several CSVs. This iteration adds a human-facing validation layer that turns the consensus queue into concrete review packets.

Changes:

- Added `scripts/dominant_rdg_validation_dossiers.R`.
- Added pipeline step `rdg_validation_dossiers` after `rdg_consensus`.
- Generates one full candidate table over all consensus rows, a selected primary-review set, top matched branch contexts, a manual review template, a gene summary, summary metrics, a Markdown review report, and PNG/PDF/HTML figures.
- The manual review template includes packet ID, transcript, design family, expected clean-CDS and branch pattern, first branch context, count/replicate gates, warnings, and blank columns for browser validation labels.
- Extended invariant tests to check dossier output existence, bounded validation priority scores, unique packet IDs, manual-template row consistency, primary rank continuity, branch-context mapping, and metric consistency.

Results:

- Candidate rows: 994.
- Primary browser-review packets: 113.
- Primary-review genes: 41.
- Matched top branch-context rows: 242.
- Viral primary-review rows: 21.
- Browser-supported primary rows: 2.
- Branch contexts passing the 30-count gate on both case and control sides: 237 of 242.
- Branch contexts with merged replicates on both sides: 242 of 242.
- First review packets are DDIT3 viral infection, PPP1R15A viral infection, ATF4 viral infection, IFIH1 viral infection, VEGFA viral infection, LDHA ISR/stress, and TFAM ISR/stress.

Interpretation:

- This is not new biological evidence; it is the practical bridge from atlas ranking to browser validation.
- The manual review template is now the preferred file for recording validated positives, weak positives, negatives, confounded rows, and not-reviewable rows.
- The next modeling improvement should consume those validation labels to recalibrate motif and branch priors instead of adding another independent score layer.

## 2026-06-12: RDG consensus atlas review layer

Motivation: the RDG atlas now has many partly overlapping evidence layers, including atlas-card confidence, branch redistribution, model agreement, robustness, context interactions, auxiliary-state warnings, browser validation, and therapeutic prioritization. This iteration adds a single ranked review queue that explains why a gene/design is prioritized and separates consensus, context-specific, compensation, therapeutic, fragile, and low-support rows.

Changes:

- Added `scripts/dominant_rdg_consensus_atlas.R`.
- Added pipeline step `rdg_consensus` after `therapeutic_perturbation_predictions`.
- Joins atlas cards with auxiliary-state atlas annotations, therapeutic gene evidence, and perturbation-prediction summaries.
- Computes bounded evidence components for atlas confidence, branch support, model agreement, leave-one-study/grouped-omission robustness, browser support, effect size, context specificity, auxiliary-state context, therapeutic relevance, sparse/on-off evidence, novelty, fragility, and warning penalties.
- Writes consensus candidates, review queue, evidence matrix, gene summary, metrics, and PNG/PDF/HTML figures.
- Extended the internal invariant tests to check consensus output existence, rank ordering, bounded scores, review queue consistency, evidence-matrix component ranges, and metric consistency.

Results:

- Atlas rows compressed: 994.
- Consensus review rows: 729.
- Genes represented: 71.
- Candidate classes: 1 browser-validated consensus, 28 replicated robust consensus, 45 multi-model branch consensus, 13 context-specific review, 9 buffering/compensation review, 93 therapeutic-translation review, 230 general atlas review, 304 fragile/confounded review, and 271 low-support background rows.
- Top consensus candidates are DDIT3 viral infection, PPP1R15A viral infection, LDHA ISR/stress, ATF4 viral infection, TFAM ISR/stress, IFIH1 viral infection, LMNB1 ISR/stress, CDK1 ISR/stress, DDIT4 ISR/stress, and ATF4 ISR/stress.
- Only DDIT3 viral infection currently exceeds a consensus score of 0.70; 24 rows are at or above 0.50.

Interpretation:

- The consensus score is a deterministic review-priority score, not a p-value and not independent biological validation.
- DDIT3 is the current browser-validated consensus anchor.
- PPP1R15A and ATF4 remain strong viral branch candidates; IFIH1 remains the important context-dependent lung viral hypothesis.
- The practical next review step is to inspect rows from the consensus queue and feed validated positives and negatives back into motif calibration and the future joint multilevel branch model.

## 2026-06-06: Auxiliary-state RDG context audit

Motivation: the auxiliary dominant-state screen identified candidate new state axes, but the atlas still did not say whether a given RDG branch effect occurred inside a strong auxiliary-state context. This iteration adds a warning/prioritization layer that joins per-run auxiliary scores back to matched RDG case/control strata, without reloading FST pages.

Changes:

- Added `scripts/dominant_rdg_auxiliary_state_audit.R`.
- Added pipeline step `auxiliary_state_rdg_audit` after `rdg_atlas_cards`.
- Reconstructs matched case/control sample sets from `study`, `CONDITION`, `AUTHOR`, `CELL_LINE`, `TISSUE`, `GENE`, `INHIBITOR`, `FRACTION`, `Cancer_type`, `Sex`, and `TIMEPOINT`.
- Computes auxiliary-state case/control z deltas for each matched RDG stratum.
- Joins those deltas to joint branch-allocation rows and atlas cards.
- Writes stratum scores, compact review branch rows, gene/design/state summaries, atlas annotations, selected-state matrices, state summaries, metrics, and PNG/PDF/HTML figures.
- Weights atlas audit priority by auxiliary-state readiness, so pipeline-ready states dominate review warnings while undercovered marker panels remain visible but discounted.
- Added invariant tests for the new registry step, output files, matched-stratum coverage, compact branch-row filtering, required annotation columns, bounded priorities, and screen-ready auxiliary-state annotations.

Results:

- RDG strata: 159.
- Strata with reconstructed auxiliary case/control scores: 86.
- Full branch/state denominator: 417,216 rows.
- Compact branch audit rows written: 55,176.
- Gene/design/state summaries: 4,592.
- Atlas rows annotated with an auxiliary-state top hit: 497 of 994.
- Atlas rows with strong auxiliary-state overlap: 112.
- Atlas rows with moderate auxiliary-state overlap: 7.
- Raw overlap winner: NRF2, but it is still marker-underpowered and needs FST expansion.
- Readiness-weighted overlap winner: mitochondrial UPR / mitonuclear stress.
- Highest weighted rows include DDIT3, MYC, ATF4, and TRIB3 viral-infection contexts, especially `PRJNA957283` BT20 infected versus WT.
- Antigen presentation is the strongest screen-ready non-viral warning layer, especially for tumor/cancer contexts.

Interpretation:

- This is a context/confounder audit, not a causal adjustment.
- Mitochondrial-stress overlap strengthens the idea that viral/post-viral RDG effects are entangled with mitonuclear/ISR/OXPHOS state.
- Antigen-presentation overlap gives a practical warning that immune composition can drive or mask tumor/infection atlas rows.
- NRF2/proteostasis/autophagy raw overlap should motivate marker/FST expansion before core promotion.
- Next useful browser checks are DDIT3, ATF4, TRIB3, and MYC viral-infection rows in BT20 infected context, followed by antigen-shifted tumor rows such as PPP1R15A, ATF4, MMP2, IGFBP5, OAS2, and TFAM.

## 2026-06-06: Auxiliary dominant-state discovery screen

Motivation: adding more dominant states could improve discovery, but adding them directly to the core clean-CDS scorer would expand FST requirements and risk creating aliases of ISR, IFN, OXPHOS, or post-viral fatigue. This iteration adds a cheap pre-screen that tests new candidate state families using the existing clean-CDS expression matrix before committing to FST-heavy marker expansion.

Changes:

- Added `scripts/dominant_state_auxiliary_module_screen.R`.
- Added pipeline step `auxiliary_state_screen` after `score_audit` and before the metadata/RDG layers.
- The screen scores eight candidate state families: proteostasis/heat shock, autophagy/lysosome, NRF2/oxidative-ferroptosis defense, DNA damage/p53 checkpoint, antigen presentation/immune composition, secretory/ER proteostasis, mitochondrial UPR/mitonuclear stress, and host-shutoff/antiviral translation repression.
- It computes per-run auxiliary scores, marker availability, redundancy against existing dominant-state scores, categorical/numeric metadata associations, RDG atlas gene links, recommendation classes, summary metrics, and PNG/PDF/HTML figures.
- Added invariant tests for the new registry step, output files, bounded novelty scores, marker availability, and expected antigen/p53 candidate calls.
- Updated the scientific note with the literature rationale and current interpretation.

Results:

- Candidate modules screened: 8.
- Candidate marker rows: 125.
- Unique candidate genes already available in the current clean-CDS matrix: 44.
- Pipeline-candidate screeners: DNA damage/p53 checkpoint, antigen presentation/immune composition, and mitochondrial UPR/mitonuclear stress.
- Modules needing targeted FST marker expansion: proteostasis/heat shock, autophagy/lysosome, and NRF2/oxidative-ferroptosis defense.
- Secretory/ER proteostasis is currently mostly an ISR-axis alias.
- Host-shutoff/antiviral translation repression is fully scoreable, but it is close to the existing post-viral fatigue/ribosome-stress score and remains exploratory.
- The strongest clean new signal is antigen presentation/immune composition: LCL samples all fall in the top decile for that auxiliary score, and blood/hematopoietic/lymphoma contexts are strongly enriched.

Interpretation:

- Yes, new dominant states can increase discovery power, but they should enter as screened confounder/context variables first.
- Antigen presentation is the most useful immediate addition because immune composition can masquerade as viral or inflammatory translation regulation.
- DNA damage/p53 may separate checkpoint arrest from ordinary low proliferation, but many top metadata hits are perturbation/study specific and need matched-control interpretation.
- Mitochondrial UPR is relevant for post-viral fatigue and energy stress, but it overlaps ISR/OXPHOS and should not yet be treated as an independent mechanism.
- The next practical step is targeted marker-page expansion for proteostasis, autophagy/lysosome, and NRF2/ferroptosis if those states become biologically important.

## 2026-06-05: Partial-pooling evidence integrated into atlas cards and browser

Motivation: the partial-pooled context model was useful as a side table, but browser review and candidate ranking still needed those evidence tiers directly on each atlas card. The atlas should distinguish shrinkage-supported context shifts from single-context and fragile raw context signals without requiring manual joins.

Changes:

- Merged `dominant_rdg_context_partial_pooling_gene_design_summary.csv` into `dominant_rdg_atlas_cards.csv`.
- Added explicit `partial_pooling_*` fields, a bounded `partial_pooling_component`, and conservative contribution to atlas confidence/review-priority scores.
- Added atlas warnings for partial-pooled context support, pooled review, single-context raw signals, and fragile-after-pooling rows.
- Added partial-pooling fields to branch-context exports and atlas summary metrics.
- Updated the pipeline registry so `rdg_atlas_cards` depends on `context_partial_pooling`.
- Updated the Shiny atlas browser with a `Partial-pooling evidence only` filter, summary counters, hover text, selected-card fields, and a dedicated `Partial Pooling` table.
- Extended internal invariants to parse the Shiny app and check partial-pooling atlas/browser integration.

Results:

- Atlas-level partial-pooling rows: 2 pooled-supported, 4 pooled-review, 9 single-context, and 3 fragile-after-pooling.
- All current atlas-level partial-pooling evidence rows are viral-infection rows.
- ATF4 and TRIB3 are the current shrinkage-supported viral context examples.
- IFIH1 remains a high-priority lung viral hypothesis, but it is now correctly labeled as `single_context_raw_signal_not_pooled`.
- Atlas priority scores remain finite, and the invariant suite passes after the integration.

Interpretation:

- The browser can now be used as the primary review surface for partial-pooling tiers.
- ATF4/TRIB3 should be treated as shrinkage-positive context examples; IFIH1 should be treated as a strong lung-focused single-context hypothesis until sibling context evidence exists.
- This keeps the atlas from overclaiming single-context effects while still making them easy to inspect.

## 2026-06-05: Partial-pooled context-interaction review

Motivation: explicit context interactions were useful, but they treated each tissue/cell-line contrast as a separate row. The next iteration needed a shrinkage-aware layer that distinguishes context effects supported across sibling context values from strong but single-context-only evidence.

Changes:

- Added `scripts/dominant_rdg_context_partial_pooling_model.R` and pipeline step `context_partial_pooling`.
- Reuses `dominant_rdg_context_interaction_branch_context_effects.csv`; no FST-heavy quantification is rerun.
- Applies empirical-Bayes shrinkage across context values within each gene, transcript, design family, branch class, and context axis.
- Classifies rows as partial-pooled supported shifts, pooled review, single-context supported/review, raw-signal fragile after pooling, pooled-emergent review, or neutral.
- Writes row, gene-design, viral-review, raw-versus-pooled, and summary-metric CSVs plus top-effect and IFIH1 lung PNG/HTML figures.
- Added pipeline registry outputs/dependencies and invariant tests for bounded q-values, bounded priorities, bounded pooling weights, supported pooled rows, and the IFIH1 lung single-context classification.

Results:

- Evaluable context rows: 908.
- Gene-designs: 72.
- Partial-pooled supported rows: 5.
- Partial-pooled review rows: 17.
- Raw supported rows: 51.
- Raw supported rows surviving partial pooling: 20.
- Raw supported rows shrunk to fragile: 12.
- Single-context raw signals: 19.
- Pooled-emergent review rows: 2.
- Strict pooled viral rows are ATF4 overlapping-uORF context shifts in Huh7/liver/lung and TRIB3 clean-CDS shifts in lung/A549.
- IFIH1 lung leader-uORF and clean-CDS remain the top review rows, but they are now labeled `single_context_supported_not_pooled`; the shrinkage-aware estimate equals the raw interaction because no sibling tissue contexts are available for shrinkage.

Interpretation:

- IFIH1 remains a strong lung-focused viral RDG browser-review hypothesis, but it should not be called a partial-pooled context effect yet.
- ATF4 and TRIB3 are the current viral examples that survive the new shrinkage-aware context tier.
- This layer is still an approximation. The next statistical step is a true multilevel branch model with study, gene, design, branch, context, and motif-prior terms rather than more independent review screens.

## 2026-06-05: Allocation-versus-output disagreement audit

Motivation: the model-agreement layer identified nine rows where net clean-CDS output and clean-CDS branch allocation point in opposite directions. These could be real buffering/total-load compensation, compensatory reinitiation, or aggregation artifacts. The next iteration needed a compact audit before treating them as biology.

Changes:

- Added `scripts/dominant_rdg_allocation_output_disagreement_audit.R` and pipeline step `allocation_output_disagreements`.
- Reuses atlas cards, model-agreement summaries, atlas branch-context rows, and study evidence; no FST pages are reloaded.
- Computes clean-CDS output change, clean-CDS allocation delta, upstream allocation delta, a total-load compensation proxy, a compensation evidence score, and a heterogeneity risk score.
- Classifies rows as strong compensation candidates, possible compensation, weak-output/strong-allocation review, aggregation/heterogeneity risk, or low-magnitude review.
- Writes audit, review, context, study-evidence, summary-metric CSVs, plus scatter and branch-heatmap PNG/HTML figures.
- Added pipeline registry outputs/dependencies and invariant tests for bounded scores plus SESN2/ATF3 row presence.

Results:

- Nine allocation-versus-output disagreement rows were found.
- SESN2 nutrient starvation is the strongest buffering/total-load compensation candidate: clean-CDS output up about 94%, clean-CDS allocation down about 0.094, upstream allocation up about 0.083.
- ATF3 nutrient starvation is a possible compensation candidate: output up about 95%, clean-CDS allocation down about 0.046, upstream allocation up about 0.056, with higher heterogeneity risk.
- SNAI2 ISR/stress is a weak-output/strong-allocation review row: output up only about 4%, but clean-CDS allocation down about 0.154 and upstream allocation up about 0.150.
- RRM1 and PCNA are aggregation/heterogeneity-risk rows and should not be treated as primary discoveries without browser support.

Interpretation:

- The audit makes the disagreement rows more useful: SESN2/ATF3 are now the top browser-review targets for compensation, while RRM1/PCNA should be down-weighted.
- A true compensation claim still requires browser validation and ideally direct total-RPF or protein-output evidence.
- The next model should connect this audit to a total-ribosome-load/branch-allocation decomposition rather than only clean-CDS output.

## 2026-06-05: Quantitative therapeutic perturbation predictions

Motivation: the previous therapeutic layer ranked intervention hypotheses but did not convert them into gene-level effect targets. For experimental planning, the useful output is a predicted CDS-buffering and branch-allocation direction, with uncertainty and assay caveats.

Changes:

- Added `scripts/dominant_rdg_therapeutic_perturbation_predictions.R` and pipeline step `therapeutic_perturbation_predictions`.
- Reads existing atlas cards and therapeutic-hypothesis scores only; no FST pages are reloaded.
- Defines disease-like design families per hypothesis and computes explicit full-normalization and half-normalization CDS-buffering targets.
- Carries branch-allocation deltas for clean CDS, leader uORF, and overlapping uORF into the perturbation target table.
- Writes gene-level predictions, hypothesis summaries, compact review rows, summary metrics, and PNG/HTML figures.
- Added pipeline registry outputs/dependencies and invariant checks for bounded confidence scores and core ISR perturbation targets.

Results:

- ISR/eIF2B-ATF4 attenuation remains the strongest perturbation hypothesis, prediction score about 0.703.
- Quantitative-review-ready ISR targets: DDIT3, ATF4, and PPP1R15A.
- Half-normalization CDS targets: DDIT3 about -46%, ATF4 about -42%, PPP1R15A about -63%.
- JAK/STAT IFN attenuation ranks second, prediction score about 0.631.
- IFIH1 is the strongest IFN-axis perturbation target, with about -20% half-normalization CDS target.
- IFIH1/MDA5 modulation ranks third as target-discovery/stratification, not as a direct treatment claim.

Interpretation:

- The half-normalization target is an explicit experimental assumption, not a measured drug efficacy estimate.
- The immediate experimental question is whether ISR or IFN-axis perturbation moves RDG branch allocation and CDS output toward these directions.
- The next statistical improvement is to learn the response fraction from actual perturbation contexts instead of assuming 30-50% normalization.

## 2026-06-05: Therapeutic-hypothesis prioritization layer

Motivation: the current RDG atlas can nominate biology, but the next practical question is which intervention concepts are worth testing ex vivo or using as biomarker-stratified follow-up hypotheses. This should be reproducible and clearly separated from clinical treatment claims.

Changes:

- Added `scripts/dominant_rdg_therapeutic_hypotheses.R` and pipeline step `therapeutic_hypotheses`.
- Uses existing atlas, model-agreement, transportability, and context-interaction summaries only; no FST pages are reloaded.
- Defines reviewable hypothesis classes: ISR/eIF2B-ATF4 attenuation, JAK/STAT IFN attenuation, IFIH1/MDA5 modulation, HIF/hypoxia modulation, mTOR/translation-capacity modulation, mitochondrial/OXPHOS support, and antigen/complement inflammation.
- Scores each hypothesis from atlas evidence, viral specificity, consensus/transport/context support, browser support, practical readiness, mechanistic alignment, and safety/risk penalty.
- Writes gene evidence, ranked score, compact review, summary-metric CSVs, plus score and evidence-matrix PNG/HTML figures.
- Added pipeline registry outputs/dependencies and invariant tests for score-table presence, bounded scores, and the two main ISR/JAK hypotheses.

Results:

- ISR/eIF2B-ATF4 attenuation ranks first, score about 0.827, as a preclinical mechanistic perturbation hypothesis. Top genes include DDIT3, ATF4, PPP1R15A, DDIT4, ATF3, HSPA5, and TRIB3.
- JAK/STAT IFN attenuation ranks second, score about 0.758, as a clinical-adjacent biomarker hypothesis. Top genes include IFIH1, IRF7, STAT1, CCL2, MX1, IFIT1, and OAS2.
- IFIH1/MDA5 axis modulation ranks third, score about 0.586, as a target-discovery/stratification hypothesis rather than a direct treatment claim.
- HIF/hypoxia and mitochondrial/OXPHOS hypotheses remain exploratory but visible in the component matrix.

Interpretation:

- The layer does not identify a trial-ready treatment. It gives a reproducible perturbation-priority table.
- The best near-term experiment is to test whether ISR or IFN-axis perturbations shift RDG branch allocation and CDS-buffering levels in the predicted direction.
- The output should eventually become quantitative: predicted branch shift, predicted CDS-buffering percentage, uncertainty, validation assay, and safety caveat per perturbation.

## 2026-06-05: Explicit viral context-interaction branch model

Motivation: grouped omission showed that IFIH1 reverses after removing all lung studies, but omission tests do not directly estimate the lung branch effect against the non-lung complement. The next iteration needed an explicit tissue/cell-line context comparison without reloading FST pages.

Changes:

- Added `scripts/dominant_rdg_context_interaction_model.R` and pipeline step `context_interactions`.
- Reuses `dominant_rdg_joint_branch_allocation_class_rows.csv`; no FST-heavy quantification is rerun.
- Collapses row-level branch deltas within study and tissue/cell-line context before comparing each context value with its non-matching complement.
- Calls `context_direction_reversal`, `context_specific_supported_shift`, `context_magnitude_modulated`, `broad_context_consistent`, weak, and non-evaluable classes.
- Writes study-context effects, branch-context effects, review rows, gene-design summaries, viral consensus summaries, summary metrics, and IFIH1/viral heatmap figures.
- Integrated context-interaction fields into atlas cards, atlas warnings, summary metrics, review priority, branch-context exports, and the Shiny atlas browser.
- Added a Shiny `Context interaction only` filter and selected-card `Context Interactions` table.
- Extended invariant tests to require the IFIH1 lung clean-CDS/leader-uORF direction-reversal rows and atlas integration.
- Fixed a context-summary display issue so the top reported context row prioritizes actual context-interaction classes over weaker rows with high inherited review priority.

Results:

- Study-context rows: 5,560.
- Branch-context effects: 4,020.
- Review rows: 422.
- Viral gene-design summaries: 77.
- Context-evaluable viral gene-designs: 72.
- Prioritized gene-level context classes: 13 direction reversals, 1 context-specific supported shift, and 4 magnitude-modulated shifts.
- Atlas rows with context-interaction classes: 18.
- Viral consensus rows with context-direction-reversal gene class: IFIH1 and ATF4.
- IFIH1 lung leader-uORF branch: context delta about -0.092, non-lung complement about +0.028, interaction delta about -0.120, `q` about 0.006.
- IFIH1 lung clean-CDS branch: context delta about +0.091, non-lung complement about -0.013, interaction delta about +0.104, `q` about 0.009.

Interpretation:

- IFIH1 is now the strongest explicit tissue-dependent viral RDG hypothesis: lung infection shifts allocation away from leader-uORF and toward clean CDS relative to non-lung infected complements.
- DDIT3 and PPP1R15A remain the stronger cross-context anchors by grouped-omission support robustness.
- ATF4 and VEGFA remain support-sensitive/context-modulated rather than universal viral branch effects.
- The context-interaction layer is a review screen, not a causal model; context values can still mix tissue, cell line, virus, time point, protocol, and publication effects.
- The full invariant suite passes after regenerating context interactions and atlas cards.

## 2026-06-04: Grouped-omission transportability and context dependence

Motivation: leave-one-study-out showed whether one accession drove a result, but related studies can share tissue, cell line, protocol, publication, or lab effects. A candidate suitable for a predictive RDG atlas should be tested against omission of those related families.

Changes:

- Added `scripts/dominant_rdg_model_agreement_transportability.R` and pipeline step `model_agreement_transport`.
- Grouped studies by publication, author/lab, coarse protocol, inhibitor, tissue, cell line, and footprint-length family.
- Refit supported branch effects after omitting groups containing at least two studies, while requiring at least two independent studies to remain.
- Separated grouped-omission support robustness, direction-only robustness, direction fragility, and non-evaluability.
- Added explicit branch and omission-axis coverage fractions to prevent sparsely testable candidates from receiving broad-transportability scores.
- Added grouped-omission evidence to atlas confidence/review priority, warnings, summary metrics, and preferred card columns.
- Added a grouped-omission filter and selected-card `Transport Evidence` table to the Shiny atlas browser.
- Fixed the viral atlas priority figure so browser-support circles are drawn only for browser-validated rows.
- Added invariant tests for omission-size gates, bounded stability scores, IFIH1 lung-context regression behavior, pipeline wiring, and atlas integration.

Results:

- Consensus gene-designs: 74.
- Grouped-omission-evaluable consensus gene-designs: 19.
- Consensus support-robust: 5.
- Consensus direction-robust/support-fragile: 10.
- Consensus direction-fragile: 4.
- Consensus non-evaluable: 55.
- All five viral consensus rows are grouped-omission evaluable.
- Viral support-robust: DDIT3 and PPP1R15A.
- Viral direction-robust/support-fragile: ATF4 and VEGFA.
- Viral direction-fragile: IFIH1.
- Median evaluable omission-axis coverage is three of seven axes; DDIT3 and PPP1R15A are each currently evaluable on three axes.
- IFIH1 reverses direction only after all four lung studies are removed; the remaining three non-lung estimates are weak and cross zero.

Interpretation:

- Grouped omission is much stricter than leave-one-study-out and exposes context dependence hidden by accession-level influence tests.
- DDIT3 and PPP1R15A are the strongest current cross-context viral branch candidates.
- IFIH1 should be treated as a lung-dependent viral RDG hypothesis, not a universal viral effect or proof of an opposite non-lung effect.
- Publication/author axes are often non-evaluable for viral candidates because independent publications/labs remain limited.
- The next model should estimate explicit gene-by-context-by-condition branch effects.

## 2026-06-04: Leave-one-study-out agreement stability and study evidence

Motivation: estimator agreement was useful but could still be driven by one influential study. The atlas also needed an inspectable per-study explanation behind each agreement tier.

Changes:

- Added `scripts/dominant_rdg_model_agreement_loo.R` and pipeline step `model_agreement_loo`.
- Preserved study-level effects from the hierarchical-joint and posterior-predictive branch models.
- Refit each evaluable branch/model after omitting one study.
- Separated support-robust, direction-robust/support-fragile, and single-study/non-evaluable calls.
- Added compact per-study evidence and LOO filters/details to the Shiny atlas browser.
- Integrated LOO stability into atlas confidence only for robust consensus; fragility/non-evaluability increase review priority instead.
- Added pipeline invariants for study effects, LOO outputs, bounded stability scores, and atlas integration.

Results:

- Consensus gene-designs: 74.
- Study-omission-evaluable consensus gene-designs: 29.
- Consensus support-robust: 14.
- Consensus direction-robust/support-fragile: 15.
- Consensus single-study/non-evaluable: 45.
- Study-evaluable consensus direction failures: 0.
- All five viral consensus candidates are study-evaluable.
- DDIT3, IFIH1, ATF4, and PPP1R15A retain direction and support after every omission.
- VEGFA retains direction after every omission but loses thresholded support after some omissions.

Interpretation:

- The viral consensus set is not driven by one accession.
- LOO remains an influence analysis, not independent validation.
- Support retention uses within-branch p-values as a practical proxy rather than recomputing global BH correction after every omission.
- Single-study calls are non-evaluable, not evidence of fragility.

## 2026-06-04: Compact RDG model-agreement layer added

Motivation: the atlas had several useful but overlapping queues: clean-CDS effects, joint allocation, study-aware hierarchical allocation, posterior-predictive branch allocation, motif alignment, and manual browser validation. Review needed one compact table that distinguishes estimator consensus from disagreement without treating every model output as independent evidence.

Changes:

- Added `scripts/dominant_rdg_model_agreement.R`.
- Added pipeline step `model_agreement` after `dirichlet_multinomial_allocation` and before `rdg_atlas_cards`.
- Reuses existing summary CSVs; no FST pages are loaded.
- Compares hierarchical-joint and posterior-predictive deltas per branch.
- Creates explicit gene/design tiers:
  - `A_browser_validated_consensus`;
  - `B_replicated_branch_consensus`;
  - `C_multi_model_consensus`;
  - `D_model_disagreement_review`;
  - `E_single_model_review`;
  - `F_low_support`.
- Separates opposite branch-model direction from opposite clean-CDS allocation versus clean-CDS output.
- Merged agreement scores, tiers, support signatures, and disagreement flags into atlas cards, warnings, review priority, summary metrics, and the Shiny browser.
- Added model-agreement filtering and a disagreement-only filter to the Shiny browser.
- Added invariant tests for agreement outputs, bounded scores, pipeline wiring, and atlas integration.

Outputs:

- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_branch_comparison.csv`
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_branch_summary.csv`
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_gene_design.csv`
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_review.csv`
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_viral_review.csv`
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_summary_metrics.csv`
- `dominant_rdg_model_agreement/figures/model_agreement_branch_comparison.*`
- `dominant_rdg_model_agreement/figures/model_agreement_viral_matrix.*`

Result:

- Gene-designs: 994 over 71 genes.
- Supported gene-designs: 320.
- Browser-validated consensus: 1.
- Replicated branch consensus: 22.
- Additional multi-model consensus: 51.
- Model-disagreement review: 9.
- Single-model review: 237.
- Low support: 674.
- Hierarchical-joint and posterior-predictive models both support 138 branch comparisons:
  - clean CDS: 69;
  - leader uORF: 36;
  - overlapping uORF: 33.
- All 138 common-supported branches agree in direction.
- Overall branch-delta Pearson correlations are above 0.99 for all three branch classes.
- Viral consensus gene-designs: DDIT3, IFIH1, VEGFA, ATF4, PPP1R15A.
- Viral model-disagreement rows: 0.
- The nine disagreement-tier rows are all clean-CDS allocation-versus-output disagreements: SNAI2, RRM1, SESN2, PCNA, ATF3, IFIH1, DDIT4, TOMM20, and ZEB2 in non-viral design families.

Interpretation:

- The high hierarchical/posterior-predictive agreement demonstrates estimator stability, not independent validation, because both models reuse the same grouped branch counts.
- PPP1R15A is no longer best described as a branch-direction disagreement. Its posterior-predictive clean-CDS delta is positive like the hierarchical delta, but it fails the posterior-predictive support threshold; it therefore remains multi-model rather than replicated branch consensus.
- Allocation-versus-output disagreement may reveal buffering or total-ribosome-load compensation and is a useful mechanistic review class.
- The new viral evidence matrix makes the current consensus set and weaker single-model candidates immediately visible.
- The model-agreement and atlas steps were rerun; full pipeline invariant tests and Shiny app sourcing pass.

## 2026-06-03: Empirical-Bayes Dirichlet-multinomial branch-allocation layer added

Motivation: the study-aware hierarchical joint layer asks whether branch deltas replicate across studies, but it does not directly ask whether a case branch count is surprising given the control counts and a background RDG prior. The next iteration needed a more likelihood-like posterior-predictive screen without reloading FST pages.

Changes:

- Added `scripts/dominant_rdg_dirichlet_multinomial_branch_model.R`.
- Added pipeline step `dirichlet_multinomial_allocation` after `hierarchical_joint_allocation` and before `rdg_atlas_cards`.
- Reuses existing joint allocation rows; no new FST-heavy quantification is required.
- Tests branch classes `leader_uORF`, `overlapping_uORF`, `clean_CDS`, and `other`.
- Builds weak branch priors from same-gene other-study rows, same-gene all-study rows, same-design other-study rows, or global control background, then combines those priors with row control counts.
- Uses beta-binomial posterior-predictive branch marginals as a tractable Dirichlet-multinomial approximation.
- Aggregates row-level calls within study and then across studies into gene/design branch effects.
- Merged Dirichlet-multinomial support into atlas cards, atlas confidence, browser priority, app table fields, effect-space hover text, selected-card details, and summary metrics.
- Added invariant tests for posterior-predictive branch rows, branch effects, gene-design summaries, priority bounds, and atlas-card integration.

Outputs:

- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_branch_rows.csv`
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_branch_effects.csv`
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_gene_design_summary.csv`
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_review.csv`
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_viral_review.csv`
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_summary_metrics.csv`
- `dominant_rdg_dirichlet_multinomial/figures/dirichlet_multinomial_top_effects.*`
- `dominant_rdg_dirichlet_multinomial/figures/dirichlet_multinomial_viral_review.*`
- `dominant_rdg_dirichlet_multinomial/figures/dirichlet_multinomial_delta_heatmap.*`

Result:

- Branch effects: 1,050.
- Gene-designs: 522.
- Genes: 82.
- Replicated posterior-predictive shift branches: 47.
- Replicated-review branches: 40.
- Single-study strong branches: 113.
- Viral branch effects: 158 over 77 gene-designs.
- Viral replicated shift branches: 7.
- Viral replicated-review branches: 3.
- Viral median absolute branch delta: about 0.0045.
- Viral median `I2`: about 0.87, so heterogeneity remains high.
- Atlas rows with Dirichlet-multinomial support: 103.
- Atlas replicated Dirichlet-multinomial rows: 43.
- Viral atlas rows with Dirichlet-multinomial support: 6, all replicated or replicated-review.
- Viral supported gene-design list: DDIT3, ATF4, VEGFA, MRPS12, IFIH1, SNAI2.

Interpretation:

- This is a posterior-predictive screen, not a full multilevel Bayesian model.
- DDIT3, ATF4, VEGFA, and IFIH1 now appear in both conservative/study-aware and posterior-predictive viral support lists.
- PPP1R15A remains supported by the study-aware hierarchical joint layer, but the posterior-predictive layer demotes it because direction stability is weaker across the current contexts. That makes PPP1R15A a useful disagreement case for browser and denominator review.
- The viral posterior-predictive review figure now fades neutral/review rows so replicated support is visually distinct.
- Full pipeline invariant tests pass with the new step.

## 2026-06-03: Study-aware hierarchical joint allocation layer added

Motivation: row-level joint branch-allocation is good for finding local branch shifts, but broad designs can mix many contexts. PABPC1 tumor/kidney is the motivating example: one kidney tumor row has a visible leader-uORF/clean-CDS allocation shift, while the broad tumor-cancer design should remain near neutral if other tumor contexts do not agree. The pipeline needed a middle layer between row-level branch contexts and atlas-wide cards.

Changes:

- Added `scripts/dominant_rdg_hierarchical_joint_allocation_model.R`.
- Added pipeline step `hierarchical_joint_allocation` after `joint_branch_allocation` and before `rdg_atlas_cards`.
- Collapses joint-allocation rows within study, then estimates study-aware random-effects branch deltas for each gene/transcript/design/branch.
- Saves branch-effect, gene-design, review, viral-review, metrics, and figure outputs under `dominant_rdg_hierarchical_joint_allocation/`.
- Merged hierarchical joint columns into atlas cards, app hover text, selected-card details, atlas confidence/review priority, and atlas summary metrics.
- Added invariant tests for the new model outputs and atlas-card integration.
- Adjusted the new heatmap so missing branch classes do not visually mean zero.

Outputs:

- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_branch_effects.csv`
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv`
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_review.csv`
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_viral_review.csv`
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_summary_metrics.csv`
- `dominant_rdg_hierarchical_joint_allocation/figures/hierarchical_joint_allocation_top_effects.*`
- `dominant_rdg_hierarchical_joint_allocation/figures/hierarchical_joint_allocation_viral_review.*`
- `dominant_rdg_hierarchical_joint_allocation/figures/hierarchical_joint_allocation_delta_heatmap.*`

Result:

- Branch effects: 2,088.
- Gene-designs: 522.
- Genes: 82.
- Replicated shift branches: 64.
- Replicated-review branches: 5.
- Single-study strong branches: 84.
- Viral branch effects: 308 over 77 gene-designs.
- Viral replicated branch shifts: 9.
- Atlas rows with hierarchical joint support: 79.
- Atlas replicated hierarchical joint rows: 36.
- Viral atlas rows with hierarchical joint support: 6, all replicated.
- Replicated viral gene-design list: PPP1R15A, DDIT3, MMP2, ATF4, VEGFA, IFIH1.

Interpretation:

- This is a useful conservative ranking layer, not a final Bayesian RDG likelihood.
- The broad PABPC1 tumor-cancer design is correctly neutral in the study-aware layer, while the PABPC1 tumor/kidney branch-context row remains available for browser review.
- Viral infection now has a smaller replicated branch-allocation list suitable for the next browser-validation pass.
- The top all-design plot is still dominated by single-study UV/heat/stress rows, so support class must be considered before biological interpretation.
- Full pipeline invariant tests pass with the new step.

## 2026-06-03: Compact atlas branch-context table added

Motivation: the branch-aware Shiny browser initially loaded the full joint branch-allocation row table to populate selected-card contexts. That works, but it is unnecessary for review and makes the app depend on a large model-internal table. The atlas step should own a compact app-facing context export.

Changes:

- Updated `scripts/dominant_rdg_atlas_cards.R`.
- Added `dominant_rdg_atlas/dominant_rdg_atlas_branch_contexts.csv`.
- Added this file to the `rdg_atlas_cards` pipeline outputs.
- Updated `apps/rdg_atlas_browser/app.R` to prefer the compact branch-context CSV, with fallback to the full joint rows for older output folders.
- Tightened the app's selected-card branch-context lookup to `gene_symbol + tx_id + design_family`.
- Added invariant checks for the compact branch-context file.

Result:

- Compact branch-context rows: 11,289.
- Genes represented: 71.
- The PABPC1 tumor/kidney row is retained as a selected-card context:
  - `PABPC1`, `ENST00000318607`, `tumor_cancer_context`;
  - context rank 4;
  - `PRJNA256316-homo_sapiens | primary | kidney | tumor | vs WT`;
  - strict joint branch shift;
  - leader-uORF allocation delta about `+0.086`;
  - clean-CDS allocation delta about `-0.086`.
- Full pipeline invariant tests pass with the new output.

Interpretation:

- The browser now has an explicit review-facing evidence table instead of scanning a model-internal table.
- This is the right pattern for future atlas layers: compute compact review tables during the pipeline, then keep Shiny lightweight.

## 2026-06-03: RDG atlas browser made branch-aware

Motivation: the atlas heatmap showed only clean-CDS buffering effect. This made rows such as PABPC1 tumor/kidney look neutral even when the RDG PDF and joint branch-allocation rows showed a clear branch redistribution. The browser needed to separate "CDS output changed" from "ribosomes changed route through the RDG."

Changes:

- Updated `apps/rdg_atlas_browser/app.R`.
- Added a `Branch Heatmap` tab that plots joint branch-allocation shift magnitude from `joint_top_l1_shift`.
- Added a `Branch Contexts` tab that shows the selected gene/design's row-level joint-allocation contexts, sorted by `joint_review_priority`.
- Added joint/grouped/count-aware branch details to the effect-space hover text.
- Split effect-space warnings onto one hover line per warning.
- Made sparse branch-heatmap aggregation robust to missing gene/design cells.

Result:

- PABPC1 tumor/kidney can now be found directly in the selected-card branch-context table:
  - `PRJNA256316-homo_sapiens`, `CONDITION=tumor`, `TISSUE=kidney`, `tumor` vs `WT`;
  - strict joint branch shift;
  - leader-uORF allocation delta about `+0.086`;
  - clean-CDS allocation delta about `-0.086`.
- The clean-CDS heatmap remains intentionally unchanged and continues to show `atlas_clean_cds_pct_change`.

Interpretation:

- A near-white clean-CDS heatmap cell is not evidence that no RDG branch point exists.
- Future review should use the clean-CDS heatmap for output buffering and the branch heatmap/context table for routing changes.

## 2026-06-02: Joint RDG branch-allocation model added

Motivation: independent beta-posterior branch screens test one feature against "everything else" and can miss whole-RDG redistribution. The next iteration needed a conservative joint branch-composition view that uses existing grouped branch counts and does not reload FST pages.

Changes:

- Added `scripts/dominant_rdg_joint_branch_allocation_model.R`.
- Added pipeline step `joint_branch_allocation` after `motif_prior_calibration` and before `rdg_atlas_cards`.
- Aggregates grouped branch counts into four classes:
  - leader uORF;
  - overlapping uORF;
  - clean CDS;
  - residual other branch mass.
- Uses a uniform Dirichlet posterior for primary case/control branch composition.
- Adds motif-prior sensitivity/alignment columns without hard-coding motifs into the primary posterior call.
- Caps non-evaluable low-count rows at low review priority.
- Merges joint allocation support into atlas cards, atlas warnings, atlas summary metrics, and browser-review priority.
- Added figures:
  - `joint_branch_allocation_top_shifts.*`
  - `joint_branch_allocation_viral_review.*`
  - `joint_branch_allocation_delta_heatmap.*`
- Added invariant tests that joint case/control posterior proportions sum to one and that atlas joint components are finite.

Outputs:

- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_rows.csv`
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_class_rows.csv`
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_gene_design_summary.csv`
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_review.csv`
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_viral_review.csv`
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_summary_metrics.csv`

Result:

- Joint allocation rows: 13,038.
- Evaluable rows: 8,623.
- Strict rows: 5,865.
- Viral-infection rows: 1,148.
- Viral-infection evaluable rows: 814.
- Viral-infection strict rows: 247.
- Row-level review classes:
  - strict joint branch shifts: 712;
  - joint branch shift review: 173;
  - motif-discordant review: 33;
  - high motif prior / weak joint evidence: 279.
- Atlas-supported joint rows: 179.
- Atlas strict joint rows: 137.
- Viral atlas-supported joint rows: 27.
- Viral atlas strict joint rows: 11.
- Motif strength versus L1 allocation shift:
  - all rows Spearman: 0.316;
  - viral-infection Spearman: 0.298.
- Top viral joint-allocation queue: ATF4, PPP1R15A, DDIT3, MYC, IFIH1, SNAI2, NDUFV1, IGFBP5, MRPS12, DDIT4, VEGFA, TBP, SESN2, TOMM20, FEN1, CDK1, TFAM, TRIB3, LMNB1.

Interpretation:

- The joint model is a better review primitive than independent feature screens because it asks whether the whole RDG allocation changes.
- It is still a screening approximation, not a final hierarchical likelihood.
- ATF4, PPP1R15A, and DDIT3 now co-occur as top viral whole-RDG redistribution candidates.
- DDIT3 remains especially interesting because it is browser-supported but motif-discordant, making it a useful model-failure/training example.

## 2026-06-02: Motif priors calibrated against count-aware branch-allocation evidence

Motivation: the constrained motif model showed weak but biologically useful viral signal. The next step was to use motif scores as priors against observed branch-allocation evidence, separating cases where the motif and data agree, disagree, or where the motif is high but branch evidence is weak.

Changes:

- Added `scripts/dominant_rdg_motif_prior_branch_calibration.R`.
- Added pipeline step `motif_prior_calibration` after `mechanistic_motif_model`.
- Joined bounded motif scores to count-aware posterior branch contrasts.
- Created branch-level prior classes:
  - motif/branch concordant;
  - motif/branch discordant;
  - high motif prior with weak branch evidence;
  - low/ambiguous prior.
- Saved gene-design review queues and viral-specific review tables.
- Added figures:
  - `motif_prior_vs_branch_evidence.*`
  - `motif_prior_viral_alignment.*`
- Added invariant tests for motif-prior calibration outputs and bounded priority scores.

Outputs:

- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_calibration_summary.csv`
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_branch_rows.csv`
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_gene_design_calibration.csv`
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_review_queue.csv`
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_viral_review.csv`

Result:

- Evaluable branch rows: 20,312.
- Gene-designs: 522.
- Global Spearman prior strength vs absolute branch delta: 0.117.
- Clean-CDS Spearman prior strength vs absolute branch delta: 0.252.
- Viral-infection Spearman prior strength vs absolute branch delta: 0.093.
- Review queue classes:
  - high prior / weak branch evidence: 64;
  - concordant branch evidence: 29;
  - discordant branch evidence: 18;
  - unresolved: 1.
- Viral actionable rows: TFAM, ATF4, PPP1R15A, LMNB1, DDIT3, STAT1, CDK1, SLC2A3, MMP2, DDIT4, SESN2, SNAI2, MKI67, MX1.
- DDIT3 is the main viral discordant row: browser-supported routing change, but simple motif prior and branch direction disagree.

Interpretation:

- Motif priors are weakly calibrated against branch evidence, especially for clean-CDS rows.
- The viral queue is useful for browser validation because it separates missing-evidence cases from model-failure cases.
- DDIT3 should be treated as a calibration failure mode to learn from, not as evidence against the viral/post-viral biology.

## 2026-06-02: Mechanistic RDG motif scores tested after raw graph geometry failed

Motivation: raw graph-geometry covariates worsened leave-one-gene-out prediction. The next iteration compressed architecture into a small set of constrained mechanistic scores that should be easier to interpret and less prone to overfitting.

Changes:

- Added `scripts/dominant_rdg_mechanistic_motif_model.R`.
- Added pipeline step `mechanistic_motif_model` after `graph_geometry_model`.
- Built bounded `[0, 1]` motif scores:
  - inhibitory-overlap burden;
  - reinitiation-distance opportunity;
  - dense-uORF collision burden;
  - strong-start leakage burden;
  - downstream-rescue opportunity;
  - derived overlap/rescue interaction scores.
- Fitted leave-one-gene-out ridge models without gene identity:
  - metadata/context only;
  - context plus motif main effects;
  - context plus motif-design interactions.
- Added figures:
  - `mechanistic_motif_model_rmse.*`
  - `mechanistic_motif_model_viral_genes.*`
  - `mechanistic_motif_model_coefficients.*`
- Added invariant tests for bounded motif scores and motif-model outputs.

Outputs:

- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_model_summary.csv`
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_model_loo.csv`
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_model_gene_design.csv`
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_model_viral_gene_design.csv`
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_scores.csv`

Result:

- All rows: context/no-gene RMSE 0.428, context + motif RMSE 0.4284, context + motif-design RMSE 0.4295.
- Viral infection: context/no-gene RMSE 0.542, context + motif RMSE 0.539, context + motif-design RMSE 0.545.
- Viral q95 absolute error improved by about 0.030 log2FC with motif main effects.
- Motif-design interactions overfit and were worse than motif main effects.
- Viral genes most improved by motifs: SESN2, IGFBP5, MCM5, ATP5MC3, LMNB1, IFIT1, DDIT3, DDIT4, TRIB3.
- Viral genes most worsened by motifs: TFAM, HLA-DRB1, NDUFS1, VWF, CDK1, MRPL12, ATF3.

Interpretation:

- Constrained motif scores are better than raw graph geometry for viral-infection prediction, but still weak.
- The current hand-weighted motifs should be treated as interpretable priors/review annotations, not as a final rule model.
- The next step is to calibrate motif weights against browser-validated RDGs or use motif scores as priors in a joint branch-allocation model.

## 2026-06-02: RDG graph-geometry generalization tested with leave-one-gene-out blocking

Motivation: the branch-covariate model showed only weak predictive value. The next question was whether explicit RDG architecture can generalize across genes, which is closer to the long-term goal of learning reusable uORF/CDS buffering rules.

Changes:

- Added `scripts/dominant_rdg_graph_geometry_generalization_model.R`.
- Added pipeline step `graph_geometry_model` after `hierarchical_branch_covariates`.
- Built architecture covariates from:
  - `uorf_structure_gene_geometry.csv`;
  - `rdg_feature_annotation.csv`.
- Fitted leave-one-gene-out ridge models without a gene identity covariate:
  - metadata/context only;
  - metadata/context plus graph-geometry covariates.
- Added automatic dropping of single-level factor covariates before model-matrix construction.
- Added figures:
  - `graph_geometry_model_rmse.*`
  - `graph_geometry_model_viral_genes.*`
  - `graph_geometry_model_coefficients.*`
- Added invariant tests for graph-geometry output tables and architecture covariates.

Outputs:

- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_summary.csv`
- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_loo.csv`
- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_gene_design.csv`
- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_viral_gene_design.csv`
- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_covariates.csv`

Result:

- All rows: context/no-gene RMSE 0.428, context + graph-geometry RMSE 0.450.
- Viral infection: context/no-gene RMSE 0.542, context + graph-geometry RMSE 0.572.
- Graph geometry worsened RMSE in every design family in this raw feature representation.
- Viral genes most improved by geometry: SESN2, NDRG1, ATP5MC3, VCAM1, MMP2.
- Viral genes most worsened by geometry: MRPL12, CDK1, ATF3, VWF, TFAM, ATF4, PPP1R15A.

Interpretation:

- The raw architecture feature dump is not a good generalizable rule model.
- This is useful negative evidence. The next geometry model should use constrained mechanistic graph scores or low-dimensional motif components, not more unstructured covariates.
- Key ISR genes getting worse under raw geometry suggests that ATF-like decision graphs require more precise branch-order and reinitiation modeling than simple counts/overlap bins.

## 2026-06-02: Leakage-controlled branch covariates tested as CDS-buffering predictors

Motivation: count-aware branch allocation was already used as evidence in the atlas and hierarchy. The next question was stricter: do branch features improve held-out CDS-buffering prediction when they are actual model covariates, not just priority labels?

Changes:

- Added `scripts/dominant_rdg_hierarchical_branch_covariate_model.R`.
- Added pipeline step `hierarchical_branch_covariates` after `hierarchical_context_model`.
- Branch covariates are built from count-aware posterior branch contrasts while excluding the held-out study.
- Fitted two leave-one-study-out ridge residual models:
  - metadata/context only;
  - metadata/context plus branch covariates.
- Saved all-scope, design-family, gene-design, viral gene-design, coefficient, lambda, and covariate audit outputs.
- Added figures:
  - `branch_covariate_model_rmse.*`
  - `branch_covariate_model_viral_genes.*`
- Added invariant tests for the new pipeline step and output metrics.

Outputs:

- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariate_model_summary.csv`
- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariate_model_loo.csv`
- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariate_model_gene_design.csv`
- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariate_model_viral_gene_design.csv`
- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariates_training_covariates.csv`

Result:

- All rows: hierarchical RMSE 0.480, metadata/context RMSE 0.4776, metadata/context + branch RMSE 0.4775.
- Viral infection: hierarchical RMSE 0.570, metadata/context RMSE 0.5633, metadata/context + branch RMSE 0.5628.
- The all-scope branch-vs-metadata RMSE gain is only 0.000023 log2FC.
- The viral branch-vs-metadata RMSE gain is 0.000492 log2FC.
- Viral genes with the largest small branch-covariate improvements are DDIT3, PPP1R15A, and ATF4.
- VCAM1 is the clearest small viral case where branch covariates worsen residual prediction.

Interpretation:

- Branch covariates carry weak directionally sensible signal in ISR-like viral genes, but the improvement is too small to treat current branch summaries as strong predictive features.
- The next model should use richer RDG graph geometry and a joint branch-allocation likelihood instead of simply adding more branch-summary columns.

## 2026-06-02: Count-aware branch evidence added to hierarchical CDS-buffering summaries

Motivation: count-aware branch allocation was visible in the atlas, but not yet in the hierarchical CDS-buffering model outputs. The next iteration needed to connect stricter branch support to hierarchical review priority without pretending branch allocation directly estimates CDS expression.

Changes:

- `scripts/dominant_rdg_hierarchical_effect_model.R` now reads `dominant_rdg_branch_allocation_gene_design_summary.csv`.
- Count-aware branch-allocation fields are merged into `dominant_rdg_hierarchical_gene_design_effects.csv`.
- Added `hierarchical_count_aware_branch_component`.
- Added `hierarchical_mechanistic_support_class`.
- Count-aware rows only contribute to hierarchical priority when they are actual branch-shift support classes.
- `count_aware_evaluable_no_strong_shift` rows are retained in outputs but contribute zero to the hierarchical count-aware branch component.
- The hierarchical heatmap now marks count-aware-supported cells with a purple diamond.
- `hierarchical_model` now declares `branch_allocation` as a pipeline dependency.
- `tests/integration/test_pipeline_invariants.R` now checks hierarchical count-aware integration and nonblank mechanistic-support labels.

Outputs:

- Updated `dominant_rdg_hierarchical/dominant_rdg_hierarchical_gene_design_effects.csv`.
- Updated `dominant_rdg_hierarchical/dominant_rdg_hierarchical_summary_metrics.csv`.
- Updated `dominant_rdg_hierarchical/figures/hierarchical_gene_design_heatmap.*`.
- Updated downstream residual/context model outputs and atlas cards through the normal pipeline chain.

Result:

- Hierarchical count-aware supported rows: 51.
- Hierarchical replicated count-aware rows: 19.
- Viral hierarchical count-aware supported rows: 3.
- Viral hierarchical replicated count-aware rows: 1.
- Leave-one-study-out MAE remains 0.291 log2FC.
- Leave-one-study-out RMSE remains 0.480 log2FC.
- Calibrated interval coverage remains about 0.95.
- Median hierarchical priority is now 0.413 after removing the small boost from count-aware-evaluable/no-strong-shift rows.

Interpretation:

- This improves ranking semantics without changing the clean-CDS posterior model itself.
- The next real model test is whether count-aware branch features improve held-out CDS-buffering prediction when used as covariates, not just labels.

## 2026-06-02: Count-aware branch-allocation evidence integrated into atlas cards

Motivation: the count-aware branch-allocation screen was useful but separate from the main atlas review queue. The next step was to make those stricter branch-count calls visible in atlas cards and browser-priority outputs.

Changes:

- `scripts/dominant_rdg_atlas_cards.R` now reads `dominant_rdg_branch_allocation_gene_design_summary.csv`.
- Count-aware fields are merged by `gene_symbol + tx_id + design_family`.
- Count-aware feature-id columns are renamed during merge to avoid collisions with older branch-usage feature-id columns.
- Added `count_aware_branch_allocation_component` to atlas confidence and browser-priority scoring.
- Added `count_aware_branch_allocation_supported` warning labels.
- Added count-aware columns near the branch evidence columns in `dominant_rdg_atlas_cards.csv`.
- Added count-aware rows to atlas summary metrics.
- Updated the viral review-priority plot: circle = browser supported; purple diamond = count-aware branch-allocation support.
- `rdg_atlas_cards` now declares `branch_allocation` as a pipeline dependency.
- `tests/integration/test_pipeline_invariants.R` now checks that atlas cards contain finite count-aware components.

Outputs:

- Updated `dominant_rdg_atlas/dominant_rdg_atlas_cards.csv`.
- Updated `dominant_rdg_atlas/dominant_rdg_atlas_browser_validation_queue.csv`.
- Updated `dominant_rdg_atlas/dominant_rdg_atlas_viral_infection_review.csv`.
- Updated `dominant_rdg_atlas/dominant_rdg_atlas_summary_metrics.csv`.
- Updated `dominant_rdg_atlas/figures/dominant_rdg_atlas_viral_review_priority.*`.

Result:

- Count-aware branch-allocation supported atlas rows: 51.
- Replicated count-aware atlas rows: 19.
- Viral count-aware-supported atlas rows: 3.
- Viral replicated count-aware atlas rows: 1.
- The top viral browser queue remains DDIT3, ATF3, PPP1R15A, RRM1, MYC, ATF4, RPLP0, ATP5MC1, IFIT1, and SESN2.
- DDIT3 and ATF3 remain browser-supported but not count-aware replicated branch calls.
- ATF4 viral infection is the clearest replicated count-aware branch-allocation row.
- PPP1R15A viral infection has single count-aware branch support.

Interpretation:

- This improves review honesty. The atlas now shows when a row is high because of browser validation, clean-CDS effect, metadata/hierarchical evidence, or stricter count-aware branch allocation.
- The viral/post-viral story remains plausible but narrow under strict branch-count evidence.

## 2026-06-02: Count-aware branch-allocation screen and invariant tests

Motivation: grouped branch usage fixed the single-sample sparsity problem, but the next model needed a stricter count-aware layer. The goal was to avoid over-reading raw relative-use shifts by adding beta-posterior shrinkage and explicit invariant checks.

Changes:

- Added `scripts/dominant_rdg_count_aware_branch_allocation.R`.
- Added pipeline step `branch_allocation` after `grouped_branch_usage`.
- The new step reads grouped same-condition replicate branch counts and computes independent beta-posterior case/control usage deltas for each feature.
- The step writes posterior contrasts, feature summaries, gene-design summaries, review rows, metrics, and PNG/PDF/HTML figures.
- Fixed `dominant_next_model_candidate_discovery.R` so `candidate_rank` is the current rank alias; the old preliminary rank is retained as `preliminary_candidate_rank`.
- Added `tests/integration/test_pipeline_invariants.R`.
- Updated `tests/integration/README.md`.

Outputs:

- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_posterior_contrasts.csv`
- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_feature_summary.csv`
- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_gene_design_summary.csv`
- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_review.csv`
- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_summary_metrics.csv`
- `dominant_rdg_branch_allocation/figures/dominant_rdg_branch_allocation_top_shifts.html`
- `dominant_rdg_branch_allocation/figures/dominant_rdg_branch_allocation_viral_review.html`

Result:

- The count-aware layer produced 46,110 posterior contrast rows, 20,312 evaluable rows, 10,851 strict rows, 2,320 feature summaries, 656 gene-design summaries, and 299 review rows.
- It called 28 replicated count-aware branch-shift rows overall.
- Viral infection has 60 count-aware review rows, but only one replicated count-aware viral branch-shift row: ATF4 clean-CDS up.
- Top all-state count-aware rows include ATF4 nutrient-starvation clean-CDS up, PPP1R15A nutrient-starvation clean-CDS up, DDIT4/LMNB1/SESN2/LDHA ISR/stress clean-CDS shifts, and ATF4 viral-infection clean-CDS up.
- The new invariant test passed after rerunning affected steps.

Interpretation:

- This is a useful conservative filter. It keeps the known ATF4 viral/ISR biology alive but narrows the strongest viral branch-allocation claim.
- DDIT3 and ATF3 remain browser-supported viral candidates, but they are not yet replicated count-aware branch-allocation calls.
- The next statistical step should be a true beta-binomial or Dirichlet-multinomial model and integration of count-aware branch evidence into the atlas/hierarchical summaries.

## 2026-06-02: Project review, atlas priority bug fix, and note rewrite

Motivation: the project had accumulated many iteration notes and needed a full sanity check against the original goals: dominant-state QC, clean-CDS scoring, RDG atlas construction, uORF-structure rules, post-viral/fatigue testing, and predictive CDS-buffering estimates.

Checks:

- Parsed all 48 R scripts under `dominant_cell_states/scripts`.
- Checked the pipeline registry: 35 steps, no missing script files, no duplicated step names.
- Audited grouped branch-usage denominators:
  - summed case feature counts equal case branch totals for every contrast group;
  - summed control feature counts equal control branch totals for every contrast group;
  - no duplicated features were found within contrast groups.
- Reviewed current atlas, branch coverage, grouped branch usage, hierarchical model, next-model, and post-viral/fatigue summary metrics.

Bug fix:

- Fixed `scripts/dominant_rdg_atlas_cards.R`.
- Missing `context_calibration_status` or `context_model_status` values previously propagated to `NA` priority-score components.
- `setorder()` then placed NA atlas browser-priority scores at the top of the review queue.
- Missing context components are now scored as zero and non-finite final priority scores are guarded.
- Reran the `rdg_atlas_cards` step.
- The browser queue now has 0 NA priority scores out of 870 review rows.
- The top review rows now start with the expected viral candidates: DDIT3, ATF3, PPP1R15A, RRM1, MYC, RPLP0, ATF4, and IFIT1.

Documentation:

- Rewrote `md/dominant_cell_state_scientific_note.md` as a coherent project review and scientific story instead of an iteration log.
- Replaced the stale handoff summary with a current 2026-06-02 version that includes branch coverage QC, grouped branch usage, hierarchical context models, atlas cards, and current caveats.

Interpretation:

- The project is on track as an RDG-atlas and candidate-generation framework.
- Current biological claims should remain review-grade until browser validation and count-aware branch models are stronger.
- DDIT3 and ATF3 viral-infection rows remain the strongest browser-supported review leads; ATF4 remains the core mechanistic positive control.

## 2026-06-02: Merged-replicate branch-usage adapter and atlas integration

Motivation: branch coverage QC showed that single samples are too sparse for general small-uORF branch inference. The next step was to build the branch-usage model around merged same-condition replicate groups instead of individual runs.

Changes:

- Added `scripts/dominant_rdg_grouped_branch_usage_model.R`.
- Added pipeline step `grouped_branch_usage` after `branch_coverage_qc`.
- The step builds a grouped adapter from `dominant_uorf_feature_expression.csv`.
- Runs are merged by `study + AUTHOR + CELL_LINE + TISSUE + GENE + CONDITION + INHIBITOR + FRACTION + Cancer_type + Sex + TIMEPOINT`.
- Only `clean_CDS`, `leader_uORF`, and `overlapping_uORF` enter the branch-allocation denominator.
- Exact case/control matching uses study, author, cell line, tissue, metadata gene, inhibitor, fraction, cancer type, sex, and time point.
- A relaxed study/cell/tissue/gene fallback is allowed only when no exact case-feature match exists.
- The initial per-term loop was too slow, so the implementation was replaced with a vectorized aggregate-and-join model using a compact `feature_key`.
- Atlas cards now ingest grouped branch-usage evidence as a separate `grouped_branch_usage_*` component.

Outputs:

- `dominant_rdg_grouped_branch_usage/dominant_rdg_grouped_branch_usage_adapter.csv`
- `dominant_rdg_grouped_branch_usage/dominant_rdg_grouped_branch_usage_contrast_table.csv`
- `dominant_rdg_grouped_branch_usage/dominant_rdg_grouped_branch_usage_feature_summary.csv`
- `dominant_rdg_grouped_branch_usage/dominant_rdg_grouped_branch_usage_gene_design_summary.csv`
- `dominant_rdg_grouped_branch_usage/dominant_rdg_grouped_branch_usage_review.csv`
- `dominant_rdg_grouped_branch_usage/dominant_rdg_grouped_branch_usage_summary_metrics.csv`
- `dominant_rdg_grouped_branch_usage/figures/dominant_rdg_grouped_branch_usage_top_shifts.html`
- `dominant_rdg_grouped_branch_usage/figures/dominant_rdg_grouped_branch_usage_viral_review.html`

Result:

- Built 546,650 grouped adapter rows from 1,885 same-condition groups, including 961 groups with at least two runs.
- Produced 46,110 grouped branch contrasts: 31,900 exact-context rows and 14,210 relaxed fallback rows.
- Summarized 2,320 gene x design-family x feature rows and 656 gene x design-family rows.
- Wrote 299 grouped branch review rows.
- Feature-level support classes include 32 `replicated_exact_group_branch_shift` rows, 40 `single_exact_group_branch_shift` rows, and 5 relaxed fallback review rows.
- Viral infection has 59 grouped review rows, but only four feature-level rows pass beyond `group_evaluable_no_strong_shift`: ATF4 clean-CDS up, ATF4 overlapping-uORF down, PPP1R15A clean-CDS up in one study, and a relaxed VEGFA uORF review row.
- Atlas cards now report 45 grouped branch-usage supported rows, including 19 replicated exact grouped rows.

Interpretation:

- This is a stricter evidence layer than the earlier matched branch-usage pass.
- ATF4 is now the clearest viral branch-redistribution example under merged-replicate count gates.
- DDIT3 and ATF3 remain browser-supported viral review leads, but they are not formal grouped branch-shift calls in this model. That distinction is useful: browser-visible routing, clean-CDS buffering, and count-based branch redistribution are related but not interchangeable.
- The adapter CSV is intentionally machine-facing and large; human review should use the feature summary, gene-design summary, review table, atlas cards, and plots.

## 2026-06-02: Branch coverage QC before sample-level branch modeling

Motivation: before building a true `Run x gene x RDG feature` branch-usage model, we needed to test whether single samples have enough coverage on small uORFs. The working assumption is now that single runs are usually too sparse, so the minimum practical unit should be merged same-condition replicates.

Changes:

- Added `scripts/dominant_rdg_branch_coverage_qc.R`.
- Added pipeline step `branch_coverage_qc` after `branch_usage`.
- The QC compares `single_run` units against merged same-condition replicate groups using `study + AUTHOR + CELL_LINE + TISSUE + GENE + CONDITION + INHIBITOR + FRACTION + Cancer_type + Sex + TIMEPOINT`.
- The allocation denominator is restricted to `clean_CDS`, `leader_uORF`, and `overlapping_uORF`.
- `iORF` and `NTE/NTT` features are excluded from clean-CDS relative allocation because they are diagnostic or alternative branch nodes, not ordinary clean-CDS competitors.

Outputs:

- `dominant_rdg_branch_coverage_qc/branch_coverage_gate_summary.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_feature_summary.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_gene_summary.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_unit_summary.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_replicate_group_summary.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_small_uorf_single_vs_merged.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_recommendation.csv`
- `dominant_rdg_branch_coverage_qc/branch_coverage_summary_metrics.csv`
- `dominant_rdg_branch_coverage_qc/figures/branch_coverage_count_gate_by_length.html`
- `dominant_rdg_branch_coverage_qc/figures/branch_coverage_small_uorf_top_features.html`

Result:

- Input feature table contains 1,241,954 rows across 3,857 runs, 82 genes, and 322 features.
- There are 1,885 same-condition merged groups, including 961 groups with at least two runs.
- For small uORFs, single-run observations reach at least 10 raw reads only 25.7% of the time and at least 30 raw reads only 14.3% of the time.
- For merged `n >= 2` same-condition groups, small-uORF observations reach at least 10 raw reads 38.6% of the time and at least 30 raw reads 24.6% of the time.
- The minimal branch-model gate improves from 37.3% of single-run allocation-feature rows to 50.5% of merged `n >= 2` rows.

Interpretation:

- Single runs should not be used as the primary unit for branch modeling.
- The next DOTSeq-like or Dirichlet-multinomial layer should use merged same-condition replicate groups as the minimum unit, keep raw-count gates, and only model covered branches.
- This does not invalidate strong high-count genes such as ATF4, DDIT3, LDHA, EEF1A1, HSPA5, or RPLP0. It means the general atlas model needs replicate aggregation before interpreting small-uORF routing.

## 2026-06-01: Study-blocked context-aware residual model

Motivation: peer-gene context calibration showed that context information reduces some large residuals, but that test used genes from the held-out context and is therefore diagnostic rather than externally predictive. The next iteration asked whether a model trained only on other studies can learn residual corrections from gene, design-family, cell-line, tissue, condition, and gene/context terms.

Changes:

- Added pipeline step `hierarchical_context_model` after `hierarchical_context_calibration` and before `rdg_atlas_cards`.
- Wired `scripts/dominant_rdg_hierarchical_context_model.R` into the single pipeline entry point.
- The model fits a sparse ridge residual learner with study-blocked prediction.
- It now evaluates both `lambda.min` and `lambda.1se` under the same leave-one-study-out loop and selects the lower held-out RMSE.
- Atlas cards now ingest context-model status, RMSE/q95/MAE deltas, worst context, and improved contexts.
- The Shiny atlas browser now has a `Context model improved only` filter and displays context-model fields.

Outputs:

- `dominant_rdg_hierarchical_context_model/dominant_rdg_hierarchical_context_model_loo.csv`
- `dominant_rdg_hierarchical_context_model/dominant_rdg_hierarchical_context_model_gene_design.csv`
- `dominant_rdg_hierarchical_context_model/dominant_rdg_hierarchical_context_model_summary.csv`
- `dominant_rdg_hierarchical_context_model/dominant_rdg_hierarchical_context_model_viral_gene_design.csv`
- `dominant_rdg_hierarchical_context_model/dominant_rdg_hierarchical_context_model_coefficients.csv`
- `dominant_rdg_hierarchical_context_model/dominant_rdg_hierarchical_context_model_lambda_selection.csv`
- `dominant_rdg_hierarchical_context_model/figures/hierarchical_context_model_rmse.html`
- `dominant_rdg_hierarchical_context_model/figures/hierarchical_context_model_viral_genes.html`

Result:

- The selected penalty is `lambda.min` (`lambda = 3.017`), not the conservative `lambda.1se` (`lambda = 83.692`).
- Across all 10,224 leave-one-study-out predictions, RMSE improves from 0.480 to 0.478 log2FC and MAE improves from 0.291 to 0.288 log2FC. The q95 residual slightly worsens from 0.966 to 0.971 log2FC.
- In viral infection, RMSE improves from 0.570 to 0.563 log2FC and q95 residual improves from 1.191 to 1.159 log2FC, but MAE worsens from 0.370 to 0.374 log2FC and direction accuracy worsens from 0.525 to 0.511.
- 15 gene x design-family rows pass the conservative context-model-improved threshold; 66 are context-model-worse.
- Viral-infection rows passing the threshold are `IFIT1` and `DDIT3`.

Interpretation:

- This is the first context-aware residual correction that does not borrow peer genes from the held-out study.
- The model improves large-error tails more than average viral predictions, which is useful for triage but not yet a reliable forecasting layer.
- The strongest biological interpretation remains local: viral `DDIT3` and `IFIT1` are contexts where gene/cell/tissue/design terms help predict residual structure. The broader lesson is that viral infection needs explicit context and branch-level terms rather than a single broad design-family average.

## 2026-06-01: Peer-gene context calibration test

Motivation: the residual audit showed that broad `viral_infection` errors cluster in specific study/cell/tissue contexts. The next question was whether a simple context-level residual offset can reduce held-out-study error. This is not a fully causal or external forecasting model: it uses peer genes in the same held-out context to estimate a context residual offset, so it is best interpreted as a within-compendium context-calibration test.

Changes:

- Added `scripts/dominant_rdg_hierarchical_context_calibration.R`.
- Added pipeline step `hierarchical_context_calibration` after `hierarchical_residual_audit` and before `rdg_atlas_cards`.
- The model computes peer-gene residual offsets while excluding the target gene.
- Offset candidates are tried in order:
  - exact design-family x study x cell-line x tissue x condition context;
  - design-family x study;
  - design-family x cell-line x tissue x condition;
  - design-family x cell-line x tissue;
  - design-family x condition;
  - design-family.
- Offsets are shrinkage-weighted by peer count and applied to held-out predictions.
- The script reports original versus context-corrected MAE, RMSE, q90/q95 residuals, and direction accuracy.
- Atlas cards now ingest context-calibration status, RMSE delta, q95 delta, peer-corrected fraction, and worst context after correction.
- The Shiny atlas browser now has a `Context calibration improved only` filter.

Outputs:

- `dominant_rdg_hierarchical_context_calibration/dominant_rdg_hierarchical_context_calibration_loo.csv`
- `dominant_rdg_hierarchical_context_calibration/dominant_rdg_hierarchical_context_calibration_gene_design.csv`
- `dominant_rdg_hierarchical_context_calibration/dominant_rdg_hierarchical_context_calibration_summary.csv`
- `dominant_rdg_hierarchical_context_calibration/dominant_rdg_hierarchical_context_calibration_offset_summary.csv`
- `dominant_rdg_hierarchical_context_calibration/dominant_rdg_hierarchical_context_calibration_viral_gene_design.csv`
- `dominant_rdg_hierarchical_context_calibration/figures/hierarchical_context_calibration_rmse.html`
- `dominant_rdg_hierarchical_context_calibration/figures/hierarchical_context_calibration_viral_genes.html`

Result:

- Across all 10,224 leave-one-study-out predictions, context correction slightly improves RMSE from 0.480 to 0.474 log2FC and q95 residual from 0.966 to 0.949 log2FC, but slightly worsens MAE from 0.291 to 0.294 log2FC.
- In viral infection, context correction improves RMSE from 0.570 to 0.562 log2FC and q95 residual from 1.191 to 1.167 log2FC, but worsens MAE from 0.370 to 0.381 log2FC and direction accuracy from 0.525 to 0.493.
- 102 gene x design-family rows are classified as context-improved by RMSE, and 118 are context-worse.
- In viral infection, 7 rows are context-improved, 17 are context-worse, and 47 are neutral.
- Viral rows improved by context calibration include `PPP1R15A`, `SESN2`, `ATF3`, `HERPUD1`, `STAT1`, `TRIB3`, and `UBC`.
- Viral rows worsened by context calibration include low-residual or context-sensitive rows such as `MS4A1`, `VWF`, `C1QB`, `C1QC`, `CDH1`, `MRPL12`, `ATP5MC3`, `PCNA`, `LMNB1`, `HLA-DRB1`, `NDUFA1`, and `MKI67`.
- Exact same-context peer offsets are usually available; median peer count is 70 and median absolute exact-context offset is about 0.046 log2FC. The offsets are small relative to the largest gene-specific residuals.

Interpretation:

- Simple additive context offsets help a little, especially for viral RMSE/q95 residuals, but they do not solve the viral prediction problem.
- The fact that MAE and direction accuracy can worsen means the residual structure is not just a shared study/cell/tissue shift. Gene x context interactions and branch-specific mechanisms are likely important.
- ATF3 is a useful example: it is high-residual but context calibration reduces its viral RMSE. DDIT3 remains high-residual but only context-neutral under this simple correction.
- The next model should include explicit context covariates and gene-context interactions in the hierarchical model itself, rather than applying a post hoc average offset.

## 2026-06-01: Hierarchical residual heterogeneity audit

Motivation: residual-calibrated hierarchical intervals fixed leave-one-study-out coverage, but they became very wide. The next question was not how to widen intervals further, but where the model fails. Broad design-family labels such as `viral_infection` likely mix distinct study, cell-line, tissue, virus, and time-point effects.

Changes:

- Added `scripts/dominant_rdg_hierarchical_residual_audit.R`.
- Added pipeline step `hierarchical_residual_audit` between `hierarchical_model` and `rdg_atlas_cards`.
- The audit attaches each leave-one-study-out residual back to matched-design context from the clean-CDS buffering table:
  - held-out study;
  - cell line;
  - tissue;
  - case condition;
  - control condition;
  - matched design signature;
  - strict/support counts.
- It summarizes residual error by `gene x transcript x design_family`.
- It summarizes residual error by design-family context dimensions: study, cell line, tissue, and condition token.
- It creates viral-infection focused residual context and gene-design tables.
- Atlas cards now ingest residual risk class, residual heterogeneity score, worst residual context, top residual studies/cell lines/tissues/conditions, and residual error metrics.
- The Shiny atlas browser now has a `High hierarchical residual only` filter and displays residual-risk fields on cards.

Outputs:

- `dominant_rdg_hierarchical_residual_audit/dominant_rdg_hierarchical_residual_audit_loo_context.csv`
- `dominant_rdg_hierarchical_residual_audit/dominant_rdg_hierarchical_residual_audit_gene_design.csv`
- `dominant_rdg_hierarchical_residual_audit/dominant_rdg_hierarchical_residual_audit_context_summary.csv`
- `dominant_rdg_hierarchical_residual_audit/dominant_rdg_hierarchical_residual_audit_viral_contexts.csv`
- `dominant_rdg_hierarchical_residual_audit/dominant_rdg_hierarchical_residual_audit_viral_gene_design.csv`
- `dominant_rdg_hierarchical_residual_audit/dominant_rdg_hierarchical_residual_audit_summary_metrics.csv`
- `dominant_rdg_hierarchical_residual_audit/figures/hierarchical_residual_gene_design_scatter.html`
- `dominant_rdg_hierarchical_residual_audit/figures/hierarchical_residual_viral_contexts.html`

Result:

- Audited 10,224 leave-one-study-out residual rows.
- Produced residual summaries for 781 gene x design-family rows.
- Flagged 191 high-residual-heterogeneity rows and 241 moderate-residual-heterogeneity rows.
- Viral infection has 71 residual-audited gene rows, with 19 high-residual-heterogeneity rows.
- Viral median residual MAE is 0.306 log2FC and median residual RMSE is 0.404 log2FC.
- Highest viral residual gene rows include `CCL2`, `IFIT1`, `DDIT3`, `ATF3`, `PPP1R15A`, `TRIB3`, `ATF4`, `SESN2`, `OAS2`, `IRF7`, `HERPUD1`, `DDIT4`, `MX1`, `SLC2A3`, `IFIH1`, `CCND1`, `STAT1`, `VEGFA`, and `CFB`.
- Strongest viral residual contexts include:
  - NPC/brain, `PRJNA854905`;
  - BT20/breast, `PRJNA957283`;
  - HFF/skin, `PRJNA533074`;
  - A549/lung, including `PRJNA1063244`, `PRJNA395384`, and `PRJNA324409`;
  - MRC5/lung, `PRJNA1062064`;
  - Huh7/liver, `PRJNA285961`.
- ATF3 and DDIT3 remain top browser-supported viral leads, but both are now explicitly marked as high residual heterogeneity rows. Their worst held-out contexts differ: ATF3 is worst in Huh7/liver `PRJNA285961`, while DDIT3 is worst in BT20/breast `PRJNA957283`.

Interpretation:

- This confirms that broad `viral_infection` is not a single transferable effect.
- The next predictive model needs explicit viral-context covariates, not only wider intervals.
- The audit gives practical browser-review and modeling targets: inspect top residual contexts and decide whether they represent virus/time/cell biology, metadata mixing, or study/protocol effects.
- High residual risk now increases browser-review priority because these are the rows most informative for improving the atlas, even when they are not calibrated effect claims.

## 2026-06-01: Residual-calibrated hierarchical intervals

Motivation: the first hierarchical clean-CDS effect model was useful for shrinkage and ranking, but its leave-one-study-out interval coverage was only about 43%. That made the `hierarchical_atlas_ready` and interval-supported labels too optimistic.

Changes:

- Added conformal-style residual calibration to `scripts/dominant_rdg_hierarchical_effect_model.R`.
- The calibration uses leave-one-study-out absolute prediction errors.
- Design-family 90% and 95% residual quantiles are used when a design family has at least 50 calibration predictions.
- Design families without enough leave-one-study-out predictions fall back to the global residual quantile.
- The model now writes calibrated posterior interval fields, calibrated interval-support flags, calibration scope, calibration residual quantiles, and uncalibrated evidence labels for audit.
- Added `dominant_rdg_hierarchical/dominant_rdg_hierarchical_interval_calibration.csv`.
- Atlas cards now ingest calibrated hierarchical intervals and mark rows where an uncalibrated hierarchical interval is lost after calibration.
- The Shiny atlas browser now exposes calibrated hierarchical interval fields and keeps the old labels only as backward-compatible filters.

Outputs:

- `dominant_rdg_hierarchical/dominant_rdg_hierarchical_interval_calibration.csv`
- updated `dominant_rdg_hierarchical/dominant_rdg_hierarchical_gene_design_effects.csv`
- updated `dominant_rdg_hierarchical/dominant_rdg_hierarchical_calibration_summary.csv`
- updated `dominant_rdg_atlas/dominant_rdg_atlas_cards.csv`
- updated `dominant_rdg_atlas/dominant_rdg_atlas_summary_metrics.csv`

Result:

- Leave-one-study-out point error is unchanged: MAE 0.291 log2FC and RMSE 0.480 log2FC.
- Original posterior interval coverage remains 0.432.
- Calibrated interval coverage is 0.950.
- The global 95% residual half-width is 0.966 log2FC.
- Viral infection is especially heterogeneous: viral-infection residual MAE is 0.370 log2FC and the 95% residual half-width is 1.191 log2FC.
- Calibrated hierarchical interval support collapses from 251 uncalibrated interval rows to 1 calibrated interval row.
- Calibrated hierarchical atlas-ready support collapses from 105 uncalibrated rows to 1 calibrated atlas-ready row: `LMNB1` in `acute_interferon_cytokine`, with posterior clean-CDS change about -38% and calibrated interval about -61% to -2%.
- No viral-infection row has a calibrated hierarchical interval away from zero. ATF3 and DDIT3 remain top browser-supported review leads, but the calibrated hierarchy treats their viral effects as replicated/no-calibrated-interval evidence.

Interpretation:

- This is a useful correction. The hierarchy is still a ranking and shrinkage layer, but it is not yet a precise effect-estimation model.
- The calibrated intervals are deliberately wide because held-out-study error is large.
- Viral infection should be framed as high-priority, heterogeneous, browser-supported review biology rather than calibrated quantitative prediction.
- The next model should reduce residual error by adding study, cell type, tissue, time point, virus, treatment dose, and branch-feature covariates instead of relying on a single broad design-family label.

## 2026-06-01: First hierarchical clean-CDS effect model

Motivation: the atlas cards were still mostly a heuristic integration layer. They ranked matched clean-CDS buffering, branch usage, sparse/on-off evidence, browser notes, and post-viral relevance, but they did not yet provide a condition-family posterior or an explicit calibration test. The next step toward a predictive RDG atlas is to learn gene x design-family clean-CDS effects with empirical shrinkage and ask how well those effects predict held-out studies.

Changes:

- Added `scripts/dominant_rdg_hierarchical_effect_model.R`.
- Added pipeline step `hierarchical_model` after `fatigue_summary` and before `rdg_atlas_cards`.
- The model reads matched clean-CDS buffering rows, estimates row uncertainty from bootstrap/error fields where possible, and summarizes effects by `gene_symbol x tx_id x design_family`.
- It builds empirical gene-design, gene, design-family, and global priors, then shrinks each direct clean-CDS effect toward the best available prior.
- It runs leave-one-study-out prediction for gene-design rows with more than one study.
- It writes static and interactive plots for the hierarchical gene-design heatmap and leave-one-study-out calibration.
- The RDG atlas cards and Shiny browser now ingest hierarchical evidence class, posterior percent effect, priority score, prior source, and calibration fields.

Outputs:

- `dominant_rdg_hierarchical/dominant_rdg_hierarchical_gene_design_effects.csv`
- `dominant_rdg_hierarchical/dominant_rdg_hierarchical_counterfactual_predictions.csv`
- `dominant_rdg_hierarchical/dominant_rdg_hierarchical_leave_one_study_calibration.csv`
- `dominant_rdg_hierarchical/dominant_rdg_hierarchical_calibration_summary.csv`
- `dominant_rdg_hierarchical/dominant_rdg_hierarchical_summary_metrics.csv`
- `dominant_rdg_hierarchical/figures/hierarchical_gene_design_heatmap.html`
- `dominant_rdg_hierarchical/figures/hierarchical_leave_one_study_calibration.html`

Result:

- Modeled 14,058 matched clean-CDS rows.
- Produced 994 gene x design-family posterior rows across 71 genes and 14 design families.
- Classified 105 rows as `hierarchical_atlas_ready`, 68 as `hierarchical_replicated_interval_review`, and 251 as interval-supported.
- Leave-one-study-out calibration produced 10,224 predictions, with MAE 0.291 log2FC and RMSE 0.480 log2FC.
- Interval coverage was only 0.432, so the current posterior intervals are under-dispersed and should not be treated as final calibrated uncertainty.
- Top viral-infection posterior rows include `ATF3`, `DDIT3`, `MYC`, `PPP1R15A`, `RPLP0`, `STAT1`, `HERPUD1`, `VEGFA`, `IFIT1`, `TNFAIP3`, `OAS2`, and `IFIH1`.

Interpretation:

- This is the first atlas-level empirical-Bayes prediction layer.
- It is useful for ranking and shrinkage because it demotes single noisy contrasts and makes gene/design evidence explicit.
- It is not yet a reliable forecasting model because held-out-study interval coverage is too low and the counterfactual grid currently has no missing gene x design-family combinations.
- The next model should move from design-family labels to sample-level covariates and branch-feature counts, so it can predict unobserved combinations such as stress in a new cell line instead of only summarizing observed gene x family rows.

## 2026-05-28: Pause frame/off-frame branch pending infected p-shift fixes

Motivation: user QC indicates that much of the infected-group off-frame signal is likely pure p-shift error rather than biology. The B2M read-length check already showed a one-nucleotide mismatch between the current BAM+shift-table reconstruction and the older aggregate FST/browser track for `SRR27450688`. That means the frame-sensitive branch is useful for QC, but it should not drive the normal biological pipeline until the infected p-shifts and FST pages are regenerated.

Changes:

- Removed frame/off-frame steps from the default `rdg`, `manual`, and `all` pipeline scopes.
- Added a separate `frame_qc` scope containing:
  - `cds_frame_bumpiness_qc`
  - `infected_off_frame_scan`
  - `infected_off_frame_frame_risk`
  - `infected_off_frame_hotspots`
  - `infected_off_frame_review`
  - `infected_off_frame_frame_diagnostics`
  - `infected_off_frame_readlength_diagnostics`
- Existing frame/off-frame outputs are retained on disk, but should be treated as a p-shift QC sandbox, not as candidate biology.

Interpretation:

- `run_dominant_cell_state_pipeline(scope = "all")` is now the normal biological pipeline again.
- `run_dominant_cell_state_pipeline(scope = "frame_qc")` is available for later, after shifts/FST pages are updated.
- Until then, prioritize the RDG atlas, clean-CDS buffering, uORF/clean-CDS allocation, matched-control branch usage, post-viral state models, and browser-validated DDIT3/ATF3-style candidates.

## 2026-05-28: B2M read-length BAM decomposition

Motivation: the aggregate FST frame-risk diagnostic kept `B2M` codon 55 as a browser-supported lead, but exact frame assignment needed raw read-length evidence. The user provided the two supporting collapsed BAMs and the two study-level shifting-table RDS files. Because these BAMs are collapsed, counts must be read through `ORFik::readBam()`, which detects collapsed qnames and stores the collapsed multiplicity in `mcols(reads)$score`.

Changes:

- Added `scripts/dominant_b2m_readlength_frame_diagnostics.R`.
- Added pipeline step `infected_off_frame_readlength_diagnostics`.
- The step reads `SRR2052924.bam` and `SRR27450688.bam` from `readlength_inputs/`.
- It extracts the run-specific shift table by matching `basename(names(shift_list))` to `<Run>.ofst`.
- It applies `ORFik::shiftFootprints()` and then profiles B2M by read length, transcript position, and CDS frame.
- It writes a BAM-default-p-shift versus aggregate-FST comparison table.

Outputs:

- `dominant_infected_off_frame_readlength_diagnostics/b2m_readlength_shift_tables_used.csv`
- `dominant_infected_off_frame_readlength_diagnostics/b2m_readlength_run_summary.csv`
- `dominant_infected_off_frame_readlength_diagnostics/b2m_readlength_psite_profile.csv`
- `dominant_infected_off_frame_readlength_diagnostics/b2m_readlength_window_frame_summary.csv`
- `dominant_infected_off_frame_readlength_diagnostics/b2m_codon55_readlength_summary.csv`
- `dominant_infected_off_frame_readlength_diagnostics/b2m_codon55_psite_offset_sensitivity.csv`
- `dominant_infected_off_frame_readlength_diagnostics/b2m_codon55_fst_vs_bam_pshift_counts.csv`
- `dominant_infected_off_frame_readlength_diagnostics/figures/`

Result:

- `SRR2052924`: codon55 has 506 total shifted BAM counts; frame1 at transcript position 194 is dominant, with 326 counts and a frame1/frame0 ratio of 25.1. The strongest read length is 32 nt at position 194. This exactly matches the aggregate FST codon55 frame1 count at position 194.
- `SRR27450688`: codon55 has 1,739 total shifted BAM counts, but the supplied shift table places the dominant signal at transcript position 195/frame2, with 1,574 counts. The aggregate FST/browser frame-risk table has the same 1,574-count peak at transcript position 194/frame1. The offset-sensitivity table shows that shifting the BAM-derived p-site coordinates one nucleotide upstream moves the peak back to 194/frame1.
- The peak itself is robust in both BAMs. The exact frame assignment in `SRR27450688` is p-site-convention-sensitive by one nucleotide.

Interpretation:

- `B2M` codon55 remains a real local peak, not a low-count artifact.
- For `SRR2052924`, the read-length decomposition supports the aggregate FST/browser frame1 call directly.
- For `SRR27450688`, the biology-vs-artifact question now splits into two parts: the local peak is real, but the exact frame label depends on a one-nt p-site convention difference between the supplied BAM+shift reconstruction and the aggregate FST/browser track.
- The next iteration should inspect how the aggregate FST pages were generated for `PRJNA1062064`, especially whether the effective shift for `SRR27450688` was -13 rather than the supplied -14, or whether there is an export/indexing convention difference.

## 2026-05-28: Hotspot start/end frame diagnostics

Motivation: browser inspection of `B2M` showed a real-looking codon 55 off-frame peak, but the run-level frame behavior was not a clean p-shift. The TIS/start region could peak off the expected coordinate while the CDS end stayed in frame. We needed a script-level diagnostic that keeps strict and frame-risk rows separate while asking whether a hotspot behaves like a simple global frame shift or a localized/mixed event.

Changes:

- Added `scripts/dominant_infected_off_frame_frame_diagnostics.R`.
- Added pipeline step `infected_off_frame_frame_diagnostics` after strict hotspot review.
- The step combines strict no-frame-risk hotspots and frame-risk hotspots.
- For each hotspot/run, it measures the hotspot 31-nt window, the CDS-start window, and the CDS-end window.
- It classifies each hotspot as `localized_or_mixed_frame_behavior`, `simple_global_frame_shift_like`, or `ambiguous_needs_read_length_decomposition`.
- Added a B2M-specific diagnostic window profile for codon 55.

Outputs:

- `dominant_infected_off_frame_frame_diagnostics/infected_off_frame_hotspot_frame_diagnostics_runs.csv`
- `dominant_infected_off_frame_frame_diagnostics/infected_off_frame_hotspot_frame_diagnostics_hotspots.csv`
- `dominant_infected_off_frame_frame_diagnostics/infected_off_frame_hotspot_frame_diagnostics_class_summary.csv`
- `dominant_infected_off_frame_frame_diagnostics/infected_off_frame_hotspot_frame_diagnostics_summary_metrics.csv`
- `dominant_infected_off_frame_frame_diagnostics/infected_off_frame_b2m_codon55_window_profiles.csv`
- `dominant_infected_off_frame_frame_diagnostics/figures/infected_off_frame_hotspot_frame_diagnostic_quadrants.html`
- `dominant_infected_off_frame_frame_diagnostics/figures/infected_off_frame_b2m_codon55_diagnostic_windows.html`

Result:

- Diagnosed 289 hotspot rows and 341 run-hotspot rows.
- Classified 160 hotspots as localized/mixed, 44 as simple-global-frame-shift-like, and 85 as ambiguous.
- Strict recurrent localized/mixed leads include `NDUFA1`, `RPL27`, `UQCRFS1`, and `IFITM3`.
- `B2M` codon 55 frame1 remains the key frame-risk stress test: it has 2 supporting runs, strong local peak concentration, 9.7-25.1x more signal than the same codon's frame0 base, mixed CDS-start top frames (`frame0;frame2`), and frame0 CDS-end top peaks.

Interpretation:

- This does not replace true read-length-specific p-site decomposition.
- It does make the B2M result more interesting: the pattern is not a simple, consistent whole-CDS frame shift.
- The next hard test is to decompose the same positions by read length and p-site model, then ask whether the localized peaks survive denoising.

## 2026-05-27: All-CDS frame metadata and frame-risk hotspot review

Motivation: the metadata table already contains global all-CDS frame usage (`Frame_usage_0`, `Frame_usage_1`, `Frame_usage_2`). The previous off-frame scan used a marker-panel global frame estimate and a housekeeping frame estimate, but did not use this all-CDS metadata directly. That omission let some possible library/protocol frame-bias cases enter the strict off-frame queue. At the same time, browser inspection of `B2M` showed that some technically flagged cases may still be biologically interesting, so they should be separated rather than silently discarded.

Changes:

- Added metadata-derived frame QC columns: `metadata_frame0_fraction`, `metadata_frame1_fraction`, `metadata_frame2_fraction`, `metadata_canonical_frame_fraction`, `metadata_dominant_frame`, and `metadata_possible_misshift`.
- Updated technical-frame gating so `technical_frame_bias_possible` is true if metadata all-CDS frame QC, marker-panel global frame QC, or housekeeping frame QC flags the run.
- Split off-frame peak evidence into strict no-frame-bias candidates and a separate frame-risk review layer.
- Added `scripts/dominant_infected_off_frame_frame_risk_hotspots.R`.
- Added pipeline step `infected_off_frame_frame_risk` after `infected_off_frame_scan`.

Outputs:

- `dominant_infected_off_frame_cds_scan/infected_off_frame_cds_frame_risk_runs.csv`
- `dominant_infected_off_frame_cds_scan/infected_off_frame_cds_frame_risk_gene_scores.csv`
- `dominant_infected_off_frame_frame_risk/infected_off_frame_frame_risk_hotspots.csv`
- `dominant_infected_off_frame_frame_risk/infected_off_frame_frame_risk_gene_summary.csv`
- `dominant_infected_off_frame_frame_risk/figures/infected_off_frame_frame_risk_hotspots.html`

Result:

- The strict no-frame-bias scan narrowed to 6 positive genes: `IFITM3`, `PPP1R15A`, `UQCRFS1`, `KRT18`, `NDUFA1`, and `RPL27`.
- The strict hotspot queue now has 71 hotspot rows across 6 genes, 0 tier-1 matched rows, and 16 tier-2 recurrent unmatched rows.
- `B2M` moved out of the strict queue because its relevant runs include abnormal all-CDS frame metadata.
- The frame-risk layer contains 27 genes, 24 runs, 218 hotspot rows, and 11 recurrent frame-risk hotspots.
- Recurrent frame-risk genes include `YWHAZ`, `ATP5MC3`, `PCNA`, `RPS6`, and `B2M`.
- `B2M` still has a recurrent frame-risk hotspot at codon 55, transcript position 194, frame1, in `ENST00000648006`, matching browser inspection. It appears in 2 of 3 B2M frame-risk runs.

Interpretation:

- All-CDS frame usage should be used as a technical-risk flag, but not as a hard biological rejection rule.
- The strict queue is now cleaner for conservative claims.
- The frame-risk queue is the right place for cases where abnormal global frame behavior may itself be part of infection-associated ribosome behavior or protocol-specific biology.
- `B2M` is no longer a strict candidate, but it remains a useful browser-validated frame-risk lead.

## 2026-05-27: Off-frame hotspot browser-review queue

Motivation: the hotspot localization table was still too large for direct browser inspection. We needed exact review rows, run maps, and local window plots that separate matched-control evidence from unmatched recurrent hotspots.

Changes:

- Added `scripts/dominant_infected_off_frame_review_queue.R`.
- Added pipeline step `infected_off_frame_review` after `infected_off_frame_hotspots`, then seeded its cache manifest.
- The script converts exact hotspot positions into a ranked browser-review queue, assigns matched/recurrent/context tiers, adds a run map with compact metadata, and plots local +/-30 codon hotspot windows.
- Tier 1 is now reserved for hotspots with matched-control evidence, so strong but unmatched recurrent rows do not get mislabeled as matched review rows.

Outputs:

- `dominant_infected_off_frame_review/infected_off_frame_hotspot_review_queue.csv`
- `dominant_infected_off_frame_review/infected_off_frame_hotspot_review_run_map.csv`
- `dominant_infected_off_frame_review/infected_off_frame_hotspot_review_summary.csv`
- `dominant_infected_off_frame_review/infected_off_frame_hotspot_review_metrics.csv`
- `dominant_infected_off_frame_review/figures/`

Result:

- Produced 159 hotspot review rows across 13 genes.
- Produced 1,135 hotspot/run mapping rows for browser inspection.
- Found 1 tier-1 matched-review row: `B2M` codon 55, frame1. It is enriched in infected samples but also present in matched controls, so it remains a shared/amplified hotspot lead rather than infection-specific proof.
- Found 51 tier-2 recurrent unmatched hotspot rows. Top review genes include `NDUFA1`, `NDUFB8`, `PCNA`, `RPL27`, `ATP5F1E`, `NPM1`, and `UQCRFS1`.

Interpretation:

- This iteration improves human validation efficiency more than biological certainty.
- The exact review queue should be used before promoting any off-frame CDS peak as infection-specific ribosome reprogramming.
- The next useful step is read-length-specific p-site/frame decomposition at the queued hotspot positions, starting with the B2M matched row and the strongest recurrent unmatched rows.

## 2026-05-27: Infected off-frame hotspot localization

Motivation: the previous scan ranked genes/runs with strong off-frame top-CDS signal, but a run-level score is not enough. A biological ribosome-routing candidate should show recurrent positions, or at least a restricted set of local hotspots, rather than arbitrary scattered peaks.

Changes:

- Added `scripts/dominant_infected_off_frame_hotspots.R`.
- The script uses the off-frame candidate gene/run table, loads cached transcript coverage for prioritized genes, localizes top-1% middle-CDS positions, and keeps off-frame hotspots.
- It compares candidate infected runs against matched controls when controls exist, and separately labels unmatched recurrent hotspots.
- Added pipeline step `infected_off_frame_hotspots` after `infected_off_frame_scan`, then seeded its cache manifest.

Outputs:

- `dominant_infected_off_frame_hotspots/infected_off_frame_hotspot_gene_summary.csv`
- `dominant_infected_off_frame_hotspots/infected_off_frame_hotspot_recurrence.csv`
- `dominant_infected_off_frame_hotspots/infected_off_frame_hotspot_run_positions.csv`
- `dominant_infected_off_frame_hotspots/infected_off_frame_hotspot_aggregate_profiles.csv`
- `dominant_infected_off_frame_hotspots/figures/`

Result:

- Localized 47 candidate infected runs across 13 prioritized genes.
- Found 159 off-frame hotspot positions.
- Found 0 matched-specific recurrent infected hotspots under the strict matched-control gate.
- Found 1 matched-enriched shared hotspot: `B2M` codon 55, frame1, present in 4/4 candidate infected runs but also 2/3 matched controls.
- Found 51 unmatched recurrent infected hotspots. The strongest unmatched genes are `NPM1`, `UQCRFS1`, `PCNA`, `NDUFB8`, `RPL27`, `IFITM3`, `NDUFA1`, and `ATP5F1E`.

Interpretation:

- The position-level pass strengthens the idea that some infected-context genes have reproducible off-frame CDS peaks, but it does not yet prove infection specificity.
- `B2M` is the only current candidate with a matched-control comparison, and its best hotspot is enriched but shared, not infection-specific.
- The many EBV/LCL recurrent hotspots are useful browser-review targets, but they remain transformed-cell or cohort-context leads until matched controls are added.
- Next useful step: for top hotspot positions, extract read-length-specific/frame-specific signal if local FST pages support it, or make a browser-review queue with exact gene/codon/run coordinates.

## 2026-05-27: General infected off-frame CDS peak scan

Motivation: browser review suggested that the infected-group signal is not only generic CDS bumpiness. The more important phenotype may be strong CDS peaks that land out of the canonical frame, sometimes with bumpiness and sometimes without it.

Changes:

- Added `scripts/dominant_infected_off_frame_cds_scan.R`.
- The scan reuses `dominant_cds_frame_bumpiness_qc/cds_frame_bumpiness_run_metrics.csv`, so it does not reload FST coverage.
- Added a local metadata sanitizer so rows with `uninfected`, `mock`, `control`, `WT`, or related text are treated as controls before the local metadata CSV is synced from the Google Sheet.
- Restricted the context aggregation to metadata contexts that contain at least one infected row, reducing the scan from about 250k infected/control-labelled rows to 61,773 rows across 22 infected contexts.
- Prioritized top-CDS peak signal in off-frame positions over bumpiness alone. Current hard filters are at least 100 CDS counts and at least 20 counts in the top-peak signal.
- Added pipeline steps `cds_frame_bumpiness_qc` and `infected_off_frame_scan`, then seeded their cache manifests so `run_dominant_cell_state_pipeline("all")` does not rerun the large QC table unless inputs or scripts change.

Outputs:

- `dominant_infected_off_frame_cds_scan/infected_off_frame_cds_candidate_gene_scores.csv`
- `dominant_infected_off_frame_cds_scan/infected_off_frame_cds_candidate_contexts.csv`
- `dominant_infected_off_frame_cds_scan/infected_off_frame_cds_timepoint_contexts.csv`
- `dominant_infected_off_frame_cds_scan/infected_off_frame_cds_top_candidate_runs.csv`
- `dominant_infected_off_frame_cds_scan/figures/`

Result:

- The scan covered 177 genes, 3,857 runs, 21,489 infected evaluable gene-run rows, and 6,467 matched-control evaluable gene-run rows.
- No gene is high-priority under the strict matched-context gate.
- Four genes are medium-priority off-frame peak candidates: `UQCRFS1`, `B2M`, `RPL27`, and `NPM1`.
- Nine additional positive-scoring low-priority candidates are `IFITM3`, `PPP1R15A`, `KRT18`, `ATP5F1E`, `PCNA`, `COX5A`, `IGFBP7`, `NDUFA1`, and `NDUFB8`.
- `RRM1`, `MYC`, and `EGLN3` are not prioritized by the off-frame-peak-specific gate. `RRM1` remains a context-specific bumpy-CDS lead from the previous pass, but it is not currently an off-frame top-peak lead.
- `RPL12` is not in the current 177-gene FST-backed metric table. `MRPL12` is present but not prioritized.

Interpretation:

- The current strongest pattern is a subset effect, not a universal infected-state signature.
- Several top rows are from EBV/LCL, HCoV-OC43/MRC5, Huh7 infection, and selected Calu3/A549/RD contexts. EBV/LCL should be treated carefully because it may represent transformed-cell context rather than acute viral infection.
- The most useful next test is read-length-specific and position-specific review of the top candidate runs, especially `UQCRFS1`, `IFITM3`, `B2M`, `RPL27`, and `NPM1`.

## 2026-05-27: Viral bumpy-CDS context and matched-control audit

Motivation: the previous MYC/RRM1 bumpiness pass showed that `RRM1` had several infected libraries with target-specific bumpy CDS coverage after same-run housekeeping QC, but it did not say whether those runs were broadly distributed across viral infection or concentrated in specific studies/cell contexts. A broad viral claim needs different evidence than a context-specific design lead.

Changes:

- Extended `scripts/dominant_cds_frame_bumpiness_qc.R` to keep an all-run target-vs-housekeeping QC table for MYC/RRM1, not only infected cases.
- Added viral design run maps from the matched-control feature table.
- Added `dominant_cds_frame_bumpiness_qc/myc_rrm1_infected_library_context_summary.csv`.
- Added `dominant_cds_frame_bumpiness_qc/myc_rrm1_viral_design_housekeeping_contrast_summary.csv`.
- Added context and matched-design review figures in `dominant_cds_frame_bumpiness_qc/figures/`.
- Updated the scientific note to distinguish broad viral effects from context-specific review leads.

Result:

- `MYC` is no longer a broad viral bumpy-CDS candidate. It has 1 of 78 candidate-bumpy infected libraries after housekeeping QC, from Huh7/PRJNA285961; that same Huh7 study also contains many technical-frame-bias flags.
- `RRM1` remains the stronger lead, with 13 of 78 candidate-bumpy infected libraries after housekeeping QC.
- `RRM1` bumpy libraries are context-structured: repeated case signals appear in HFF/PRJNA533074, primary fibroblast/PRJNA287112, RD/PRJNA401882, MRC5/PRJNA1062064, and Huh7/PRJNA285961.
- The matched case-control audit is narrower: HFF/PRJNA533074 and MRC5/PRJNA1062064 are the clearest case-enriched `RRM1` designs, while Huh7/PRJNA285961 and CD4/PRJNA667051 are control-enriched.

Interpretation:

- This does not support a universal infected-state `RRM1` CDS-routing claim.
- It does support a context-specific review path: check HFF and MRC5 matched infected/control browser views first, then test read-length-specific p-site/frame behavior around the RRM1 bumpy regions.
- The next model should score viral candidates by cross-design reproducibility, matched-control direction, housekeeping cleanliness, and hotspot recurrence together rather than treating all infected libraries as exchangeable.

## 2026-05-26: Viral CDS frame and bumpiness QC

Motivation: browser review suggested that `RRM1`, `MYC`, and `EGLN3` have bumpy CDS coverage with many apparent out-of-frame peaks in infected samples. The main risk was mistaking library-level p-shift/protocol bias or low-depth sampling noise for viral ribosome reprogramming. The user also noted that `DDIT3` looks correct in the browser, so it can act as a same-run negative-control gene only when it is sufficiently expressed in the same libraries.

Changes:

- Added `scripts/dominant_cds_frame_bumpiness_qc.R`.
- Added per-run middle-CDS frame and bumpiness metrics for all modeled marker CDSs.
- Added run-level global CDS-frame QC across marker genes.
- Added target-gene summaries for `DDIT3`, `EGLN3`, `MYC`, and `RRM1`.
- Added viral design summaries for those target genes.
- Added same-run `DDIT3` control tables requiring at least 30 middle-CDS counts for both `DDIT3` and the candidate.
- Added a housekeeping baseline from always-on genes already available in local FST pages.
- Added explicit per-infected-library `MYC`/`RRM1` frame and bumpiness metrics with housekeeping deltas, DDIT3 control availability, and coarse QC classes.
- Added MYC/RRM1 hotspot localization: per-library top CDS hotspots, recurrent hotspot summaries, and group-normalized CDS profiles by infected-library QC class.
- Added QC figures for bumpiness versus depth, gene frame versus global frame, aggregate CDS-frame profiles, and candidate versus `DDIT3` same-run bumpiness.

Result:

- `MYC` and `RRM1` have positive matched-depth top-1% bumpiness z-scores overall; `EGLN3` and `DDIT3` are weak or negative.
- In same-run viral-case comparisons with sufficient `DDIT3`, `MYC` and `RRM1` each have 3 runs where candidate bumpiness exceeds `DDIT3` by at least one matched-depth z-score unit.
- No same-run `DDIT3` control rows are globally flagged as misshifted.
- `EGLN3` does not survive the same-run `DDIT3` control.
- The housekeeping baseline covers 3,821 runs with median housekeeping canonical frame around 0.59 and 234 housekeeping-frame misshift candidates.
- In the explicit infected-library table, `RRM1` has 13 of 78 libraries classified as candidate-bumpy versus housekeeping without a housekeeping/global misshift flag, while `MYC` has 1 of 78.
- Hotspot recurrence is weak. `MYC` has too few target-specific bumpy libraries to infer recurrence. `RRM1` has a maximum recurrent hotspot frequency of 3 of 13 target-specific bumpy libraries, or about 23%.

Interpretation:

- The data do not support a single simple p-shift artifact explanation for every bumpy viral CDS observation.
- Some candidate behavior still tracks run-level global frame, so library/protocol bias remains an important covariate.
- `RRM1` is now the stronger viral bumpy-CDS review lead, but the weak hotspot recurrence argues against a simple common viral-routing event. `MYC` remains interesting in a small number of libraries but is weakened by the housekeeping and recurrence checks. `EGLN3` should be deprioritized until browser review or a stronger control supports it.
- The next test should use read-length-specific p-site/frame decomposition around the bumpy regions.

## 2026-05-26: DOTSeq-inspired branch-usage model

Motivation: the DOTSeq preprint suggested that the atlas needs a formal branch-usage statistic, not only relative-use deltas and clean-CDS buffering summaries.

Changes:

- Added `scripts/dominant_rdg_branch_usage_model.R`.
- Added the cache-aware pipeline step `branch_usage`.
- Added `dominant_rdg_branch_usage/dominant_rdg_branch_usage_contrast_table.csv`.
- Added `dominant_rdg_branch_usage/dominant_rdg_branch_usage_feature_summary.csv`.
- Added `dominant_rdg_branch_usage/dominant_rdg_branch_usage_gene_design_summary.csv`.
- Added `dominant_rdg_branch_usage/dominant_rdg_branch_usage_review.csv`.
- Added branch-usage review figures in `dominant_rdg_branch_usage/figures/`.
- Merged branch-usage evidence into `dominant_rdg_atlas_cards.csv`.
- Added branch-usage fields and a `Branch-usage shift only` filter to the Shiny atlas browser.
- Updated the scientific note with a "DOTSeq-Inspired Branch-Usage Pass" section.

Result:

- Matched feature contrast rows: 55,044.
- Feature summaries: 3,892.
- Gene x design-family summaries: 994.
- Branch-usage review rows: 452.
- Replicated branch-usage shift rows: 37.
- Single-design branch-usage shift rows: 50.
- Viral-infection branch-usage review rows: 40.
- Atlas-card rows with branch-usage shift support: 58.

Interpretation:

- The strongest viral replicated branch-usage rows are ATF4 clean-CDS up, PPP1R15A uORF5 down, and TRIB3 clean-CDS up.
- DDIT3 and ATF3 remain browser-supported viral candidates, but their viral branch-usage rows are currently `evaluable_no_strong_shift`, which separates visual/browser support from formal branch redistribution evidence.
- The next upgrade should be a real `Run x gene x RDG feature` count table instead of reusing already-aggregated matched contrasts.

## 2026-05-26: Atlas browser full-candidate default

Motivation: the atlas browser should show all review candidates by default, not only the manually verified DDIT3 and ATF3 rows.

Changes:

- Changed the Shiny gene selector default to empty, meaning all genes.
- Changed the default confidence floor to 0 so all review rows are visible unless the user raises the slider.
- Renamed the manual validation filter to `Manually verified only`.
- Kept manual verification as a filter rather than a default restriction.
- Sorted the table by browser-review priority by default.

Interpretation:

- The app now supports broad candidate exploration first, with DDIT3/ATF3 available as a focused manually verified subset.

## 2026-05-26: DOTSeq relevance note

Motivation: a new DOTSeq preprint appeared that directly addresses differential ORF usage from Ribo-seq.

Changes:

- Added a "DOTSeq Relevance" section to the scientific note.
- Interpreted DOTSeq as a useful statistical template for branch-usage testing, not as a replacement for the RDG atlas.
- Added a next-step recommendation to build a `Run x gene x RDG feature` count table and test beta-binomial or Dirichlet-multinomial branch redistribution.

Interpretation:

- DOTSeq's Differential ORF Usage model matches our need for a statistical test of condition-dependent branch allocation.
- DOTSeq still depends on reliable ORF annotation and matched experimental designs, so our manual translons, lab translon predictor, browser validation, and study-aware compendium model remain necessary.

## 2026-05-26: Browser-supported viral RDG validation lane

Motivation: the user inspected DDIT3 and ATF3 in the online RiboCrypt all-vs-infected view and saw interesting translational changes. Those browser findings should become structured atlas evidence, not informal memory.

Changes:

- Added `dominant_rdg_atlas/browser_validation_notes.csv` as a manual browser-validation register.
- Integrated browser-validation notes into `scripts/dominant_rdg_atlas_cards.R`.
- Added `atlas_browser_review_priority_score`, separate from `atlas_confidence_score`.
- Added `dominant_rdg_atlas/dominant_rdg_atlas_browser_validation_queue.csv`.
- Added `dominant_rdg_atlas/dominant_rdg_atlas_viral_infection_review.csv`.
- Added `dominant_rdg_atlas/figures/dominant_rdg_atlas_viral_review_priority.{png,pdf,html}`.
- Updated the Shiny atlas browser with a browser-supported filter and validation fields.
- Updated the cache-aware pipeline dependency list so the atlas-card step reruns when `browser_validation_notes.csv` changes.

Result:

- Browser-supported rows: 2.
- Viral-infection review rows: 68.
- ATF3 and DDIT3 viral-infection cards now rank first and second in the browser-validation queue.

Interpretation:

- This is the first formal human-in-the-loop atlas feedback cycle.
- Browser validation now raises review priority without changing the statistical confidence score.
- The next browser targets should come from the viral-infection queue, especially RRM1, PPP1R15A, MYC, EGLN3, SESN2, IFIT1, RPLP0, ATF4, HERPUD1, VEGFA, DDIT4, HSPA5, STAT1, and MRPS12.

## 2026-05-26: RDG atlas cards and Shiny browser prototype

Motivation: the pipeline had many useful tables, but no single queryable product layer for "gene x condition-family -> predicted clean-CDS buffering plus RDG evidence".

Changes:

- Added `scripts/dominant_rdg_atlas_cards.R`.
- Added the cache-aware pipeline step `rdg_atlas_cards`.
- Added `dominant_rdg_atlas/dominant_rdg_atlas_cards.csv`.
- Added `dominant_rdg_atlas/dominant_rdg_atlas_cards_review.csv`.
- Added `dominant_rdg_atlas/dominant_rdg_atlas_summary_metrics.csv`.
- Added atlas heatmap and confidence/effect figures in `dominant_rdg_atlas/figures/`.
- Added a standalone Shiny prototype in `apps/rdg_atlas_browser/app.R`.
- Added `scripts/run_rdg_atlas_browser_app.R` as a launcher.
- Updated the scientific note with an "RDG Atlas Cards" section.

Result:

- Atlas rows: 994.
- Review rows: 870.
- Genes: 71.
- Design families: 14.
- Atlas-ready replicated interval rows: 137.
- Replicated strict review rows: 513.
- Single-design strict review rows: 220.

Interpretation:

- This is the first compact atlas product layer: matched clean-CDS buffering, top shifted RDG feature, sparse on/off evidence, post-viral evidence, QC, browser-review pointers, and warnings are now in one table.
- Viral infection is currently the strongest card-level context, but high between-design heterogeneity is common and remains biologically important rather than ignorable.
- The Shiny app is intentionally simple and should become a RiboCrypt module after the card schema survives browser validation.

## 2026-05-26: AlphaGenome-inspired RDG modelling note

Motivation: the user asked whether the new AlphaGenome paper is relevant to the dominant-state/RDG atlas, and what ideas should be borrowed.

Changes:

- Added an "AlphaGenome Relevance" section to the scientific note.
- Summarized AlphaGenome as a unified sequence-to-function and variant-effect model, useful as a design template but not as a direct ribosome-routing solution.
- Added a proposed RDG analogue: shared transcript/RDG representation with output heads for clean-CDS usage, ouORF/iORF/NTE usage, branch-switch probability, dominant-state score, and CDS-buffering effect.
- Added AlphaGenome-style counterfactual scoring as a recommended next model direction.
- Noted Translatomer as the closer translation-specific precedent because it predicts cell-type-specific Ribo-seq from mRNA expression and sequence.

Interpretation:

- The key reusable ideas are multimodal shared representation, reference-versus-alternative scoring, model distillation, and in silico mutagenesis/contribution scanning.
- The key gap is that AlphaGenome models DNA regulatory and RNA/splicing readouts, while our main problem is ribosome decision flow through transcript architecture under condition-specific perturbations.

## 2026-05-26: Forecasting-inspired atlas model note

Motivation: the user asked whether fields such as weather prediction offer useful math, statistics, or algorithm ideas for the RDG atlas.

Changes:

- Added an "Algorithmic Inspiration From Adjacent Fields" section to the scientific note.
- Reviewed ensemble weather forecasting, data assimilation/ensemble Kalman filtering, multi-task Gaussian processes, causal forests, and dynamic causal modelling.
- Reframed the RDG atlas as a calibrated probabilistic forecasting system: predictive distribution plus uncertainty, nearest evidence, and RDG-route explanation.
- Added ensemble/hierarchical forecasting as a recommended next model direction.

Interpretation:

- Weather and data-assimilation methods are the strongest analogy because they update latent state estimates from noisy, incomplete, heterogeneous observations and verify predictions on held-out cases.
- Causal forests and related nonparametric heterogeneity screens are useful for triage, but should not be treated as causal proof unless a candidate is supported by matched designs and browser-level evidence.

## 2026-05-26: Shared design-family semantics

Motivation: sparse on/off candidates were still grouped mostly as `other_non_control`, which made replication hard to interpret biologically.

Changes:

- Added `scripts/dominant_condition_families.R` as a shared condition classifier.
- Added finer design families: viral infection, interferon/cytokine, ISR/ER/translation stress, nutrient starvation, hypoxia/oxidative stress, heat/osmotic/recovery, DNA damage/senescence, tumor/cancer context, drug/stimulus, loss-of-function, gain-of-function, genetic variant/mutant, immune/infection-model context, protocol/fraction, and other.
- Added `design_family` to next-model matched feature contrasts and internal-ORF matched bump evidence.
- Added `dominant_next_model_condition_design_family_audit.csv`.
- Updated sparse ouORF/CDS coupling to use `design_family` before the older broad `perturbation_class`.
- Updated post-viral/fatigue summaries to use the shared classifier and compatible post-viral family sets.

Result:

- Condition-family audit: 16 design families.
- Largest family remains `other_non_control`, so exact study-level curation is still needed.
- Sparse on/off design summaries increased to 324 because coarse classes were split into more specific families.
- Replicated sparse candidates remain 25, but now resolve into tumor/cancer context (19), viral infection (2), nutrient starvation (2), and drug/stimulus (2).
- Sparse branch-switch review candidates increased from 0 to 1: MKI67 uORF9 in tumor/cancer context, with only 3 rows across 2 studies and 2 cell lines.
- Post-viral matched viral/IFN review rows narrowed from 3,124 to 1,917, with 40 strict supported rows.

Interpretation:

- This is a better biological grouping layer, not yet final metadata curation.
- The MKI67 branch-switch row is review-grade only.
- The next iteration should distinguish real tumor biology from tumor-labeled study design and continue shrinking `other_non_control`.

## 2026-05-26: Sparse on/off design-family model

Motivation: one-sided count-supported rows are useful, but row-level triage is still too noisy. We need a compact layer that asks whether sparse on/off events replicate across studies, cell lines, and design families, while staying separate from strict matched log2FC evidence.

Changes:

- Added design-family aggregation for one-sided matched ouORF/CDS rows.
- Added sign-consistency, count-strength, replication, specificity, and pairing components to a bounded sparse on/off score.
- Added sparse on/off candidate classes:
  - `replicated_specific_ouorf_onoff_review`;
  - `broad_replicated_ouorf_onoff_review`;
  - `replicated_sparse_shared_onoff_review`;
  - `sparse_branch_switch_review`.
- Added `dominant_ouorf_cds_coupling_sparse_onoff_design_summary.csv`.
- Added `dominant_ouorf_cds_coupling_sparse_onoff_candidates.csv`.
- Added sparse candidate-score and design-family heatmap figures.

Result:

- Sparse design-family summaries: 193.
- Replicated sparse candidates: 25.
- Broad replicated ouORF on/off candidates: 15.
- Specific replicated ouORF on/off candidates: 7.
- Replicated shared ouORF + clean-CDS on/off candidates: 3.
- Sparse branch-switch candidates: 0.

Interpretation:

- ATF3 uORF2, SESN2 uORF3/uORF4, IFIT1 uORF1, CDK1 uORF5, MKI67 uORF9, and SNAI2 uORF2 are the clearest sparse review rows.
- The current sparse evidence supports on/off behavior and co-loading more than true ouORF-versus-clean-CDS routing switches.
- The next model should improve design-family semantics because many rows still collapse into `other_non_control`.

## 2026-05-26: Count-aware on/off triage for matched ouORF/CDS rows

Motivation: the strict both-side count gate correctly removes sparse rows from default matched log2FC review, but some rows have strong raw-count support on one side and near-absent support on the other. Those can be real induction/repression events, so they need a separate review lane rather than being discarded or mixed with ordinary fold-change rows.

Changes:

- Added raw-count log2(case/control), approximate count-ratio SE, z, p, and q values for one-sided matched ouORF/CDS rows.
- Added `onoff_pattern`, `onoff_review_class`, `onoff_evidence_score`, and interpretation columns.
- Added feature-level on/off summaries.
- Added `ouorf_cds_coupling_onoff_review` PNG/PDF/HTML plots.
- Added the new on/off outputs to the pipeline output map.

Result:

- One-sided matched ouORF/CDS review rows: 1,092.
- Features represented in the on/off summary: 27.
- Paired ouORF + clean-CDS on/off rows: 142.
- The PRJNA256316 CDK1 uORF5 case-supported/control-low example now appears in this on/off triage layer, not in default strict matched review.

## 2026-05-26: Low-count gates for matched ouORF/CDS review

Motivation: manual inspection of PRJNA256316 CDK1 uORF5 WT_tumor versus WT showed only 1-2 reads over the reviewed uORF, but the row could still appear as an enriched matched contrast.

Changes:

- Required at least 2 samples in both case and control for matched contrasts.
- Raised the core matched feature raw-count floor to 30.
- Added any-side and both-side count-gate columns to matched feature outputs.
- Changed the default ouORF/CDS matched summary and quadrant plot to require 30 raw reads on both case and control sides for both ouORF and clean-CDS.
- Kept low-count rows in CSV outputs with explicit labels: `low_ouorf_counts`, `low_clean_cds_counts`, `one_sided_low_ouorf_counts`, and `one_sided_low_clean_cds_counts`.
- Raised clean-CDS raw-count floors in RDG, structure-rule, and score-semantics scripts to 30.

Result:

- Matched ouORF/CDS rows: 11,704 total.
- Sample gate: 7,056.
- Any-side count gate: 4,713.
- Strict both-side gate used by default review plots: 3,621.
- PRJNA256316 CDK1 uORF5 low-read rows are no longer default-review positives.

## Downstream-TIS rescue audit

Motivation: positive residual ouORF/CDS coupling might reflect downstream TIS usage or truncated CDS products rather than ordinary canonical CDS translation.

Changes:

- Added `dominant_downstream_tis_rescue_audit.R`.
- Reused cached transcript coverage and sequence instead of reopening FST pages.
- Scanned downstream clean-CDS sequence for ATG and near-cognate starts.
- Annotated Kozak context, frame, predicted downstream ORF length, and total post/pre coverage shape.

Result:

- Eight residual-positive ouORF/CDS review features were scanned.
- 1,838 candidate downstream starts were found.
- IGFBP5 uORF10 is currently the only sequence-plus-total-coverage downstream-TIS review lead.
- ATF4, CDK1, HSPA5, LMNB1, MKI67, MMP2, and SNAI2 remain shape-only or weaker review rows.

## ouORF/CDS coupling audit

Motivation: many overlapping-uORF and clean-CDS absolute levels correlated positively, which conflicts with a simple competition model.

Changes:

- Added `dominant_ouorf_cds_coupling_analysis.R`.
- Separated raw absolute correlation from total-route residual coupling.
- Added matched-control ouORF/clean-CDS log2FC and relative-use quadrant outputs.
- Added browser-review candidate prioritization for residual-positive and possible truncated-CDS rescue rows.

Result:

- 28 modeled ouORFs show raw positive ouORF/CDS correlation.
- 14 look like shared-route effects after residualization.
- 14 remain positive residual-coupling review rows.
- Eight are possible downstream-start/truncated-CDS rescue review rows.

## Internal-ORF evidence gates

Motivation: internal ORFs can share ordinary CDS elongation coverage, so metadata or body-coverage shifts are not enough to claim independent initiation.

Changes:

- Added matched case/control internal-ORF feature/CDS bump evidence.
- Added internal-body frame and local start-window diagnostic summaries from cached coverage.
- Downgraded canonical-frame-like internal ORFs in review outputs.

Result:

- 2,971 internal-ORF matched condition rows.
- 235 current review rows after stricter gates.
- 0 BH-supported matched-condition internal-ORF bump rows.
- 18 internal-ORF initiation-evidence rows: 16 canonical-frame-like, 2 feature-frame review rows.

## Exon-support and canonical isoform checks

Motivation: wrong representative isoforms can create false clean-CDS or uORF interpretations.

Changes:

- Switched to `canonical_isoforms(df)` after local canonical definitions became available.
- Added aggregate CDS exon support checks on the `human_all_merged_l50` track.
- Demoted transcript candidates with severely under-supported CDS exons when better candidates exist.

Result:

- ATF3 and other genes were rechecked under the corrected canonical isoform layer.
- Exon-level QC is now written to `human_dominant_cell_states_cds_exon_coverage_qc.csv`.

## Post-viral fatigue / ribosome-stress module

Motivation: test whether viral infection and post-viral biology show host ribosome routing changes involving uORFs and iORFs.

Changes:

- Added a post-viral fatigue/ribosome-stress dominant-state module.
- Added submodules for interferon/antiviral, ISR/ribosome stress, senescence/SASP, hypoxia/inflammation, suppressed growth/translation, and suppressed OXPHOS.
- Added focused post-viral matched design tables, leave-one-design and leave-one-study diagnostics, and strict random-effects proxy intervals.
- Added a validation queue combining candidate rankings, counterfactuals, matched design evidence, and iORF diagnostics.

Result:

- The fatigue axis is coherent but not disease-specific.
- Viral infection is the strongest current matched-design test case.
- RPLP0, CDK1, ATF3, and PPP1R15A currently have replicated strict viral-infection matched support.
- No post-viral gene x design-family row is atlas-ready under the strict interval gate.

## Score semantics and OXPHOS audit

Motivation: raw winner labels were dominated by OXPHOS, making the branch labels hard to interpret.

Changes:

- Added raw signed-score versus within-state percentile winner comparisons.
- Added OXPHOS arm/component summaries.
- Added design-aware score-semantics fixed-effect summaries.

Result:

- Raw OXPHOS winners remain common, but percentile labels are more balanced.
- Raw signed OXPHOS should be treated as a quantitative module balance rather than a simple dominant-state winner.
- Percentile branches are better for compendium-relative QC labels.

## Expanded metadata and condition modeling

Motivation: CELL_LINE, TISSUE, CONDITION, and GENE were too narrow for the compendium.

Changes:

- Added broad metadata modeling across selected columns.
- Added n-wise metadata terms.
- Added glmnet and random-forest summaries.
- Added PCA/latent-program summaries.
- Added control-aware condition contrasts where control samples were available.

Result:

- Metadata explains a large fraction of state-score variance.
- PCA shows the dominant-state matrix is highly compressible.
- The same result is also a warning: strong design and cell-context effects must be controlled before calling biology.

## Matched-control and CDS-buffering model

Motivation: broad enrichment cannot distinguish gene/cell/study confounding from condition effects.

Changes:

- Added matched control detection for WT/control/mock/vehicle/untreated and related labels.
- Added matched clean-CDS effect sizes, bootstrap uncertainty, empirical shrinkage, and design-diversity scoring.
- Added generic gene x perturbation-class evidence outside the post-viral special case.

Result:

- Current matched-control feature contrast rows: 53,376.
- Current clean-CDS matched buffering rows: 13,632.
- Current generic replicated strict interval review rows: 6.
- These are review-level summaries, not calibrated causal effects.

## RDG graph construction and relative-use matrices

Motivation: move from dominant-state QC scores to transcript-level ribosome decision graph inspection.

Changes:

- Implemented RDG nodes, edges, feature annotations, branch contrasts, and state-specific RDG plots.
- Added relative-use matrices by state.
- Added plot refinements for state labels, legends, and clean-CDS/uORF color interpretation.

Result:

- RDG plots are now available for the modeled gene panel and state-specific branches.
- Relative-use matrices became the easier human-review object for many genes.

## T + TC + M translon superset

Motivation: ORFik-only predicted uORFs had false negatives for several biologically important genes.

Changes:

- Added Transcode/TC translons.
- Added manual M translons from `manual_translons/translons_to_be_added_to_manual.txt`.
- Deduplicated identical start/stop features across sources.
- Restricted clean-CDS masking to true upstream ORF classes.

Result:

- Manual overlapping-uORF additions materially changed SNAI2, MMP2, LMNB1, MKI67, IGFBP5, HSPA5, DDIT4, RPLP0, LDHA, and other candidates.
- Manual curation is now a temporary bridge until the lab translon predictor improves overlapping-uORF detection.

## FST clean-CDS rewrite

Motivation: count-table CDS totals were inflated by uORFs overlapping CDS, especially for genes such as ATF4.

Changes:

- Replaced count-table CDS scoring with transcript FST coverage.
- Measured leader + CDS coverage per selected transcript.
- Masked CDS bases overlapped by upstream ORFs.
- Wrote cleaned dominant-state score tables and gene diagnostics.

Result:

- The dominant-state score table became more biologically interpretable for uORF-rich genes.
- The method became dependent on correct FST pages and correct transcript selection.

## Initial dominant-state QC prototype

Motivation: proliferation, stress, interferon, hypoxia, EMT, senescence, OXPHOS, and translation-capacity programs can dominate differential RNA/Ribo analysis.

Changes:

- Defined marker modules.
- Loaded `all_samples-Homo_sapiens`.
- Computed sample-level module scores.
- Added `Run` from `runIDs(df)` for metadata matching.
- Added metadata enrichment by cell line, tissue, condition, gene, and later broader metadata fields.

Result:

- Produced `human_dominant_cell_states.csv` and later FST-cleaned replacements.
- Established the core idea: dominant-state QC annotations are useful, but the raw winner label is too crude for the final RDG atlas.

## Current Next-Iteration Themes

- Scale from marker genes to all confidently expressed uORF/uoORF/iORF/NTE/NTT genes.
- Replace heuristic candidate ranking with a hierarchical gene x branch x study x condition model.
- Add count-aware on/off models for one-sided low-count but potentially real feature induction.
- Integrate lab translon prediction and denoised p-shifted coverage before promoting iORF or downstream-TIS claims.
- Make counterfactual CDS-buffering predictions queryable with uncertainty and nearest observed designs.
- Keep post-viral/viral infection as the main biological demonstration, but avoid disease-specific claims until matched cohorts or stronger metadata support exist.
