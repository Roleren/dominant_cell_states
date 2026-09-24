# Dominant Cell States and Ribosome Decision Graphs

Updated: 2026-08-19

## 2026-08-18 Focused Discovery Result

The first genuinely high-contrast, artifact-resistant biological candidate is
now the 51-nt ATG uORF1 in the EIF4G2/DAP5 leader. The most novel observation is
that tRNA-Glu(UUC) overexpression, but not matched metastasis-promoting
tRNA-Arg(CCG), redistributes ribosome occupancy from the EIF4G2 clean CDS into
this exact uORF. The tRNA-Glu result has two biological replicates per arm after
three sequencing lanes per replicate are correctly collapsed, complete
replicate-range separation, phase-0-dominant uORF signal, coordinate and peak
deletion robustness, stable matched total RNA/isoform composition, and the top
Glu-minus-Arg specificity rank among strict count-supported shared leader
uORFs. A direct transcriptome-wide Glu-codon-dose model is null, so the two GAG
codons are a reporter hypothesis rather than an established cause.

The final raw tRNA artifact gate now passes. All 12 tRNA-Glu/control and four
tRNA-Arg/control FASTQs were integrity checked, length preserved, uniquely
aligned and calibrated independently by run and footprint length. In the
predeclared lane-matched pool, the tRNA-Glu union and exact-uORF1 allocation
effects were -0.136 and -0.135, with non-overlapping true biological-unit
ranges; the matched tRNA-Arg effect was +0.012 and the Glu-minus-Arg
specificity contrast was -0.147. The result also passes in study-wide shared
lengths (-0.141), canonical 28--32-nt footprints (-0.145), all three matched
lanes, coordinate perturbations, every single-position deletion and greedy
three-position deletion. The sparse tRNA 20--23-nt stratum failed its own
measurement and exact-uORF1 phase gates and is not counted as evidence.

The strongest independent recurrence is UV irradiation in HeLa. All four raw
FASTQs were reprocessed with read length retained, unique GRCh38 alignment,
per-run/per-length ORFik offset detection, independent CDS-start phase QC, and
exact EIF4G2 mapping. The clean-CDS allocation shift survives separately in
shared well-phased canonical 28--32-nt footprints (UV minus control = -0.341)
and short 20--23-nt footprints (-0.274); both replicate ranges are fully
separated. The all-shared effect is -0.343, remains -0.313 under the worst
coordinate perturbation, -0.330 after the worst single-position deletion, and
-0.311 after greedy removal of the three positions that most weaken it. This
closes the principal short-footprint/P-site/pooling artifact gate for the UV
recurrence. ABCE1 loss supplies a second, independent recycling perturbation.

The biological endpoint is unusually concrete: a 2025 Nature Communications
study independently demonstrated that the same EIF4G2 uORF produces an
HLA-presented peptide and is inducible in a leader/start-dependent reporter
under mitotic arrest. Current exact literature searches found no prior link
between tRNA-Glu(UUC) and translation of this EIF4G2 immunogenic uORF. The
candidate paper-level novelty is therefore a tRNA-pool-specific switch into an
experimentally proven immunogenic uORF, with collision/recycling recurrence—not
discovery of the uORF or generic tRNA/codon regulation.

This is a strong computationally validated, experiment-ready hypothesis, not a
demonstrated mechanism. It still needs tRNA-Glu perturbation reporters,
uORF-start and GAG/GAA mutations, mature-tRNA versus tRNA-fragment separation,
DAP5 protein measurements, and HLA immunopeptidomics. A later independent
MDA-LM2 dataset does not reproduce a general metastatic-state direction, so no
broad metastasis claim is made. The complete evidence, exact numbers,
falsifications, literature boundary, runtime recovery paths, and remaining
experiments are recorded in `md/eif4g2_recycling_stress_working_note.md`.

## Project Review Verdict

The project is still aligned with the original goal: use thousands of Ribo-seq samples plus metadata to compress translation behavior into interpretable dominant-state scores, ribosome decision graphs (RDGs), and CDS-buffering estimates that a biologist can inspect. The current pipeline is not yet a final discovery engine, but it is now a useful candidate generator and review system.

The strongest current result is the framework itself: marker-state scoring, clean-CDS quantification, uORF/NTE/iORF-aware branch summaries, matched metadata contrasts, browser validation queues, atlas cards, the consensus atlas, validation dossiers, evidence-gap triage, active-review batches, direct review links, the Shiny Review Workbench and Inspection Assets pages, review-feedback ingestion, and the validation-label bridge form a coherent loop. The strongest biology-supported leads remain viral/post-viral translation effects around ISR and antiviral genes. DDIT3 is the current browser-supported consensus anchor, ATF4 remains the best positive-control RDG, and PPP1R15A, IFIH1, VEGFA, BBC3, and ATF3 are prioritized follow-up genes with different evidence classes. The newest manual-translon iteration also pulls PELO and HK2 into the review space, which is useful because they connect ribosome rescue/quality control and glycolytic stress to the same atlas framework.

The current claims should be framed as review-grade evidence, not as proof of a new disease mechanism. The newest iterations tested whether count-aware branch features, raw RDG graph geometry, constrained mechanistic motif scores, joint branch allocation, a lightweight study-aware hierarchical joint-allocation layer, an empirical-Bayes Dirichlet-multinomial posterior-predictive layer, explicit viral context interactions, a first partial-pooled context layer, model-disagreement audits, therapeutic-hypothesis prioritization, quantitative therapeutic perturbation predictions, auxiliary dominant-state discovery, auxiliary-state/RDG context auditing, a consensus atlas, browser-validation dossiers, and manual translon refreshes improve the model. Branch features add only a tiny predictive improvement beyond metadata/context. Raw graph-geometry covariates still make leave-one-gene-out prediction worse than context-only. Constrained motif scores are interpretable and correlate weakly with observed branch shifts, especially for clean-CDS rows, but they are not yet predictive enough for counterfactual claims. The strongest statistical layers are now branch allocation, model agreement, leave-one-study-out stability, grouped-omission transportability, explicit context interactions, partial pooling, and validation dossiers. These separate estimator consensus, single-study influence, cross-context robustness, tissue/cell-line dependence, single-context signals, and possible buffering/total-load compensation. The therapeutic layers compress atlas evidence into preclinical/biomarker hypotheses and half-normalization perturbation targets; they are not treatment recommendations. The consensus atlas and validation dossiers are interpretation and review aids, not independent validation. Therefore branch, motif, and geometry evidence should currently be used for prioritization, mechanistic annotation, and hypothesis generation, not as strong standalone predictive layers.

## Standalone Repository Handoff

The analysis is being split out of the RiboCrypt source tree. The old path is:

```text
/home/roler/Desktop/forks/RiboCrypt/dominant_cell_states
```

The intended standalone repo path is:

```text
/home/roler/Desktop/forks/dominant_cell_states
```

Future agents starting from the standalone path should first read `AGENTS.md`, `README.md`, this scientific note, `md/dominant_rdg_atlas_browser_tutorial.md`, and `NEWS.md`. The standalone repo should track source-like files, docs, tests, app code, and the manual translon text input. It should not blindly commit the generated output tree: the current local directory is about 4.7 GB and includes multiple generated CSVs over 100 MB.

The standalone checkout is now an R package scaffold. Package/source files live
in `R/`, `tests/testthat/`, `scripts/`, `inst/shiny/`, docs, and
`manual_translons/translons_to_be_added_to_manual.txt`. Generated results and
app-derived curated CSVs live in ignored `results/`; root-level generated output
folders or compatibility symlinks should not be kept. FST pages, ORFik
annotations, BAMs, reference indexes, and other Bio_data artifacts are external
runtime inputs, not package data.

The analysis still depends on local RiboCrypt code. Use `devtools::load_all()` against the local RiboCrypt repo, not an installed package. In standalone mode set:

```bash
export RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt
```

Then run R commands as:

```bash
env -u LC_ALL R --vanilla -q -e 'devtools::load_all(Sys.getenv("RIBOCRYPT_REPO")); ...'
```

The pipeline entry point now separates the RiboCrypt package root from the analysis directory. From the standalone path, `repo_root` should resolve from `RIBOCRYPT_REPO`, while `analysis_dir` should resolve to the current `dominant_cell_states` checkout.

## 2026-07-10 Current State Summary

The project now has three separable products.

First, it has a dominant-state scoring layer. Per-run scores are computed from marker gene sets after FST-backed clean-CDS quantification. The active core states are proliferation/cell cycle, mTOR/TOP translation capacity, ISR/ER stress, interferon/antiviral, hypoxia/HIF1A, EMT/mesenchymal shift, senescence/SASP, OXPHOS/mitochondrial state, post-viral fatigue/ribosome stress, and baseline controls. EMT, senescence, OXPHOS, and post-viral fatigue are signed modules; the other biological states are mostly up-only marker programs. OXPHOS has already been expanded to 60 markers because the original compact set made mitochondrial state too dominant and too easy to confuse with housekeeping or study effects.

Second, it has an RDG quantification layer. For each selected canonical transcript, the pipeline loads FST coverage, masks clean CDS against predicted and manual overlapping uORFs, treats NTE/NTT/iORF candidates as branch features rather than clean-CDS denominators, and aggregates grouped branch counts before testing. This grouping is essential: single samples usually do not have enough reads on short uORFs. The branch model currently includes independent count-aware beta contrasts, a joint branch-allocation approximation over leader-uORF/overlapping-uORF/clean-CDS/other mass, a study-aware hierarchical branch layer, and an empirical-Bayes Dirichlet-multinomial posterior-predictive screen. These are ranking and review layers, not a final generative likelihood.

Third, it has an interpretation layer. Metadata enrichment, matched controls, control-aware contrasts, broad metadata models, glmnet/latent summaries, leave-one-study-out tests, grouped-omission transportability, explicit viral context interactions, partial pooling, auxiliary-state/RDG overlap warnings, therapeutic hypothesis tables, consensus atlas rows, validation dossiers, evidence-gap priorities, active-review batches, review links, review-feedback labels, validation-label bridge tables, and the Shiny atlas app compress thousands of samples into reviewable packets. The highest-value output is now the validation dossier plus evidence-gap queue plus active-review loop: the dossier says what to inspect in RiboCrypt, the evidence-gap queue says which evidence layer is missing, the batch/link files turn that into a practical review plan, the feedback ingest turns completed manual rows into calibration-ready tables, and the label bridge merges curated browser notes with feedback labels before ranking the next label gaps.

The pipeline now separates RDG table regeneration from full RDG figure regeneration. `run_dominant_cell_state_pipeline(..., rdg_plot_mode = "tables")` or the CLI flag `--skip-rdg-plots` reruns the RDG table layer while preserving existing plot manifests from already written PNG/PDF files. Full plotting remains the default through `rdg_plot_mode = "all"`. The plot mode is written into `.runtime/rdg_plot_mode.txt` and included in the `rdg_graphs` cache fingerprint, so table-only and full-plot runs do not reuse the wrong cache state.

Current result scale:

- Atlas: 1,680 gene-design rows over 120 genes and 14 design families.
- Consensus review queue: 1,293 rows, including 1 browser-validated consensus row, 48 replicated robust rows, 82 multi-model branch rows, 29 context-specific rows, 23 buffering/compensation rows, 103 therapeutic-translation rows, 387 low-support background rows, and 606 fragile/confounded rows.
- Validation dossiers: 184 primary review packets across 72 genes, with 380 matched branch contexts; 372 of those contexts pass the 30-count gate on both sides and all 380 have merged replicates on both sides. Among primary packets, 145 have exact strict pooling proof and 53 carry a pooling/protocol hazard flag.
- Evidence-gap triage: 1,680 rows ranked across 120 genes; 105 high-value acquisition rows; 5 near-counterfactual candidates; 129 browser-validation rows; 10 translon/measurement-review rows; 435 context/metadata-audit rows; 449 rows with protocol/pooling need >= 0.55, of which 203 now get `protocol_pooling_audit` as their primary action; 901 rows whose primary need is more replication or robustness; and 250 suggested negative-control rows.
- Active-review batches: 154 rows over 69 genes in 7 batches: 22 candidate-positive browser checks, 11 counterfactual-readiness checks, 20 protocol/pooling audits, 28 context/confounder audits, 15 measurement/translon reviews, 28 replication/transport follow-ups, and 30 negative controls. The counterfactual overlay intentionally reuses ATF4, PPP1R15A, DDIT3, ATF3, VEGFA, and IFIH1 viral rows from other batches because it asks a different question.
- Review links: all 154 active-review rows now have RiboCrypt Observatory URLs. With the metadata drive mounted, 134 rows use exact context case/control subsets and 20 use study+condition fallback because some context fields did not match both sides. Zero rows fell back to gene-only links. Links include an empty third `All merged` selection as a temporary Observatory restore guard.
- Review feedback ingest: the feedback step scans the 154 linked review rows, overlays any draft labels saved from the atlas browser, and writes normalized label, training, browser-validation-note candidate, summary, and batch-summary outputs. The current linked template has no filled manual rows and no production draft labels yet, so the label, training, and browser-note candidate CSVs are intentionally empty and annotated as expected in the result-sanity inventory.
- Validation-label bridge: the new bridge combines the two curated browser-supported viral labels, DDIT3 and ATF3, with future active-review feedback labels. It currently has 2 positive seed labels over 2 gene-designs, 970 motif-overlay rows, and a 154-row label-gap queue with 152 rows still unlabeled. The label-balance table shows the expected current weakness: positives exist, but weak positives, negatives, confounded rows, measurement-problem rows, and negative controls are all missing.
- Merged-replicate translon-shift gate: 22,896 joint rows, 15,840 exact same-condition merged rows, 8,938 clean-CDS/uORF opposition rows, and 6,434 exact merged opposition rows. The inspectable row queue has 95 rows, collapsed to 60 candidate-direction rows. Eight collapsed candidates have multi-study exact recurrence, one has multi-condition exact recurrence, and 18 have direction conflicts elsewhere in the same gene/transcript/design. No row is claim-ready; the strongest viral rows are review leads requiring protocol/start-site/pooling audit.
- Viral model-agreement set: DDIT3, IFIH1, VEGFA, ATF4, BBC3, and PPP1R15A. DDIT3 is the current browser-validated consensus anchor. ATF3 remains manually interesting, but the formal viral model-agreement layer currently ranks it below the replicated branch-consensus set.
- New manual-translon candidates: PELO is now rank 7 in the next-model candidate table and HK2 is rank 11. PELO viral infection is context-specific rather than robust consensus; HK2 has stronger ISR/hypoxia-style evidence than viral evidence.
- Strongest biological interpretation: viral/post-viral contexts repeatedly touch ISR and antiviral RDGs, but the evidence is still candidate-grade. DDIT3 and PPP1R15A remain cross-context anchors, ATF4 and VEGFA are directionally robust but support-sensitive, IFIH1 is context-dependent, and PELO/HK2 are new biology-expansion leads needing browser validation.
- Strongest negative lesson: raw RDG geometry and hand-built motif covariates do not yet generalize well enough across genes for prediction. The useful direction is a calibrated, constrained, branch-aware model trained against validated graph motifs and explicit branch-allocation evidence.

## 2026-06-15 Literature And Gap Analysis

The literature supports the main direction but argues for more biological separation. Viral infection can suppress host translation, degrade host mRNAs, block nuclear export, and impair translation of induced immune transcripts; SARS-CoV-2 is a concrete example where host protein synthesis is hit through multiple mechanisms ([Finkel et al., Nature 2021](https://www.nature.com/articles/s41586-021-03610-3)). Nsp1 can bind the 40S mRNA channel and inhibit mRNA binding while viral leader sequence helps viral RNA translation, which is directly relevant to start choice and branch routing ([Schubert et al., Nature Structural & Molecular Biology 2020](https://www.nature.com/articles/s41594-020-0511-8)). The eIF2-alpha/ISR literature also fits the RDG model: eIF2-alpha phosphorylation lowers global initiation while allowing selected uORF-controlled transcripts such as ATF4 and related stress genes to be preferentially translated ([Wek, Cold Spring Harbor Perspectives in Biology 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6028073/)). Long COVID and ME/CFS reviews point to heterogeneous mixtures of viral persistence/reactivation, immune activation, autonomic/vascular symptoms, and energy-metabolism abnormalities, which means post-viral RDG work should be framed as a mechanistic sub-hypothesis, not as a complete disease model ([Davis et al., Nature Reviews Microbiology 2023](https://www.nature.com/articles/s41579-022-00846-2)). AlphaGenome-like sequence-to-function models are relevant as inspiration for counterfactual and sequence-prior modeling, but they do not replace Ribo-seq branch measurements because they predict genome tracks and variant effects rather than RDG branch allocation ([AlphaGenome, Nature 2026](https://www.nature.com/articles/s41586-025-10014-0)). The Tierney ribosome decision graph paper remains the conceptual basis for treating a transcript as a flow graph rather than a single CDS value ([PMC10659439](https://pmc.ncbi.nlm.nih.gov/articles/PMC10659439/)).

Biology gaps to address next:

- Split the ISR output score into upstream cause axes: PKR/EIF2AK2 antiviral dsRNA, PERK/EIF2AK3 ER stress, GCN2/EIF2AK4 nutrient/ribosome-collision stress, and HRI/EIF2AK1 heme/mitochondrial stress. Current ISR markers mostly see ATF4-axis consequences, not the initiating kinase.
- Add RNA sensing, RNase L, stress-granule, and host-shutoff markers. The current antiviral module is too close to ISG expression and underrepresents DDX58/MAVS/TBK1/IRF3, OAS2/OAS3/RNASEL, ZC3HAV1, DDX3X, G3BP1/2, TIA1/TIAL1, NUP98/RAE1/NXF1, and EIF4F remodeling.
- Expand proteostasis/heat-shock, autophagy/lysosome, and NRF2/ferroptosis before making those auxiliary states core. They are biologically plausible, but current FST marker coverage is too thin.
- Keep antigen presentation, endothelial/coagulation, complement, sex, cancer status, and immune composition primarily as context/confounder layers. They can dominate metadata enrichment but are not automatically RDG mechanisms.
- Treat ribosome quality control and collision biology as a missing bridge. ZNF598/GIGYF2/PELO/HBS1L/ASCC3-like pathways could explain branch-specific accumulation or apparent stalled coverage without requiring a classical uORF rule.
- Keep DecodeME and other ME/CFS genetics as a prior layer, not a replacement for translation-state analysis. Common variants can change risk modestly while infection/stress still changes RDG state in ordinary public Ribo-seq samples.

Algorithm gaps:

- Replace the current collection of branch screens with a true multilevel branch likelihood, ideally beta-binomial or logistic-normal/Dirichlet-multinomial, with gene, branch, study, design, tissue/cell line, protocol, and condition effects.
- Separate causal/context questions explicitly. Grouped omission and partial pooling are useful stress tests, but they do not identify whether infection, tissue, protocol, or cell line caused a branch shift.
- Add a measurement model for isoform choice, exon support, branch overlap, p-shift/read-length bias, and branch-specific count uncertainty. Frame/iORF claims should stay diagnostic until this layer is stronger.
- Use browser/manual validation as active-learning labels. ATF4, DDIT3, ATF3, PPP1R15A, PABPC1 tumor/kidney, and negative controls should calibrate motif priors.
- Build a counterfactual prediction target around protein-output effect size: "condition X is predicted to change functional clean-CDS output by Y percent with uncertainty." That is the biologist-facing endpoint; RDG branch diagrams explain how the endpoint arises.
- Compress data into R-callable sufficient statistics, not token-heavy tables: per-gene RDG branch matrices, context summaries, model components, and validation packets should be the primary artifacts.

FST/page product:

The FST-page discovery layer has been collapsed into one canonical product. The script `scripts/dominant_cell_state_required_fst_pages.R` resolves the full atlas target set against the mounted `coverage_index.fst`: base dominant-state markers, expanded OXPHOS markers, post-viral/fatigue markers, curated additional prior genes, and genes with manual translon annotations. The output is the only page list needed for preparing refreshed FST pages.

Canonical generated page-manifest files:

- `results/required_fst_pages.csv`: unique FST page table with indexed paths and gene usage.
- `results/required_fst_pages.txt`: full indexed paths for the required pages.
- `results/required_fst_page_basenames.txt`: basename-only list for server-side matching/export.
- `results/required_fst_pages_by_gene.csv`: audit table mapping genes/transcripts/regions to pages.
- `results/required_fst_selected_transcripts.csv`: transcript-selection audit.
- `results/required_fst_pages_summary.txt`: compact run summary.
- `results/curated_inputs/curated_additional_gene_priorities.csv`: curated auxiliary-prior gene input used by the target resolver.

The current regenerated product contains 290 target genes, 290 selected transcripts, and 252 unique required FST pages: 129 forward pages and 123 reverse pages. It covers chromosomes 1-22 and X. The selected transcript audit currently reports 283 canonical transcript selections and 44 selected transcripts with at least one CDS-overlapping uORF. The server-side tar/export helper now reads `required_fst_pages.csv` by default through `scripts/dominant_cell_state_pack_fst_pages.R`.

## Original Goals

The project started with a practical QC question: can dominant biological states explain large global differences in RNA-seq and Ribo-seq differential expression? Proliferation was the first example, but the project expanded to include ISR/ER stress, interferon/antiviral response, hypoxia, EMT, senescence, OXPHOS, mTOR/translation capacity, and post-viral fatigue/ribosome stress.

The second goal was to make the scores biologically honest for translation. CDS counts alone are misleading when overlapping uORFs, N-terminal extensions, truncations, or internal ORFs contribute coverage inside the canonical CDS. The model therefore moved from count-table CDS values to FST-backed transcript coverage and clean-CDS regions, with predicted and manual translons defining branch points.

The long-term goal is now clearer:

1. Build an RDG atlas for all relevant uORF genes, and later add NTE and NTT effects as explicit branch types.
2. Use all Ribo-seq samples and metadata to learn how conditions change branch usage and CDS output.
3. Generalize from observed condition effects to predictive CDS-buffering levels. For example, estimate how much arsenite stress in LCL-like samples would increase ATF4 protein output even if that exact experiment is not present.
4. Use post-viral fatigue/ribosome stress as the high-value test case: first prove the atlas works on known dominant states such as ISR, then ask whether viral infection produces a broader and possibly novel reprogramming of ouORFs, iORFs, and CDS buffering.

The practical biological output should not only be complicated RDG branch diagrams. The most useful quantity for many biologists will be the downstream effect size on protein production: for example, "this condition is predicted to increase functional ATF4 CDS output by about X percent."

The 2026-06-15 literature and gap section above supersedes the older quick literature checks. The conclusion is unchanged: add candidate states cautiously, prioritize targeted FST expansion for undercovered biology, and keep RDG branch measurements separate from sequence-to-function priors.

## Current Pipeline

The historical nested entry point is:

```r
source("dominant_cell_states/scripts/run_dominant_cell_state_pipeline.R")
run_dominant_cell_state_pipeline("all")
```

From the future standalone analysis repo, start in
`/home/roler/Desktop/forks/dominant_cell_states`, load local RiboCrypt first,
then source the same pipeline script:

```bash
export RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt
env -u LC_ALL R --vanilla -q -e 'source("scripts/run_dominant_cell_state_pipeline.R"); run_dominant_cell_state_pipeline(steps = "merged_replicate_translon_shift", dry_run = TRUE, rdg_plot_mode = "tables")'
```

Use `rdg_plot_mode = "tables"` for most model/table iterations and reserve full
plotting for figure or visual-asset refreshes.

The pipeline currently checks manual translons first, then runs clean-CDS scoring, score auditing, dominant-state activity metadata export, auxiliary dominant-state screening, enrichment, control-aware metadata contrasts, broad metadata modeling, glmnet/latent summaries, known-uORF sheets, uORF-structure rules, RDG figures, branch usage, grouped branch usage, count-aware branch allocation, coupling/rescue analyses, fatigue summaries, hierarchical models, context calibration, a branch-covariate residual prediction test, graph-geometry and mechanistic-motif generalization tests, motif-prior branch calibration, joint branch allocation, study-aware hierarchical joint allocation, empirical-Bayes Dirichlet-multinomial branch allocation, model-agreement compression, leave-one-study-out agreement testing, grouped study-family/protocol/context omission, explicit viral context-interaction testing, partial-pooled context review, atlas cards, pooling/protocol-bias auditing, auxiliary-state/RDG context auditing, allocation-output disagreement auditing, therapeutic-hypothesis prioritization, therapeutic perturbation predictions, RDG consensus compression, RDG validation dossier generation, evidence-gap prioritization, active-review batch generation, active-review link/card generation, active-review feedback ingestion, and validation-label bridging. Frame QC scripts still exist, but the project interpretation currently treats frame/iORF calls as diagnostic rather than as primary claims because p-shift issues are still being repaired upstream.

The RDG atlas browser now separates clean-CDS buffering from branch redistribution. The original heatmap remains a clean-CDS effect view based on `atlas_clean_cds_pct_change`. A branch heatmap uses joint-allocation shift magnitude, a branch-context table exposes the study/condition rows behind the selected atlas card, a study-evidence table exposes the per-study hierarchical and posterior-predictive branch deltas behind model-agreement and LOO labels, a transport-evidence table shows which grouped omissions preserve or break the signal, a context-interaction table shows tissue/cell-line branch effects against their viral complements, and a partial-pooling table shows which context rows survive empirical-Bayes shrinkage versus remain single-context evidence. The atlas cards now include a `partial_pooling_component`, warnings, filter, hover text, and selected-card fields. This matters because a gene can show a strong RDG branch point without a large aggregate clean-CDS effect. PABPC1 tumor/kidney is the current example: the tumor versus WT kidney row has a strict leader-uORF allocation increase and matching clean-CDS allocation decrease, while the broad tumor-cancer clean-CDS atlas cell is near neutral. The study-aware hierarchical joint layer also calls the broad PABPC1 tumor-cancer design neutral, because the tumor/kidney row is diluted by many other tumor contexts. That is the desired behavior: the branch-context table keeps the local signal visible, while the hierarchical aggregate refuses to overstate it.

Manual translons are intentionally part of the model. ORFik/Transcode predictions still miss important overlapping uORFs for several genes, so the manual list is a temporary bridge until the lab's next-generation translon predictor and p-shift denoising model are ready.

## Data Compression Strategy

The project is too large to inspect sample by sample. The useful compression layers are:

- Dominant-state marker scores per sample.
- Auxiliary dominant-state screens that test candidate state families before committing to FST-heavy marker expansion.
- Auxiliary-state/RDG overlap audits that flag atlas rows whose branch effects occur in strong immune, p53, mitochondrial-stress, host-shutoff, or undercovered stress-marker contexts.
- Metadata enrichments and contrasts per biological level.
- Clean-CDS, branch, and relative-usage summaries per gene and condition.
- Matched-control contrasts when WT/control samples exist in the same study/design.
- Grouped branch usage to avoid single-sample sparsity for short uORFs.
- Hierarchical, leave-one-study-out, and grouped-omission models to identify effects that generalize beyond one study or context family.
- Atlas cards and browser queues to turn large tables into a small ranked review list.
- Consensus atlas candidates that combine atlas output, branch redistribution, model agreement, robustness, browser support, context specificity, auxiliary-state warnings, therapeutic relevance, and fragility penalties into one review queue.
- Validation dossiers that turn the consensus queue into concrete browser packets, expected clean-CDS/branch patterns, matched branch-context rows, a manual review template, and figures for the first-pass human review workload.
- Evidence-gap priorities that assign each atlas row one next action and separate acquisition value from counterfactual readiness.
- Active-review batches that turn the evidence-gap queue into a small balanced label-acquisition plan with positive candidates, context audits, measurement reviews, replication follow-up, and negative controls.
- Review-feedback tables that convert completed linked-template rows into label classes, calibration weights, candidate browser-validation notes, and training rows without overwriting curated validation inputs.
- Validation-label bridge tables that combine curated browser labels and active-review feedback, overlay labels on motif priors, and rank the next label gaps.

This is the right direction. The key project risk is not lack of data; it is false precision from sparse branches, wrong isoforms, missing translons, and metadata confounding. The pipeline should continue to prefer conservative gates and review-priority outputs over very long unfiltered candidate tables.

## Pooling, Bias, And Claim Language

The atlas should now use explicit evidence language. The strongest compendium proof is a project-specific matched contrast: same study, same cell or tissue context, same protocol/fraction/inhibitor, same time point where available, and condition versus WT/mock/control with enough merged replicates and branch counts on both sides. This is the only current tier that should be described as a condition-specific RDG effect. It still needs browser inspection and ideally an orthogonal assay, but it preserves the local experimental design.

Same-condition replicate pooling is a measurement step, not a biological claim. It sums biological or technical replicates within the same metadata stratum so short uORFs can pass count gates. It gains count stability and reduces single-run sparsity, but it loses replicate heterogeneity, single-run dominance, and outlier information. The audit language should therefore distinguish "merged-replicate support" from "replicated biological effect".

Study/cell/gene fallback pooling and design-family pooling are recurrence screens. They ask whether an RDG branch shift appears across related contexts such as viral infection, nutrient starvation, ISR/translation stress, tumor context, or cytokine/interferon treatment. They gain power and reveal broad patterns, but they can mix virus type, tissue, inhibitor, footprint length, lab, time point, and sample composition. These rows should be called review-grade or context-modulated unless leave-one-study, grouped-omission, partial-pooling, and browser evidence all agree.

Normalization reduces depth and scale artifacts; it does not remove protocol or state bias. Branch allocation mostly compares counts inside the same gene/transcript denominator, so library size is less central than for clean-CDS output. However, inhibitor and footprint protocols can still bias where reads pile up. LTM and harringtonine-like start-site enrichment can create stronger initiation peaks, CHX can preserve elongating footprints differently, and fractionation/polysome choices can change the sampled ribosome population. Clean-CDS fpkm-like/log2FC normalization similarly cannot by itself distinguish true protein-output buffering from global translation-load changes.

Biological condition families can also skew the whole distribution. Viral infection can produce host shutoff, immune induction, RNA sensing, and tissue-composition shifts, so an infected-versus-control branch delta is not automatically a universal viral rule. Amino-acid or nutrient starvation can increase stalling, ribosome collisions, and GCN2/ISR effects, making branch-specific peaks more likely to reflect a mixture of initiation control and elongation/collision biology. Broad stress or drug families therefore need protocol/context covariates and negative controls, not only pooled effect sizes.

Partial pooling is a shrinkage diagnostic. A pooled-supported row means a context effect survives borrowing information across sibling contexts within the same gene/design/branch/axis. A single-context row can still be biologically real, but it should be named single-context evidence. A raw signal that shrinks away is a fragility warning. The all-merged profile is only a visual background for transcript shape and should not be cited as condition-specific evidence.

This iteration adds a pooling/protocol-bias audit to quantify these distinctions. It reports coverage gains from merged replicates, rows supported only by exact local contrasts versus relaxed or design-family pooling, partial-pooling survival versus shrinkage, inhibitor/protocol signatures such as start-peak enrichment, global frame/misshift risk, and condition-family hazards such as viral distribution shifts or nutrient-starvation stalling. The output should guide which claims are ready for browser validation, which need a project-specific WT/control contrast, and which are broad pooled hypotheses only.

The current audit numbers support the merged-replicate-first policy. Small-uORF support at the 30-read gate rises from 12.9% in single runs to 23.2% in merged n>=2 groups. The allocation minimal gate rises from 38.5% to 53.2%. The grouped branch table contains 82,521 contrast rows: 18,925 are exact strict project-matched rows, while 12,916 are relaxed evaluable pooled rows. Partial pooling separates 11 pooled-supported context rows from 28 single-context rows and 27 raw-supported rows that shrink to fragile. The protocol audit covers 3,857 unique runs and finds 114 LTM/harringtonine-like start-site-enriched runs, so inhibitor/start-peak effects must remain explicit covariates rather than hidden inside design-family labels.

The newest algorithm iteration moves that language into the atlas itself. `dominant_rdg_atlas_cards.csv` now carries `pooling_proof_tier`, `pooling_claim_language`, `pooling_proof_component`, `pooling_hazard_class`, `pooling_hazard_score`, `pooling_hazard_sources`, and compact protocol/context summaries derived from the joint branch-context rows. Exact strict WT/mock/control-like local contrasts score as the strongest proof tier. Pooled recurrence, context-specific rows, and partial-pooling rows remain visible but are labeled as review leads rather than condition-specific claims. The consensus atlas now treats pooling proof as a positive evidence component and pooling hazard as a fragility/penalty component. The atlas browser exposes these fields as visual encodings, card-table columns, selected-card details, and a `Pooling/protocol risk` preset. This makes the first question for each card explicit: is this a local matched-control claim candidate, a pooled recurrence screen, or a potentially protocol-confounded review lead?

The downstream review loop now uses the same language. Validation dossiers carry pooling proof tier and hazard class into the primary manual-review template. Evidence-gap scoring has a separate `protocol_pooling_need` and `protocol_pooling_audit` action instead of hiding protocol risk inside generic context audits. Active-review batches now include B03, a 20-row pooling/protocol confounder audit focused on inhibitor/fraction/start-site enrichment, relaxed context pooling, and local same-protocol controls. This is important because a row can be biologically interesting and still be the wrong kind of claim: exact local proof supports condition-specific wording, while high pooling hazard supports inspection and qualification, not a broad mechanism statement.

## Code Review Findings

The current codebase under `dominant_cell_states/scripts` parses cleanly. The analysis folder now contains 75 R scripts, and the pipeline registry contains 61 steps, with no missing script files and no duplicate step names.

One real bug was found and fixed in `dominant_rdg_atlas_cards.R`. The browser-priority score used context-calibration and context-model components that became `NA` when the corresponding status was missing. Since `setorder()` placed those NA priority rows at the top, the review queue began with low-information heat/osmotic rows instead of the strongest viral candidates. The fix makes missing context components score as zero and guards against non-finite final priority values. After rerunning the atlas-card step, there are zero NA browser-priority scores and the queue now starts with DDIT3, ATF3, PPP1R15A, RRM1, MYC, RPLP0, ATF4, and IFIT1 viral-infection rows.

The candidate-ranking table was also cleaned up. `candidate_rank_v4` was already correct and score-sorted, but the legacy `candidate_rank` column still contained a preliminary rank. `candidate_rank` is now a current-rank alias, while the old value is retained as `preliminary_candidate_rank`.

Grouped branch-usage denominators were checked directly. For every contrast group, summed case feature counts equal the case branch total, summed control feature counts equal the control branch total, and no duplicated features were found within a contrast group. This is important because earlier relative-CDS artifacts came from denominator mistakes involving NTE or non-allocation features. The current grouped-branch implementation appears internally consistent.

A new internal invariant test now checks script parsing, Shiny app parsing, pipeline step wiring, manual translon output, candidate rank ordering, auxiliary dominant-state outputs, auxiliary-state/RDG audit outputs, atlas priority scores, consensus atlas ranking, validation dossier consistency, evidence-gap pooling/protocol scoring, active-review pooling/protocol batches, grouped branch denominators, count-aware branch-allocation outputs, hierarchical count-aware integration, branch-covariate residual-model outputs, graph-geometry generalization outputs, mechanistic motif-score outputs, motif-prior calibration outputs, joint branch-allocation posterior sums, study-aware and posterior-predictive branch outputs, model agreement, leave-one-study-out stability, grouped-omission transportability, explicit viral context interactions, partial-pooled context outputs, pooling/protocol-bias outputs, atlas integration, and partial-pooling atlas/browser integration. This should catch several of the failure modes that previously required manual inspection.

Known code and analysis risks remain:

- The largest grouped-branch adapter CSV is machine-facing and large. It is useful for auditability, but future versions should also save compact RDS/FST summaries.
- Many scripts are robust analysis scripts, not formal package tests. The project needs small regression tests for manual translon import, clean-CDS exclusion, NTE unique-region handling, grouped denominators, and atlas-priority ordering.
- Cache invalidation is gene/file based and useful, but it can hide internal algorithm changes unless force flags or version keys are used carefully.
- Metadata terms are still noisy. Broad labels like `infected`, `WT`, or study-specific condition names can dominate results unless matched controls and hierarchical residual checks are applied.
- The branch-usage model now has independent beta screens, a joint Dirichlet composition approximation, study-aware branch aggregation, and an empirical-Bayes posterior-predictive screen. The next statistical upgrade is a true hierarchical beta-binomial or Dirichlet-multinomial model with study, gene, design, and branch random effects.
- The auxiliary-state/RDG audit reconstructs matched case/control sample sets from metadata. It is useful for warning and prioritization, but it is not a causal correction because study, tissue, condition, and auxiliary state can still be entangled.

## Main Results

The atlas currently contains 1,680 gene-design rows over 120 genes and 14 design families. The broader atlas review table has 1,532 rows. Browser-supported rows are still rare, which is expected because only a few candidates have been manually checked in RiboCrypt. The current browser-supported candidates are DDIT3 and ATF3 in viral infection, with DDIT3 now the only browser-supported row that also sits in the top model-agreement consensus tier.

The consensus atlas compresses those 1,680 rows into a ranked 1,293-row review queue. The consensus score combines clean-CDS atlas confidence, branch redistribution, model agreement, leave-one-study/grouped-omission robustness, browser support, effect size, context specificity, auxiliary-state warnings, therapeutic relevance, sparse/on-off evidence, novelty, pooling proof, pooling/protocol hazard, and fragility penalties. It creates one browser-validated consensus row, 48 replicated robust consensus rows, 82 multi-model branch consensus rows, 29 context-specific review rows, 23 buffering/compensation review rows, 103 therapeutic-translation review rows, 401 general atlas review rows, 387 low-support background rows, and 606 fragile/confounded review rows. The top consensus rows are DDIT3 viral infection, PPP1R15A viral infection, LDHA ISR/stress, ATF4 viral infection, TFAM ISR/stress, EIF4G1 ISR/stress, IFIH1 viral infection, DDIT4 ISR/stress, LMNB1 ISR/stress, and CDK1 ISR/stress. This ranking is useful for deciding what to inspect next; it should not be read as a new statistical p-value.

The validation-dossier layer turns the consensus atlas into a practical browser-review set. It keeps all 1,680 consensus rows in a candidate table, selects 184 primary review packets across 72 genes, and attaches 380 top matched branch contexts. Of those branch contexts, 372 pass the 30-count gate on both case and control sides, and all 380 have merged replicates on both sides. Among primary review packets, 145 have exact pooling proof and 53 carry a pooling/protocol hazard flag, so the template now distinguishes strong local claim candidates from protocol-confounded review leads. The manual review template is the preferred human-facing file: it records the packet ID, transcript, design, first branch context, expected clean-CDS and branch pattern, count/replicate gates, pooling proof/hazard fields, warnings, and blank columns for validation status. The first packets remain centered on DDIT3 viral infection, PPP1R15A viral infection, ATF4 viral infection, IFIH1 viral infection, VEGFA viral infection, LDHA ISR/stress, and TFAM ISR/stress, now with a broader second tier that includes PELO/HK2 contexts. This is not new evidence; it is the bridge between the atlas and browser/manual validation.

The auxiliary dominant-state screen tested eight candidate state families using only genes already present in the clean-CDS expression matrix, so it did not reload FST coverage. The screened families were proteostasis/heat shock, autophagy/lysosome, NRF2/oxidative-ferroptosis defense, DNA damage/p53 checkpoint, antigen presentation/immune composition, secretory/ER proteostasis, mitochondrial UPR/mitonuclear stress, and host-shutoff/antiviral translation repression. Three modules currently look pipeline-ready as screeners: DNA damage/p53 checkpoint, antigen presentation/immune composition, and mitochondrial UPR/mitonuclear stress. Antigen presentation is the cleanest new dominant-state candidate because it is not just another ISR/OXPHOS alias and it strongly identifies LCL, blood, hematopoietic, lymphoma, and antigen-presentation contexts. DNA damage/p53 separates checkpoint-like conditions from simple low proliferation, but the top metadata hits are study/perturbation-heavy and need matched-control interpretation. Mitochondrial UPR is marker-rich but partially overlaps ISR/OXPHOS and should be treated as a bridge module, not a new mechanistic claim. Host-shutoff/antiviral repression is fully scoreable, but it is close to the existing post-viral fatigue/ribosome-stress score, so it should remain exploratory until we decide whether the post-viral module should be split.

The auxiliary-state/RDG audit connects these per-sample auxiliary scores back to matched RDG branch contexts. It reconstructed auxiliary case/control deltas for 86 of 159 RDG strata, generated 13,104 gene/design/state summaries, annotated 840 of 1,680 atlas rows, and wrote a compact 145,352-row branch audit subset from a 1,190,592-row denominator. After the literature-gap expansion, antigen presentation is now both the raw and readiness-weighted top auxiliary overlap state. This is useful but also a warning: antigen/immune composition can dominate tumor and infection contexts and should be treated as a context/confounder layer before it becomes a mechanistic RDG claim. Interpretation: new dominant states are useful, but their first role should be context/confounder annotation and browser-review triage, not automatic addition to the core dominant-state winner table.

Grouped branch usage is the most important recent improvement. It aggregates replicates and condition groups before testing branch shifts, which is much more realistic than expecting individual samples to quantify short uORFs. The grouped branch tables contain 1,885 groups, 961 groups with at least two samples, and 299 review rows after filtering. This matches the user's concern that single samples are usually insufficient for small uORFs.

The count-aware branch-allocation layer reuses the grouped replicate counts and computes beta-posterior case/control usage deltas for each branch feature. It produced 46,110 posterior contrast rows, 20,312 evaluable rows, 10,851 strict rows, 2,320 feature summaries, 656 gene-design summaries, and 299 review rows. It calls 28 replicated count-aware branch shifts overall, but only one replicated viral-infection row: ATF4 clean-CDS up. That is an important conservative result. It says the viral branch signal is not absent, but after count gates and count-aware shrinkage, the replicated viral evidence is currently narrow.

The count-aware layer is integrated into atlas cards, browser review priority, and the hierarchical CDS-buffering summaries. After the expanded/manual run, 95 atlas rows have count-aware branch-allocation support and 36 are replicated. Viral infection has five count-aware-supported rows and one replicated row. The viral review plot marks browser-supported rows with circles and count-aware branch-allocation support with purple diamonds, making it easier to separate "manually visible in the browser" from "survives stricter branch-count modeling."

The hierarchical integration is deliberately conservative. Count-aware rows only boost hierarchical mechanistic support when the support class is a real branch-shift call. `count_aware_evaluable_no_strong_shift` rows are retained in the CSVs but do not increase hierarchical priority. This keeps the hierarchy from treating "we had enough counts to look" as evidence for a branch mechanism.

The branch-covariate prediction test asks a stricter question: if branch-allocation evidence is converted into covariates and the held-out study is excluded from those covariates, does it improve held-out CDS-buffering prediction? Overall, the metadata/context residual model improves hierarchical RMSE from 0.480 to 0.4776 log2FC, and adding branch covariates changes this only to 0.4775 log2FC. In viral infection, RMSE improves from 0.570 to 0.5633 with metadata/context and to 0.5628 with branch covariates. This is a real but very small signal. DDIT3, PPP1R15A, and ATF4 are the viral genes where branch covariates most reduce residual error, but all gene-level effects remain below the practical support threshold. Interpretation: branch features are directionally useful in the expected ISR genes, but the current representation is too weak/noisy for strong prediction.

The branch coverage QC confirms that single-sample support is limited. Only about 13% of single-sample small-uORF features reach 30 counts, while merged groups with at least two samples reach about 23%. For allocation features overall, the minimal gate improves from about 39% in single samples to about 53% in merged replicate groups. This supports a merged-replicate-first strategy.

The hierarchical CDS-buffering model is useful but not yet publication-grade. It has about 0.29 log2FC leave-one-out MAE and 0.48 log2FC RMSE, but uncalibrated intervals under-cover. Calibrated intervals are more reliable but wide. This means the model can prioritize candidates and estimate broad effect sizes, but precise condition-specific forecasts still need better features and validation.

Adding count-aware branch evidence did not materially change the predictive calibration metrics, which is expected because it currently changes prioritization and mechanistic labels rather than the posterior CDS effect itself. The leave-one-study-out MAE remains about 0.291 log2FC, RMSE about 0.480 log2FC, and calibrated interval coverage about 0.95.

Adding count-aware branch covariates as a residual model also does not materially change predictive accuracy. This is useful negative evidence: the next improvement should not simply add more branch-summary columns to the same model. It should use better graph features, better translon calls, and a joint allocation model that can distinguish total ribosome-load changes from true branch redistribution.

The graph-geometry generalization model asked a different question: if gene identity is withheld, do uORF/RDG architecture features generalize across genes? In the expanded run, the context-only leave-one-gene-out model reached 0.401 RMSE overall and 0.508 RMSE in viral infection. Adding many explicit geometry covariates worsened this to 0.412 overall and 0.525 in viral infection. This is a useful failure. It means raw feature counts, overlap classes, spacing bins, Kozak summaries, and simple interactions are too crude or too correlated to act as a general rule model. The next geometry model should not be a larger unconstrained feature dump; it should use constrained mechanistic scores, validated motif labels, or dimensionality reduction over branch-allocation evidence.

The constrained mechanistic motif model compressed architecture into interpretable scores: inhibitory-overlap burden, reinitiation-distance opportunity, dense-uORF collision burden, strong-start leakage burden, downstream-rescue opportunity, and derived overlap/rescue interactions. This is more interpretable than the raw graph-geometry dump but is still not a strong predictor. Overall context/no-gene RMSE was 0.401 and motif main effects were 0.401, essentially neutral/slightly worse. In viral infection, motif main effects worsened RMSE from 0.508 to 0.509, and motif-design interactions worsened it further to 0.513. This updates the earlier interpretation: hand-built motifs are useful priors and review annotations, but they are not yet a predictive layer.

The motif-prior calibration step compares these bounded motif scores against observed count-aware branch-allocation contrasts. It now spans 37,597 branch rows over 970 gene-designs. The global Spearman correlation between motif prior strength and observed absolute branch delta is weak but positive, about 0.149; for clean-CDS rows it is stronger, about 0.340. Viral infection is also weakly positive, about 0.115. This means the hand-built motifs are not random, but they are far from calibrated. The next use of motifs should be label-driven: browser-validated positives and negatives should teach the motif weights rather than relying on fixed hand weights.

The joint branch-allocation model aggregates grouped branch counts into one composition per contrast: leader uORF, overlapping uORF, clean CDS, and residual other branch mass. The primary posterior uses a uniform Dirichlet prior; motif scores are used as sensitivity/alignment annotations rather than hard-coded into the primary call. This produced 22,896 allocation rows, 16,548 evaluable rows, and 11,265 strict rows in the expanded run. Viral infection has 2,016 allocation rows, 1,534 evaluable rows, and 466 strict rows. After atlas support gating, 323 atlas rows have joint branch-allocation support, 251 are strict, 47 are viral-supported, and 21 are viral strict. The motif-strength versus L1 allocation-shift correlation is now moderate for a hand-built prior: about 0.396 overall and 0.370 in viral infection.

The top viral joint-allocation queue starts with ATF4, PPP1R15A, DDIT3, MYC, IFIH1, SNAI2, NDUFV1, IGFBP5, MRPS12, DDIT4, VEGFA, TBP, SESN2, TOMM20, FEN1, CDK1, TFAM, TRIB3, and LMNB1. This is more biologically interesting than the independent beta screen because it asks whether the whole RDG branch composition changes, not only whether one feature differs from the rest. It also confirms the direction of several user-reviewed viral candidates: ATF4 and PPP1R15A show viral clean-CDS allocation increases, while DDIT3 shows a viral leader-uORF allocation decrease and remains motif-discordant. That discordance is useful: it marks DDIT3 as a candidate where the current motif prior is probably missing the real decision rule.

The study-aware hierarchical joint-allocation layer is a conservative aggregation on top of the row-level joint model. It first collapses contexts within each study, then estimates random-effects branch deltas across studies for each gene, transcript, design, and branch class. This produced 3,880 branch effects over 970 gene-designs and 144 genes. Overall, 102 branches are replicated shifts and 15 are replicated-review shifts. Viral infection has 556 branch effects over 139 gene-designs, with 11 replicated branch shifts. After atlas-card integration, 145 rows have hierarchical joint-allocation support, 61 are replicated, and viral infection has seven supported rows, all replicated.

The replicated viral hierarchical joint-allocation list is now compact and biologically coherent: PPP1R15A, DDIT3, MMP2, ATF4, VEGFA, and IFIH1. PPP1R15A, DDIT3, MMP2, ATF4, and IFIH1 show clean-CDS allocation increases in infected contexts, while VEGFA shows a leader-uORF allocation increase. Most top viral contexts come from infected versus mock/control designs such as `PRJNA401882-homo_sapiens | RD | muscle | infected | vs mock`. This is one of the strongest current outputs because it survives study-aware aggregation instead of relying on a single browser-visible row. At the same time, the top all-design plot is still dominated by single-study UV/heat/stress rows, so those should remain review leads rather than replicated claims.

The empirical-Bayes Dirichlet-multinomial layer is the newest count-aware branch-allocation screen. It uses the row control counts plus weak priors learned from same-gene, same-design, or global background rows, then asks whether the case branch count is surprising under a posterior-predictive beta-binomial marginal. It produced 1,891 branch effects in the expanded run, with 82 replicated-review branches overall. Viral infection has 280 branch effects, with 8 replicated shifts and 6 replicated-review shifts. After atlas integration, 205 rows have Dirichlet-multinomial support, 87 are replicated, and viral infection has 10 supported rows, all replicated.

The viral posterior-predictive queue is short: DDIT3 clean-CDS up, ATF4 clean-CDS up, VEGFA leader-uORF up, MRPS12 leader-uORF up, IFIH1 clean-CDS up, and SNAI2 overlapping-uORF up. This partly agrees with the study-aware hierarchical list but is not identical. PPP1R15A remains a strong branch-effect row in the older hierarchical joint layer, but the posterior-predictive layer demotes it because the branch direction is not stable enough across the current contexts. This is a support-threshold disagreement, not an opposite-direction call: candidates supported by both layers are stronger review targets, while threshold disagreements mark places where the prior, study mix, or manual browser evidence should be inspected.

The model-agreement layer turns the overlapping queues into one compact evidence table over 1,680 gene-designs. It assigns 1 browser-validated consensus row, 42 replicated branch-consensus rows, 88 additional multi-model consensus rows, 23 disagreement-review rows, and 1,100 low-support rows. The hierarchical-joint and posterior-predictive layers still agree in direction for supported shared branches; there are currently no strong hierarchical-versus-Dirichlet-multinomial direction disagreements. This does not constitute independent validation because the estimators reuse the same grouped branch counts, but it shows that the branch-effect direction is stable to the two current aggregation/shrinkage strategies.

The six viral model-agreement consensus gene-designs are DDIT3, IFIH1, VEGFA, ATF4, BBC3, and PPP1R15A. DDIT3 is the only browser-validated consensus row. IFIH1, VEGFA, ATF4, and BBC3 are replicated branch-consensus rows. PPP1R15A is multi-model consensus rather than replicated posterior-predictive consensus: its posterior-predictive clean-CDS delta remains positive, but it does not pass the posterior-predictive support threshold. ATF3 remains important because it is manually validated, but it currently has only clean-CDS model support and therefore stays below the formal viral consensus tier.

The leave-one-study-out agreement layer asks whether each supported branch retains direction and thresholded support after omitting one study. It preserves per-study effects from both branch models, refits every evaluable gene/design/branch/model combination, and writes a compact per-study evidence table for the atlas browser. Of 131 consensus gene-designs, 49 are genuinely evaluable by study omission. Twenty retain both direction and thresholded support after all relevant omissions, 29 retain direction but sometimes fall below support thresholds, and 82 are non-evaluable by this criterion. None of the study-evaluable consensus rows reverses direction under study omission. This is stronger than estimator agreement alone, but it remains an influence analysis rather than independent validation. The LOO support gate uses within-branch p-values as a robustness proxy instead of recomputing global BH correction in every omission.

All six viral model-agreement consensus rows are study-evaluable. Four are support-robust under leave-one-study-out and two are direction-robust but support-sensitive. This makes the current viral consensus set materially stronger than a queue driven by one influential study, while preserving the distinction between robustness and biological validation.

The grouped-omission transportability layer asks a harder question: does the effect survive removal of an entire related family of studies rather than one accession? It groups studies by publication, author/lab, coarse protocol, inhibitor, tissue, cell line, and footprint-length family; a test is only evaluable when at least two studies are removed and at least two independent studies remain. Of 131 consensus gene-designs, 25 are evaluable under at least one grouped omission. Eight are support-robust across all evaluable axes, 10 retain direction but sometimes lose thresholded support, one is partially evaluable, six are fragile, and 106 are not evaluable because the current compendium lacks enough independent context families. This is intentionally much stricter than leave-one-study-out. Median axis coverage among evaluable rows is only three of seven axes, so the stability score explicitly includes axis coverage and prevents robustness on a narrow testable subset from scoring like broad transportability.

The grouped-omission transportability layer is stricter. Among the six viral consensus rows, two are support-robust, two are direction-robust but support-sensitive, one is not evaluable across grouped omissions, and one is fragile. DDIT3 and PPP1R15A remain the strongest cross-context viral candidates. ATF4 and VEGFA retain useful directional evidence but are more support-sensitive. IFIH1 remains context-dependent in the present compendium, especially around lung-like viral contexts. The correct interpretation is not that all infection has one universal RDG effect; it is that viral RDG rewiring is recurrent but context-modulated.

The explicit viral context-interaction layer now tests that interpretation directly instead of inferring it only from grouped omission. It reuses the joint branch-allocation class rows, collapses rows within each study and tissue/cell-line context, and compares each context-specific viral branch delta to the remaining non-matching viral contexts. No FST pages are reloaded. The expanded run produces 10,496 study-context rows, 7,488 branch-context effects, 819 review rows, and 139 viral gene-design summaries; 133 gene-designs are context-evaluable. The prioritized gene-level classes contain 27 context direction reversals, 2 context-specific supported shifts, and 12 magnitude-modulated shifts.

IFIH1 is the cleanest result from this layer. In lung viral-infection contexts, IFIH1 leader-uORF allocation decreases by about 0.092 while the non-lung complement trends up by about 0.028, giving an interaction delta of about -0.120 (`q` about 0.006). The clean-CDS branch moves in the opposite direction: lung is about +0.091 and non-lung is about -0.013, with interaction delta about +0.104 (`q` about 0.009). This supports the earlier grouped-omission conclusion: IFIH1 is currently a lung-focused viral RDG branch hypothesis. It is not yet a universal viral-infection rule.

Among the viral consensus genes, the context layer adds nuance rather than replacing the consensus layer. ATF4 is context-direction-reversal at the gene level, with a lung overlapping-uORF effect as the top interaction. PPP1R15A and DDIT3 are magnitude-modulated rather than direction-reversal calls in their top context rows. VEGFA is magnitude-modulated, with the current top atlas row driven by A549 leader-uORF context. DDIT3 and PPP1R15A therefore remain the strongest cross-context anchors by grouped omission, while IFIH1 remains the clearest tissue-dependent viral hypothesis.

The partial-pooled context layer is the first shrinkage pass over those explicit context rows. It estimates empirical-Bayes shrinkage across context values within each gene, design, branch, and context axis, then separates true pooled context support from single-context-only evidence. It now has 1,756 evaluable context rows over 133 gene-designs, with 11 pooled-supported rows, 40 pooled-review rows, 28 single-context raw signal rows, 27 raw signals shrunk to fragile, and 5 pooled-emergent rows. After atlas integration, the current atlas has 5 partial-pooling supported rows, 11 partial-pooling review rows, 13 single-context rows, and 5 fragile-after-pooling rows. The top current pooled context is a TFE3 lung leader-uORF viral-infection effect, which shows that the context layer is broadening beyond the original IFIH1 example.

The allocation-versus-output disagreement audit now separates these nine rows into interpretable review classes. All nine have the same core signature: clean-CDS output increases while clean-CDS branch allocation decreases. SESN2 nutrient starvation is the strongest buffering/total-load compensation candidate: clean-CDS output is up about 94%, clean-CDS allocation is down about 0.094, and upstream branch allocation is up about 0.083. ATF3 nutrient starvation is the second-best lead: output is up about 95%, clean-CDS allocation is down about 0.046, and overlapping-uORF allocation is up about 0.051, but heterogeneity risk is higher. SNAI2 ISR/stress remains a special mechanistic review target because output is only weakly up while allocation strongly shifts away from clean CDS toward upstream branches. TOMM20 cytokine/interferon is a possible compensation row. RRM1 and PCNA are now best treated as aggregation/heterogeneity-risk rows rather than primary biological candidates. DDIT4, ZEB2, and IFIH1 are lower-magnitude review rows.

The random forest/metadata models are useful as variance-compression tools, not as mechanistic explanations. The previous ISR value around 0.76 means metadata/features could explain roughly 76% of held-out variation in the ISR dominant-state score under that model setup. It does not mean ISR is 76% caused by one metadata term.

## Biological Interpretation

ATF4 remains the canonical positive-control RDG. It demonstrates the model type we want: uORF architecture controls whether ribosomes terminate upstream, reinitiate, bypass inhibitory branches, or reach the CDS. The project should use ATF4 to validate the logic, but not treat every candidate as ATF4-like.

DDIT3 and ATF3 are currently the strongest viral/post-viral review leads because the user manually inspected them online in RiboCrypt and saw convincing translational changes in infected versus all/uninfected contexts. These are not yet formal mechanistic discoveries, but they are exactly the kind of candidates the atlas should surface.

ATF4 viral infection, ATF4 nutrient/stress-like contexts, and PPP1R15A nutrient/stress-like contexts appear in formal grouped-branch, count-aware, joint-allocation, or matched-control outputs and should be prioritized for browser validation. They are especially useful because they connect known ISR biology to the branch-allocation framework.

The joint-allocation viral queue makes the post-viral story more concrete. ATF4, PPP1R15A, and DDIT3 appear together as whole-RDG redistribution candidates in infected-versus-control contexts. The expanded model-agreement layer identifies DDIT3, IFIH1, VEGFA, ATF4, BBC3, and PPP1R15A as the current viral consensus set. Study omission confirms that the viral consensus is not driven by one accession, while grouped omission separates robust anchors from context-modulated rows. DDIT3 and PPP1R15A remain the strongest cross-context anchors, ATF4 and VEGFA are directionally robust but support-sensitive, IFIH1 remains context-dependent, and BBC3 is a newer replicated branch-consensus entry that needs browser review before biological interpretation. ATF3 remains important because the user has manually seen convincing viral translation changes, but it still needs stronger branch-model support to match DDIT3 in the formal consensus stack.

The manual PELO and HK2 additions expand the biology rather than replacing the viral ISR story. PELO now ranks 7th in the next-model table and has a viral context-specific row with about +10% clean-CDS effect, but it is not a robust viral consensus row. That makes it a ribosome-rescue/RQC hypothesis to inspect, not a claim. HK2 ranks 11th and has a stronger ISR/stress consensus row with about -14% clean-CDS effect; its viral row is weak and therapeutically annotated rather than consensus-supported. These two genes are useful because they test whether the atlas can move beyond classical ATF4/DDIT3 ISR biology into ribosome rescue and glycolytic stress.

The earlier apparent broad frame-shift/bumpiness signal in infected samples should be quarantined for now. User inspection and QC suggested a substantial p-shifting component in some infected groups. This does not invalidate the RDG work, but it means iORF/off-frame claims should wait for the lab's p-shift denoising model.

OXPHOS was an unexpectedly dominant score winner in many samples. This is probably a mixture of real mitochondrial/energy-state biology, marker-set breadth, housekeeping expression, library composition, and metadata/study effects. The OXPHOS module should be expanded and audited before interpreting it biologically.

The auxiliary screen argues against blindly adding many more dominant states to the core table. Antigen presentation is worth adding as a confounder-aware screen because immune composition can dominate global expression and may otherwise masquerade as viral or inflammatory translation regulation. DNA damage/p53 is worth keeping because it can separate checkpoint arrest from ordinary low proliferation. Mitochondrial UPR is biologically relevant for post-viral fatigue and energy stress, but because it reuses ISR-up/OXPHOS-down structure, it should be interpreted as a decomposition of the broad fatigue/OXPHOS signal rather than an independent discovery. The new auxiliary-state/RDG audit sharpens this: DDIT3/ATF4/TRIB3 viral rows often sit in mitochondrial-stress-shifted contexts, while tumor/cancer rows often carry antigen-presentation shifts. Proteostasis, autophagy/lysosome, and NRF2/ferroptosis are biologically important, but the current FST gene set undersamples them; the right next step is targeted FST-page expansion before making them core dominant states or interpreting their raw overlap counts.

## Therapeutic Translation

The atlas now includes a small therapeutic-hypothesis prioritization layer. This is deliberately conservative: it ranks intervention concepts for ex vivo perturbation, biomarker stratification, and literature-guided follow-up, not patient treatment. The score combines RDG atlas evidence, viral specificity, model agreement, transport/context robustness, browser support, practical readiness, mechanistic alignment, and a safety/risk penalty.

The current ranking is:

1. ISR/eIF2B-ATF4 attenuation, score about 0.827. This is the best mechanistic preclinical hypothesis because DDIT3, ATF4, PPP1R15A, DDIT4, ATF3, HSPA5, and TRIB3 carry strong atlas evidence. The right experiment is not a medical trial; it is a perturbation assay asking whether ISR attenuation changes RDG allocation and functional CDS buffering in the predicted direction.
2. JAK/STAT IFN attenuation, score about 0.758. This is the best clinical-adjacent biomarker hypothesis, because it connects viral/IFN genes such as IFIH1, IRF7, STAT1, CCL2, MX1, IFIT1, and OAS2 to an intervention class that already exists pharmacologically. The claim should remain narrow: persistent post-acute IFN-high samples might be stratified by RDG/marker response. Immune suppression during uncontrolled acute infection is an obvious safety caveat.
3. IFIH1/MDA5 axis modulation, score about 0.586. This is a strong target-discovery and stratification idea because IFIH1 has explicit lung-dependent viral RDG context evidence, but direct modulation of viral sensing has high safety risk and is not trial-ready.
4. HIF/hypoxia and mitochondrial/OXPHOS hypotheses score in the exploratory range. They may become useful disease-state modifiers, but they are currently too broad for direct therapeutic interpretation.

The most useful near-term therapeutic output is a ranked perturbation plan: test whether ISR or IFN-axis perturbations shift RDG branch allocation and CDS-buffering levels in browser-validated genes. A clinically useful endpoint would be a predicted change in functional protein output, not only a branch diagram. For example, the future atlas should report whether a condition or perturbation is predicted to increase or decrease ATF4/DDIT3/IFIH1 CDS output by a calibrated percentage with uncertainty.

The new perturbation-prediction layer makes the first simple version of that endpoint. It does not estimate true drug efficacy. Instead, it reports a half-normalization target: if a perturbation reversed half of the disease-like atlas effect, what CDS-buffering change would be expected? This gives a concrete assay target while keeping the assumption explicit.

Current perturbation-prediction priorities are:

1. ISR/eIF2B-ATF4 attenuation, prediction score about 0.703. Three genes are quantitative-review-ready: DDIT3, ATF4, and PPP1R15A. The half-normalization targets are approximately DDIT3 -46% CDS, ATF4 -42% CDS, and PPP1R15A -63% CDS relative to the disease-like increase. ATF3, DDIT4, and HERPUD1 are directional follow-ups.
2. JAK/STAT IFN attenuation, prediction score about 0.631. IFIH1 is quantitative-review-ready with an approximate -20% half-normalization CDS target. IRF7, CCL2, and IFIT1 are directional follow-ups.
3. IFIH1/MDA5 axis modulation, prediction score about 0.509. This remains a stratification or target-discovery hypothesis. IFIH1 and IRF7 have directional predictions, but direct intervention on viral sensing remains high-risk.

These values are best read as experimental design targets: if an ISR or IFN-axis perturbation does not move the nominated RDG branches and CDS output toward these half-normalization directions, then the therapeutic hypothesis is probably weaker than the atlas ranking suggests.

## Post-Viral Fatigue and Infection

The post-viral fatigue/ribosome-stress module is conceptually important but currently immature. No samples are dominated by a fatigue score as a primary winner, which is not surprising because the module is new and the compendium is not a patient cohort. The useful near-term question is narrower: do viral or immune contexts repeatedly redirect ribosomes on ISR, antiviral, mitochondrial, and proteostasis genes in ways that alter CDS output?

This is still a good test case for the final paper. The proposed story is:

1. Introduce the RDG atlas and CDS-buffering model using known dominant-state effects such as ISR and ATF4.
2. Apply the same framework to viral infection and post-viral-like contexts.
3. Test whether viral infection causes recurrent branch reallocation into ouORFs or iORFs, reducing normal CDS output or changing isoform output.
4. Use the atlas to nominate specific genes and branch points for experimental validation.

The ME/CFS genetic literature does not make this project irrelevant. Common risk variants usually change probability modestly and do not imply that translation-state mechanisms require a specific patient genotype to be discoverable. Genetic predisposition, viral trigger, immune state, and ribosome programming can all interact. The project should avoid claiming ME/CFS causality, but it can sensibly search for post-viral translation programs that could contribute to fatigue biology.

## Methodological Lessons

Relative CDS usage must only compare features that belong in the same branch-allocation denominator. NTEs are branch points, but an NTE covering the full CDS should be quantified from its unique portion for expression support and should not be allowed to dilute clean-CDS relative usage. iORFs are diagnostic branch candidates, but their frame usage is hard to interpret until p-shifts are fixed.

The number, type, and spacing of uORFs matter. Clean-CDS relative usage tends to fall as predicted/manual uORF burden increases, especially for overlapping uORFs and dense uORF clusters near the CDS TIS. However, false-negative uORF predictions can invert this interpretation. SNAI2, MMP2, LMNB1, ATF3, and MKI67 examples showed why manual translons are currently necessary.

Positive correlation between ouORF and CDS levels is not automatically contradictory. It can reflect total transcript/ribosome-load changes, study-level expression differences, incomplete normalization, reinitiation after a permissive uORF, downstream alternative TIS usage, or technical branch sharing. The meaningful test is not raw ouORF-CDS correlation; it is whether branch allocation changes after controlling for gene expression, study/design, and matched controls.

Single-sample branch decisions are usually too sparse for small uORFs. The project should prefer merged replicate groups, condition-level aggregation, and count-aware uncertainty models. Individual samples are still useful for dominant-state scores and broad CDS estimates, but not for most small-branch claims.

Pooling must be named by what it preserves. Exact matched WT/control contrasts preserve the local experimental design and are the best current proof tier. Merged replicate groups preserve the condition stratum but compress replicate heterogeneity. Relaxed study/cell/gene pools and design-family pools preserve recurrence but lose causal specificity. Partial pooling preserves context-aware shrinkage information but can demote real single-context biology. All-merged profiles preserve transcript shape only. Future branch models should carry these tiers directly as model structure: study, protocol/inhibitor, tissue/cell line, time point, condition family, gene, branch, and replicate effects should be estimated together rather than audited after the fact.

Merged replicates are now the default screen for small-branch biology, but they are not automatically proof. They improve uORF count depth and reduce single-run spike artifacts, which is essential for short leader and overlapping uORFs. They also hide replicate heterogeneity and can mix protocol effects inside a condition stratum. This matters for inhibitor and treatment interpretation: LTM-like start-site enrichment can inflate apparent start/uORF peaks, viral infection can move the whole ribosome-load distribution rather than a single branch, and amino-acid starvation can increase stalling-like signal that looks like branch redistribution. Therefore the current algorithm asks two separate questions: does clean-CDS allocation move relative to uORF allocation in a same-condition merged contrast, and what evidence is still missing before that pattern can be interpreted?

Joint allocation is a better review primitive than independent feature-vs-other tests because it separates branch redistribution from total ribosome-load changes. It is not yet a final likelihood model: branch counts are aggregated from existing grouped feature counts, overlapping branch counts are guarded by denominator checks, and the current chi-square-like vector test is approximate. The useful biological output is the ranked candidate list and delta matrix, not a final p-value claim. The study-aware hierarchical joint layer improves review specificity by asking whether the branch effect survives after collapsing within studies and estimating across-study effects. It should be treated as a conservative ranking layer, not as a final Bayesian RDG model.

The empirical-Bayes Dirichlet-multinomial layer is a useful intermediate step toward a real RDG likelihood. It improves on the uniform joint allocation screen by using control counts and weak background priors, then testing posterior-predictive branch outliers. It still uses beta-binomial branch marginals rather than fitting a full multilevel Dirichlet-multinomial over all branches, and the high viral heterogeneity shows that study/context structure remains a central problem. For now, agreement between the study-aware hierarchical layer, the posterior-predictive layer, and browser inspection should be treated as the strongest evidence tier.

Model agreement is estimator stability, not independent replication. The hierarchical-joint and posterior-predictive layers share the same grouped branch-count foundation, so their very high branch-delta correlations mainly show that shrinkage and aggregation choices preserve direction. Independent evidence still requires held-out studies, matched biological controls, browser inspection, orthogonal assays, or experimental perturbation. The useful new distinction is between branch-model disagreement and allocation-versus-output disagreement. No supported branch has opposite direction between the two branch models, while nine gene-designs show opposite clean-CDS allocation and net clean-CDS output. The latter is a plausible signature of buffering or total-load compensation and should be investigated mechanistically.

Leave-one-study-out stability is a stricter influence test, but it is still not independent validation because every refit reuses the remaining compendium and the same branch-count construction. A support-robust call means no single study is required to retain the current direction and thresholded support. A direction-robust/support-fragile call means the biological direction is stable but the evidence margin is thin. A single-study row is labeled non-evaluable, not fragile.

Grouped omission is closer to a transportability stress test because it removes related studies together. It can reveal context dependence that accession-level omission hides, as seen for IFIH1 after removing all lung studies. It is still not causal identification or independent validation: metadata groups can be noisy, publication/author families are often non-evaluable, and a grouped omission can simultaneously remove biology, protocol, and study-family effects. The next model should explicitly estimate gene-by-context-by-condition effects rather than forcing every viral signal into one universal effect.

The context-interaction and partial-pooling layers are the first simple version of an explicit context model. They are useful because they turn "this breaks when lung is omitted" into "lung branch effect differs from the non-lung complement", and then ask whether that context effect survives shrinkage across sibling contexts. They are still not causal models. A context value with two or three studies can still combine tissue, cell line, protocol, virus, time point, and publication effects, and a single-context-only row can be strong without being partial-pooled. Therefore context interactions should drive browser review and model design, not final claims. The next statistical step is a true multilevel branch model where context effects are estimated jointly with infection, tissue, cell line, study family, gene, and branch terms.

Atlas review must now be explicit about which layer is being interpreted. Clean-CDS heatmap enrichment means a net CDS-buffering/output effect. Branch heatmap enrichment means ribosomes are being reallocated among RDG branches. The two can disagree for good biological reasons, especially when total ribosome load, reinitiation, or alternate branch choices compensate at the CDS level.

## Relationship To Lab Methods

Two parallel lab efforts are directly relevant.

First, the lab is developing improved Ribo-seq simulation software extending `coverageSim` and a full-scale translon predictor intended to work better for overlapping uORFs. Current models miss many overlapping uORFs, which is why the manual list exists.

Second, the lab is developing a denoising model for p-shifted Ribo-seq. The idea is to use global CDS frame ratios and protocol-specific smear behavior to infer what the signal should look like after removing frame contamination. This will be essential before making strong iORF or frame-reprogramming claims.

These tools should eventually feed the RDG atlas directly: better translon calls define better graph nodes, and denoised p-shifted signal gives better branch weights.

## Useful Outputs

Key files to inspect:

- `human_dominant_cell_states_clean_cds.csv`: per-sample dominant-state scores from clean CDS.
- `dominant_state_score_audit.csv`: marker availability and score sanity checks.
- `manual_translons/manual_translons_manifest.csv`: selected manual uORF/uoORF/NTE features, including the current PELO and HK2 additions.
- `dominant_state_auxiliary_modules/dominant_state_auxiliary_module_summary.csv`: novelty, marker availability, metadata enrichment, and recommendation table for candidate new dominant states.
- `dominant_state_auxiliary_modules/dominant_state_auxiliary_module_scores.csv`: per-run auxiliary module scores without changing the main dominant-state table.
- `dominant_state_auxiliary_modules/dominant_state_auxiliary_metadata_enrichment.csv`: metadata levels enriched for each auxiliary state.
- `dominant_state_auxiliary_modules/figures/auxiliary_state_novelty_and_availability.png`: decision plot showing novelty versus current marker availability.
- `dominant_rdg_auxiliary_state_audit/dominant_rdg_auxiliary_state_atlas_annotations.csv`: atlas rows annotated with top auxiliary-state context overlap, readiness-weighted priority, and warning text.
- `dominant_rdg_auxiliary_state_audit/dominant_rdg_auxiliary_state_gene_design_summary.csv`: gene/design/state summaries linking auxiliary-state deltas to matched RDG branch effects.
- `dominant_rdg_auxiliary_state_audit/dominant_rdg_auxiliary_state_summary_by_state.csv`: compact count of supported RDG rows occurring in each auxiliary-state context.
- `dominant_rdg_auxiliary_state_audit/figures/auxiliary_state_rdg_overlap_heatmap.png`: heatmap of top RDG gene/design rows versus auxiliary-state shifts.
- `dominant_rdg_auxiliary_state_audit/figures/auxiliary_state_vs_atlas_effect.png`: scatter showing auxiliary-state shift size versus atlas clean-CDS effect.
- `dominant_rdg_atlas_cards.csv`: main atlas table.
- `dominant_rdg_atlas_browser_validation_queue.csv`: prioritized browser-review queue.
- `dominant_rdg_consensus/dominant_rdg_consensus_candidates.csv`: one-row-per-atlas-card consensus evidence table with score components, evidence tags, classes, warnings, and review questions.
- `dominant_rdg_consensus/dominant_rdg_consensus_review_queue.csv`: consensus-ranked review queue for browser inspection.
- `dominant_rdg_consensus/dominant_rdg_consensus_evidence_matrix.csv`: long-form component matrix behind the consensus ranking.
- `dominant_rdg_consensus/figures/dominant_rdg_consensus_top_candidates.png`: top consensus candidate plot.
- `dominant_rdg_consensus/figures/dominant_rdg_consensus_evidence_heatmap.png`: component heatmap showing why each top candidate ranked highly.
- `dominant_rdg_validation_dossiers/dominant_rdg_validation_dossier_manual_review_template.csv`: preferred manual browser-review sheet with packet IDs, expected patterns, top context, count gates, warnings, and blank validation columns.
- `dominant_rdg_validation_dossiers/dominant_rdg_validation_dossier_candidates.csv`: full consensus-to-validation candidate table with primary-review selection and priority tiers.
- `dominant_rdg_validation_dossiers/dominant_rdg_validation_dossier_branch_contexts.csv`: top matched branch contexts for primary review packets.
- `dominant_rdg_validation_dossiers/dominant_rdg_validation_dossier_report.md`: compact human-readable review protocol and top packet summary.
- `dominant_rdg_validation_dossiers/figures/rdg_validation_dossier_priority.png`: prioritized browser packet plot.
- `dominant_rdg_validation_dossiers/figures/rdg_validation_dossier_evidence_heatmap.png`: evidence-component heatmap for the selected review packets.
- `dominant_rdg_validation_dossiers/figures/rdg_validation_dossier_branch_contexts.png`: branch-context matrix for the selected review packets.
- `dominant_rdg_evidence_gaps/dominant_rdg_evidence_gap_priorities.csv`: ranked acquisition queue with next action, readiness, gap score, biological value, and review reason.
- `dominant_rdg_evidence_gaps/dominant_rdg_evidence_gap_gene_summary.csv`: gene-level action and readiness summary.
- `dominant_rdg_evidence_gaps/dominant_rdg_evidence_gap_negative_controls.csv`: candidate falsification rows with low clean-CDS effect but adequate measurement/replication support.
- `dominant_rdg_evidence_gaps/figures/evidence_gap_priority_map.png`: acquisition-priority versus counterfactual-readiness map.
- `dominant_rdg_evidence_gaps/figures/evidence_gap_gene_action_heatmap.png`: top-gene action matrix.
- `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_rows.csv`: all joint-allocation rows with strongest uORF contrast, merged-replicate support, claim-readiness class, and missing-evidence flags.
- `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_review_queue.csv`: inspectable merged-replicate clean-CDS/uORF opposition queue.
- `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv`: collapsed candidate queue with one row per gene/transcript/design/shift direction, recurrence class, pooling benefit/loss class, output-allocation context, and next review action.
- `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_pooling_summary.csv`: design/pooling-scope summary showing where reviewable clean-CDS/uORF shifts come from and which rows are limited by protocol, measurement, or coverage assets.
- `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_gene_design_summary.csv`: gene/transcript/design summary of merged-translon shift evidence.
- `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_report.md`: compact interpretation and top-row summary.
- `dominant_rdg_merged_replicate_translon_shift/figures/merged_replicate_translon_shift_clean_uorf_scatter.png`: clean-CDS delta versus strongest uORF delta.
- `dominant_rdg_review_batches/dominant_rdg_review_batches.csv`: active-review batch table with one row per review task.
- `dominant_rdg_review_batches/dominant_rdg_review_batch_template.csv`: preferred manual label sheet for the next browser-validation pass.
- `dominant_rdg_review_batches/dominant_rdg_review_batch_summary.csv`: compact summary of batch size, goal, and top gene/design rows.
- `dominant_rdg_review_batches/figures/review_batch_priority_readiness.png`: review-batch effect-space map.
- `dominant_rdg_review_batches/figures/review_batch_gene_matrix.png`: gene coverage across review batches.
- `dominant_rdg_review_links/dominant_rdg_review_batch_template_with_links.csv`: active-review manual label sheet with direct Observatory URLs.
- `dominant_rdg_review_links/dominant_rdg_review_batch_links.csv`: link audit table with link mode, case/control run counts, and URL notes.
- `dominant_rdg_review_links/dominant_rdg_review_cards.html`: compact clickable browser-review card page.
- `dominant_rdg_review_feedback/dominant_rdg_review_label_drafts.csv`: optional atlas-browser draft label file used as a writeback buffer before feedback ingestion.
- `dominant_rdg_review_feedback/dominant_rdg_review_feedback_labels.csv`: normalized manual feedback labels; currently header-only until linked-template rows or app draft labels are filled.
- `dominant_rdg_review_feedback/dominant_rdg_review_feedback_training_set.csv`: calibration-ready labeled rows with model covariates; currently header-only until manual or draft labels exist.
- `dominant_rdg_review_feedback/dominant_rdg_review_feedback_browser_validation_notes_candidate.csv`: candidate rows compatible with atlas browser-validation notes, for manual inspection before merging.
- `dominant_rdg_review_feedback/dominant_rdg_review_feedback_summary.csv`: feedback-ingest counts by label and calibration use.
- `dominant_rdg_review_feedback/dominant_rdg_review_feedback_batch_summary.csv`: per-batch review and label counts.
- `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_set.csv`: unified validation labels from curated browser notes and active-review feedback.
- `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_gene_design.csv`: one-row-per-labeled-gene summary with motif-prior context.
- `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_motif_overlay.csv`: all motif-prior gene-design rows annotated with validation-label status.
- `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_gap_queue.csv`: active-review rows reranked by label-balance need and motif/context priority.
- `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_balance.csv`: current label counts versus target seed counts.
- `dominant_rdg_qualified_inspection/qualified_inspection_manifest.csv`: one row per static inspection packet with source row, matched context, count gates, readiness flags, and PNG paths.
- `dominant_rdg_qualified_inspection/qualified_review_priority.csv`: ranked inspection queue used by the browser Review Workbench.
- `dominant_rdg_qualified_inspection/qualified_profile_metrics.csv`: coverage breadth and peak-concentration metrics for each selected profile group.
- `dominant_rdg_qualified_inspection/qualified_feature_contrasts.csv`: exact case/control feature-level branch contrasts used by the Inspection Assets page.
- `dominant_rdg_qualified_inspection/qualified_replicate_summary.csv`: per-condition replicate support and single-run-dominance diagnostics.
- `dominant_rdg_qualified_inspection/qualified_sequence_features.csv`: transcript/RDG feature coordinates, start codons, frames, Kozak context, and mRNA-sequence availability.
- `dominant_rdg_sdrive_readiness/sdrive_experiment_status.csv`: mounted S-drive, ORFik all-merged library, reference FASTA/index, count-table, bigWig, and indexed FST-page status checks.
- `dominant_rdg_sdrive_readiness/sdrive_fst_page_readiness.csv`: one row per required indexed coverage page with current local existence and file-size metadata.
- `dominant_rdg_sdrive_readiness/sdrive_transcript_readiness.csv`: selected qualified-inspection rows annotated with backing-data readiness, required FST-page coverage, sequence/reference availability, and review-support class.
- `dominant_rdg_sdrive_readiness/sdrive_readiness_summary.md`: human-readable audit summary separating data availability from biological claim quality.
- `dominant_rdg_grouped_branch_usage_review.csv`: grouped branch-shift review rows.
- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_review.csv`: count-aware beta-posterior branch-allocation review rows.
- `dominant_rdg_branch_allocation/dominant_rdg_branch_allocation_gene_design_summary.csv`: count-aware gene-design branch summary.
- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariate_model_summary.csv`: leakage-controlled test of whether branch covariates improve held-out CDS-buffering prediction.
- `dominant_rdg_hierarchical_branch_covariates/dominant_rdg_hierarchical_branch_covariate_model_gene_design.csv`: gene-design level branch-covariate residual-model summary.
- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_summary.csv`: leave-one-gene-out test of whether RDG architecture generalizes beyond gene identity.
- `dominant_rdg_graph_geometry_model/dominant_rdg_graph_geometry_model_gene_design.csv`: gene-design level graph-geometry generalization summary.
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_model_summary.csv`: leave-one-gene-out test of constrained RDG motif scores.
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_scores.csv`: bounded mechanistic scores per gene/transcript.
- `dominant_rdg_mechanistic_motif_model/dominant_rdg_mechanistic_motif_model_viral_gene_design.csv`: viral gene-design summary for motif-model effects.
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_calibration_summary.csv`: motif-prior versus count-aware branch-evidence calibration.
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_review_queue.csv`: gene-designs where motif priors and branch evidence agree, disagree, or expose missing evidence.
- `dominant_rdg_motif_prior_calibration/dominant_rdg_motif_prior_viral_review.csv`: viral-infection motif-prior review queue.
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_gene_design_summary.csv`: joint Dirichlet branch-composition summary by gene/design.
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_review.csv`: joint allocation review queue.
- `dominant_rdg_joint_branch_allocation/dominant_rdg_joint_branch_allocation_viral_review.csv`: viral-infection joint allocation queue.
- `dominant_rdg_joint_branch_allocation/figures/joint_branch_allocation_delta_heatmap.png`: compact matrix of leader-uORF, overlapping-uORF, and clean-CDS allocation deltas.
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_branch_effects.csv`: study-aware random-effects branch deltas.
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv`: compact gene-design summary of the study-aware joint layer.
- `dominant_rdg_hierarchical_joint_allocation/dominant_rdg_hierarchical_joint_allocation_viral_review.csv`: viral-infection replicated/study-aware branch review queue.
- `dominant_rdg_hierarchical_joint_allocation/figures/hierarchical_joint_allocation_delta_heatmap.png`: study-aware branch delta matrix.
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_gene_design_summary.csv`: empirical-Bayes posterior-predictive branch-allocation summary by gene/design.
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_review.csv`: posterior-predictive branch review queue.
- `dominant_rdg_dirichlet_multinomial/dominant_rdg_dirichlet_multinomial_viral_review.csv`: viral-infection posterior-predictive branch review queue.
- `dominant_rdg_dirichlet_multinomial/figures/dirichlet_multinomial_viral_review.png`: viral posterior-predictive priority plot, with replicated rows emphasized and neutral rows faded.
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_gene_design.csv`: compact consensus/disagreement tier per gene/design.
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_branch_comparison.csv`: direct branch-level comparison of study-aware hierarchical and posterior-predictive deltas.
- `dominant_rdg_model_agreement/dominant_rdg_model_agreement_viral_review.csv`: viral model-agreement validation queue.
- `dominant_rdg_model_agreement/figures/model_agreement_viral_matrix.png`: compact viral evidence matrix across clean-CDS, joint, hierarchical, posterior-predictive, browser, and motif evidence.
- `dominant_rdg_model_agreement_loo/dominant_rdg_model_agreement_loo_gene_design.csv`: agreement tiers annotated with study-omission robustness and non-evaluability.
- `dominant_rdg_model_agreement_loo/dominant_rdg_model_agreement_loo_branch_summary.csv`: branch/model-level direction and support retention after study omission.
- `dominant_rdg_model_agreement_loo/dominant_rdg_model_agreement_study_evidence.csv`: compact per-study evidence behind atlas cards.
- `dominant_rdg_model_agreement_loo/figures/model_agreement_loo_viral_stability.png`: viral consensus stability under study omission.
- `dominant_rdg_model_agreement_transport/dominant_rdg_model_agreement_transport_gene_design.csv`: grouped-omission transportability class and stability per gene/design.
- `dominant_rdg_model_agreement_transport/dominant_rdg_model_agreement_transport_gene_axis_summary.csv`: compact context/protocol/publication/lab omission evidence used by the atlas browser.
- `dominant_rdg_model_agreement_transport/dominant_rdg_model_agreement_transport_viral_review.csv`: viral consensus transportability review table.
- `dominant_rdg_model_agreement_transport/figures/model_agreement_transport_viral_axes.png`: viral grouped-omission support/direction matrix.
- `dominant_rdg_context_interactions/dominant_rdg_context_interaction_branch_context_effects.csv`: tissue/cell-line viral branch effects compared with non-matching complements.
- `dominant_rdg_context_interactions/dominant_rdg_context_interaction_gene_design_summary.csv`: compact gene-design summary of context-direction reversals, context-specific shifts, and magnitude-modulated effects.
- `dominant_rdg_context_interactions/figures/context_interaction_ifih1_tissue.png`: IFIH1 lung versus non-lung branch-effect figure.
- `dominant_rdg_context_interactions/figures/context_interaction_viral_consensus_tissue_heatmap.png`: context-interaction matrix for the viral consensus genes.
- `dominant_rdg_context_partial_pooling/dominant_rdg_context_partial_pooling_rows.csv`: shrinkage-aware context rows separating pooled support, single-context-only evidence, fragile raw calls, and pooled-emergent review rows.
- `dominant_rdg_context_partial_pooling/dominant_rdg_context_partial_pooling_viral_review.csv`: compact viral gene-design review table from the partial-pooling layer.
- `dominant_rdg_context_partial_pooling/figures/context_partial_pooling_viral_top_effects.png`: top viral context effects after shrinkage-aware classification.
- `dominant_rdg_context_partial_pooling/figures/context_partial_pooling_ifih1_lung.png`: IFIH1 lung raw versus shrinkage-aware estimate, explicitly labeled as single-context evidence.
- `dominant_rdg_pooling_bias_audit/pooling_bias_report.md`: compact report defining the evidence tiers and current pooling/protocol-bias results.
- `dominant_rdg_pooling_bias_audit/pooling_claim_language.csv`: approved claim vocabulary for exact matched contrasts, replicate pooling, design-family pooling, single-context evidence, partial pooling, and all-merged backgrounds.
- `dominant_rdg_pooling_bias_audit/normalization_limits.csv`: what each normalized quantity controls and which protocol or biology biases remain.
- `dominant_rdg_pooling_bias_audit/pooling_contrast_tier_summary.csv`: exact strict project-matched rows versus exact evaluable, relaxed pooled, and low-count branch contrasts.
- `dominant_rdg_pooling_bias_audit/design_family_pooling_bias_summary.csv`: design-family confounding hazard scores from branch-context diversity in protocol, tissue, cell line, and condition.
- `dominant_rdg_pooling_bias_audit/protocol_bias_run_summary.csv`: inhibitor/protocol categories, TIS peak strength, global frame QC, and condition-family summaries across unique runs.
- `dominant_rdg_pooling_bias_audit/partial_pooling_tradeoff_summary.csv`: context rows that survive pooling, remain single-context only, emerge after pooling, or shrink to fragile.
- `dominant_rdg_allocation_output_disagreements/dominant_rdg_allocation_output_disagreement_review.csv`: compact audit of rows where clean-CDS output and clean-CDS branch allocation point in opposite directions.
- `dominant_rdg_allocation_output_disagreements/dominant_rdg_allocation_output_disagreement_contexts.csv`: branch-context rows behind the disagreement audit.
- `dominant_rdg_allocation_output_disagreements/figures/allocation_output_disagreement_scatter.png`: output-versus-allocation disagreement map.
- `dominant_rdg_allocation_output_disagreements/figures/allocation_output_disagreement_branch_heatmap.png`: clean-CDS, leader-uORF, overlapping-uORF, and other branch deltas for disagreement rows.
- `dominant_rdg_therapeutic_hypotheses/dominant_rdg_therapeutic_hypothesis_scores.csv`: preclinical/biomarker intervention-hypothesis ranking from atlas evidence.
- `dominant_rdg_therapeutic_hypotheses/dominant_rdg_therapeutic_hypothesis_review.csv`: compact review table with examples, caveats, top genes, and score components.
- `dominant_rdg_therapeutic_hypotheses/figures/therapeutic_hypothesis_scores.png`: therapeutic-hypothesis ranking plot.
- `dominant_rdg_therapeutic_hypotheses/figures/therapeutic_hypothesis_evidence_matrix.png`: component matrix behind the therapeutic ranking.
- `dominant_rdg_therapeutic_perturbation_predictions/dominant_rdg_therapeutic_perturbation_gene_predictions.csv`: gene-level disease-like effects and full/half-normalization perturbation targets.
- `dominant_rdg_therapeutic_perturbation_predictions/dominant_rdg_therapeutic_perturbation_hypothesis_summary.csv`: perturbation-hypothesis summary with prediction priority scores.
- `dominant_rdg_therapeutic_perturbation_predictions/dominant_rdg_therapeutic_perturbation_review.csv`: compact review table for top gene-level perturbation targets.
- `dominant_rdg_therapeutic_perturbation_predictions/figures/therapeutic_perturbation_top_gene_effects.png`: top half-normalization CDS-buffering targets.
- `dominant_rdg_therapeutic_perturbation_predictions/figures/therapeutic_perturbation_hypothesis_summary.png`: perturbation-prediction priority plot.
- `dominant_rdg_atlas/dominant_rdg_atlas_branch_contexts.csv`: compact row-level branch-context table used by the atlas browser.
- `R/run-app.R`: exported `run_dcs_app()` launcher for the packaged Shiny atlas browser.
- `inst/shiny/rdg_atlas_browser/app.R`: Shiny atlas browser with clean-CDS heatmap, branch heatmap, effect-space plot, selected-card details, row-level branch contexts, per-study agreement evidence, grouped-omission transport evidence, context-interaction evidence, and a partial-pooling tab/filter for shrinkage-supported versus single-context rows.
- `branch_coverage_summary_metrics.csv`: evidence that grouped replicates are needed for short branches.
- `dominant_next_model_candidate_rankings_v4.csv`: candidate ranking table.
- `md/postviral_fatigue_translation_note.md`: separate literature summary on viral translation and fatigue.
- `results/curated_inputs/known_uorf_effects_with_citations.csv`: curated known-uORF effects.

Figures and HTML reports are distributed under:

- `dominant_rdg_atlas/`
- `dominant_rdg_figures/`
- `dominant_rdg_grouped_branch_usage/`
- `dominant_rdg_joint_branch_allocation/`
- `dominant_rdg_hierarchical_joint_allocation/`
- `dominant_rdg_context_partial_pooling/`
- `dominant_rdg_pooling_bias_audit/`
- `dominant_rdg_consensus/`
- `dominant_rdg_validation_dossiers/`
- `dominant_rdg_evidence_gaps/`
- `dominant_rdg_review_batches/`
- `dominant_rdg_review_links/`
- `dominant_rdg_auxiliary_state_audit/`
- `dominant_rdg_allocation_output_disagreements/`
- `dominant_rdg_therapeutic_hypotheses/`
- `dominant_rdg_therapeutic_perturbation_predictions/`
- `dominant_state_auxiliary_modules/`
- `dominant_uorf_structure_rules/`
- `relative_usage_matrices/`
- `postviral_fatigue/`
- `inst/shiny/rdg_atlas_browser/`

## 2026-06-16 Literature-Gap FST Expansion Run

The 90 newly downloaded FST pages completed the literature-gap expansion. The readiness table now has 118 required literature-gap page rows and all of them point to existing local FST files. The clean-CDS layer measured 290 marker genes, including the new viral/fatigue, ISR-kinase, RNA-sensing, host-shutoff, collision/RQC, DNA-damage, antigen-presentation, complement/endothelial, and DecodeME-prior genes.

One implementation bug was fixed during the run: a few translons mapped outside the displayed leader+CDS matrix and produced invalid zero/negative row indices during coverage summation. `dominant_uorf_regulation_inference.R` now sanitizes mapped positions to finite 1-based rows inside the coverage matrix before computing feature counts and measured bases.

The expanded run first completed with cached upstream steps and downstream RDG/model/atlas reruns. After that, the manual translon list was updated again and the manual-scope pipeline was rerun with 47 successful steps and 2 expected cache-skipped FST-readiness checks. No pipeline step failed. R printed many warnings, so warnings should still be audited before treating this as a frozen release, but the registered outputs were regenerated successfully.

Key run outputs:

- `human_dominant_cell_states_clean_cds.csv`: expanded clean-CDS dominant-state scores.
- `dominant_uorf_feature_expression.csv`: per-run branch feature counts for all modeled genes.
- `dominant_rdg_outputs/`: RDG figures for the expanded modeled panel.
- `dominant_next_model/dominant_next_model_candidate_rankings_v4.csv`: current candidate ranking.
- `dominant_rdg_atlas/dominant_rdg_atlas_viral_infection_review.csv`: focused viral-infection atlas review table.
- `dominant_rdg_consensus/dominant_rdg_consensus_review_queue.csv`: cross-model consensus review queue.
- `dominant_rdg_validation_dossiers/`: browser-validation dossier tables and manual-review template.
- `dominant_rdg_evidence_gaps/`: action-priority queue, gene/action summaries, negative controls, and evidence-gap figures.
- `dominant_rdg_review_batches/`: balanced active-review batches and manual label template.
- `dominant_rdg_review_links/`: direct Observatory links and compact review cards for active-review rows.
- `dominant_rdg_review_feedback/`: normalized review labels, calibration/training rows, candidate browser-validation notes, and feedback summaries.
- `dominant_rdg_validation_label_bridge/`: unified validation labels, motif overlay, label-balance table, and label-gap queue.
- `dominant_rdg_therapeutic_hypotheses/` and `dominant_rdg_therapeutic_perturbation_predictions/`: preclinical perturbation hypothesis layer.

The result is not a final biological claim yet, but it is a stronger atlas than before. The top next-model candidates are now IFIT1, DDIT4, EGLN3, MKI67, ATF3, CDK1, PELO, LMNB1, RPLP0, MMP2, HK2, and HSPA5. These are candidate-discovery rankings, not the same as the strongest consensus biology. The strongest consensus biology is still centered on viral ISR and antiviral sensing:

- `DDIT3 | viral_infection`: current top consensus row and browser-supported anchor, clean-CDS up about 87%.
- `PPP1R15A | viral_infection`: replicated and transport-support robust, clean-CDS up about 128%.
- `ATF4 | viral_infection`: replicated robust viral branch candidate, clean-CDS up about 81%, but transport support is direction-robust rather than support-robust.
- `IFIH1 | viral_infection`: replicated branch candidate with context dependence, especially lung-like contexts, clean-CDS up about 38%.
- `BBC3 | viral_infection`: replicated branch-consensus row from the expanded model-agreement layer, needing browser review.
- `ATF3 | viral_infection`: atlas-ready clean-CDS lead and manually interesting browser case, but currently weaker in cross-model consensus than DDIT3/PPP1R15A/ATF4/IFIH1.
- `MYC`, `VEGFA`, `DDIT4`, `MRPS12`, `RRM1`, `G3BP1`, `MAVS`, `PELO`, and `HK2` are useful review candidates, but should be treated as mechanistic or context-specific until browser validation is stronger.

## 2026-06-19 Creative Iteration: Evidence Gaps

The missing layer was not another FST-heavy model. The atlas already has many effect estimators; what it lacked was an acquisition policy. `scripts/dominant_rdg_evidence_gap_prioritization.R` now reads the atlas, consensus, validation dossier, feature-confidence, robustness, context, auxiliary-state, and therapeutic outputs and converts each gene/design row into one primary next action.

The new product is deliberately pragmatic. It computes an evidence-gap score, a biological-value score, an acquisition-priority score, and a counterfactual-readiness score. A high acquisition-priority row is not automatically a stronger biological claim; it means the row is valuable and missing a specific evidence layer. A high counterfactual-readiness row is closer to being useful for prediction, because it already has stronger replication, model agreement, measurement quality, context handling, and validation support.

The first run supports the current project direction. ATF4 and PPP1R15A viral infection are high-priority browser-validation rows and near-counterfactual candidates. DDIT3 viral infection is highly ready but classified as context/metadata audit because it already has browser support and now needs careful context interpretation rather than another generic validation label. ATF3 viral infection is flagged as perturbation follow-up rather than as the strongest consensus row, which matches the current interpretation: it is biologically interesting and browser-supported enough to care about, but still needs stronger model robustness. The queue also makes the negative lesson explicit: hundreds of rows are not ready for mechanistic interpretation because their primary missing layer is replication/transport robustness, not a new scoring formula.

## 2026-06-23 Active-Review Batches

The next missing layer was a practical review design. A 1,680-row acquisition queue is useful for scripts but still too large for manual RiboCrypt inspection. `scripts/dominant_rdg_active_review_batches.R` now turns the evidence-gap queue into six small review batches with a shared manual-review template.

The batches are:

1. Browser validation for viral/stress anchors.
2. Counterfactual-readiness and perturbation follow-up.
3. Context and metadata confounder audit.
4. Translon and measurement review.
5. Replication and transport follow-up.
6. Negative controls and falsification.

This is an active-learning design. It deliberately mixes positive candidates, borderline predictive candidates, measurement-risk rows, confounded rows, replication rows, and negative controls. That is more useful than validating only the top positive-looking rows, because the next model needs labels for positives, weak positives, negatives, confounded cases, and not-reviewable rows. The key output is `dominant_rdg_review_batches/dominant_rdg_review_batch_template.csv`. Fill the blank review columns there and keep the batch/order columns unchanged; later scripts can use those labels to recalibrate motif priors, context filters, and counterfactual-readiness scores.

## 2026-06-23 Review Links

The active-review batches are now operational rather than just tabular. `scripts/dominant_rdg_review_batch_links.R` adds RiboCrypt Observatory URLs to every review row and writes compact HTML/Markdown review cards. When the metadata CSV is available, the script reconstructs the best case/control run selections from the branch-context metadata. If exact context fields fail to recover both sides, it falls back to study+condition matching and marks the row as `study_condition`.

The preferred manual file is now `dominant_rdg_review_links/dominant_rdg_review_batch_template_with_links.csv`. It contains the same blank review columns as the batch template plus a direct Observatory URL, link mode, case/control run counts, and link notes. This is the right file to use for browser validation because it preserves the active-learning batch design while removing most manual subset reconstruction.

## 2026-06-24 Review Feedback Ingest

The review loop now has a readback layer. `scripts/dominant_rdg_review_feedback_ingest.R` reads the linked active-review template, overlays draft labels saved from `Review -> Review Labels` in the atlas browser, detects rows where manual review columns have been filled, normalizes the labels, and writes feedback products under `dominant_rdg_review_feedback/`. It deliberately does not overwrite `dominant_rdg_atlas/browser_validation_notes.csv`; instead it writes a candidate browser-validation-notes CSV that should be inspected before curated rows are merged into the atlas input.

Draft labels saved from the browser now carry the Workbench acquisition context: target label, label-gap class, label-acquisition score, protocol/pooling need, pooling proof tier, pooling hazard class/score, and pooling-hazard sources. The feedback ingest preserves those fields in normalized label and training outputs, so later calibration rows can be traced back to the reason they were selected rather than only to the manual label assigned after inspection.

The first real run found zero completed manual-review rows, which is expected because the linked template is still blank and no production draft label file has been saved. It still wrote stable outputs: header-only label, training, and browser-note candidate files; a 15-row feedback summary; a six-row batch summary; and a short Markdown report. Synthetic filled-template and draft-overlay checks confirmed the non-empty path for positive, negative-control-pass, and absent-browser-signal labels.

## 2026-06-24 Validation-Label Bridge

The next layer makes the review labels usable before a full manual pass exists. `scripts/dominant_rdg_validation_label_bridge.R` combines the curated atlas browser-validation notes with active-review feedback labels, summarizes labels by gene/transcript/design, overlays label status onto motif-prior calibration rows, and writes a label-gap queue that accounts for which label classes are missing.

The current seed set has two curated browser-supported positive labels: `DDIT3 | viral_infection` and `ATF3 | viral_infection`. There are no active-review feedback labels yet. The bridge therefore reports a deliberately imbalanced label set: two positive labels, zero weak positives, zero negatives, zero confounded rows, zero measurement-problem rows, and zero negative-control labels. This is useful because the next review pass should not only add more positive anchors; it should also collect weak/negative/confounded/measurement and negative-control labels.

The first gap-queue rerank contains 129 active-review rows, of which 126 are still unlabeled after matching the curated browser notes. The top rows are currently ATF4 and PPP1R15A viral-infection counterfactual/browser checks, followed by IFIH1/VEGFA and early negative-control rows such as BECN1, TFE3, and PELO. The queue now carries motif-prior review class and motif strength for all review rows where motif evidence exists.

The viral/fatigue layer is now broad enough to act as the main test case. It produced 55 strict matched viral/IFN rows and four matched design cards with strict support, but zero rows passed the strict interval-supported atlas-ready gate in the focused fatigue summary. The broader atlas, which combines more evidence layers, has six viral model-agreement consensus rows, four viral LOO support-robust consensus rows, and two viral grouped-omission transport-support-robust rows. This distinction matters: the signal exists, but the strictest uncertainty gates still prevent overclaiming.

The therapeutic layer should be interpreted conservatively. It ranks `ISR/eIF2B-ATF4 attenuation` highest, with `DDIT3`, `ATF4`, and `PPP1R15A` as quantitative review-ready viral-infection predictions. The best current statement is a preclinical perturbation hypothesis: if chronic ISR-like viral translation rewiring is causal, partial normalization should reduce the disease-like clean-CDS output shifts in these genes. It is not yet a medical-trial candidate. `JAK/STAT IFN attenuation` is more clinical-adjacent as a biomarker/stratification hypothesis, but immune suppression is unsafe to infer without disease stage, infection status, and patient-level data.

## 2026-07-03 Atlas Browser Review Workbench

The atlas browser now has four pages aimed at faster qualified inspection. `Review -> Label Dashboard` is the label-acquisition control surface: it shows the validation-label balance, pending app draft labels, current bridge labels, feedback-ingest outputs, candidate browser-note rows, and a filtered label-gap queue. Selecting a label-gap row focuses the selected card and exact active-review task, so the reviewer can jump directly to `Review Labels`, `Inspection Assets`, or the Workbench. `Start Label Session` turns the current target-label filter into a focused labeling pass by opening the first open label-gap task, setting the Workbench to the same target, and landing on `Review Labels`.

`Review -> Review Workbench` combines four review queues into one selectable table: the 154-row active review batch, the 154-row validation-label gap queue, the 1,680-row evidence-gap acquisition table, and the qualified-inspection priority table. Label-gap rows now recover their active-review order, so selecting either a batch row or a label-gap row updates the selected atlas card and the review-label task. The Workbench also has target-label and gene filters, displays review destination/gene/transcript/design before queue metadata, and exposes target deficit, label-acquisition score, `protocol_pooling_need`, pooling proof, and pooling hazard fields. The left sidebar gene filter applies when Workbench scope is `Current filters`; the Workbench-local gene selector provides an explicit start point, and the scope note states whether Workbench genes, sidebar genes, selected-card scope, or all queued rows are active. Workbench-local gates now restrict the current table to active label tasks, static inspection packets, or protocol/pooling audit rows. The dashboard above the table reports destination counts, protocol audits, pooling-risk rows, and exact-proof rows for the current scope, while the selected-row card shows the destination, transcript, queue rank, active task, inspection source, label target, pooling proof/hazard, review reason, and Observatory link. `Use Sidebar Genes`, `Use Selected Gene`, and selected-row `Gene Session` convert the current visual context into a Workbench gene session by filling the Workbench gene selector and switching scope to all queued rows. `Open Selected` routes the clicked row to the strongest available surface: active label tasks open `Review Labels`, inspection-only rows open `Inspection Assets`, and rows without either open `Selected Card`. Rows with static inspection packets update the qualified-inspection source row. `Start Session` opens the first open active-review task under the current Workbench queue, target, gene, and scope filters. This makes the browser the primary triage surface for both label-balance acquisition and protocol/pooling audits.

`Review -> Review Labels` is the new lightweight writeback surface for completed inspection decisions. It mirrors the manual review columns from the linked CSV, previews the normalized label class, and saves draft rows to `dominant_rdg_review_feedback/dominant_rdg_review_label_drafts.csv`. It now has a task action row that opens Observatory, Inspection Assets, Selected Card, the Workbench, or Files directly from the label page, plus quick-fill buttons for the common positive, weak-positive, negative, confounded, measurement-problem, not-reviewable, and negative-control patterns. It writes the target-label and pooling/protocol context from the Workbench into each draft row. `Save & Next` turns a filtered Workbench into a review session by saving the current draft and advancing to the next open active-review task in the current queue, target, and scope while skipping tasks that already have draft or ingested feedback labels. The linked CSV remains valid for spreadsheet-style editing, but the app draft file is faster for Codex/human review because the selected card, review task, inspection packet, and acquisition reason stay synchronized.

`Evidence -> Inspection Assets` now exposes the selected qualified-inspection packet directly in the app. It renders the coverage profile, feature-allocation plot, RDG-flow annotation, and replicate-diagnostics plot, plus compact tables for coverage metrics, exact case/control feature contrasts, replicate summaries, and transcript/RDG sequence features. This matters for the current claim gate: Codex or a human reviewer should be able to distinguish a broad branch shift from a single-position spike, see whether the top feature is carried by multiple runs, check whether mRNA sequence and RDG annotation were available, and verify whether translon coordinates match the plotted branch interpretation before opening an interactive Observatory session.

`Evidence -> Data Readiness` adds the mounted-data audit needed for faster Codex/human inspection. With the S drive mounted, the current audit finds the ORFik all-merged library, reference FASTA/index, count tables, bigWigs, the all-sample FST index, and all 252 required indexed coverage pages. It also marks all 129 qualified static-inspection source rows as `sdrive_ready`, while only 29 are `coverage_reviewable` by the existing visual-support class. That split is important: data readiness means the row can be reopened and inspected against real backing data; it does not mean the static coverage supports a biological claim. The page therefore prevents wasted review time on missing files while keeping the claim gate conservative.

The new inspection workflow does not replace RiboCrypt. It decides which rows are worth interactive inspection and which rows should be marked not-claimable, confounded, measurement-risk, or negative-control candidates. The next review pass should use the workbench to acquire balanced labels, then rerun the review-feedback and validation-label bridge steps so motif calibration and context filters learn from positives, weak positives, negatives, confounded rows, and not-reviewable rows.

## 2026-07-09 Merged-Replicate Translon-Shift Gate

The algorithmic focus is now sharper: find cases where a condition-level merge changes clean-CDS allocation relative to other translons, especially leader and overlapping uORFs. `scripts/dominant_rdg_merged_replicate_translon_shift.R` reads the joint branch-allocation rows, selects the strongest uORF branch per contrast, scores clean-CDS/uORF opposition, and merges atlas, evidence-gap, S-drive readiness, and qualified-inspection context.

The first run found 22,896 joint rows, 15,840 exact same-condition merged rows, 8,938 clean-CDS/uORF opposition rows, and 6,434 exact merged opposition rows. The inspectable review queue is deliberately much smaller: 95 rows after requiring merged replicates, nontrivial clean-CDS/uORF contrast, and some available S-drive/static inspection support. That queue contains 44 viral-infection rows, 16 nutrient-starvation rows, and 9 tumor/cancer-context rows.

The next iteration collapses the 95 context rows into 60 candidate-direction rows so review starts from biology rather than duplicated contexts. Eight collapsed candidates have multi-study exact recurrence, one has multi-condition exact recurrence, and 18 have opposite-direction rows elsewhere in the same gene/transcript/design. This is the right split: recurrence is the main benefit of pooling, while direction conflicts and protocol/start-site bias are the main loss of interpretability.

No row is currently promoted to `claim_candidate_exact_merged_review_ready`. This is the correct conservative outcome. The top exact rows still lack browser labels and model-robustness support, while 39 collapsed candidates are marked for protocol/start-site/pooling audit before interpretation. Examples include `PPP1R15A`, `TFE3`, `IFIH1`, `FBXL4`, `DDIT3`, and `ATF4` viral-infection candidates. These are good review targets, not final claims.

The practical interpretation rule is now: a large merged-replicate branch shift only starts the review. The candidate becomes biologically interpretable only after the missing-evidence flags are cleared or explicitly accepted. The critical fields are `recurrence_class`, `pooling_benefit_class`, `pooling_loss_class`, `candidate_interpretation_class`, `next_review_action`, missing browser label, protocol/pooling audit needed, weak transport or leave-one-study-out support, missing static coverage, single-run dominance in the static asset, missing mRNA sequence, and missing RDG annotation track. This turns the model from "top score wins" into "what evidence layer blocks this claim?"

## 2026-06-16 Manual PELO/HK2 Translon Iteration

The manual translon refresh selected 27 manual features across 18 genes. The newest biology-relevant additions are PELO and HK2. PELO contributes two leader uORFs and one overlapping uORF on `ENST00000274311`; HK2 contributes four manual features on `ENST00000290573`, including leader/overlapping candidates and one manual NTE-like feature. All PELO/HK2 manual features were missing from the current ORFik/Transcode union, which confirms why the manual bridge is still needed.

The clean-CDS/RDG pipeline picked up the manual features without requiring additional FST pages. The feature-expression run reported translon source counts of T = 60,784, TC = 28,359, and M = 27 before downstream deduplication. After deduplication in the marker transcript set, the branch-feature source composition included 26 manual-only features and one TC+manual duplicate. PELO and HK2 RDG figures were written in the general output and in state-specific folders.

PELO is now a top-10 next-model candidate: rank 7 with next-iteration score about 7.50. Its viral-infection row is classified as context-specific review, not robust consensus, with about +10% clean-CDS effect and a leader-uORF allocation signal. This is interesting for ribosome rescue/RQC biology, but it needs browser validation before interpretation.

HK2 is now rank 11 with next-iteration score about 7.29. Its strongest consensus row is ISR/ER translation stress, classified as replicated robust consensus with about -14% clean-CDS effect. HK2 also appears in hypoxia/therapeutic-style rows, which fits glycolytic stress biology. The viral HK2 row is weak and should not be treated as a viral discovery candidate yet.

The model lesson from this iteration is sharper than the biology. Manual translon recovery can materially change candidate rankings and branch interpretation, while raw geometry and hand-built motif models still do not generalize well enough for prediction. PELO/HK2 should be used as new browser-validation and motif-calibration cases, not as final discoveries.

## Are We On Track?

Yes. The project remains aligned with the main goal: build a predictive RDG atlas that compresses thousands of Ribo-seq samples into interpretable branch, CDS-buffering, and context effects.

The project is now past a simple dominant-state QC prototype. It has become an atlas workflow:

1. Clean CDS measurement removes known uORF overlap artifacts.
2. RDG feature counts describe ribosome-routing options.
3. Grouped replicate branch models avoid the single-sample small-uORF coverage trap.
4. Consensus, LOO, and grouped-omission transport filters separate browser-priority candidates from broad low-support background.
5. The therapeutic layer translates atlas rows into testable preclinical effect-size hypotheses.

The strongest current claim is methodological. RDG-aware clean-CDS and branch-allocation models expose translation changes that simple CDS count summaries would miss. The strongest biological hypothesis is that viral infection produces reproducible ISR/antiviral RDG rewiring, with DDIT3 as the current anchor and PPP1R15A, ATF4, IFIH1, ATF3, and VEGFA as prioritized follow-up candidates. That is interesting and plausibly novel at the atlas/compendium level, but it is not yet proven disease mechanism.

The near-term paper shape still looks good:

1. Introduce a compendium-scale RDG atlas and CDS-buffering framework.
2. Validate the framework on known ISR/uORF biology, especially ATF4 and DDIT3/CHOP-axis behavior.
3. Use viral/post-viral translation as the real discovery test, with clear uncertainty gates.
4. End with preclinical perturbation hypotheses, not clinical treatment claims.

## Next Model

Immediate priorities:

1. Use `Review -> Review Workbench` in the atlas browser as the immediate triage surface. Start with the target-label filter when the goal is balanced label acquisition, then use the active review batch and label-gap queues, including B03 protocol/pooling audit rows. Switch to evidence-gap or qualified-inspection rows when the goal is negative controls, measurement/translon review, context replication, or robustness follow-up.
2. Start each biological interpretation from `pooling_proof_tier` and `pooling_hazard_class`. Exact strict project-matched control rows can be called claim candidates after browser review. Pooled recurrence, context-specific, partial-pooled, and moderate/high hazard rows should be labeled as review leads, confounded, or context-limited until local controls and robustness agree.
3. Use `Evidence -> Data Readiness` first to confirm the row is backed by mounted S-drive coverage/reference data, then use `Evidence -> Inspection Assets` before making claims from a candidate. Check coverage breadth, feature allocation, RDG-flow annotation, replicate diagnostics, exact case/control feature contrasts, and sequence-feature coordinates. A row with low top-feature counts, single-run dominance, missing mRNA sequence, or nonmatching RDG annotation should be labeled as measurement-risk or not reviewable before it becomes a biological claim.
4. Use `Review -> Label Dashboard` before each labeling session to choose the needed label class. The current bridge is positive-heavy, so negative, weak-positive, confounded, measurement-problem, and negative-control rows should often outrank more positive-looking candidates.
5. Use `Review -> Review Labels` or `dominant_rdg_review_links/dominant_rdg_review_batch_template_with_links.csv` as the writeback surface for completed manual review. From the label page, open Observatory or Inspection Assets for the selected task, then use the quick-fill buttons to set the common label pattern and fill review status, browser signal grade, branch/CDS direction checks, context match, translon status, metadata-confounder status, negative-control result, and notes. For a focused labeling pass, start from `Label Dashboard -> Start Label Session` or set the Workbench queue/target/gene/scope and click `Start Session`; use the label-task, inspection-packet, and protocol-audit gates when the next action should be constrained by review surface rather than queue source. Use `Use Sidebar Genes`, `Use Selected Gene`, or selected-row `Gene Session` when starting from the current visual context, and check the Workbench note/dashboard before starting so ignored sidebar filters and selected-row destinations are visible. Click a specific Workbench row and use `Open Selected` when the review should begin at a named gene such as `GABARAPL1`; it will open Review Labels, Inspection Assets, or Selected Card depending on the row's available context. Then use `Save & Next` so completed draft or ingested-feedback rows drop out of the open task set. After each review pass, rerun `run_dominant_cell_state_pipeline(steps = c("rdg_review_feedback", "rdg_validation_labels"), force = TRUE)` and inspect `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_gap_queue.csv`, `dominant_rdg_validation_label_bridge/dominant_rdg_validation_label_balance.csv`, `dominant_rdg_review_feedback/dominant_rdg_review_feedback_labels.csv`, and `dominant_rdg_review_feedback/dominant_rdg_review_feedback_browser_validation_notes_candidate.csv`.
6. Use `dominant_rdg_evidence_gaps/dominant_rdg_evidence_gap_priorities.csv` as the full acquisition queue when you need more rows. The first question should not be "what is the top score?" but "what evidence layer is missing?" Browser labels, translon fixes, protocol/pooling audits, context replication, metadata audits, perturbation follow-up, and negative controls are different jobs.
7. Browser-validate the top batch rows first: ATF4, PPP1R15A, DDIT3, ATF3, VEGFA, IFIH1, DDIT4, TRIB3, MRPS12, CCL2, IRF7, MYC, TFAM, LDHA, HSPA5, and the batch-specific negative controls. Mark each row as positive, weak positive, negative, confounded, or not reviewable.
8. Feed those manual labels back into motif calibration. The current motif-prior correlation with branch evidence is positive but weak, and the newest predictive motif model is neutral/slightly worse than context-only. Validated positives/negatives are more valuable than another unsupervised scoring layer.
9. Audit PELO and HK2 in RiboCrypt. PELO is a ribosome-rescue/RQC hypothesis with viral context-specific evidence; HK2 is a glycolytic-stress/hypoxia/ISR hypothesis with weak viral evidence.
10. Use `rdg_plot_mode = "tables"` for most table/model/app-data iterations and reserve full RDG plotting for manual-translon, figure-layout, and visual-review refreshes.
11. Treat grouped same-condition replicate models as the default branch evidence layer. Single runs are not generally sufficient for short uORFs.
12. Keep the viral/fatigue test focused on consensus plus transport evidence: DDIT3 and PPP1R15A are the strongest robust anchors; ATF4 and VEGFA are direction-robust but thinner; IFIH1 is context-dependent; ATF3 is a strong clean-CDS/browser lead that needs stronger branch-model support; BBC3 is a new replicated branch-consensus candidate needing browser review.
13. Audit the manual translon list against browser validation. The manual list is useful, but it is a bridge until the new translon predictor handles overlapping uORFs better.
14. Keep iORF/frame claims diagnostic until the p-shift denoising and improved translon prediction projects are available.
15. Extend the atlas from marker genes to all confident uORF genes only after a stricter transcript/translon QC gate is in place. More genes will help discovery, but only if false uORF/iORF structure does not dominate the model.
16. Replace independent evidence screens with a true joint multilevel branch model: gene, study, design, branch class, context, motif prior, count depth, protocol/inhibitor, tissue/cell line, and condition effects should be in one model rather than combined after the fact.
17. Learn perturbation response fractions from data. The current full/half-normalization predictions are explicit review targets, but the response fraction should eventually be estimated from matched perturbation assays or validated perturbation contexts.
18. Use `dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv` when the review goal is specifically merged-replicate clean-CDS/uORF redistribution. Start with `next_review_action`, `pooling_benefit_class`, `pooling_loss_class`, and `missing_evidence_flags`, not the score alone. A candidate with strong opposition but protocol/start-site audit need, direction conflicts, missing static coverage, missing RDG annotation, or single-run dominance is a review lead, not a claim. Use `merged_replicate_translon_shift_review_queue.csv` only after selecting the collapsed candidate to inspect its context rows.

Medium-term priorities:

1. Use the lab translon predictor to replace most manual overlapping-uORF annotations.
2. Use denoised p-shifted signal to revisit iORF and off-frame viral reprogramming.
3. Add dimensionality reduction over RDG motifs and state scores so the model becomes less token-heavy and more R-callable.
4. Build a RiboCrypt-facing atlas viewer around gene cards, browser links, branch matrices, validation labels, and CDS-buffering estimates.
5. Move from candidate discovery to calibrated counterfactual prediction: given gene, cell type, condition, and RDG structure, estimate expected CDS buffering and uncertainty for unobserved combinations.

Bottom line: the project is on track. The expanded FST run strengthens the atlas, gives concrete viral ISR candidates, and produces a usable validation and acquisition queue. The next decisive step is not another broad screen; it is browser validation, negative controls, context replication, and label-driven recalibration.

## 2026-09-24 First Full-Server Run and Iteration Review

This repository was developed and validated against a subset of `.fst` files on
a local machine. This date marks the first attempt to run the full pipeline
against the complete backing dataset on the analysis server. Getting there
required fixing several real environment bugs that were invisible on the local
subset (a non-existent `ORFik::txdbFile()` call, `.qs`-versus-`.rds` count
table extensions, a stale ORFik config cached under isolated XDG/BiocFileCache
state, and a project-specific `collection_tables_indexed` path convention that
seven scripts had silently missed) plus one real circular bootstrap dependency
between `manual_translons` and `clean_cds` on a fresh checkout. The EIF4G2 raw
tRNA/UV reprocessing was also independently reproduced on the server by reusing
massiveNGSpipe's already-trimmed reads instead of re-downloading, reproducing
the manuscript's qualitative result.

A full forward-looking review -- what is good, what is not optimal, what is
missing, and prioritized next steps for both the EIF4G2 and post-viral-fatigue
stories -- is in `dominant_cell_state_iteration_review_2026-09-24.md`. The
short version: EIF4G2 is ready for its first lab experiment today; the
fatigue/Long-COVID story has comparably strong computational candidates but
still has no lab go/no-go plan, and writing one is the single highest-priority
next step.
