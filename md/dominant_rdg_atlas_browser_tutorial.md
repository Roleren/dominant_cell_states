# RDG Atlas Browser Tutorial

Updated: 2026-07-10

## Purpose

The RDG Atlas Browser is a review surface for dominant-state ribosome
decision graph candidates. It is not a replacement for RiboCrypt coverage
inspection. Use it to choose candidates, understand which evidence layers
support them, and open the right RiboCrypt Observatory subset for manual
validation.

The browser answers four practical questions:

- Which gene-design rows are highest priority?
- Is the signal broad, context-specific, fragile, or model-discordant?
- Which evidence table explains the priority label?
- Which coverage view should be opened next in the Observatory?

## Starting From The Standalone Repo

This analysis is being moved from:

```text
/home/roler/Desktop/forks/RiboCrypt/dominant_cell_states
```

to:

```text
/home/roler/Desktop/forks/dominant_cell_states
```

When starting a new coding/review session from the standalone location, read
these files first:

- `AGENTS.md`
- `README.md`
- `md/dominant_cell_state_scientific_note.md`
- `md/dominant_rdg_atlas_browser_tutorial.md`
- `NEWS.md`

The analysis still depends on local RiboCrypt. Set the package repo path before
running R code:

```bash
export RIBOCRYPT_REPO=/home/roler/Desktop/forks/RiboCrypt
```

Use `devtools::load_all(Sys.getenv("RIBOCRYPT_REPO"))`, not an installed
RiboCrypt package. For pipeline runs, source
`scripts/run_dominant_cell_state_pipeline.R` and let
`run_dominant_cell_state_pipeline()` load RiboCrypt after runtime setup. For
most table/model iterations, use `rdg_plot_mode = "tables"` so existing plot
assets are preserved.

Do not blindly commit generated outputs to GitHub. The local analysis folder is
large and contains generated CSVs, caches, and static inspection assets. The
standalone repo should track source-like files, docs, tests, app code, and small
curated inputs; generated outputs should stay local under ignored `results/`
unless an explicit artifact strategy is chosen. FST pages, ORFik annotations,
BAMs, reference indexes, and other Bio_data artifacts are external runtime
inputs and should not be copied into the package.

## First Pass Review

Start in `Review -> Cards`.

1. Keep the default `Review queue` preset for first-pass inspection.
2. Use the left sidebar to narrow genes, design families, evidence classes, or
   stability classes.
3. Watch the sidebar summary. It shows how many cards, genes, and design
   families remain after filtering.
4. Use the default `Review summary` table view to pick rows quickly.
5. Click a row, then open `Selected Card`.
6. Use the `Open in Observatory` button for the case/control coverage view.

The default table intentionally hides many columns. Switch to `Model and
stability`, `Context evidence`, or `Full scan columns` only when the compact
view is not enough.

For the current merged-replicate clean-CDS/uORF redistribution goal, start from
the collapsed candidate table before opening individual context rows:

```text
dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv
```

The canonical local result path is:

```text
results/dominant_rdg_merged_replicate_translon_shift/merged_replicate_translon_shift_candidate_queue.csv
```

Read `next_review_action`, `pooling_benefit_class`, `pooling_loss_class`, and
`missing_evidence_flags` before trusting the rank. Strong viral candidates such
as `PPP1R15A`, `DDIT3`, `ATF4`, `TFE3`, `IFIH1`, and `FBXL4` are currently
review leads, not claims, because many still require protocol/start-site or
pooling audit.

Use `Review -> Label Dashboard` when the question is "what label does the
model need next?" It shows current bridge labels, pending draft labels,
feedback-ingest outputs, label deficits, and the active label-gap queue. Filter
by target label class when you need negatives, confounded rows,
measurement-problem rows, or negative controls rather than another
positive-looking candidate. Selecting a gap row focuses the selected card and
active-review task, and the page has direct jumps to `Review Labels`,
`Inspection Assets`, and the Workbench. `Start Label Session` opens the first
open label-gap task for the current target filter, configures the Workbench to
the same target, and lands on `Review Labels`.

Use `Review -> Review Workbench` when the question is "what should be reviewed
next?" rather than "what is in the atlas card table?" The workbench combines the
active review batch, validation-label gap queue, full evidence-gap acquisition
queue, and qualified-inspection priority table. Selecting a workbench row
updates the selected card and, when available, the exact qualified-inspection
source row. Use the target-label filter when the review goal is a missing label
class rather than a queue type. The left sidebar `Genes` selector affects the
Workbench when `Scope` is `Current filters`; `Scope = All queued rows` ignores
the left sidebar filters. The Workbench note reports which gene scope is active,
including when sidebar genes are ignored. For an explicit gene start, use the
Workbench-local `Workbench genes` selector, type a gene such as `GABARAPL1`,
and click `Start Session`. `Use Sidebar Genes` copies the left sidebar genes
into the Workbench gene selector and switches scope to all queued rows; `Use
Selected Gene` does the same for the selected card. Use the Workbench-local
gates when you want only active label tasks, only rows with static inspection
packets, or only protocol/pooling audit rows. The Workbench dashboard reports
how many rows in the current scope open to `Review Labels`, `Inspection Assets`,
or `Selected Card`, plus protocol-audit, pooling-risk, and exact-proof counts.
To open a specific visible row, click the row and use `Open Selected`: active
label tasks open `Review Labels`, inspection-only rows open `Inspection Assets`,
and rows without either fall back to `Selected Card`. The selected-row card
shows the destination, transcript, queue rank, active task, inspection source,
target label, pooling proof/hazard, review reason, and Observatory link; `Gene
Session` turns that selected row's gene into the current Workbench gene session.
The table shows `review_destination`, `gene_symbol`, `tx_id`, and
`design_family` before the queue metadata, then carries target-label deficit,
label-acquisition score, `protocol_pooling_need`, pooling proof, and pooling
hazard fields, so B03 protocol/pooling audit rows can be reviewed without
leaving the Workbench. This is the fastest route for label acquisition, negative
controls, context audits, and measurement/translon review rows. `Start Session`
opens the first open active-review task under the current Workbench filters.

Use `Review -> Review Labels` after inspecting a workbench row. It records the
active-review task, browser signal grade, branch/CDS/context direction checks,
translon/isoform status, metadata confounders, negative-control result, reviewer
notes, and next action. The action row opens the selected task in Observatory,
`Inspection Assets`, `Selected Card`, the Workbench, or `Files` without leaving
the review context behind. The quick-fill buttons set the common positive,
weak-positive, negative, confounded, measurement-problem, not-reviewable, and
negative-control patterns before you add row-specific notes. Use `Save Draft
Label` when you want to stay on the same task, or `Save & Next` to advance to
the next open active-review task in the current Workbench queue, target, and
scope. For target-balanced review, start from `Label Dashboard -> Start Label
Session`; for protocol or evidence audits, start from `Review Workbench -> Start
Session`. Saving a draft writes
`dominant_rdg_review_feedback/dominant_rdg_review_label_drafts.csv`; it does not
edit the linked review CSV or curated browser-validation notes. Draft rows also
retain the Workbench target-label and protocol/pooling context so the feedback
outputs can be audited back to the acquisition reason. Rerun the review-feedback
step after a review pass so these drafts become normalized label rows and
calibration candidates.

## Selected Card

The `Selected Card` page is the main single-candidate summary. The metric strip
at the top gives the fast triage view:

- `confidence`: atlas confidence score.
- `clean-CDS effect`: estimated clean-CDS output change.
- `model agreement`: agreement tier across branch models and browser labels.
- `study omission`: leave-one-study-out robustness class.
- `grouped omission`: transportability under grouped omission.
- `top branch`: dominant branch class from joint allocation.
- `context interaction`: whether a tissue or cell-line context changes the
  branch effect.
- `partial pooling`: whether a context effect survives shrinkage.

Use `Evidence summary` for normal review. Use the model/context detail views
when you need to understand why a row was ranked. Use `Full field dump` only
when debugging a pipeline column or checking a rare warning.

## Observatory Links

The `Links` column and `Open in Observatory` button encode the best available
case/control subset.

The intended inspection pattern is:

1. Confirm the transcript isoform and branch annotations.
2. Compare case and control coverage around the predicted branch point.
3. Check whether the signal is a broad region shift or a single-position spike.
4. Check whether clean-CDS, leader-uORF, overlapping-uORF, and internal-ORF
   patterns match the expected branch interpretation.
5. Record browser validation only when the coverage supports the model label.

If the link note says the subset fell back from exact context to study and
condition, treat the view as weaker evidence. If the row is gene/transcript
only, rebuild the subset manually before making a claim.

## Evidence Pages

Use the `Evidence` group after selecting a card.

`Inspection Assets` is the first page to open for qualified review. It shows
the selected review source row, readiness flags, count gates, coverage profile,
feature allocation, RDG-flow annotation, replicate diagnostics, feature
contrasts, replicate summaries, and transcript/RDG feature coordinates. Use this
page to decide whether the signal is broad enough to claim, whether the top
feature is supported by multiple runs, and whether the visible coverage matches
the translon/RDG annotation.

`Data Readiness` is the mounted-data audit. It answers a different question:
whether the selected card and global atlas page set can be reopened from the
current S drive. It checks the ORFik all-merged library, reference FASTA/index,
required all-sample FST pages, static inspection PNGs, and sequence flags. A
row can be `sdrive_ready` and still not be claimable from coverage; use this
page to separate data-access problems from biological/coverage-support
problems.

`Branch Contexts` shows the matched study/condition rows behind the selected
card. This is the best page for seeing whether a row is driven by one local
context or repeated across several contexts.

`Study Evidence` shows per-study support behind model agreement and
leave-one-study-out labels. A strong row should keep the same direction when
important studies are omitted.

`Transport Evidence` shows grouped omission. This is stricter than
leave-one-study-out because it removes related study families, protocols,
tissues, cell lines, or condition groups.

`Context Interactions` compares one viral tissue or cell-line context against
the remaining viral contexts. This page is important for genes such as IFIH1,
where lung-like contexts can behave differently from the broad viral design.

`Partial Pooling` shows whether context-specific effects survive empirical
Bayes shrinkage. Pooled-supported context rows are stronger than raw
single-context signals.

## Visual Pages

The `Visuals` group gives global filtered views.

`Heatmap` shows clean-CDS buffering effect by gene and design family. This is a
protein-output-oriented view.

`Branch Heatmap` shows joint branch-allocation shift magnitude. This can reveal
branch redistribution even when the aggregate clean-CDS effect is modest.

`Effect Space` shows confidence versus clean-CDS effect, colored by evidence
class. Use it to see whether the current filter is dominated by large effects,
high-confidence rows, or low-support exploratory rows.

## Reports

The `Reports` page embeds generated global HTML reports. Use it when the
question is broader than one selected card:

- dominant RDG review queues and model diagnostics;
- uORF and ORF-start statistics;
- global state, frame, infected/off-frame, and metadata QC.

Reports are static files, so switching reports is fast. Open a report in a new
tab when you need more screen width.

## Claim Gates

Use conservative language until the evidence is strong.

Review-grade evidence usually needs:

- enough branch counts in merged case and control groups;
- a coverage pattern that is not a single-position spike;
- agreement between the selected-card expectation and the Observatory view;
- study omission that preserves direction;
- grouped omission or explicit context evidence when making broad claims;
- translon coordinates that match the visible coverage and branch annotation.

Context-specific claims need explicit context support. A row can be real in
lung, A549, RD, LCL, or another context without being a universal viral
infection rule.

Counterfactual or therapeutic claims need more than the atlas browser. Treat
therapeutic pages as hypothesis prioritization, not treatment evidence.

## Practical Review Order

For the current viral/post-viral review goal, use this order:

1. `Review -> Label Dashboard`, open rows only.
2. Pick the target label class that matches the job: positive/negative,
   confounded, measurement problem, or negative control.
3. Select a gap row, then jump to `Inspection Assets` or `Review Labels`.
4. `Review -> Review Workbench`, scope `Current filters`, when you need to mix
   active review, evidence gaps, and qualified-inspection priorities.
5. Pick the target-label class first if the job is label-balance driven; pick
   the queue first if the job is active review, evidence gap, or qualified
   inspection triage.
6. Select a candidate and read `Selected Card -> Evidence summary`.
7. Open `Evidence -> Data Readiness` to confirm the selected row is backed by
   mounted S-drive coverage/reference data.
8. Open `Evidence -> Inspection Assets` and check coverage, feature allocation,
   RDG flow, replicate diagnostics, and sequence-feature coordinates.
9. Open the Observatory link and inspect coverage interactively when the static
   inspection assets look plausible or ambiguous.
10. Record the decision in `Review -> Review Labels`.
11. Check `Branch Contexts` and `Study Evidence`.
12. Check `Transport Evidence` for broad claims.
13. Check `Context Interactions` and `Partial Pooling` for lung or other
   context-specific hypotheses.
14. Use `Reports` for the global picture before promoting a candidate.

DDIT3 is the current browser-supported consensus anchor. ATF4, PPP1R15A,
IFIH1, VEGFA, BBC3, ATF3, PELO, and HK2 are useful review examples, but each
has different stability and context-dependence labels.
