# EIF4G2 manuscript workspace

Updated: 2026-08-19

This directory is the editable source for a possible Nature Article about the
tRNA-Glu-sensitive EIF4G2 upstream open reading frame (uORF). It is deliberately
split into files that are easy to review in Git but can be combined into one
Word document for Google Docs.

## Scientific status

The existing analysis supports a large, tRNA-selective redistribution of
ribosome occupancy into the exact 51-nt immunogenic EIF4G2 ATG uORF1. The result
survives biological-unit collapse, matched tRNA-Arg controls and a complete
16-library raw reconstruction with run-by-length P-site calibration,
study/lane-shared and canonical footprint strata, exact-uORF1 phase,
coordinate and position-deletion gates, plus a matched total-RNA/isoform audit.
The decisive raw tRNA-Glu and tRNA-Arg effects are -0.136 and +0.012. UV irradiation independently reproduces
the redistribution in raw, separately calibrated short and canonical
footprints, and ABCE1 loss supplies a second recycling-related recurrence.

The existing analysis does **not** establish that the two GAG codons cause the
effect, distinguish mature tRNA decoding from a tRNA-derived fragment route,
show altered DAP5 protein, or show increased EIF4G2-derived HLA peptide after
tRNA-Glu perturbation. A simple transcriptome-wide Glu-codon-dose model is
null, and a later metastatic-state dataset does not reproduce the original
direction. Those negative results remain in the manuscript.

The paper is therefore a rigorous working draft, not a submission-ready Nature
paper. The experiments that separate those states are specified in
`eif4g2_claims_evidence_and_gaps.md`; the smallest economical first decision is
laid out in `eif4g2_lab_go_no_go_plan.md`.

## Files

- `eif4g2_article_draft.md`: Nature-style title, summary paragraph, main text,
  and numbered main references.
- `eif4g2_figure_plan_and_legends.md`: four-main-figure story, current assets,
  missing panels, and draft legends.
- `eif4g2_methods.md`: reproducible computational Methods, statistics, and
  draft availability statements.
- `eif4g2_extended_data.md`: Extended Data plan and draft legends.
- `eif4g2_claims_evidence_and_gaps.md`: claim ledger, fatal alternatives,
  submission gates, and exact experiment priorities.
- `eif4g2_lab_go_no_go_plan.md`: two-screen reporter/UV triage with explicit
  go, stop and spending criteria.
- `eif4g2_cover_letter.md`: an editor-facing pitch that is intentionally held
  below the final causal claim.

The detailed scientific audit behind the manuscript remains in
`../md/eif4g2_recycling_stress_working_note.md`. Machine-readable local outputs
are under ignored `../results/` directories and are rebuilt by the source
scripts in `../scripts/`.

## Google Docs build

Run from the repository root:

```bash
scripts/build_eif4g2_manuscript.sh
```

The build first runs a focused manuscript audit against the saved evidence
tables, then writes `manuscript/generated/eif4g2_working_manuscript.docx` and
`manuscript/generated/eif4g2_manuscript_audit.txt`. Upload the Word file to
Google Drive and open it with Google Docs. The generated directory is ignored
because the Markdown files are canonical and Word files are binary.

The build is intentionally text-first. Existing analysis figures are listed in
the figure plan but are not silently presented as submission-quality main
figures. Once the main panels are finalized, they can be inserted into the
Google Doc or into a review PDF.

## Current Nature format target

The structure follows the Nature author pages accessed on 2026-08-19:

- [Formatting guide](https://www.nature.com/nature/for-authors/formatting-guide)
- [Initial submission](https://www.nature.com/nature/for-authors/initial-submission)
- [Reporting standards](https://www.nature.com/nature/editorial-policies/reporting-standards)

The working target is a typical six-page Article: a fully referenced summary
paragraph of about 200 words, about 2,500 words of main text, four modest main
display items, no more than about 50 main references, figure legends below 250
words each, and online Methods containing sufficient detail to reproduce the
analysis. Nature currently requires exact sample sizes, a statistics section,
a life-sciences reporting summary, and separate Data and Code Availability
statements. The initial submission should normally be one Word or PDF file with
figures and legends together and line numbers. Current Nature guidance also
requires disclosure of large-language-model use; a draft disclosure is kept in
the Methods and must not be removed unless journal policy changes.

## Editing rules

1. Do not convert a hypothesis into past-tense experimental prose before the
   corresponding result exists.
2. Preserve the distinction between biological replicates and sequencing
   lanes. The tRNA-Glu comparison is `n = 2` versus `n = 2`, not six versus six.
3. Call the primary quantity a ribosome-occupancy allocation, not an initiation
   rate, protein output, or antigen abundance.
4. Do not describe raw-MAD outlier scores as z-scores or rank fractions as
   inferential P values.
5. Keep the null codon-dose test and metastatic-state non-replication visible.
6. Replace placeholders for authors, affiliations, funding and repository DOI
   only with verified information.
