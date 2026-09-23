# Main figure plan and draft legends

This is a four-display-item Nature Article plan. It separates completed
public-data evidence from the experiments required to make the mechanistic and
immune claims. A panel labelled **required** must not appear as a result until
the experiment exists. All final figures should show individual biological
replicates, use a colour-blind-safe palette, and retain the same condition
colours across figures.

Suggested palette:

- tRNA-Glu(UUC): blue (`#0072B2`)
- tRNA-Arg(CCG): orange (`#E69F00`)
- UV or recycling stress: purple (`#7B3294`)
- matched control: dark grey (`#4D4D4D`)
- uORF occupancy: vermilion (`#D55E00`)
- clean-CDS occupancy: green-blue (`#009E73`)

## Figure 1 | tRNA-Glu selectively redistributes EIF4G2 ribosomes into uORF1

### Story

This figure earns the word *selective*. It should establish the exact branch,
the real biological sample size, the matched tRNA control and the unbiased
feature-wide rank before technical details appear.

### Panels

- **a, discovery and validation flow.** Compendium contrast → 95 review rows →
  60 candidate-direction rows → exact biological-unit and artifact gates →
  EIF4G2 uORF1. Show that no row was accepted directly from the initial queue.
- **b, transcript architecture.** ENST00000339995: 308-nt leader; uORF1 at
  23–73 (ATG, 51 nt); uORF2 at 128–220 (TTG, 93 nt); clean CDS at 309–3032
  (GUG start). Place the uORF1 peptide below the feature and indicate the two
  GAG codons and phase-0 position-35 pause.
- **c, biological-unit allocation.** Four points for tRNA-Glu versus control
  and four points for tRNA-Arg versus its control, joined only when the
  biological pairing is defensible. The tRNA-Glu values are 0.700, 0.706,
  0.826 and 0.851; tRNA-Arg values are 0.863, 0.842, 0.832 and 0.840.
- **d, feature decomposition.** uORF1 and uORF2 counts or densities per
  biological unit, making clear that uORF1 carries the phased signal.
- **e, matched feature-wide specificity.** Rank/outlier plot for all eligible
  non-overlapping leader uORFs at the strict 20-count gate; label EIF4G2 rank
  1 of 18 and show the same rank at 5- and 10-count thresholds in an inset.
- **f, claim boundary.** Small schematic: observed tRNA-selective ribosome
  allocation → candidate effect on downstream DAP5 and peptide presentation;
  arrows beyond occupancy should be dashed until experimentally measured.

### Current source data and assets

- `results/dominant_rdg_eif4g2_evidence_packet/eif4g2_biological_unit_evidence.csv`
- `results/dominant_rdg_eif4g2_evidence_packet/eif4g2_uorf1_biological_unit_evidence.csv`
- `results/dominant_rdg_trna_codon_branch_audit/trna_codon_feature_audit.csv`
- `results/dominant_rdg_trna_codon_branch_audit/eif4g2_specificity_rank_sensitivity.csv`
- Existing specificity panel:
  `results/dominant_rdg_trna_codon_branch_audit/trna_glu_specificity_outlier.png`
- Existing overview:
  `results/dominant_rdg_eif4g2_evidence_packet/eif4g2_evidence_overview.pdf`

### Draft legend

**Figure 1 | tRNA-Glu selectively redistributes EIF4G2 ribosome occupancy into
an ATG uORF.** **a**, Focused candidate-discovery and validation workflow.
Initial context rows were collapsed by gene, transcript, design and direction
before exact biological-unit, count, phase, shape, coordinate and position-
deletion tests. **b**, Architecture of EIF4G2 ENST00000339995. Transcript
positions are one-based; the main CDS begins at the native GUG codon. The
uORF1 peptide and its two GAG codons are shown. **c**, Clean-CDS allocation for
tRNA-Glu(UUC) and tRNA-Arg(CCG) overexpression and their matched controls. Each
point is one biological unit; three technical lanes were summed within each
tRNA-Glu unit (`n = 2` biological units per arm for both perturbations).
Horizontal bars denote aggregate count-based allocation. **d**, Individual-uORF
occupancy and phase support for the same biological units. **e**,
tRNA-Glu-specific uORF-allocation effects among non-overlapping leader uORFs
meeting the strict count gate in both tRNA perturbations. EIF4G2 is labelled;
the robust distance uses unscaled median absolute deviation and is descriptive,
not a z-score. **f**, Completed occupancy result and untested downstream model.
Source data are provided as a Source Data file.

## Figure 2 | Orthogonal falsification of the tRNA-Glu branch shift

### Story

This figure earns the phrase *artifact-resistant*. Raw tRNA-Glu and tRNA-Arg
reads, not the old pooled compendium tracks, should carry its main panels. The
original lane-linkage audit remains visible in Extended Data as an independent
reconstruction check.

### Panels

- **a, raw reconstruction.** Sixteen FASTQs → adapter trimming with 15--40-nt
  inserts retained → unique GRCh38 mapping → 4,670 CDS-start windows →
  run-by-length offsets → study- and lane-shared phase gates → exact EIF4G2
  mapping.
- **b, run and calibration matrix.** All 12 tRNA-Glu/control and four
  tRNA-Arg/control runs, showing passing lengths, exact-uORF1 phase and target
  counts. Keep lane and biological-unit labels explicit.
- **c, biological-unit allocation.** Exact-uORF1 and union-uORF allocation in
  the all-lane-matched, study-shared, canonical and short strata. Show all true
  biological units; tRNA-Arg is the specificity control. Mark the sparse short
  stratum as a failed measurement/target-phase gate rather than supportive
  evidence.
- **d, exact-length and lane consistency.** Allocation contrast by exact
  footprint length and by the three matched tRNA-Glu lanes. Indicate which
  lengths pass every required run rather than selecting by effect direction.
- **e, positional falsification.** Raw uORF1 profile with position 35 labelled,
  beside observed, coordinate-minimum, worst-single-position and greedy-three-
  position-deletion effects.
- **f, RNA and failed codon alternatives.** Compact total-RNA/isoform panel and
  the null 18-feature Glu-codon-dose test. Label GAG sensing as a reporter-level
  hypothesis.

### Current source data and assets

- `results/dominant_rdg_eif4g2_trna_length_resolved/psite_offset_qc.csv`
- `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_unit_length_metrics.csv`
- `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_length_contrast_summary.csv`
- `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_length_specificity_summary.csv`
- `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_position_profile.csv`
- `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_position_leaveout_sensitivity.csv`
- `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_trna_length_figure_source_data.csv`
- Current raw overview:
  `results/dominant_rdg_eif4g2_trna_length_resolved/eif4g2_trna_length_falsification.pdf`
- `results/dominant_rdg_eif4g2_trna_lane_subset_audit/eif4g2_trna_glu_lane_subset_audit.csv`
- `results/dominant_rdg_eif4g2_evidence_packet/eif4g2_position_leaveout_evidence.csv`
- `results/dominant_rdg_eif4g2_goodarzi_rna_audit/eif4g2_isoform_rna_audit_summary.csv`
- `results/dominant_rdg_trna_codon_branch_audit/trna_codon_analysis_set.csv`
- Existing audit plot:
  `results/dominant_rdg_trna_codon_branch_audit/trna_glu_codon_branch_audit.png`

### Draft legend

**Figure 2 | Raw length-resolved reads validate the tRNA-selective EIF4G2
branch shift.** **a**, Raw-read processing and P-site calibration for the 12
tRNA-Glu/control and four matched tRNA-Arg/control libraries. Offsets were
inferred independently for every run and read length against 4,670
non-overlapping CDS-start windows. **b**, Per-run calibration and exact-uORF1
phase support. Sequencing lanes remain technical observations and are not
inferential replicates. **c**, EIF4G2 clean-CDS allocation for exact uORF1 and
the uORF union in study-shared and lane-matched footprint strata (`n = 2`
biological units per arm). **d**, Effect direction by exact footprint length
and matched technical lane; filled symbols denote lengths that passed every
required run. **e**, Raw EIF4G2 position profile and effect remaining after
coordinate perturbation, worst single-position deletion and greedy deletion
of three positions. **f**, Source-provided NM_001418 RNA/isoform stability and
the null 18-feature glutamate-codon-dose analysis. In the decisive lane-matched
pool, the tRNA-Glu union and exact-uORF1 allocation effects were −0.136 and
−0.135, compared with +0.012 for the matched tRNA-Arg union. Study-shared and
canonical tRNA-Glu effects were −0.141 and −0.145; the sparse short stratum
failed its predeclared measurement and exact-uORF1 phase gates. Coordinate,
single-position and greedy-three-position lower bounds were −0.116, −0.100 and
−0.060 in the observed direction. Source data are provided as a Source Data
file.

## Figure 3 | Raw short and canonical footprints reproduce the switch after UV

### Story

This is an orthogonal confirmation that the EIF4G2 branch responds under an
unrelated ribosome-stress condition. It validates locus responsiveness and the
raw workflow; it is not evidence that UV and tRNA-Glu share a mechanism and is
not the paper's primary claim.

### Panels

- **a, raw workflow.** Four FASTQs → library-specific trimming → unique GRCh38
  mapping → 4,670 CDS-start windows → run-by-length offsets → shared phase
  gate → exact EIF4G2 mapping.
- **b, per-run strata.** Clean-CDS allocation for all shared, canonical
  28–32-nt and short 20–23-nt footprints. Show all four runs in all strata.
- **c, exact-length effects.** Delta for each 15–34-nt class with target counts;
  identify shared well-phased lengths and low-count non-gated lengths.
- **d, positional profile.** Shared well-phased UV and control profiles across
  the leader and beginning of the CDS; label uORF1, uORF2 and GUG CDS start.
- **e, robustness.** Observed, coordinate-minimum, worst-single-deletion and
  greedy-three-deletion deltas for all shared, canonical and short strata.
- **f, calibration QC.** Compact phase/offset heatmap for all run-length pairs,
  with the rejected impossible 16-nt control offset visibly excluded.

### Current source data and assets

- `results/dominant_rdg_eif4g2_uv_length_resolved/psite_offset_qc.csv`
- `results/dominant_rdg_eif4g2_uv_length_resolved/eif4g2_run_length_metrics.csv`
- `results/dominant_rdg_eif4g2_uv_length_resolved/eif4g2_length_contrast_summary.csv`
- `results/dominant_rdg_eif4g2_uv_length_resolved/eif4g2_position_profile.csv`
- Current four-panel figure:
  `results/dominant_rdg_eif4g2_evidence_packet/eif4g2_uv_length_falsification.pdf`

### Draft legend

**Figure 3 | UV-induced EIF4G2 redistribution recurs independently in short and
canonical footprints.** **a**, Raw-read processing and P-site calibration.
Offsets were inferred separately for each run and read length against 4,670
non-overlapping CDS-start windows; decisive lengths passed the phase gate in
all four libraries. **b**, EIF4G2 clean-CDS allocation in untreated and
UV-irradiated HeLa biological replicates (`n = 2` per arm) for all shared
well-phased, canonical 28–32-nt and short 20–23-nt footprints. **c**, Allocation
contrast by exact footprint length. Point size denotes target feature counts;
filled points passed the shared phase gate. **d**, Transcript-position profiles
for shared well-phased footprints. **e**, Observed effects and lower-bound
robustness after coordinate perturbation, worst single-position deletion and
greedy deletion of three positions. **f**, independent CDS-start phase quality
control by run and length. The physically impossible 16-nt offset in one
control was rejected before target mapping. Source data are provided as a
Source Data file.

## Figure 4 | Causal routing connects the EIF4G2 branch to DAP5 and HLA display

### Status

**This figure is required and currently missing.** It is the difference between
an unusually strong reanalysis and a Nature-level biological paper. No panel
below may be described in the Results until completed.

### Minimal panels

- **a, factorial reporter design.** Native EIF4G2 leader and GUG main start;
  WT, uORF1 ATG mutant, GAG→GAA synonymous mutant, codon-5 pause mutant and a
  peptide-disrupting/non-Glu control. Cross each with vector, tRNA-Glu and
  tRNA-Arg in at least three independently performed experiments.
- **b, branch-specific translation.** Dual readout that separately reports
  uORF1 translation and downstream GUG output. The decisive statistic is a
  construct × tRNA perturbation interaction, not a within-construct P value.
- **c, mature tRNA versus fragment route.** Northern/small-RNA measurements plus
  decoding-competent, anticodon-altered and cleavage-altered tRNA constructs.
- **d, endogenous validation.** Targeted Ribo-seq, toeprinting or equivalent
  showing the same endogenous uORF1/CDS redistribution in newly generated
  perturbations, with DAP5 immunoblot or quantitative proteomics.
- **e, immune output.** Targeted immunopeptidomics or peptide–HLA assay for the
  endogenous EIF4G2-uORF product after tRNA-Glu, tRNA-Arg and vector controls.
  Include the uORF-start mutant where technically possible.
- **f, integrated model.** Solid arrows only for measured links. Separate
  direct GAG decoding, tRNA-fragment and recycling-state models if they remain
  experimentally indistinguishable.

### Draft legend template

**Figure 4 | [Revise only after experiments: tRNA-Glu causally routes EIF4G2
translation into an immunogenic uORF].** **a**, EIF4G2 leader constructs and
tRNA perturbations. **b**, uORF1 and downstream GUG reporter output. **c**,
mature tRNA-Glu and tRNA-derived fragment abundance. **d**, endogenous EIF4G2
ribosome allocation and DAP5 protein. **e**, EIF4G2-uORF-derived peptide–HLA
abundance. **f**, experimentally supported model. For every panel, specify the
exact biological `n`, centre and error summaries, test, sidedness, adjustment,
effect size and definition of independent replication. Do not retain this
template text in a submission.

## Main-figure production checklist

- Use 180–183 mm double-column width for the multi-panel figures.
- Use Arial or Helvetica consistently, 5–7 pt at final size.
- Use bold upright lower-case panel letters.
- Avoid red–green contrast and rainbow colour scales.
- Export editable vector PDF/SVG for line art; keep text editable.
- Define every symbol, line, error bar, normalization and exact `n` in the
  legend; legends should remain below 250 words.
- Provide a machine-readable source-data table for every plotted value.
- Do not reuse panels from source publications; redraw public-data results from
  the reanalysis outputs.
