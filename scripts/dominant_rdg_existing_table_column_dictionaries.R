#!/usr/bin/env Rscript
#
# Column dictionaries for two EXISTING output tables that were not
# self-documenting and directly caused confusion this session:
#
#   - dominant_rdg_consensus_candidates.csv (dominant_rdg_consensus_atlas.R)
#   - dominant_rdg_context_interaction_branch_context_effects.csv
#     (dominant_rdg_context_interaction_model.R) -- its `branch_class`
#     column required reading the producing script directly to understand,
#     which delayed (and initially caused a wrong conclusion in) the VEGFA
#     investigation. See
#     papers/to_be_explored_more_later/VEGFA_uorf_viral_infection/README.md
#     and md/dominant_state_pipeline_lessons_2026-09-28.md gap #4.
#
# This is deliberately a separate, documentation-only script rather than a
# modification of the two heavy, already-validated producing scripts --
# it reads their already-written output tables and writes a column
# dictionary next to each, using the same shared helper new scripts use
# going forward (dominant_column_dictionary.R). Column descriptions below
# reflect this project's actual working understanding of these tables as
# of 2026-09-28; columns not confidently understood are left undocumented
# on purpose rather than guessed at.

suppressPackageStartupMessages(library(data.table))

source(file.path(Sys.getenv("DOMINANT_CELL_STATES_DIR", "."), "scripts",
                 "dominant_script_common.R"), local = TRUE)

analysis_dir <- find_analysis_dir()
results_dir <- Sys.getenv("DOMINANT_CELL_STATES_RESULTS_DIR",
                          unset = file.path(analysis_dir, "results"))
result_file <- function(...) file.path(results_dir, ...)
source(file.path(analysis_dir, "scripts", "dominant_column_dictionary.R"))

# --- context_interaction_branch_context_effects.csv --------------------
context_file <- result_file("dominant_rdg_context_interactions",
                            "dominant_rdg_context_interaction_branch_context_effects.csv")
write_column_dictionary(
  context_file,
  c(
    gene_symbol = "Gene this row is about.",
    tx_id = "Transcript ID scored.",
    design_family = "Design family this row's context comparison is scored within. As of 2026-09-28 this table only covers 'viral_infection' -- other design families are not yet evaluable here (see dominant_rdg_consensus_context_consistency_audit.R's own summary output).",
    branch_class = "Which structural feature this row's delta refers to: leader_uORF, overlapping_uORF, clean_CDS, or other (see scripts/dominant_rdg_context_interaction_model.R's branch_levels). Rows for the same gene/context but different branch_class are NOT contradictory measurements of the same thing -- they are different features, and are expected to move in opposite directions when ribosomes are being reallocated (e.g. more uORF engagement, less clean-CDS usage). This is the column that caused the VEGFA investigation's initial confusion.",
    context_axis = "Which metadata field defines the 'context' being tested (e.g. cell_line_context, tissue_context).",
    context_value = "The specific level of context_axis this row's context group is (e.g. 'Calu3' for cell_line_context).",
    context_studies = "Number of studies contributing to this specific context group.",
    complement_studies = "Number of studies contributing to the complement (everything NOT in this context group, within the same design_family).",
    total_axis_studies = "Total studies across all levels of this context_axis.",
    context_delta = "This specific context's own effect size for branch_class (NOT the pooled/global estimate -- see global_axis_delta for that).",
    context_delta_se = "Standard error of context_delta.",
    context_delta_p = "P-value for context_delta.",
    context_direction_fraction = "Fraction of studies within this context whose individual effect agrees in direction with context_delta.",
    context_ci_excludes_zero = "TRUE if context_delta's confidence interval excludes zero.",
    context_support_class = "This context's own support/robustness tier.",
    complement_delta = "The complement group's (everything NOT in this context) own effect size for branch_class.",
    global_axis_delta = "The POOLED, non-context-specific estimate for branch_class across the whole design_family -- constant across all context_value rows for a given gene/design_family/branch_class. This is the number to trust over any single context_delta or over the separate consensus table's own headline number when they disagree (see dominant_rdg_consensus_context_consistency_audit.R).",
    global_axis_delta_se = "Standard error of global_axis_delta.",
    global_axis_delta_p = "P-value for global_axis_delta.",
    interaction_delta = "context_delta minus complement_delta -- how much this specific context's effect differs from the rest of the design family.",
    interaction_delta_p = "P-value for interaction_delta.",
    context_opposite_complement_direction = "TRUE if this context's effect direction is opposite to its complement's.",
    context_interaction_evaluable = "Whether this row had enough data to compute an interaction test at all.",
    context_interaction_class = "Overall interaction-evidence tier for this row (e.g. context_magnitude_modulated, weak_context_difference, not_evaluable).",
    context_study_names = "Semicolon-separated list of studies contributing to this context group.",
    complement_study_names = "Semicolon-separated list of studies contributing to the complement group.",
    model_agreement_score = "Cross-model agreement score for this row.",
    model_agreement_tier = "Cross-model agreement tier label.",
    transport_stability_score = "Leave-one-axis-out transport robustness score for this row.",
    transport_agreement_class = "Transport robustness tier label.",
    transport_consensus_class = "Whether this row's direction/support survives the transport (leave-one-axis-out) robustness check.",
    context_review_priority = "This table's own suggested priority for manual review of this row."
  ),
  title = "Per-context, per-branch effect sizes and pooled global estimates for the viral_infection design family (currently the only design family this deeper context-interaction model covers)."
)

# --- consensus_candidates.csv -------------------------------------------
consensus_file <- result_file("dominant_rdg_consensus", "dominant_rdg_consensus_candidates.csv")
write_column_dictionary(
  consensus_file,
  c(
    consensus_rank = "This row's rank by consensus_score within the whole table.",
    gene_symbol = "Gene this row is about.",
    tx_id = "Transcript ID scored.",
    design_family = "Design family (condition category) this consensus claim is scored under.",
    consensus_class = "Overall tier label combining the many *_component scores below (e.g. B_replicated_robust_consensus, C_multi_model_branch_consensus).",
    consensus_score = "This table's single headline score (higher = stronger consensus). A weighted combination of the *_component columns near the end of this table.",
    evidence_tags = "Semicolon-separated tags summarizing this row's evidence profile (e.g. branch_supported, pooling_hazard, fragile, single_context) -- read these before trusting consensus_score alone.",
    atlas_evidence_class = "Evidence tier from the base RDG atlas (independent of the consensus-specific scoring below).",
    atlas_clean_cds_pct_change = "This table's own headline percent-change number for clean-CDS usage. IMPORTANT: this can disagree in sign with the more carefully pooled global_axis_delta in dominant_rdg_context_interaction_branch_context_effects.csv for the same gene/design_family (found for VEGFA, MMP2, LMNB1 under viral_infection as of 2026-09-28) -- cross-check with dominant_rdg_consensus_context_consistency_audit.csv before quoting this number.",
    atlas_direction = "This table's own direction label for atlas_clean_cds_pct_change.",
    n_studies = "Number of studies supporting this row anywhere in the atlas's evidence.",
    strict_supported_studies = "Number of those studies meeting a stricter support bar.",
    total_case_samples = "Total case-condition samples across all supporting studies.",
    total_control_samples = "Total control-condition samples across all supporting studies.",
    best_case_condition = "The single condition label contributing this row's strongest evidence.",
    best_study = "The single study contributing this row's strongest evidence.",
    top_shifted_feature_id = "The specific structural feature (e.g. a specific uORF) driving this row's headline effect.",
    top_shifted_feature_type = "Feature type of top_shifted_feature_id (e.g. leader_uorf, clean_cds).",
    model_agreement_tier = "Cross-model agreement tier (see dominant_rdg_model_agreement.R).",
    loo_consensus_class = "Whether this row's direction/support survives leave-one-out (across studies) robustness checking.",
    transport_consensus_class = "Whether this row's direction/support survives leave-one-axis-out (transport) robustness checking.",
    pooling_hazard_class = "This design_family's own pooling-confound hazard tier (see dominant_rdg_pooling_bias_audit's pooling_bias_report.md) -- a property of the design_family in general, not specific to this gene.",
    pooling_hazard_score = "Numeric version of pooling_hazard_class.",
    context_interaction_gene_class = "Whether this gene showed a meaningful context-specific interaction effect (see the context-interaction table above) -- only evaluable for design families that table covers.",
    therapeutic_gene_hypothesis = "Whether this gene is part of a registered therapeutic hypothesis (see dominant_rdg_therapeutic_hypotheses.R).",
    validation_assay = "Suggested follow-up validation assay for this row, if any."
  ),
  notes = c(
    "The many `*_component` columns near the end of this table",
    "(`atlas_component`, `branch_component`, `model_component`,",
    "`robustness_component`, `browser_component`, `output_effect_component`,",
    "`context_component`, `aux_context_component`, `therapeutic_component`,",
    "`pooling_component`, `pooling_hazard_component`, `sparse_component2`,",
    "`novelty_component`, `fragility_component`, `warning_penalty`) are the",
    "individual weighted sub-scores that sum into `consensus_score`. Their",
    "exact formulas live in `scripts/dominant_rdg_consensus_atlas.R` and are",
    "not restated here -- treat them as internal scoring detail rather than",
    "independently-interpretable evidence unless you've read that script.",
    ""
  ),
  title = "This atlas's single flagship gene x design_family consensus table -- 82 columns; only the ones most commonly relied on this session are documented above, not all of them."
)
