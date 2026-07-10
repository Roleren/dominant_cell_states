#!/usr/bin/env Rscript

# Rank what the RDG atlas needs next: browser validation, stronger context
# replication, translon review, perturbation follow-up, or negative controls.

dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
output_dir <- file.path(analysis_dir, "dominant_rdg_evidence_gaps")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

clip01 <- function(x) pmin(1, pmax(0, x))

safe_fread <- function(path) {
  if (!file.exists(path)) return(data.table())
  data.table::fread(path, showProgress = FALSE, nThread = 1)
}

has_text <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(trimws(x))
}

num_col <- function(dt, column, default = 0) {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- suppressWarnings(as.numeric(dt[[column]]))
  x[!is.finite(x)] <- default
  x
}

chr_col <- function(dt, column, default = "") {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- as.character(dt[[column]])
  x[is.na(x)] <- default
  x
}

bool_col <- function(dt, column, default = FALSE) {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[column]]
  if (is.logical(x)) {
    x[is.na(x)] <- default
    return(x)
  }
  y <- tolower(trimws(as.character(x)))
  fifelse(y %chin% c("true", "t", "1", "yes", "y"), TRUE,
          fifelse(y %chin% c("false", "f", "0", "no", "n"), FALSE, default))
}

text_has <- function(x, pattern) {
  grepl(pattern, as.character(x), ignore.case = TRUE, perl = TRUE)
}

warning_count <- function(x) {
  vapply(strsplit(as.character(x), "\\s*[|;]\\s*"), function(parts) {
    sum(nzchar(trimws(parts)))
  }, integer(1))
}

merge_new_columns <- function(x, y, by) {
  if (nrow(y) == 0) return(x)
  add <- setdiff(names(y), names(x))
  if (!length(add)) return(x)
  merge(x, y[, c(by, add), with = FALSE], by = by, all.x = TRUE,
        sort = FALSE)
}

candidate_file <- file.path(
  analysis_dir,
  "dominant_rdg_validation_dossiers",
  "dominant_rdg_validation_dossier_candidates.csv"
)
template_file <- file.path(
  analysis_dir,
  "dominant_rdg_validation_dossiers",
  "dominant_rdg_validation_dossier_manual_review_template.csv"
)
atlas_file <- file.path(
  analysis_dir,
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_cards.csv"
)
feature_file <- file.path(
  analysis_dir,
  "dominant_next_model",
  "dominant_next_model_feature_confidence.csv"
)

priority_out <- file.path(output_dir, "dominant_rdg_evidence_gap_priorities.csv")
gene_out <- file.path(output_dir, "dominant_rdg_evidence_gap_gene_summary.csv")
action_out <- file.path(output_dir, "dominant_rdg_evidence_gap_action_summary.csv")
negative_out <- file.path(output_dir, "dominant_rdg_evidence_gap_negative_controls.csv")
metrics_out <- file.path(output_dir, "dominant_rdg_evidence_gap_summary_metrics.csv")
report_out <- file.path(output_dir, "dominant_rdg_evidence_gap_report.md")

candidates <- safe_fread(candidate_file)
if (nrow(candidates) == 0) {
  stop("Missing or empty validation dossier candidates: ", candidate_file,
       call. = FALSE)
}

keys3 <- c("gene_symbol", "tx_id", "design_family")
if (!all(keys3 %in% names(candidates))) {
  stop("Candidate table is missing one of: ", paste(keys3, collapse = ", "),
       call. = FALSE)
}

atlas <- safe_fread(atlas_file)
if (nrow(atlas) > 0 && all(keys3 %in% names(atlas))) {
  atlas_keep <- c(
    keys3,
    "atlas_confidence_score", "atlas_browser_review_priority_score",
    "abs_atlas_clean_cds_pct_change", "weighted_clean_cds_pct_change",
    "n_designs", "n_cell_lines", "n_tissues", "supported_rows",
    "supported_studies", "postviral_atlas_condition_evidence_gate",
    "browser_validation_status", "browser_validation_score",
    "n_uorfs", "n_overlapping_uorfs", "has_overlapping_uorf",
    "n_measured_features", "n_measured_uorf_ouorf",
    "n_measured_internal_orf", "n_measured_nte_ntt",
    "transcript_qc_class", "cds_exon_qc_status", "cds_exon_qc_min_ratio",
    "primary_clean_cds_usage", "n_rdg_features", "manual_m_source_gene",
    "hidden_uorf_candidate", "warnings", "n_design_rows",
    "n_case_conditions", "sparse_n_onoff_rows",
    "n_features_with_signal_bump", "browser_review_priority_score",
    "pooling_proof_tier", "pooling_claim_language",
    "pooling_proof_component", "pooling_hazard_class",
    "pooling_hazard_score", "pooling_hazard_sources",
    "pooling_top_protocol_classes", "pooling_exact_context_fraction",
    "pooling_joint_strict_fraction",
    "pooling_any_start_site_enriched_protocol",
    "pooling_any_relaxed_context"
  )
  atlas_keep <- intersect(atlas_keep, names(atlas))
  candidates <- merge_new_columns(
    candidates,
    unique(atlas[, atlas_keep, with = FALSE], by = keys3),
    by = keys3
  )
}

features <- safe_fread(feature_file)
if (nrow(features) > 0 && all(c("gene_symbol", "tx_id") %in% names(features))) {
  features[, has_signal_counts := num_col(.SD, "total_feature_raw_counts") >= 30]
  features[, manual_feature := text_has(
    paste(chr_col(.SD, "translon_source"),
          chr_col(.SD, "translon_sources_merged")),
    "\\bM\\b|manual"
  )]
  features[, overlapping_feature := text_has(
    paste(chr_col(.SD, "category"), chr_col(.SD, "feature_type")),
    "overlapping|ouorf"
  )]
  features[, internal_feature := text_has(
    paste(chr_col(.SD, "feature_type"), chr_col(.SD, "category")),
    "internal|iorf"
  )]
  feature_summary <- features[, .(
    n_feature_confidence_rows = .N,
    n_signal_supported_features = sum(has_signal_counts, na.rm = TRUE),
    n_low_count_features = sum(!has_signal_counts, na.rm = TRUE),
    n_manual_features = sum(manual_feature, na.rm = TRUE),
    n_overlapping_features = sum(overlapping_feature, na.rm = TRUE),
    n_internal_features = sum(internal_feature, na.rm = TRUE),
    max_feature_raw_counts = max(num_col(.SD, "total_feature_raw_counts"),
                                 na.rm = TRUE)
  ), by = .(gene_symbol, tx_id)]
  feature_summary[!is.finite(max_feature_raw_counts),
                  max_feature_raw_counts := 0]
  candidates <- merge_new_columns(
    candidates,
    feature_summary,
    by = c("gene_symbol", "tx_id")
  )
}

template <- safe_fread(template_file)
if (nrow(template) > 0 && "review_packet_id" %in% names(template)) {
  template_keep <- intersect(
    c("review_packet_id", "validation_status", "browser_signal_grade",
      "branch_direction_verified", "clean_cds_direction_verified",
      "context_match_verified", "likely_failure_mode", "reviewer_notes",
      "reviewer", "review_date", "next_action"),
    names(template)
  )
  template <- template[, template_keep, with = FALSE]
  candidates <- merge_new_columns(candidates, template, by = "review_packet_id")
}

dt <- copy(candidates)

warning_text <- paste(
  chr_col(dt, "warning_items"),
  chr_col(dt, "consensus_warning_items"),
  chr_col(dt, "warnings"),
  chr_col(dt, "aux_audit_warning"),
  sep = "; "
)
dt[, warning_item_count := warning_count(warning_text)]
dt[, severe_warning := text_has(
  warning_text,
  "fragile|confound|heterogeneity|direction_not_stable|not_evaluable|low_support|lost_after_calibration|single_context|suspect_low_cds_exon|not_measurable"
)]

effect_pct <- pmax(
  abs(num_col(dt, "atlas_clean_cds_pct_change")),
  abs(num_col(dt, "abs_atlas_clean_cds_pct_change")),
  abs(num_col(dt, "half_normalization_cds_pct_change")),
  na.rm = TRUE
)
effect_pct[!is.finite(effect_pct)] <- 0
branch_abs <- pmax(
  abs(num_col(dt, "main_branch_delta")),
  abs(num_col(dt, "joint_top_branch_delta")),
  abs(num_col(dt, "hierarchical_joint_top_delta")),
  abs(num_col(dt, "dm_top_delta")),
  abs(num_col(dt, "top_shifted_weighted_relative_use_delta")),
  na.rm = TRUE
)
branch_abs[!is.finite(branch_abs)] <- 0

selected_primary <- bool_col(dt, "selected_for_primary_review")
manual_review_done <- has_text(chr_col(dt, "validation_status")) |
  has_text(chr_col(dt, "reviewer_notes")) |
  has_text(chr_col(dt, "browser_signal_grade"))
browser_supported <- text_has(
  paste(chr_col(dt, "browser_validation_label"),
        chr_col(dt, "browser_validation_status"),
        chr_col(dt, "browser_signal_grade")),
  "interesting|positive|supported|browser"
)

model_good <- text_has(chr_col(dt, "model_agreement_tier"),
                       "A_|B_|C_|browser|replicated|multi")
loo_good <- text_has(chr_col(dt, "loo_consensus_class"),
                     "support_robust|direction_robust")
transport_good <- text_has(chr_col(dt, "transport_consensus_class"),
                           "support_robust|direction_robust")
transport_bad <- text_has(chr_col(dt, "transport_consensus_class"),
                          "fragile|not_evaluable|disagreement")
loo_bad <- text_has(chr_col(dt, "loo_consensus_class"),
                    "fragile|not_evaluable|disagreement|non_consensus")

context_signal <- text_has(chr_col(dt, "context_interaction_gene_class"),
                           "direction_reversal|magnitude_modulated|specific_supported")
context_weak <- text_has(chr_col(dt, "context_interaction_gene_class"),
                         "context_evaluable_weak")
partial_signal <- text_has(chr_col(dt, "partial_pooling_gene_class"),
                           "supported|review|single_context|raw_signal|pooled")
aux_signal <- text_has(chr_col(dt, "aux_top_overlap_class"),
                       "supported|shift") |
  has_text(chr_col(dt, "aux_top_state"))
therapeutic_signal <- has_text(chr_col(dt, "therapeutic_gene_label")) |
  has_text(chr_col(dt, "perturbation_label"))
therapeutic_ready <- text_has(chr_col(dt, "perturbation_support_class"),
                              "quantitative|directional")
pooling_exact_proof <- chr_col(dt, "pooling_proof_tier") ==
  "A_exact_strict_project_matched_control"
pooling_review_lead <- chr_col(dt, "pooling_proof_tier") %chin% c(
  "C_pooled_recurrence_or_transport",
  "D_context_specific_or_partial_pooled"
)
pooling_hazard_score <- num_col(dt, "pooling_hazard_score")
pooling_protocol_hazard <- pooling_hazard_score >= 0.40 |
  chr_col(dt, "pooling_hazard_class") %chin% c(
    "high_pooling_confounding_hazard",
    "moderate_pooling_confounding_hazard"
  )
pooling_high_hazard <- pooling_hazard_score >= 0.65 |
  chr_col(dt, "pooling_hazard_class") == "high_pooling_confounding_hazard"
pooling_relaxed_context <- bool_col(dt, "pooling_any_relaxed_context") |
  num_col(dt, "pooling_exact_context_fraction", 1) < 0.75 |
  text_has(chr_col(dt, "pooling_hazard_sources"), "relaxed_context")
start_site_protocol <- bool_col(
  dt, "pooling_any_start_site_enriched_protocol"
) |
  text_has(
    paste(chr_col(dt, "pooling_hazard_sources"),
          chr_col(dt, "pooling_top_protocol_classes"),
          warning_text),
    "start_site|ltm|harringtonine"
  )

replication_ready <- clip01(
  0.35 * pmin(num_col(dt, "n_studies"), 4) / 4 +
    0.25 * pmin(num_col(dt, "strict_supported_studies"), 3) / 3 +
    0.20 * pmin(num_col(dt, "total_case_samples"), 150) / 150 +
    0.20 * pmin(num_col(dt, "total_control_samples"), 150) / 150
)
robustness_ready <- clip01(
  0.35 * as.numeric(model_good) +
    0.25 * as.numeric(loo_good) +
    0.25 * as.numeric(transport_good) +
    0.15 * clip01(num_col(dt, "consensus_score"))
)
measurement_ready <- clip01(
  0.20 * as.numeric(chr_col(dt, "transcript_qc_class") %chin%
                      c("canonical_qc_priority", "noncanonical_exon_qc_pass", "")) +
    0.20 * as.numeric(chr_col(dt, "cds_exon_qc_status") %chin%
                        c("pass", "")) +
    0.20 * pmin(num_col(dt, "n_measured_features", 1), 5) / 5 +
    0.15 * pmin(num_col(dt, "n_signal_supported_features", 1),
                pmax(1, num_col(dt, "n_feature_confidence_rows", 1))) /
      pmax(1, num_col(dt, "n_feature_confidence_rows", 1)) +
    0.15 * pmin(num_col(dt, "max_feature_raw_counts"), 300) / 300 +
    0.10 * as.numeric(!dt$severe_warning)
)
context_ready <- clip01(
  0.45 * as.numeric(!context_signal & !partial_signal) +
    0.10 * as.numeric(pooling_exact_proof & !pooling_protocol_hazard) +
    0.25 * as.numeric(text_has(chr_col(dt, "partial_pooling_gene_class"),
                               "supported|pooled")) +
    0.20 * as.numeric(transport_good) +
    0.10 * pmin(num_col(dt, "n_cell_lines") + num_col(dt, "n_tissues"), 8) / 8 -
    0.10 * as.numeric(pooling_high_hazard)
)

dt[, browser_need := clip01(
  (0.55 * as.numeric(selected_primary) +
     0.25 * clip01(num_col(dt, "validation_priority_score")) +
     0.20 * as.numeric(text_has(chr_col(dt, "review_priority_tier"),
                                "immediate|context|branch|buffering|anchor"))) *
    as.numeric(!manual_review_done & !browser_supported)
)]
dt[, measurement_need := clip01(1 - measurement_ready +
                                  0.20 * as.numeric(severe_warning))]
dt[, replication_need := clip01(1 - replication_ready)]
dt[, robustness_need := clip01(1 - robustness_ready +
                                 0.20 * as.numeric(loo_bad | transport_bad))]
dt[, context_need := clip01(
  0.55 * as.numeric(context_signal) +
    0.35 * as.numeric(partial_signal) +
    0.15 * as.numeric(context_weak) +
    0.15 * as.numeric(aux_signal) -
    0.25 * as.numeric(context_ready > 0.65)
)]
dt[, metadata_confounder_need := clip01(
  0.40 * as.numeric(aux_signal) +
    0.30 * as.numeric(text_has(chr_col(dt, "aux_audit_warning"),
                               "state|context|confound|driver")) +
    0.20 * as.numeric(context_signal) +
    0.20 * as.numeric(pooling_protocol_hazard) +
    0.10 * as.numeric(chr_col(dt, "design_family") %chin%
                        c("tumor_cancer_context", "viral_infection"))
)]
dt[, protocol_pooling_need := clip01(
  0.42 * as.numeric(pooling_protocol_hazard) +
    0.24 * as.numeric(start_site_protocol) +
    0.16 * as.numeric(pooling_relaxed_context) +
    0.12 * as.numeric(pooling_review_lead) +
    0.10 * clip01(pooling_hazard_score) +
    0.08 * as.numeric(text_has(chr_col(.SD, "pooling_hazard_sources"),
                                "multi_protocol")) -
    0.14 * as.numeric(pooling_exact_proof & !pooling_protocol_hazard)
)]
dt[, therapeutic_need := clip01(
  0.75 * as.numeric(therapeutic_signal & !therapeutic_ready) +
    0.25 * as.numeric(therapeutic_signal) -
    0.20 * as.numeric(manual_review_done)
)]

dt[, biological_value_score := clip01(
  0.26 * clip01(num_col(dt, "validation_priority_score")) +
    0.22 * clip01(effect_pct / 100) +
    0.18 * clip01(branch_abs / 0.35) +
    0.12 * as.numeric(chr_col(dt, "design_family") %chin%
                        c("viral_infection", "postviral_fatigue",
                          "ISR_ER_stress", "ISR / ER stress (ATF4 axis)")) +
    0.10 * as.numeric(therapeutic_signal) +
    0.08 * as.numeric(browser_supported) +
    0.04 * clip01(num_col(dt, "pooling_proof_component")) +
    0.04 * clip01(num_col(dt, "atlas_browser_review_priority_score"))
)]
dt[, evidence_gap_score := clip01(
  0.18 * browser_need +
    0.16 * measurement_need +
    0.16 * replication_need +
    0.14 * robustness_need +
    0.12 * context_need +
    0.10 * metadata_confounder_need +
    0.09 * protocol_pooling_need +
    0.05 * therapeutic_need
)]
dt[, counterfactual_readiness_score := round(100 * clip01(
  0.20 * clip01(num_col(dt, "consensus_score")) +
    0.18 * replication_ready +
    0.18 * robustness_ready +
    0.16 * measurement_ready +
    0.12 * context_ready +
    0.05 * as.numeric(pooling_exact_proof) -
    0.07 * clip01(pooling_hazard_score) +
    0.08 * as.numeric(browser_supported | manual_review_done) +
    0.08 * as.numeric(therapeutic_ready | !therapeutic_signal)
), 2)]
dt[, acquisition_priority_score := round(100 * clip01(
  (0.58 * biological_value_score + 0.42 * evidence_gap_score) *
    fifelse(manual_review_done, 0.55, 1)
), 2)]

dt[, negative_control_score := clip01(
  0.35 * (1 - clip01(effect_pct / 30)) +
    0.25 * measurement_ready +
    0.20 * replication_ready +
    0.10 * as.numeric(!severe_warning) +
    0.10 * as.numeric(!browser_supported & !therapeutic_signal)
)]

dt[, action_label := fifelse(
  measurement_need >= 0.62 & biological_value_score >= 0.25,
  "add_or_verify_translon",
  fifelse(
    browser_need >= 0.55 & biological_value_score >= 0.30,
    "browser_validate_now",
      fifelse(
        therapeutic_need >= 0.55 & counterfactual_readiness_score >= 45,
        "perturbation_followup",
        fifelse(
          protocol_pooling_need >= 0.55 & biological_value_score >= 0.25,
          "protocol_pooling_audit",
          fifelse(
        metadata_confounder_need >= 0.50 | context_need >= 0.55,
        "metadata_context_audit",
        fifelse(
          replication_need >= 0.55 | robustness_need >= 0.60,
          "add_context_replication",
          fifelse(
            negative_control_score >= 0.70 & effect_pct <= 15,
            "negative_control_review",
            "defer_low_information"
          )
        )
      )
    )
  )
  )
)]

dt[, claim_stage := fifelse(
  counterfactual_readiness_score >= 70 & biological_value_score >= 0.45,
  "near_counterfactual_candidate",
  fifelse(
    acquisition_priority_score >= 45,
    "high_value_acquisition",
    fifelse(
      counterfactual_readiness_score >= 55,
      "stable_review_grade",
      fifelse(evidence_gap_score >= 0.55, "needs_more_evidence", "background")
    )
  )
)]

reason_flags <- list(
  "needs browser label" = dt$browser_need >= 0.55,
  "translon/measurement risk" = dt$measurement_need >= 0.62,
  "needs more replicated context" = dt$replication_need >= 0.55,
  "LOO/transport robustness gap" = dt$robustness_need >= 0.60,
  "context-specific branch behavior" = dt$context_need >= 0.55,
  "metadata/aux-state confounder" = dt$metadata_confounder_need >= 0.50,
  "pooling/protocol confounder" = dt$protocol_pooling_need >= 0.55,
  "therapeutic estimate needs validation" = dt$therapeutic_need >= 0.55
)
dt[, review_reason := vapply(seq_len(.N), function(i) {
  reasons <- names(reason_flags)[vapply(reason_flags, function(flag) {
    isTRUE(flag[[i]])
  }, logical(1))]
  if (!length(reasons)) "low current acquisition value" else
    paste(reasons, collapse = "; ")
}, character(1))]

priority_cols <- intersect(c(
  "gene_symbol", "tx_id", "design_family", "action_label", "claim_stage",
  "acquisition_priority_score", "counterfactual_readiness_score",
  "evidence_gap_score", "biological_value_score",
  "browser_need", "measurement_need", "replication_need", "robustness_need",
  "context_need", "metadata_confounder_need", "protocol_pooling_need",
  "therapeutic_need",
  "negative_control_score", "review_reason",
  "review_packet_id", "selected_for_primary_review", "review_priority_tier",
  "validation_priority_score", "consensus_rank", "consensus_class",
  "consensus_score", "atlas_evidence_class", "atlas_clean_cds_pct_change",
  "main_branch_class", "main_branch_delta", "top_shifted_feature_id",
  "top_shifted_feature_type", "top_shifted_feature_family",
  "model_agreement_tier", "loo_consensus_class",
  "transport_consensus_class", "context_interaction_gene_class",
  "partial_pooling_gene_class", "browser_validation_label",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_proof_component", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "pooling_top_protocol_classes", "pooling_exact_context_fraction",
  "pooling_joint_strict_fraction",
  "pooling_any_start_site_enriched_protocol",
  "pooling_any_relaxed_context",
  "browser_signal_grade", "validation_status", "aux_top_state",
  "aux_top_overlap_class", "therapeutic_gene_label", "perturbation_label",
  "perturbation_support_class", "n_studies", "strict_supported_studies",
  "total_case_samples", "total_control_samples", "best_case_condition",
  "best_control_conditions", "best_study", "transcript_qc_class",
  "cds_exon_qc_status", "n_measured_features",
  "n_signal_supported_features", "n_low_count_features", "n_manual_features",
  "n_overlapping_features", "n_internal_features", "warning_item_count",
  "warning_items", "warnings"
), names(dt))

priority <- dt[, priority_cols, with = FALSE]
setorder(priority, -acquisition_priority_score, -counterfactual_readiness_score,
         gene_symbol, design_family)
priority[, acquisition_rank := .I]
setcolorder(priority, c("acquisition_rank",
                        setdiff(names(priority), "acquisition_rank")))
fwrite(priority, priority_out)

gene_summary <- priority[, {
  best <- which.max(acquisition_priority_score)
  action_tab <- sort(table(action_label), decreasing = TRUE)
  .(
    n_designs = .N,
    max_acquisition_priority_score = max(acquisition_priority_score,
                                         na.rm = TRUE),
    mean_counterfactual_readiness_score = round(
      mean(counterfactual_readiness_score, na.rm = TRUE), 2
    ),
    max_evidence_gap_score = round(max(evidence_gap_score, na.rm = TRUE), 3),
    best_design_family = design_family[best],
    best_action_label = action_label[best],
    best_claim_stage = claim_stage[best],
    n_browser_validate = sum(action_label == "browser_validate_now"),
    n_translon_review = sum(action_label == "add_or_verify_translon"),
    n_context_audit = sum(action_label == "metadata_context_audit"),
    n_protocol_pooling_audit = sum(action_label == "protocol_pooling_audit"),
    n_replication_needed = sum(action_label == "add_context_replication"),
    n_negative_controls = sum(action_label == "negative_control_review"),
    dominant_actions = paste(names(action_tab)[seq_len(min(3, length(action_tab)))],
                             collapse = "; ")
  )
}, by = gene_symbol]
setorder(gene_summary, -max_acquisition_priority_score, gene_symbol)
fwrite(gene_summary, gene_out)

action_summary <- priority[, {
  top_rows <- .SD[order(-acquisition_priority_score)][seq_len(min(.N, 12))]
  .(
    n_rows = .N,
    n_genes = uniqueN(gene_symbol),
    median_acquisition_priority_score = round(
      median(acquisition_priority_score, na.rm = TRUE), 2
    ),
    max_acquisition_priority_score = max(acquisition_priority_score,
                                         na.rm = TRUE),
    median_counterfactual_readiness_score = round(
      median(counterfactual_readiness_score, na.rm = TRUE), 2
    ),
    top_gene_designs = paste(
      paste(top_rows$gene_symbol, top_rows$design_family, sep = " | "),
      collapse = "; "
    )
  )
}, by = action_label]
setorder(action_summary, -max_acquisition_priority_score)
fwrite(action_summary, action_out)

negative_controls <- priority[
  negative_control_score >= 0.70 &
    abs(atlas_clean_cds_pct_change) <= 15 &
    n_studies >= 2 &
    total_case_samples >= 30 &
    total_control_samples >= 30
]
setorder(negative_controls, -negative_control_score,
         -counterfactual_readiness_score, gene_symbol)
negative_controls <- negative_controls[seq_len(min(.N, 250))]
fwrite(negative_controls, negative_out)

metrics <- data.table(
  metric = c(
    "candidate_rows", "genes", "design_families",
    "high_value_acquisition_rows", "near_counterfactual_candidate_rows",
    "browser_validate_rows", "translon_review_rows",
    "context_audit_rows", "protocol_pooling_audit_rows",
    "protocol_pooling_need_rows", "pooling_exact_proof_rows",
    "pooling_hazard_rows", "pooling_start_site_protocol_rows",
    "replication_needed_rows", "perturbation_followup_rows",
    "negative_control_rows",
    "median_acquisition_priority", "median_counterfactual_readiness"
  ),
  value = c(
    nrow(priority),
    uniqueN(priority$gene_symbol),
    uniqueN(priority$design_family),
    sum(priority$claim_stage == "high_value_acquisition"),
    sum(priority$claim_stage == "near_counterfactual_candidate"),
    sum(priority$action_label == "browser_validate_now"),
    sum(priority$action_label == "add_or_verify_translon"),
    sum(priority$action_label == "metadata_context_audit"),
    sum(priority$action_label == "protocol_pooling_audit"),
    sum(priority$protocol_pooling_need >= 0.55, na.rm = TRUE),
    sum(pooling_exact_proof, na.rm = TRUE),
    sum(pooling_protocol_hazard, na.rm = TRUE),
    sum(start_site_protocol, na.rm = TRUE),
    sum(priority$action_label == "add_context_replication"),
    sum(priority$action_label == "perturbation_followup"),
    nrow(negative_controls),
    round(median(priority$acquisition_priority_score, na.rm = TRUE), 2),
    round(median(priority$counterfactual_readiness_score, na.rm = TRUE), 2)
  )
)
fwrite(metrics, metrics_out)

priority[, plot_label := fifelse(
  acquisition_rank <= 24,
  paste(gene_symbol, design_family, sep = "\n"),
  ""
)]
plot_dt <- priority[acquisition_rank <= 120]
priority_plot <- ggplot(
  plot_dt,
  aes(
    x = counterfactual_readiness_score,
    y = acquisition_priority_score,
    color = action_label,
    size = biological_value_score,
    label = plot_label
  )
) +
  geom_point(alpha = 0.82) +
  geom_text(check_overlap = TRUE, size = 2.4, vjust = -0.7,
            show.legend = FALSE) +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(
    title = "RDG Atlas Evidence Gaps",
    subtitle = "High priority means the row is biologically valuable and missing a concrete evidence layer; readiness estimates current counterfactual usability.",
    x = "Counterfactual readiness score",
    y = "Evidence acquisition priority",
    color = "Next action"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(figure_dir, "evidence_gap_priority_map.png"),
       priority_plot, width = 10, height = 7, dpi = 160)
ggsave(file.path(figure_dir, "evidence_gap_priority_map.pdf"),
       priority_plot, width = 10, height = 7)
if (exists("dominant_save_ggplotly") && requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    priority_plot,
    file.path(figure_dir, "evidence_gap_priority_map.html"),
    title = "evidence_gap_priority_map"
  )
}

action_plot <- ggplot(
  action_summary,
  aes(x = reorder(action_label, n_rows), y = n_rows, fill = action_label)
) +
  geom_col(width = 0.72, show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "RDG Evidence-Gap Workload",
    subtitle = "Rows are assigned one primary next action for triage.",
    x = NULL,
    y = "Atlas rows"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(figure_dir, "evidence_gap_action_counts.png"),
       action_plot, width = 8, height = 4.8, dpi = 160)
ggsave(file.path(figure_dir, "evidence_gap_action_counts.pdf"),
       action_plot, width = 8, height = 4.8)

top_genes <- head(gene_summary$gene_symbol, 36)
heat <- priority[gene_symbol %chin% top_genes, .(
  max_priority = max(acquisition_priority_score, na.rm = TRUE),
  n_rows = .N
), by = .(gene_symbol, action_label)]
heat[, gene_symbol := factor(gene_symbol, levels = rev(top_genes))]
heat_plot <- ggplot(
  heat,
  aes(x = action_label, y = gene_symbol, fill = max_priority)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = n_rows), size = 2.2) +
  scale_fill_gradient(low = "#f2f2f2", high = "#8b0000",
                      name = "Max priority") +
  labs(
    title = "Action Map For Top RDG Genes",
    subtitle = "Numbers are rows per action; color is the strongest acquisition priority for that gene/action.",
    x = "Next action",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
ggsave(file.path(figure_dir, "evidence_gap_gene_action_heatmap.png"),
       heat_plot, width = 10.5, height = 8.5, dpi = 160)
ggsave(file.path(figure_dir, "evidence_gap_gene_action_heatmap.pdf"),
       heat_plot, width = 10.5, height = 8.5)

top_priority <- priority[seq_len(min(.N, 15))]
report_lines <- c(
  "# RDG Evidence-Gap Prioritization",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Purpose",
  "",
  paste(
    "This layer asks what the atlas is missing before a row becomes",
    "actionable: a browser label, a translon/measurement correction,",
    "more replicated context, a metadata-confounder audit, a protocol/pooling",
    "audit, perturbation validation, or a negative control."
  ),
  "",
  "## Summary",
  "",
  paste0("- Rows ranked: ", nrow(priority)),
  paste0("- Genes ranked: ", uniqueN(priority$gene_symbol)),
  paste0("- High-value acquisition rows: ",
         sum(priority$claim_stage == "high_value_acquisition")),
  paste0("- Near-counterfactual candidates: ",
         sum(priority$claim_stage == "near_counterfactual_candidate")),
  paste0("- Browser-validation rows: ",
         sum(priority$action_label == "browser_validate_now")),
  paste0("- Translon/measurement-review rows: ",
         sum(priority$action_label == "add_or_verify_translon")),
  paste0("- Protocol/pooling-audit rows: ",
         sum(priority$action_label == "protocol_pooling_audit")),
  paste0("- Rows with protocol/pooling need >= 0.55: ",
         sum(priority$protocol_pooling_need >= 0.55, na.rm = TRUE)),
  paste0("- Negative-control rows: ", nrow(negative_controls)),
  "",
  "## Top Acquisition Rows",
  "",
  paste0(
    seq_len(nrow(top_priority)), ". ",
    top_priority$gene_symbol, " | ", top_priority$design_family,
    " | action=", top_priority$action_label,
    " | priority=", top_priority$acquisition_priority_score,
    " | readiness=", top_priority$counterfactual_readiness_score,
    " | reason=", top_priority$review_reason
  ),
  "",
  "## Interpretation",
  "",
  paste(
    "The atlas now has many evidence tables, but its bottleneck is not",
    "another broad screen. The useful next data are labels and falsification",
    "points: validated positives, validated negatives, translon corrections,",
    "and context-replicated rows. This file turns that into a ranked queue."
  )
)
writeLines(report_lines, report_out, useBytes = TRUE)

message("Wrote evidence-gap priorities: ", priority_out)
message("Wrote evidence-gap report: ", report_out)
