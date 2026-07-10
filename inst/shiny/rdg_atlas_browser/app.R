#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(shiny)
  library(DT)
  library(plotly)
})

find_analysis_dir <- function(start = getwd()) {
  env_root <- Sys.getenv("DOMINANT_CELL_STATES_DIR", unset = "")
  if (nzchar(env_root)) {
    env_root <- normalizePath(env_root, mustWork = FALSE)
    if (file.exists(file.path(env_root, "results", "dominant_rdg_atlas",
                              "dominant_rdg_atlas_cards.csv"))) {
      return(env_root)
    }
  }

  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        (file.exists(file.path(here, "dominant_rdg_atlas",
                               "dominant_rdg_atlas_cards.csv")) ||
           file.exists(file.path(here, "results", "dominant_rdg_atlas",
                                 "dominant_rdg_atlas_cards.csv")))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (file.exists(file.path(candidate, "dominant_rdg_atlas",
                              "dominant_rdg_atlas_cards.csv")) ||
        file.exists(file.path(candidate, "results", "dominant_rdg_atlas",
                              "dominant_rdg_atlas_cards.csv"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states/dominant_rdg_atlas cards from: ",
       start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
analysis_dir_norm <- normalizePath(analysis_dir, mustWork = TRUE)
analysis_results_dir <- file.path(analysis_dir_norm, "results")
path_exists_or_link <- function(path) {
  link <- Sys.readlink(path)
  file.exists(path) || (!is.na(link) && nzchar(link))
}
analysis_path <- function(..., create_parent = FALSE) {
  rel <- file.path(...)
  legacy <- file.path(analysis_dir_norm, rel)
  result <- file.path(analysis_results_dir, rel)
  if (path_exists_or_link(legacy)) {
    return(legacy)
  }
  if (path_exists_or_link(result) || isTRUE(create_parent)) {
    if (isTRUE(create_parent)) {
      dir.create(dirname(result), recursive = TRUE, showWarnings = FALSE)
    }
    return(result)
  }
  legacy
}
rdg_file_resource_prefix <- "rdg_atlas_files"
shiny::addResourcePath(rdg_file_resource_prefix, analysis_dir_norm)
cards_file <- analysis_path("dominant_rdg_atlas",
                            "dominant_rdg_atlas_cards.csv")
cards <- fread(cards_file, showProgress = FALSE)
cards[, card_key := paste(gene_symbol, tx_id, design_family, sep = "\r")]
branch_contexts_file <- analysis_path(
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_branch_contexts.csv"
)
joint_rows_file <- analysis_path(
  "dominant_rdg_joint_branch_allocation",
  "dominant_rdg_joint_branch_allocation_rows.csv"
)
branch_contexts <- if (file.exists(branch_contexts_file)) {
  fread(branch_contexts_file, showProgress = FALSE)
} else if (file.exists(joint_rows_file)) {
  fread(joint_rows_file, showProgress = FALSE)
} else {
  data.table()
}
study_evidence_file <- file.path(
  analysis_path("dominant_rdg_model_agreement_loo"),
  "dominant_rdg_model_agreement_study_evidence.csv"
)
study_evidence <- if (file.exists(study_evidence_file)) {
  fread(study_evidence_file, showProgress = FALSE)
} else {
  data.table()
}
transport_evidence_file <- file.path(
  analysis_path("dominant_rdg_model_agreement_transport"),
  "dominant_rdg_model_agreement_transport_gene_axis_summary.csv"
)
transport_evidence <- if (file.exists(transport_evidence_file)) {
  fread(transport_evidence_file, showProgress = FALSE)
} else {
  data.table()
}
context_interaction_file <- file.path(
  analysis_path("dominant_rdg_context_interactions"),
  "dominant_rdg_context_interaction_branch_context_effects.csv"
)
context_interactions <- if (file.exists(context_interaction_file)) {
  fread(context_interaction_file, showProgress = FALSE)
} else {
  data.table()
}
partial_pooling_file <- file.path(
  analysis_path("dominant_rdg_context_partial_pooling"),
  "dominant_rdg_context_partial_pooling_rows.csv"
)
partial_pooling <- if (file.exists(partial_pooling_file)) {
  fread(partial_pooling_file, showProgress = FALSE)
} else {
  data.table()
}
read_analysis_csv <- function(...) {
  path <- analysis_path(...)
  if (file.exists(path)) fread(path, showProgress = FALSE) else data.table()
}
review_batch_template <- read_analysis_csv(
  "dominant_rdg_review_links",
  "dominant_rdg_review_batch_template_with_links.csv"
)
validation_label_gap_queue <- read_analysis_csv(
  "dominant_rdg_validation_label_bridge",
  "dominant_rdg_validation_label_gap_queue.csv"
)
validation_label_balance <- read_analysis_csv(
  "dominant_rdg_validation_label_bridge",
  "dominant_rdg_validation_label_balance.csv"
)
validation_label_set <- read_analysis_csv(
  "dominant_rdg_validation_label_bridge",
  "dominant_rdg_validation_label_set.csv"
)
review_feedback_labels <- read_analysis_csv(
  "dominant_rdg_review_feedback",
  "dominant_rdg_review_feedback_labels.csv"
)
review_feedback_candidates <- read_analysis_csv(
  "dominant_rdg_review_feedback",
  "dominant_rdg_review_feedback_browser_validation_notes_candidate.csv"
)
review_feedback_summary <- read_analysis_csv(
  "dominant_rdg_review_feedback",
  "dominant_rdg_review_feedback_summary.csv"
)
review_feedback_batch_summary <- read_analysis_csv(
  "dominant_rdg_review_feedback",
  "dominant_rdg_review_feedback_batch_summary.csv"
)
evidence_gap_priorities <- read_analysis_csv(
  "dominant_rdg_evidence_gaps",
  "dominant_rdg_evidence_gap_priorities.csv"
)
qualified_inspection_manifest <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_inspection_manifest.csv"
)
qualified_review_priority <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_review_priority.csv"
)
qualified_profile_metrics <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_profile_metrics.csv"
)
qualified_review_flags <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_review_flags.csv"
)
qualified_replicate_summary <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_replicate_summary.csv"
)
qualified_replicate_effects <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_replicate_effects.csv"
)
qualified_sequence_features <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_sequence_features.csv"
)
qualified_feature_counts <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_feature_counts.csv"
)
qualified_feature_contrasts <- read_analysis_csv(
  "dominant_rdg_qualified_inspection",
  "qualified_feature_contrasts.csv"
)
sdrive_experiment_status <- read_analysis_csv(
  "dominant_rdg_sdrive_readiness",
  "sdrive_experiment_status.csv"
)
sdrive_fst_page_readiness <- read_analysis_csv(
  "dominant_rdg_sdrive_readiness",
  "sdrive_fst_page_readiness.csv"
)
sdrive_transcript_readiness <- read_analysis_csv(
  "dominant_rdg_sdrive_readiness",
  "sdrive_transcript_readiness.csv"
)
sdrive_readiness_summary_metrics <- read_analysis_csv(
  "dominant_rdg_sdrive_readiness",
  "sdrive_readiness_summary_metrics.csv"
)
sdrive_readiness_class_counts <- read_analysis_csv(
  "dominant_rdg_sdrive_readiness",
  "sdrive_readiness_class_counts.csv"
)
review_label_drafts_file <- file.path(
  analysis_path(
    "dominant_rdg_review_feedback",
    "dominant_rdg_review_label_drafts.csv",
    create_parent = TRUE
  )
)
review_label_manual_cols <- c(
  "manual_review_status",
  "browser_signal_grade",
  "branch_direction_verified",
  "clean_cds_direction_verified",
  "context_match_verified",
  "translon_structure_status",
  "metadata_confounder_status",
  "negative_control_result",
  "reviewer_notes",
  "reviewer",
  "review_date",
  "next_action"
)
review_label_draft_cols <- c(
  "global_review_order", "batch_id", "batch_name", "batch_review_order",
  "batch_role", "gene_symbol", "tx_id", "design_family",
  "review_packet_id", "action_label", "claim_stage",
  "review_question", "review_reason", "observatory_link_note",
  "target_label_class", "label_gap_class", "label_acquisition_score",
  "protocol_pooling_need", "pooling_proof_tier", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  review_label_manual_cols,
  "review_label_class_preview", "draft_saved_at"
)
empty_review_label_drafts <- function() {
  as.data.table(setNames(
    replicate(length(review_label_draft_cols), character(0), simplify = FALSE),
    review_label_draft_cols
  ))
}
read_review_label_drafts <- function() {
  if (!file.exists(review_label_drafts_file)) return(empty_review_label_drafts())
  out <- fread(review_label_drafts_file, showProgress = FALSE)
  for (col in setdiff(review_label_draft_cols, names(out))) out[, (col) := ""]
  out <- out[, ..review_label_draft_cols]
  out[, global_review_order := suppressWarnings(as.integer(global_review_order))]
  out
}
for (column in c("atlas_browser_review_priority_score",
                 "browser_validation_score",
                 "browser_validation_status",
                 "browser_validation_label",
                 "browser_validation_notes",
                 "branch_usage_top_score",
                 "branch_usage_top_support_class",
                 "branch_usage_top_feature_id",
                 "branch_usage_top_feature_class",
                 "branch_usage_top_context",
                 "grouped_branch_usage_top_feature_id",
                 "grouped_branch_usage_top_feature_class",
                 "grouped_branch_usage_top_support_class",
                 "grouped_branch_usage_top_weighted_usage_delta",
                 "grouped_branch_usage_top_context",
                 "count_aware_top_feature_id",
                 "count_aware_top_feature_class",
                 "count_aware_top_support_class",
                 "count_aware_top_weighted_usage_delta",
                 "count_aware_top_context",
                 "joint_top_review_class",
                 "joint_top_branch_class",
                 "joint_top_branch_delta",
                 "joint_top_context",
                 "hierarchical_joint_top_support_class",
                 "hierarchical_joint_top_branch_class",
                 "hierarchical_joint_top_delta",
                 "hierarchical_joint_top_priority",
                 "hierarchical_joint_top_context",
                 "dm_top_support_class",
                 "dm_top_branch_class",
                 "dm_top_delta",
                 "dm_top_priority",
                 "dm_top_context",
                 "model_agreement_class",
                 "model_agreement_tier",
                 "model_agreement_score",
                 "model_agreement_review_priority",
                 "model_support_signature",
                 "model_disagreement_flag",
                 "model_disagreement_flags",
                 "model_agreement_top_context",
                 "loo_agreement_class",
                 "loo_consensus_class",
                 "loo_stability_score",
                 "loo_review_priority",
                 "loo_relevant_supported_branches",
                 "loo_support_robust_branches",
                 "loo_direction_robust_support_fragile_branches",
                 "loo_direction_fragile_branches",
                 "loo_not_evaluable_branches",
                 "loo_min_direction_fraction",
                 "loo_min_support_retention",
                 "loo_min_model_studies",
                 "loo_support_robust_branch_names",
                 "loo_fragile_branch_names",
                 "loo_worst_omitted_studies",
                 "transport_agreement_class",
                 "transport_consensus_class",
                 "transport_stability_score",
                 "transport_review_priority",
                 "transport_min_direction_fraction",
                 "transport_min_support_retention",
                 "transport_min_remaining_studies",
                 "transport_coverage_fraction",
                 "transport_axis_coverage_fraction",
                 "transport_evaluable_axis_names",
                 "transport_support_robust_axis_names",
                 "transport_fragile_axis_names",
                 "transport_worst_omitted_groups",
                 "transport_worst_omitted_studies",
                 "context_interaction_gene_class",
                 "context_interaction_top_class",
                 "context_interaction_top_score",
                 "context_interaction_review_priority",
                 "context_interaction_top_axis",
                 "context_interaction_top_value",
                 "context_interaction_top_branch_class",
                 "context_interaction_context_delta",
                 "context_interaction_complement_delta",
                 "context_interaction_delta",
                 "context_interaction_delta_q",
                 "context_interaction_context_studies",
                 "context_interaction_complement_studies",
                 "context_interaction_context_study_names",
                 "context_interaction_complement_study_names",
                 "context_interaction_summary",
                 "partial_pooling_gene_class",
                 "partial_pooling_top_class",
                 "partial_pooling_top_survival_class",
                 "partial_pooling_review_priority",
                 "partial_pooling_component",
                 "partial_pooling_top_axis",
                 "partial_pooling_top_value",
                 "partial_pooling_top_branch_class",
                 "partial_pooling_top_delta",
                 "partial_pooling_top_q",
                 "partial_pooling_top_weight",
                 "partial_pooling_top_shrinkage_fraction",
                 "partial_pooling_rows",
                 "partial_pooling_pooled_supported_rows",
                 "partial_pooling_pooled_review_rows",
                 "partial_pooling_raw_supported_rows",
                 "partial_pooling_raw_survived_rows",
                 "partial_pooling_raw_shrunk_rows",
                 "partial_pooling_single_context_raw_signal_rows",
                 "partial_pooling_pooled_emergent_rows",
                 "partial_pooling_top_context_summary",
                 "partial_pooling_context_values",
                 "pooling_proof_tier",
                 "pooling_claim_language",
                 "pooling_proof_component",
                 "pooling_hazard_class",
                 "pooling_hazard_score",
                 "pooling_hazard_sources",
                 "pooling_context_rows",
                 "pooling_exact_context_rows",
                 "pooling_relaxed_context_rows",
                 "pooling_joint_strict_rows",
                 "pooling_joint_evaluable_rows",
                 "pooling_exact_context_fraction",
                 "pooling_joint_strict_fraction",
                 "pooling_joint_evaluable_fraction",
                 "pooling_context_wt_like_control_rows",
                 "pooling_context_strong_signal_rows",
                 "pooling_top_protocol_classes",
                 "pooling_top_inhibitors",
                 "pooling_top_tissues",
                 "pooling_top_cell_lines",
                 "pooling_top_case_conditions",
                 "pooling_top_controls",
                 "pooling_any_start_site_enriched_protocol",
                 "pooling_any_relaxed_context",
                 "hierarchical_evidence_class",
                 "hierarchical_priority_score",
                 "hierarchical_clean_cds_pct_change",
                 "hierarchical_prior_source",
                 "hierarchical_evidence_weight",
                 "context_calibration_status",
                 "context_model_status",
                 "context_model_rmse_delta_log2fc",
                 "context_model_fraction_improved",
                 "context_model_worst_context")) {
  if (!column %in% names(cards)) cards[, (column) := NA]
}
cards[!is.finite(atlas_browser_review_priority_score),
      atlas_browser_review_priority_score := atlas_confidence_score]
setorder(cards, -atlas_browser_review_priority_score, -atlas_confidence_score)

`%||%` <- function(x, y) if (is.null(x)) y else x

rdg_missing_tokens <- c(
  "", "NA", "N/A", "MISSING", "UNKNOWN", "NONE", "NULL", "NOT_AVAILABLE"
)

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

is_informative_text <- function(x) {
  x <- clean_text(x)
  nzchar(x) & !(toupper(x) %chin% rdg_missing_tokens)
}

is_scalar_informative <- function(x) {
  length(x) > 0 && isTRUE(is_informative_text(x[1]))
}

scalar_text <- function(x) {
  if (!length(x)) return("")
  clean_text(x[1])
}

split_condition_values <- function(x) {
  x <- scalar_text(x)
  if (!is_scalar_informative(x)) return(character())
  vals <- unlist(strsplit(x, "\\s*(;|,|\\|)\\s*", perl = TRUE),
                 use.names = FALSE)
  vals <- clean_text(vals)
  unique(vals[is_informative_text(vals)])
}

safe_choices <- function(x) sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
fmt_pct <- function(x) {
  ifelse(is.na(x), "NA", paste0(round(x, 1), "%"))
}
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
}
fmt_warning_hover <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- "none"
  gsub("\\s*;\\s*", "<br>", x)
}
fmt_pp <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(round(100 * x, digits), " pp"))
}
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}
safe_first_text <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  x[1]
}
available_card_choices <- function(choices) {
  choices[unname(choices) %chin% names(cards)]
}
choice_label <- function(choices, value) {
  idx <- match(value, unname(choices))
  if (is.na(idx)) return(value)
  names(choices)[idx]
}
safe_choice_value <- function(value, choices, default) {
  if (!default %chin% unname(choices)) default <- unname(choices)[1]
  value <- as.character(value %||% default)[1]
  if (!is_scalar_informative(value) || !(value %chin% unname(choices))) {
    return(default)
  }
  value
}
scale_marker_size <- function(x, fixed = 8, min_size = 6, max_size = 18) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (!any(is.finite(x))) return(rep(fixed, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) {
    out <- rep((min_size + max_size) / 2, length(x))
  } else {
    out <- min_size + (x - rng[1]) / diff(rng) * (max_size - min_size)
  }
  out[!is.finite(out)] <- min_size
  out
}
effect_y_choices <- available_card_choices(c(
  "Atlas confidence" = "atlas_confidence_score",
  "Review priority" = "atlas_browser_review_priority_score",
  "Strict studies" = "strict_supported_studies",
  "Strict designs" = "strict_supported_designs",
  "Heterogeneity I2" = "atlas_i2",
  "Sparse on/off score" = "sparse_onoff_score",
  "Joint branch L1 shift" = "joint_top_l1_shift",
  "Joint branch delta" = "joint_top_branch_delta",
  "Model agreement score" = "model_agreement_score",
  "LOO stability" = "loo_stability_score",
  "Transport stability" = "transport_stability_score",
  "Context interaction delta" = "context_interaction_delta",
  "Partial-pooling delta" = "partial_pooling_top_delta",
  "Partial-pooling q" = "partial_pooling_top_q",
  "Pooling proof" = "pooling_proof_component",
  "Pooling hazard" = "pooling_hazard_score",
  "Hierarchical priority" = "hierarchical_priority_score",
  "Residual heterogeneity" = "residual_heterogeneity_score"
))
effect_color_choices <- available_card_choices(c(
  "Atlas evidence" = "atlas_evidence_class",
  "Model agreement" = "model_agreement_tier",
  "Study omission" = "loo_consensus_class",
  "Grouped omission" = "transport_consensus_class",
  "Context interaction" = "context_interaction_gene_class",
  "Partial pooling" = "partial_pooling_gene_class",
  "Pooling proof tier" = "pooling_proof_tier",
  "Pooling hazard" = "pooling_hazard_class",
  "Hierarchy" = "hierarchical_evidence_class",
  "Residual risk" = "residual_risk_class",
  "Browser validation" = "browser_validation_status",
  "Design family" = "design_family"
))
effect_size_choices <- c(
  "Fixed" = "none",
  available_card_choices(c(
    "Review priority" = "atlas_browser_review_priority_score",
    "Strict studies" = "strict_supported_studies",
    "Strict designs" = "strict_supported_designs",
    "Model agreement score" = "model_agreement_score",
    "LOO stability" = "loo_stability_score",
    "Transport stability" = "transport_stability_score",
    "Joint branch L1 shift" = "joint_top_l1_shift",
    "Joint branch delta magnitude" = "joint_top_branch_delta",
    "Context interaction magnitude" = "context_interaction_delta",
    "Partial-pooling magnitude" = "partial_pooling_top_delta",
    "Pooling proof" = "pooling_proof_component",
    "Pooling hazard" = "pooling_hazard_score",
    "Hierarchical priority" = "hierarchical_priority_score"
  ))
)
effect_palette <- c(
  "#0b5cab", "#b45309", "#15803d", "#be123c", "#0891b2", "#7c3aed",
  "#475569", "#c2410c", "#047857", "#a21caf", "#4d7c0f", "#334155"
)

url_encode_path <- function(path) {
  path <- gsub("\\\\", "/", path)
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  paste(vapply(parts, utils::URLencode, character(1), reserved = TRUE),
        collapse = "/")
}

make_selected_file_link <- function(path_value) {
  raw_path <- scalar_text(path_value)
  if (!is_scalar_informative(raw_path)) return(NULL)
  raw_path <- gsub("\\\\", "/", raw_path)
  abs_path <- if (grepl("^/", raw_path)) {
    normalizePath(raw_path, mustWork = FALSE)
  } else {
    normalizePath(analysis_path(raw_path), mustWork = FALSE)
  }
  root_prefix <- paste0(analysis_dir_norm, .Platform$file.sep)
  inside_analysis <- identical(abs_path, analysis_dir_norm) ||
    startsWith(abs_path, root_prefix)
  if (!inside_analysis) {
    return(list(abs_path = abs_path, href = NA_character_,
                exists = file.exists(abs_path)))
  }
  rel_path <- if (identical(abs_path, analysis_dir_norm)) {
    ""
  } else {
    substring(abs_path, nchar(root_prefix) + 1L)
  }
  list(
    abs_path = abs_path,
    href = paste0(rdg_file_resource_prefix, "/", url_encode_path(rel_path)),
    exists = file.exists(abs_path)
  )
}

dt_col <- function(dt, column, default = NA) {
  if (column %in% names(dt)) return(dt[[column]])
  if (length(default) == nrow(dt)) return(default)
  rep(default, nrow(dt))
}

add_card_key <- function(dt) {
  if (nrow(dt) &&
      all(c("gene_symbol", "tx_id", "design_family") %chin% names(dt))) {
    dt[, card_key := paste(gene_symbol, tx_id, design_family, sep = "\r")]
  }
  dt
}

review_task_key_cols <- c(
  "gene_symbol", "tx_id", "design_family", "batch_id", "batch_role"
)

make_review_task_lookup <- function() {
  if (!nrow(review_batch_template)) {
    return(data.table(
      gene_symbol = character(),
      tx_id = character(),
      design_family = character(),
      batch_id = character(),
      batch_role = character(),
      active_review_order = integer()
    ))
  }
  lookup <- copy(review_batch_template)
  for (column in review_task_key_cols) {
    if (!column %in% names(lookup)) lookup[, (column) := ""]
    lookup[, (column) := clean_text(get(column))]
  }
  keep <- intersect(c(review_task_key_cols, "global_review_order"),
                    names(lookup))
  lookup <- lookup[, ..keep]
  lookup[, active_review_order :=
           suppressWarnings(as.integer(global_review_order))]
  lookup[, global_review_order := NULL]
  unique(lookup, by = review_task_key_cols)
}

review_task_lookup <- make_review_task_lookup()

attach_review_task_lookup <- function(dt) {
  if (!nrow(dt)) return(dt)
  dt <- copy(dt)
  for (column in review_task_key_cols) {
    if (!column %in% names(dt)) dt[, (column) := ""]
    dt[, (column) := clean_text(get(column))]
  }
  if (!nrow(review_task_lookup)) {
    if (!"active_review_order" %in% names(dt)) {
      dt[, active_review_order := NA_integer_]
    }
    return(dt)
  }
  if ("active_review_order" %in% names(dt)) dt[, active_review_order := NULL]
  merge(dt, review_task_lookup, by = review_task_key_cols,
        all.x = TRUE, sort = FALSE)
}

round_numeric_columns <- function(dt, digits = 4) {
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  for (col in num_cols) set(dt, j = col, value = round(dt[[col]], digits))
  dt
}

make_queue_observatory_links <- function(dt, label_col = "gene_symbol") {
  if (!nrow(dt)) return(dt)
  if (!"observatory_url" %in% names(dt)) dt[, observatory_url := ""]
  if (!"observatory_link_note" %in% names(dt)) dt[, observatory_link_note := ""]
  labels <- clean_text(dt[[label_col]])
  labels[!is_informative_text(labels)] <- clean_text(dt$gene_symbol)
  links <- as.character(htmltools::htmlEscape(labels))
  has_url <- is_informative_text(dt$observatory_url)
  if (any(has_url)) {
    urls <- as.character(htmltools::htmlEscape(dt$observatory_url[has_url],
                                              attribute = TRUE))
    notes <- as.character(htmltools::htmlEscape(dt$observatory_link_note[has_url],
                                               attribute = TRUE))
    link_labels <- as.character(htmltools::htmlEscape(labels[has_url]))
    links[has_url] <- sprintf(
      "<a class=\"observatory-link\" href=\"%s\" target=\"_blank\" rel=\"noopener noreferrer\" title=\"%s\">%s</a>",
      urls, notes, link_labels
    )
  }
  dt[, Links := links]
  dt
}

norm_review_text <- function(x) {
  x <- tolower(clean_text(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

infer_review_label_class <- function(manual_status, browser_grade, branch_dir,
                                     clean_cds_dir, context_match,
                                     translon_status, confounder_status,
                                     negative_control, batch_role) {
  manual <- norm_review_text(manual_status)
  browser <- norm_review_text(browser_grade)
  translon <- norm_review_text(translon_status)
  confounder <- norm_review_text(confounder_status)
  neg <- norm_review_text(negative_control)
  role <- norm_review_text(batch_role)
  direction <- norm_review_text(paste(branch_dir, clean_cds_dir, context_match))
  combined <- norm_review_text(paste(
    manual, browser, translon, confounder, neg, direction, role
  ))
  combined_without_confounder <- norm_review_text(paste(
    manual, browser, translon, neg, direction, role
  ))
  confounder_positive <- grepl(
    "^(possible_confounder|confounded|wrong_context)$|metadata|context_wrong",
    confounder
  )

  if (nzchar(neg) || grepl("negative_control", role)) {
    if (grepl("fail|failed|unexpected|positive|present|strong_signal|signal_present",
              neg)) {
      return("negative_control_fail")
    }
    if (grepl("pass|passed|absent|no_signal|none|clean|negative", neg)) {
      return("negative_control_pass")
    }
  }
  if (confounder_positive ||
      grepl("confound|confounded|metadata|context_wrong|wrong_context",
            combined_without_confounder)) {
    return("confounded")
  }
  if (grepl("measurement|translon|isoform|wrong_isoform|wrong_transcript|p_shift|frameshift|frame_problem|overlap_problem",
            combined)) {
    return("measurement_problem")
  }
  if (grepl("not_reviewable|unreviewable|unclear|ambiguous|too_sparse|low_count|insufficient|mixed",
            manual)) {
    return("not_reviewable")
  }
  if (grepl("negative|reject|rejected|refute|refuted|absent|no_signal|none_visible|not_supported|unsupported|no_support|false|(^|_)no($|_)",
            manual) ||
      grepl("negative|absent|no_signal|none|not_supported|unsupported|no_support",
            browser)) {
    return("negative")
  }
  if (grepl("weak|partial|suggestive|maybe", manual) ||
      grepl("weak|partial|suggestive", browser)) {
    return("weak_positive")
  }
  if (grepl("positive|support|supported|valid|validated|confirm|confirmed|yes|true|interesting",
            manual) ||
      grepl("strong|clear|support|supported|positive|interesting", browser)) {
    return("positive")
  }
  if (grepl("opposite", direction)) return("negative")
  ""
}

write_review_label_draft <- function(draft_row) {
  dir.create(dirname(review_label_drafts_file), recursive = TRUE,
             showWarnings = FALSE)
  drafts <- read_review_label_drafts()
  draft_row <- copy(draft_row)
  for (col in setdiff(review_label_draft_cols, names(draft_row))) {
    draft_row[, (col) := ""]
  }
  draft_row <- draft_row[, ..review_label_draft_cols]
  draft_row[, global_review_order := suppressWarnings(as.integer(global_review_order))]
  drafts <- drafts[global_review_order != draft_row$global_review_order[1]]
  drafts <- rbindlist(list(drafts, draft_row), fill = TRUE)
  drafts[, global_review_order := suppressWarnings(as.integer(global_review_order))]
  setorder(drafts, global_review_order)
  fwrite(drafts, review_label_drafts_file)
  drafts
}

delete_review_label_draft <- function(global_review_order) {
  drafts <- read_review_label_drafts()
  if (!nrow(drafts)) return(drafts)
  order_id <- suppressWarnings(as.integer(global_review_order))
  drafts <- drafts[global_review_order != order_id]
  dir.create(dirname(review_label_drafts_file), recursive = TRUE,
             showWarnings = FALSE)
  fwrite(drafts, review_label_drafts_file)
  drafts
}

label_balance_status <- function(drafts = data.table()) {
  balance <- copy(validation_label_balance)
  if (!nrow(balance)) {
    return(data.table(
      label_class = character(),
      target_minimum = numeric(),
      current_count = numeric(),
      pending_draft_count = integer(),
      label_deficit = numeric(),
      remaining_after_drafts = numeric(),
      deficit_fraction = numeric(),
      acquisition_priority = numeric()
    ))
  }
  if (nrow(drafts) &&
      all(c("review_label_class_preview") %chin% names(drafts))) {
    draft_counts <- drafts[
      is_informative_text(review_label_class_preview),
      .(pending_draft_count = .N),
      by = .(label_class = review_label_class_preview)
    ]
    balance <- merge(balance, draft_counts, by = "label_class",
                     all.x = TRUE, sort = FALSE)
  } else {
    balance[, pending_draft_count := 0L]
  }
  if (!"pending_draft_count" %in% names(balance)) {
    balance[, pending_draft_count := 0L]
  }
  balance[is.na(pending_draft_count), pending_draft_count := 0L]
  if (!"label_deficit" %in% names(balance)) balance[, label_deficit := 0]
  balance[, remaining_after_drafts := pmax(
    0,
    suppressWarnings(as.numeric(label_deficit)) -
      suppressWarnings(as.numeric(pending_draft_count))
  )]
  balance
}

review_label_open_rows <- function(dt) {
  if (!nrow(dt)) return(dt)
  out <- copy(dt)
  if (!"active_review_order" %in% names(out)) {
    out[, active_review_order := NA_integer_]
  }
  out[, active_review_order := suppressWarnings(as.integer(active_review_order))]
  out <- out[is.finite(active_review_order)]
  if (!"queue_source" %in% names(out)) out[, queue_source := "Label gap"]
  if (!"queue_rank" %in% names(out)) {
    out[, queue_rank := suppressWarnings(as.numeric(
      dt_col(out, "label_gap_rank", seq_len(nrow(out)))
    ))]
  }
  if (!"queue_score" %in% names(out)) {
    out[, queue_score := suppressWarnings(as.numeric(
      dt_col(out, "label_acquisition_score", NA_real_)
    ))]
  }
  for (column in c("draft_label_class", "feedback_label_class")) {
    if (!column %in% names(out)) out[, (column) := ""]
  }
  if ("is_open" %in% names(out)) {
    out <- out[is.na(is_open) | is_open == TRUE]
  }
  out <- out[
    !is_informative_text(draft_label_class) &
      !is_informative_text(feedback_label_class)
  ]
  if (!nrow(out)) return(out)
  setorder(out, queue_source, queue_rank, -queue_score,
           active_review_order)
  unique(out, by = "active_review_order")
}

choose_next_review_label_row <- function(dt, current_order = NA_integer_) {
  open <- review_label_open_rows(dt)
  if (!nrow(open)) return(open)
  current_order <- suppressWarnings(as.integer(current_order))
  if (length(current_order) == 1L && is.finite(current_order) &&
      current_order %in% open$active_review_order) {
    idx <- which(open$active_review_order == current_order)[1]
    if (idx < nrow(open)) return(open[idx + 1L])
    if (nrow(open) > 1L) return(open[1L])
    return(open[0])
  }
  open[1L]
}

review_workbench_selected_destination <- function(row) {
  if (is.null(row) || !nrow(row)) return("none")
  review_workbench_row_destination(row)[1]
}

workbench_finite_integer_flag <- function(dt, column) {
  if (!nrow(dt) || !column %in% names(dt)) return(rep(FALSE, nrow(dt)))
  value <- suppressWarnings(as.integer(dt[[column]]))
  is.finite(value) & !is.na(value)
}

review_workbench_row_destination <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(character())
  out <- rep("selected_card", nrow(dt))
  out[workbench_finite_integer_flag(dt, "inspection_source_row")] <-
    "inspection_assets"
  out[workbench_finite_integer_flag(dt, "active_review_order")] <-
    "review_labels"
  out
}

review_workbench_destination_counts <- function(dt) {
  levels <- c("review_labels", "inspection_assets", "selected_card")
  if (is.null(dt) || !nrow(dt)) {
    return(data.table(destination = levels, rows = integer(length(levels))))
  }
  destinations <- review_workbench_row_destination(dt)
  tab <- table(factor(destinations, levels = levels), useNA = "no")
  data.table(destination = levels, rows = as.integer(tab))
}

review_workbench_destination_label <- function(destination) {
  destination <- as.character(destination)[1]
  switch(
    destination,
    review_labels = "Review Labels",
    inspection_assets = "Inspection Assets",
    selected_card = "Selected Card",
    none = "No row selected",
    destination
  )
}

review_workbench_protocol_audit_mask <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(logical())
  out <- rep(FALSE, nrow(dt))
  if ("protocol_pooling_need" %in% names(dt)) {
    need <- suppressWarnings(as.numeric(dt$protocol_pooling_need))
    out <- out | (is.finite(need) & need >= 0.55)
  }
  for (column in intersect(
    c("batch_role", "recommended_action", "review_reason", "flags"),
    names(dt)
  )) {
    value <- clean_text(dt[[column]])
    out <- out | grepl(
      "protocol_pooling|protocol/pooling|protocol audit|pooling audit",
      value,
      ignore.case = TRUE
    )
  }
  out
}

inspection_rows_for_card <- function(row) {
  if (!nrow(qualified_inspection_manifest) || is.null(row)) return(data.table())
  out <- qualified_inspection_manifest[
    gene_symbol == row$gene_symbol[1] &
      tx_id == row$tx_id[1] &
      design_family == row$design_family[1]
  ]
  if (!nrow(out)) return(out)
  setorder(out, source_row)
  out
}

make_review_workbench <- function() {
  empty <- data.table()
  card_links <- unique(cards[, .(
    card_key,
    card_observatory_url = observatory_url,
    card_observatory_link_note = observatory_link_note
  )], by = "card_key")
  manifest_map <- if (nrow(qualified_inspection_manifest)) {
    tmp <- copy(qualified_inspection_manifest)
    tmp <- add_card_key(tmp)
    tmp <- tmp[is_informative_text(card_key)]
    setorder(tmp, source_row)
    unique(tmp[, .(card_key, manifest_source_row = source_row)],
           by = "card_key")
  } else {
    data.table(card_key = character(), manifest_source_row = integer())
  }
  gap_target_lookup <- if (nrow(validation_label_gap_queue)) {
    tmp <- copy(validation_label_gap_queue)
    for (column in review_task_key_cols) {
      if (!column %in% names(tmp)) tmp[, (column) := ""]
      tmp[, (column) := clean_text(get(column))]
    }
    keep <- intersect(c(
      review_task_key_cols, "label_gap_rank", "label_gap_class",
      "target_label_class", "label_acquisition_score", "already_labeled"
    ), names(tmp))
    tmp <- tmp[, ..keep]
    setorder(tmp, label_gap_rank, -label_acquisition_score)
    unique(tmp, by = review_task_key_cols)
  } else {
    data.table(
      gene_symbol = character(), tx_id = character(),
      design_family = character(), batch_id = character(),
      batch_role = character(), label_gap_rank = numeric(),
      label_gap_class = character(), target_label_class = character(),
      label_acquisition_score = numeric(), already_labeled = logical()
    )
  }

  review_batch <- if (nrow(review_batch_template)) {
    out <- data.table(
      queue_source = "Review batch",
      queue_rank = suppressWarnings(as.numeric(
        dt_col(review_batch_template, "global_review_order", seq_len(nrow(review_batch_template)))
      )),
      queue_score = suppressWarnings(as.numeric(
        dt_col(review_batch_template, "acquisition_priority_score", NA_real_)
      )),
      recommended_action = clean_text(
        dt_col(review_batch_template, "action_label", "browser_validate_now")
      ),
      status_label = clean_text(
        dt_col(review_batch_template, "manual_review_status", "unreviewed")
      ),
      gene_symbol = clean_text(review_batch_template$gene_symbol),
      tx_id = clean_text(review_batch_template$tx_id),
      design_family = clean_text(review_batch_template$design_family),
      active_review_order = suppressWarnings(as.integer(
        dt_col(review_batch_template, "global_review_order", NA_integer_)
      )),
      inspection_source_row = NA_integer_,
      batch_id = clean_text(dt_col(review_batch_template, "batch_id", "")),
      batch_name = clean_text(dt_col(review_batch_template, "batch_name", "")),
      batch_role = clean_text(dt_col(review_batch_template, "batch_role", "")),
      claim_stage = clean_text(dt_col(review_batch_template, "claim_stage", "")),
      readiness = clean_text(dt_col(review_batch_template, "review_priority_tier", "")),
      flags = clean_text(dt_col(review_batch_template, "warnings", "")),
      target_label_class = "",
      label_gap_class = "",
      label_acquisition_score = NA_real_,
      protocol_pooling_need = suppressWarnings(as.numeric(
        dt_col(review_batch_template, "protocol_pooling_need", NA_real_)
      )),
      pooling_proof_tier = clean_text(
        dt_col(review_batch_template, "pooling_proof_tier", "")
      ),
      pooling_hazard_class = clean_text(
        dt_col(review_batch_template, "pooling_hazard_class", "")
      ),
      pooling_hazard_score = suppressWarnings(as.numeric(
        dt_col(review_batch_template, "pooling_hazard_score", NA_real_)
      )),
      pooling_hazard_sources = clean_text(
        dt_col(review_batch_template, "pooling_hazard_sources", "")
      ),
      review_question = clean_text(dt_col(review_batch_template, "review_question", "")),
      review_reason = clean_text(dt_col(review_batch_template, "review_reason", "")),
      observatory_url = clean_text(dt_col(review_batch_template, "observatory_url", "")),
      observatory_link_note = clean_text(
        dt_col(review_batch_template, "observatory_link_note", "")
      ),
      is_open = !is_informative_text(
        dt_col(review_batch_template, "manual_review_status", "")
      )
    )
    if (nrow(gap_target_lookup)) {
      out <- merge(
        out,
        gap_target_lookup,
        by = review_task_key_cols,
        all.x = TRUE,
        sort = FALSE,
        suffixes = c("", "_gap")
      )
      for (column in c("target_label_class", "label_gap_class",
                       "label_acquisition_score")) {
        gap_column <- paste0(column, "_gap")
        if (gap_column %in% names(out)) {
          if (is.numeric(out[[column]])) {
            out[!is.finite(get(column)) & is.finite(get(gap_column)),
                (column) := get(gap_column)]
          } else {
            out[!is_informative_text(get(column)) &
                  is_informative_text(get(gap_column)),
                (column) := get(gap_column)]
          }
          out[, (gap_column) := NULL]
        }
      }
    }
    out
  } else empty

  label_gap <- if (nrow(validation_label_gap_queue)) {
    data.table(
      queue_source = "Label gap",
      queue_rank = suppressWarnings(as.numeric(
        dt_col(validation_label_gap_queue, "label_gap_rank", seq_len(nrow(validation_label_gap_queue)))
      )),
      queue_score = suppressWarnings(as.numeric(
        dt_col(validation_label_gap_queue, "label_acquisition_score", NA_real_)
      )),
      recommended_action = paste(
        clean_text(dt_col(validation_label_gap_queue, "label_gap_class", "")),
        clean_text(dt_col(validation_label_gap_queue, "target_label_class", "")),
        sep = " / "
      ),
      status_label = clean_text(
        dt_col(validation_label_gap_queue, "validation_label_class", "unlabeled")
      ),
      gene_symbol = clean_text(validation_label_gap_queue$gene_symbol),
      tx_id = clean_text(validation_label_gap_queue$tx_id),
      design_family = clean_text(validation_label_gap_queue$design_family),
      active_review_order = NA_integer_,
      inspection_source_row = NA_integer_,
      batch_id = clean_text(dt_col(validation_label_gap_queue, "batch_id", "")),
      batch_name = clean_text(dt_col(validation_label_gap_queue, "batch_name", "")),
      batch_role = clean_text(dt_col(validation_label_gap_queue, "batch_role", "")),
      claim_stage = clean_text(dt_col(validation_label_gap_queue, "label_gap_class", "")),
      readiness = clean_text(dt_col(validation_label_gap_queue, "target_label_class", "")),
      target_label_class = clean_text(
        dt_col(validation_label_gap_queue, "target_label_class", "")
      ),
      label_gap_class = clean_text(
        dt_col(validation_label_gap_queue, "label_gap_class", "")
      ),
      label_acquisition_score = suppressWarnings(as.numeric(
        dt_col(validation_label_gap_queue, "label_acquisition_score", NA_real_)
      )),
      protocol_pooling_need = NA_real_,
      pooling_proof_tier = "",
      pooling_hazard_class = "",
      pooling_hazard_score = NA_real_,
      pooling_hazard_sources = "",
      flags = clean_text(dt_col(validation_label_gap_queue, "motif_prior_review_class", "")),
      review_question = clean_text(dt_col(validation_label_gap_queue, "review_question", "")),
      review_reason = clean_text(dt_col(validation_label_gap_queue, "review_reason", "")),
      observatory_url = clean_text(dt_col(validation_label_gap_queue, "observatory_url", "")),
      observatory_link_note = clean_text(
        dt_col(validation_label_gap_queue, "observatory_link_note", "")
      ),
      is_open = !as.logical(dt_col(validation_label_gap_queue,
                                   "already_labeled", FALSE))
    )
  } else empty

  evidence_gap <- if (nrow(evidence_gap_priorities)) {
    data.table(
      queue_source = "Evidence gap",
      queue_rank = suppressWarnings(as.numeric(
        dt_col(evidence_gap_priorities, "acquisition_rank", seq_len(nrow(evidence_gap_priorities)))
      )),
      queue_score = suppressWarnings(as.numeric(
        dt_col(evidence_gap_priorities, "acquisition_priority_score", NA_real_)
      )),
      recommended_action = clean_text(
        dt_col(evidence_gap_priorities, "action_label", "")
      ),
      status_label = clean_text(
        dt_col(evidence_gap_priorities, "validation_status", "unreviewed")
      ),
      gene_symbol = clean_text(evidence_gap_priorities$gene_symbol),
      tx_id = clean_text(evidence_gap_priorities$tx_id),
      design_family = clean_text(evidence_gap_priorities$design_family),
      active_review_order = NA_integer_,
      inspection_source_row = NA_integer_,
      batch_id = "",
      batch_name = "",
      batch_role = "",
      claim_stage = clean_text(dt_col(evidence_gap_priorities, "claim_stage", "")),
      readiness = clean_text(dt_col(evidence_gap_priorities, "review_priority_tier", "")),
      target_label_class = "",
      label_gap_class = "",
      label_acquisition_score = NA_real_,
      protocol_pooling_need = suppressWarnings(as.numeric(
        dt_col(evidence_gap_priorities, "protocol_pooling_need", NA_real_)
      )),
      pooling_proof_tier = clean_text(
        dt_col(evidence_gap_priorities, "pooling_proof_tier", "")
      ),
      pooling_hazard_class = clean_text(
        dt_col(evidence_gap_priorities, "pooling_hazard_class", "")
      ),
      pooling_hazard_score = suppressWarnings(as.numeric(
        dt_col(evidence_gap_priorities, "pooling_hazard_score", NA_real_)
      )),
      pooling_hazard_sources = clean_text(
        dt_col(evidence_gap_priorities, "pooling_hazard_sources", "")
      ),
      flags = clean_text(dt_col(evidence_gap_priorities, "warnings", "")),
      review_question = "",
      review_reason = clean_text(dt_col(evidence_gap_priorities, "review_reason", "")),
      observatory_url = "",
      observatory_link_note = "",
      is_open = !is_informative_text(
        dt_col(evidence_gap_priorities, "validation_status", "")
      )
    )
  } else empty

  inspection_priority <- if (nrow(qualified_review_priority)) {
    data.table(
      queue_source = "Qualified inspection",
      queue_rank = suppressWarnings(as.numeric(
        dt_col(qualified_review_priority, "review_priority_rank", seq_len(nrow(qualified_review_priority)))
      )),
      queue_score = suppressWarnings(as.numeric(
        dt_col(qualified_review_priority, "display_min_exact_counts", NA_real_)
      )),
      recommended_action = clean_text(
        dt_col(qualified_review_priority, "recommended_action", "")
      ),
      status_label = clean_text(
        dt_col(qualified_review_priority, "inspection_readiness", "")
      ),
      gene_symbol = clean_text(qualified_review_priority$gene_symbol),
      tx_id = clean_text(qualified_review_priority$tx_id),
      design_family = clean_text(qualified_review_priority$design_family),
      active_review_order = NA_integer_,
      inspection_source_row = suppressWarnings(as.integer(
        dt_col(qualified_review_priority, "source_row", NA_integer_)
      )),
      batch_id = clean_text(dt_col(qualified_review_priority, "batch_id", "")),
      batch_name = "",
      batch_role = clean_text(dt_col(qualified_review_priority, "batch_role", "")),
      claim_stage = "",
      readiness = clean_text(
        dt_col(qualified_review_priority, "inspection_readiness", "")
      ),
      target_label_class = "",
      label_gap_class = "",
      label_acquisition_score = NA_real_,
      protocol_pooling_need = NA_real_,
      pooling_proof_tier = "",
      pooling_hazard_class = "",
      pooling_hazard_score = NA_real_,
      pooling_hazard_sources = "",
      flags = clean_text(dt_col(qualified_review_priority, "flags", "")),
      review_question = "",
      review_reason = paste(
        "top branch", clean_text(dt_col(qualified_review_priority, "top_branch_class", "")),
        "delta", clean_text(dt_col(qualified_review_priority, "top_branch_delta", "")),
        "feature", clean_text(dt_col(qualified_review_priority, "top_feature_id", ""))
      ),
      observatory_url = "",
      observatory_link_note = "",
      is_open = TRUE
    )
  } else empty

  out <- rbindlist(
    list(review_batch, label_gap, evidence_gap, inspection_priority),
    fill = TRUE
  )
  if (!nrow(out)) return(out)
  out <- add_card_key(out)
  out <- attach_review_task_lookup(out)
  out <- merge(out, manifest_map, by = "card_key", all.x = TRUE, sort = FALSE)
  out[!is.finite(inspection_source_row) & is.finite(manifest_source_row),
      inspection_source_row := manifest_source_row]
  out[, manifest_source_row := NULL]
  out <- merge(out, card_links, by = "card_key", all.x = TRUE, sort = FALSE)
  out[!is_informative_text(observatory_url),
      observatory_url := card_observatory_url]
  out[!is_informative_text(observatory_link_note),
      observatory_link_note := card_observatory_link_note]
  out[, c("card_observatory_url", "card_observatory_link_note") := NULL]
  out[!is_informative_text(status_label), status_label := "open"]
  out[!is_informative_text(recommended_action),
      recommended_action := "review_candidate"]
  setorder(out, queue_source, queue_rank)
  out
}

make_global_report_manifest <- function() {
  reports <- data.table::rbindlist(list(
    data.table(
      category = "dominant_rdgs",
      label = c(
        "Qualified RDG inspection index",
        "RDG review cards",
        "Atlas viral review priority",
        "Atlas confidence vs effect",
        "Atlas buffering heatmap",
        "Branch allocation viral review",
        "Branch usage viral review",
        "Grouped branch usage viral review",
        "Consensus evidence heatmap",
        "Consensus top candidates",
        "Joint branch allocation viral review",
        "Hierarchical joint allocation viral review",
        "Dirichlet-multinomial viral review",
        "Model agreement viral matrix",
        "LOO consensus stability",
        "Transport consensus axes",
        "Context interaction viral contexts",
        "Context partial pooling viral effects",
        "Validation dossier priority",
        "Validation dossier evidence heatmap",
        "Therapeutic hypothesis evidence matrix",
        "Therapeutic perturbation summary",
        "Evidence gap priority map",
        "Allocation output disagreement scatter",
        "Allocation output disagreement branch heatmap",
        "Hierarchical residual viral contexts",
        "Hierarchical context model viral genes",
        "Motif prior viral alignment",
        "Graph geometry viral genes",
        "Review batch readiness"
      ),
      rel_path = c(
        "dominant_rdg_qualified_inspection/qualified_review_index.html",
        "dominant_rdg_review_links/dominant_rdg_review_cards.html",
        "dominant_rdg_atlas/figures/dominant_rdg_atlas_viral_review_priority.html",
        "dominant_rdg_atlas/figures/dominant_rdg_atlas_confidence_vs_effect.html",
        "dominant_rdg_atlas/figures/dominant_rdg_atlas_buffering_heatmap.html",
        "dominant_rdg_branch_allocation/figures/dominant_rdg_branch_allocation_viral_review.html",
        "dominant_rdg_branch_usage/figures/dominant_rdg_branch_usage_viral_review.html",
        "dominant_rdg_grouped_branch_usage/figures/dominant_rdg_grouped_branch_usage_viral_review.html",
        "dominant_rdg_consensus/figures/dominant_rdg_consensus_evidence_heatmap.html",
        "dominant_rdg_consensus/figures/dominant_rdg_consensus_top_candidates.html",
        "dominant_rdg_joint_branch_allocation/figures/joint_branch_allocation_viral_review.html",
        "dominant_rdg_hierarchical_joint_allocation/figures/hierarchical_joint_allocation_viral_review.html",
        "dominant_rdg_dirichlet_multinomial/figures/dirichlet_multinomial_viral_review.html",
        "dominant_rdg_model_agreement/figures/model_agreement_viral_matrix.html",
        "dominant_rdg_model_agreement_loo/figures/model_agreement_loo_consensus_stability.html",
        "dominant_rdg_model_agreement_transport/figures/model_agreement_transport_viral_axes.html",
        "dominant_rdg_context_interactions/figures/context_interaction_viral_top_contexts.html",
        "dominant_rdg_context_partial_pooling/figures/context_partial_pooling_viral_top_effects.html",
        "dominant_rdg_validation_dossiers/figures/rdg_validation_dossier_priority.html",
        "dominant_rdg_validation_dossiers/figures/rdg_validation_dossier_evidence_heatmap.html",
        "dominant_rdg_therapeutic_hypotheses/figures/therapeutic_hypothesis_evidence_matrix.html",
        "dominant_rdg_therapeutic_perturbation_predictions/figures/therapeutic_perturbation_hypothesis_summary.html",
        "dominant_rdg_evidence_gaps/figures/evidence_gap_priority_map.html",
        "dominant_rdg_allocation_output_disagreements/figures/allocation_output_disagreement_scatter.html",
        "dominant_rdg_allocation_output_disagreements/figures/allocation_output_disagreement_branch_heatmap.html",
        "dominant_rdg_hierarchical_residual_audit/figures/hierarchical_residual_viral_contexts.html",
        "dominant_rdg_hierarchical_context_model/figures/hierarchical_context_model_viral_genes.html",
        "dominant_rdg_motif_prior_calibration/figures/motif_prior_viral_alignment.html",
        "dominant_rdg_graph_geometry_model/figures/graph_geometry_model_viral_genes.html",
        "dominant_rdg_review_batches/figures/review_batch_priority_readiness.html"
      )
    ),
    data.table(
      category = "uorf_statistics",
      label = c(
        "uORF structure rule heatmap",
        "Clean-CDS usage by uORF count",
        "Clean-CDS usage by uORF count bin",
        "Clean-CDS usage by overlap class",
        "Clean-CDS usage vs last leader gap",
        "Clean-CDS usage vs minimum inter-uORF gap",
        "Hidden uORF candidate residuals",
        "Overlapping uORF/CDS matched control summary",
        "Overlapping uORF/CDS shift quadrants",
        "Overlapping uORF/CDS on-off review",
        "Sparse on-off candidate scores",
        "Sparse on-off design heatmap",
        "Small-uORF branch coverage gates",
        "Downstream TIS rescue candidates",
        "Downstream TIS gene maps",
        "Next-model candidate relative use matrix",
        "Next-model CDS buffering condition effects",
        "Next-model RDG evidence space",
        "Next-model internal ORF initiation evidence",
        "Next-model internal ORF matched bump evidence"
      ),
      rel_path = c(
        "dominant_uorf_structure_rules/figures/uorf_structure_rule_heatmap.html",
        "dominant_uorf_structure_rules/figures/clean_cds_usage_vs_uorf_count.html",
        "dominant_uorf_structure_rules/figures/clean_cds_usage_by_uorf_count_bin.html",
        "dominant_uorf_structure_rules/figures/clean_cds_usage_by_overlap_class.html",
        "dominant_uorf_structure_rules/figures/clean_cds_usage_vs_last_leader_gap.html",
        "dominant_uorf_structure_rules/figures/clean_cds_usage_vs_min_inter_uorf_gap.html",
        "dominant_uorf_structure_rules/figures/hidden_uorf_candidate_residuals.html",
        "dominant_ouorf_cds_coupling/figures/ouorf_cds_coupling_matched_control_summary.html",
        "dominant_ouorf_cds_coupling/figures/ouorf_cds_coupling_matched_shift_quadrants.html",
        "dominant_ouorf_cds_coupling/figures/ouorf_cds_coupling_onoff_review.html",
        "dominant_ouorf_cds_coupling/figures/ouorf_cds_sparse_onoff_candidate_scores.html",
        "dominant_ouorf_cds_coupling/figures/ouorf_cds_sparse_onoff_design_heatmap.html",
        "dominant_rdg_branch_coverage_qc/figures/branch_coverage_small_uorf_top_features.html",
        "dominant_downstream_tis_rescue/figures/downstream_tis_rescue_best_candidates.html",
        "dominant_downstream_tis_rescue/figures/downstream_tis_rescue_gene_maps.html",
        "dominant_next_model/figures/next_model_candidate_relative_use_matrix.html",
        "dominant_next_model/figures/next_model_cds_buffering_condition_effects.html",
        "dominant_next_model/figures/next_model_rdg_evidence_space.html",
        "dominant_next_model/figures/next_model_internal_orf_initiation_evidence.html",
        "dominant_next_model/figures/next_model_internal_orf_matched_bump_evidence.html"
      )
    ),
    data.table(
      category = "global_qc",
      label = c(
        "Clean-CDS gene enrichment heatmap",
        "Clean-CDS condition enrichment heatmap",
        "Clean-CDS cell-line enrichment heatmap",
        "Clean-CDS tissue enrichment heatmap",
        "Clean-CDS interaction enrichment heatmap",
        "Dominant-state condition enrichment heatmap",
        "Dominant-state cell-line enrichment heatmap",
        "Dominant-state tissue enrichment heatmap",
        "Control-aware condition summary",
        "Broad metadata n-wise heatmap",
        "Broad metadata random-forest importance",
        "GLMNET metadata coefficients",
        "Latent program PCA scatter",
        "Latent program metadata heatmap",
        "Latent program state correlations",
        "Auxiliary state novelty and availability",
        "Auxiliary state existing correlations",
        "Auxiliary state top metadata effects",
        "Score audit winner counts",
        "Score audit winner transitions",
        "Score semantics RDG raw vs percentile",
        "Score semantics design diversity",
        "Target gene frame vs global frame",
        "CDS frame profiles",
        "Viral target vs housekeeping frame",
        "Infected off-frame hotspot aggregate profiles",
        "Infected off-frame review priority",
        "Postviral matched viral/IFN review",
        "Postviral matched gene-design evidence",
        "Postviral matched stability",
        "Postviral top RDG feature shifts",
        "Postviral state correlation heatmap"
      ),
      rel_path = c(
        "dominant_state_clean_cds_enrichment_gene_heatmap.html",
        "dominant_state_clean_cds_enrichment_condition_heatmap.html",
        "dominant_state_clean_cds_enrichment_cell_line_heatmap.html",
        "dominant_state_clean_cds_enrichment_tissue_heatmap.html",
        "dominant_state_clean_cds_enrichment_interactions_heatmap.html",
        "dominant_state_enrichment_condition_heatmap.html",
        "dominant_state_enrichment_cell_line_heatmap.html",
        "dominant_state_enrichment_tissue_heatmap.html",
        "dominant_state_control_aware_condition_summary_heatmap.html",
        "dominant_state_broad_metadata_nwise_heatmap.html",
        "dominant_state_broad_metadata_rf_importance_heatmap.html",
        "dominant_state_glmnet_metadata_coefficient_heatmap.html",
        "dominant_state_latent_program_pca_scatter.html",
        "dominant_state_latent_program_pc_metadata_heatmap.html",
        "dominant_state_latent_program_state_correlation_heatmap.html",
        "dominant_state_auxiliary_modules/figures/auxiliary_state_novelty_and_availability.html",
        "dominant_state_auxiliary_modules/figures/auxiliary_state_existing_state_correlations.html",
        "dominant_state_auxiliary_modules/figures/auxiliary_state_top_metadata_effects.html",
        "dominant_state_score_audit/figures/dominant_state_score_audit_winner_counts.html",
        "dominant_state_score_audit/figures/dominant_state_score_audit_winner_transitions.html",
        "dominant_state_score_semantics/figures/dominant_state_score_semantics_rdg_raw_vs_percentile.html",
        "dominant_state_score_semantics_design_aware/figures/dominant_state_score_semantics_design_diversity.html",
        "dominant_cds_frame_bumpiness_qc/figures/target_gene_frame_vs_global_frame.html",
        "dominant_cds_frame_bumpiness_qc/figures/target_gene_cds_frame_profiles.html",
        "dominant_cds_frame_bumpiness_qc/figures/viral_target_vs_housekeeping_frame.html",
        "dominant_infected_off_frame_hotspots/figures/infected_off_frame_hotspot_aggregate_profiles.html",
        "dominant_infected_off_frame_review/figures/infected_off_frame_hotspot_review_priority.html",
        "postviral_fatigue_outputs/figures/postviral_fatigue_matched_viral_ifn_review.html",
        "postviral_fatigue_outputs/figures/postviral_fatigue_matched_gene_design_evidence.html",
        "postviral_fatigue_outputs/figures/postviral_fatigue_matched_gene_design_loo_stability.html",
        "postviral_fatigue_outputs/figures/postviral_fatigue_top_rdg_feature_shifts.html",
        "postviral_fatigue_outputs/figures/postviral_fatigue_state_correlation_heatmap.html"
      )
    )
  ), use.names = TRUE)
  reports[, group := "Other"]
  reports[
    category == "dominant_rdgs" &
      grepl("Qualified|review cards|viral review priority|Atlas confidence|Evidence gap|Review batch",
            label, ignore.case = TRUE),
    group := "Review queues"
  ]
  reports[
    category == "dominant_rdgs" &
      grepl("Consensus|Model agreement|LOO|Transport", label),
    group := "Agreement and stability"
  ]
  reports[
    category == "dominant_rdgs" &
      grepl("Branch|allocation|buffering|Dirichlet", label,
            ignore.case = TRUE),
    group := "Branch allocation"
  ]
  reports[
    category == "dominant_rdgs" &
      grepl("Context|pooling|Validation", label, ignore.case = TRUE),
    group := "Context and validation"
  ]
  reports[
    category == "dominant_rdgs" &
      grepl("Therapeutic", label, ignore.case = TRUE),
    group := "Therapeutic hypotheses"
  ]
  reports[
    category == "dominant_rdgs" &
      grepl("residual|model viral|Motif|Graph|disagreement", label,
            ignore.case = TRUE),
    group := "Diagnostics and model internals"
  ]
  reports[
    category == "uorf_statistics" &
      grepl("structure|usage|overlap class|leader gap|inter-uORF",
            label, ignore.case = TRUE),
    group := "uORF structure rules"
  ]
  reports[
    category == "uorf_statistics" &
      grepl("coupling|matched control|shift quadrants|on-off|Sparse|Hidden",
            label, ignore.case = TRUE),
    group := "uORF/CDS coupling"
  ]
  reports[
    category == "uorf_statistics" &
      grepl("TIS|Next-model|coverage|Downstream", label, ignore.case = TRUE),
    group := "ORF starts and coverage"
  ]
  reports[
    category == "global_qc" &
      grepl("enrichment|Control-aware", label, ignore.case = TRUE),
    group := "State enrichment"
  ]
  reports[
    category == "global_qc" &
      grepl("metadata|Latent|Auxiliary", label, ignore.case = TRUE),
    group := "Metadata and latent state"
  ]
  reports[
    category == "global_qc" &
      grepl("Score|frame|off-frame", label, ignore.case = TRUE),
    group := "Score and frame QC"
  ]
  reports[
    category == "global_qc" &
      grepl("Postviral", label, ignore.case = TRUE),
    group := "Postviral fatigue"
  ]
  reports[, abs_path := normalizePath(vapply(
    rel_path,
    analysis_path,
    character(1)
  ), mustWork = FALSE)]
  reports[, exists := file.exists(abs_path)]
  reports[exists == TRUE]
}

global_reports <- make_global_report_manifest()

report_category_choices <- c(
  "Dominant RDG reports" = "dominant_rdgs",
  "uORF and ORF-start reports" = "uorf_statistics",
  "Global state and QC reports" = "global_qc"
)

report_category_note <- function(report_category) {
  switch(
    report_category,
    dominant_rdgs = paste(
      "RDG inspection, branch allocation, model agreement, validation,",
      "and therapeutic-hypothesis outputs."
    ),
    uorf_statistics = paste(
      "uORF structure, overlapping-uORF/CDS coupling, coverage gates,",
      "and downstream ORF-start outputs."
    ),
    global_qc = paste(
      "Global dominant-state, metadata, frame, infected/off-frame,",
      "and postviral fatigue outputs."
    ),
    "Generated HTML reports."
  )
}

report_choices <- function(report_category) {
  dt <- global_reports[global_reports[["category"]] == report_category]
  if (!nrow(dt)) return(stats::setNames("", "No generated HTML reports found"))
  stats::setNames(dt$rel_path, dt$label)
}

report_choice_groups <- function(report_category) {
  dt <- copy(global_reports[global_reports[["category"]] == report_category])
  if (!nrow(dt)) return(stats::setNames("", "No generated HTML reports found"))
  setorder(dt, group, label)
  split(stats::setNames(dt$rel_path, dt$label), dt$group)
}

first_report_choice <- function(report_category) {
  choices <- report_choice_groups(report_category)
  if (is.list(choices)) {
    return(unname(unlist(choices, use.names = FALSE))[1])
  }
  unname(choices[[1]])
}

report_browser_ui <- function() {
  div(
    class = "report-browser",
    div(
      class = "report-controls panel",
      div(
        class = "report-controls-grid",
        selectInput(
          "report_category",
          "Report family",
          choices = report_category_choices,
          selected = "dominant_rdgs"
        ),
        uiOutput("report_select_ui")
      ),
      div(class = "small-note", uiOutput("report_category_note"))
    ),
    uiOutput("report_frame")
  )
}

render_static_html_frame <- function(rel_path,
                                     source_rel_path = NULL,
                                     empty_message = "No document selected.",
                                     missing_message = "Selected document file is missing.",
                                     open_label = "Open document") {
  file_info <- make_selected_file_link(rel_path)
  validate(
    need(!is.null(file_info), empty_message),
    need(isTRUE(file_info$exists), missing_message)
  )
  source_row <- NULL
  if (!is.null(source_rel_path)) {
    source_info <- make_selected_file_link(source_rel_path)
    if (!is.null(source_info) && isTRUE(source_info$exists) &&
        is_scalar_informative(source_info$href)) {
      source_newer <- isTRUE(file.info(source_info$abs_path)$mtime >
                               file.info(file_info$abs_path)$mtime)
      source_row <- tagList(
        tags$br(),
        tags$a(
          class = "file-link",
          href = source_info$href,
          target = "_blank",
          rel = "noopener noreferrer",
          title = source_info$abs_path,
          "Source Markdown"
        ),
        tags$span(
          class = if (source_newer) "doc-status stale" else "doc-status",
          if (source_newer) {
            "source newer than rendered HTML"
          } else {
            "rendered HTML is current"
          }
        )
      )
    }
  }
  div(
    class = "report-frame-shell panel",
    div(
      class = "report-open-row",
      tags$a(
        class = "primary-link",
        href = file_info$href,
        target = "_blank",
        rel = "noopener noreferrer",
        open_label
      ),
      tags$span(class = "link-note", file_info$abs_path),
      source_row
    ),
    tags$iframe(
      class = "report-frame",
      src = file_info$href,
      loading = "lazy"
    )
  )
}

render_report_frame <- function(rel_path) {
  render_static_html_frame(
    rel_path,
    empty_message = "No report selected.",
    missing_message = "Selected report file is missing.",
    open_label = "Open selected report"
  )
}

load_rdg_metadata <- function() {
  candidates <- c(
    analysis_path("metadata_done_samples_extended_qc.csv"),
    "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(list(data = data.table(), path = NA_character_))

  path <- candidates[1]
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  wanted <- c(
    "Run", "study", "CONDITION", "CELL_LINE", "TISSUE", "GENE",
    "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT", "AUTHOR"
  )
  cols <- intersect(wanted, header)
  if (!all(c("Run", "CONDITION") %chin% cols)) {
    return(list(data = data.table(), path = path))
  }
  dt <- fread(path, select = cols, showProgress = FALSE)
  for (col in names(dt)) set(dt, j = col, value = clean_text(dt[[col]]))
  list(data = dt, path = path)
}

make_gene_label_map <- function(card_dt) {
  out <- unique(card_dt[, .(tx_id = clean_text(tx_id),
                            gene_symbol = clean_text(gene_symbol))])
  out[, observatory_gene_label := gene_symbol]

  if (!requireNamespace("ORFik", quietly = TRUE)) return(out)

  labels <- tryCatch({
    df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
    sy <- ORFik::symbols(df)
    setDT(sy)
    if (!all(c("ensembl_tx_name", "ensembl_gene_id") %chin% names(sy))) {
      stop("Missing expected symbol columns", call. = FALSE)
    }
    sy[, gene_label_symbol := clean_text(external_gene_name)]
    sy[!is_informative_text(gene_label_symbol),
       gene_label_symbol := clean_text(ensembl_gene_id)]
    sy[, observatory_gene_label :=
         paste(gene_label_symbol, clean_text(ensembl_gene_id))]
    sy[, observatory_gene_label :=
         sub(" ", "-", observatory_gene_label, fixed = TRUE)]
    sy[, observatory_gene_label :=
         sub("(^-)|(^NA-)", "", observatory_gene_label, perl = TRUE)]
    unique(sy[
      is_informative_text(ensembl_tx_name) &
        is_informative_text(observatory_gene_label),
      .(tx_id = clean_text(ensembl_tx_name), observatory_gene_label)
    ])
  }, error = function(e) data.table())

  if (!nrow(labels)) return(out)
  out <- merge(out[, !"observatory_gene_label"], labels,
               by = "tx_id", all.x = TRUE, sort = FALSE)
  out[!is_informative_text(observatory_gene_label),
      observatory_gene_label := gene_symbol]
  out
}

rdg_observatory_compact_selections <- function(selections) {
  if (is.null(selections) || !is.list(selections)) {
    return(list(active = "1", order = "1", labels = list("1" = ""),
                runs = character(), p = list("1" = integer()),
                d = list("1" = integer())))
  }

  ids <- as.character(selections$index %||%
                        names(selections$plot_selections %||% list("1" = NULL)))
  if (length(ids) == 0) ids <- "1"
  plot_sel <- selections$plot_selections %||% list()
  data_sel <- selections$data_table_selections %||% list()
  labels <- selections$labels %||% list()
  active <- as.character(selections$active_selection_id %||% ids[1])
  if (!(active %in% ids)) active <- ids[1]

  all_runs <- unique(unlist(c(plot_sel[ids], data_sel[ids]), use.names = FALSE))
  all_runs <- all_runs[!is.na(all_runs) & nzchar(all_runs)]
  run_to_idx <- stats::setNames(seq_along(all_runs), all_runs)

  encode_idx <- function(x) {
    if (is.null(x) || length(x) == 0) return(integer())
    as.integer(run_to_idx[as.character(x)])
  }

  p <- setNames(lapply(ids, function(selection_id) {
    encode_idx(plot_sel[[selection_id]])
  }), ids)
  d <- setNames(lapply(ids, function(selection_id) {
    encode_idx(data_sel[[selection_id]])
  }), ids)
  lb <- setNames(lapply(ids, function(selection_id) {
    as.character(labels[[selection_id]] %||% "")
  }), ids)

  list(active = active, order = ids, labels = lb, runs = all_runs, p = p, d = d)
}

rdg_observatory_state_from_inputs <- function(
  selected_experiment,
  color_by,
  view = c("umap", "browser")[1],
  browser = list(),
  selections = list()
) {
  view <- tolower(as.character(view)[1])
  if (!view %in% c("umap", "browser")) view <- "umap"
  list(
    v = 1L,
    exp = as.character(selected_experiment %||% ""),
    color_by = as.character(color_by %||% character()),
    view = view,
    browser = browser,
    selections = rdg_observatory_compact_selections(selections)
  )
}

rdg_make_observatory_url_state_param <- function(state) {
  json <- jsonlite::toJSON(state, auto_unbox = TRUE, null = "null")
  compressed <- memCompress(charToRaw(json), type = "gzip")
  encoded <- jsonlite::base64_enc(compressed)
  encoded <- gsub("[\r\n]", "", encoded)
  encoded <- chartr("+/", "-_", encoded)
  sub("=+$", "", encoded)
}

best_branch_context <- function(row) {
  if (!nrow(branch_contexts)) return(NULL)
  out <- branch_contexts[
    gene_symbol == row$gene_symbol[1] &
      design_family == row$design_family[1]
  ]
  if ("tx_id" %in% names(out) && "tx_id" %in% names(row)) {
    out <- out[tx_id == row$tx_id[1]]
  }
  if (!nrow(out)) return(NULL)
  out <- copy(out)
  if ("branch_context_rank" %in% names(out) &&
      "joint_review_priority" %in% names(out)) {
    setorder(out, branch_context_rank, -joint_review_priority)
  } else if ("branch_context_rank" %in% names(out)) {
    setorder(out, branch_context_rank)
  } else if ("joint_review_priority" %in% names(out)) {
    setorder(out, -joint_review_priority)
  }
  out[1]
}

filter_metadata_by_context <- function(meta, ctx) {
  if (!nrow(meta) || is.null(ctx)) return(meta)
  for (col in c("AUTHOR", "CELL_LINE", "TISSUE", "GENE", "INHIBITOR",
                "FRACTION", "Cancer_type", "Sex", "TIMEPOINT")) {
    if (col %in% names(meta) && col %in% names(ctx) &&
        is_scalar_informative(ctx[[col]][1])) {
      meta <- meta[get(col) == scalar_text(ctx[[col]][1])]
    }
  }
  meta
}

condition_runs <- function(meta, condition_values) {
  if (!nrow(meta) || !"Run" %in% names(meta) || !"CONDITION" %in% names(meta) ||
      !length(condition_values)) {
    return(character())
  }
  unique(meta[CONDITION %chin% condition_values, Run])
}

get_observatory_subset <- function(row, metadata) {
  ctx <- best_branch_context(row)
  study <- if (!is.null(ctx) && "study" %in% names(ctx) &&
               is_scalar_informative(ctx$study)) {
    scalar_text(ctx$study)
  } else {
    scalar_text(row$best_study)
  }
  case_condition <- if (!is.null(ctx) && "case_condition" %in% names(ctx) &&
                        is_scalar_informative(ctx$case_condition)) {
    scalar_text(ctx$case_condition)
  } else {
    scalar_text(row$best_case_condition)
  }
  control_conditions <- if (!is.null(ctx) &&
                            "control_conditions" %in% names(ctx) &&
                            is_scalar_informative(ctx$control_conditions)) {
    split_condition_values(ctx$control_conditions)
  } else {
    split_condition_values(row$best_control_conditions)
  }

  if (!nrow(metadata) || !"study" %in% names(metadata)) {
    return(list(
      case_runs = character(), control_runs = character(),
      note = "No local metadata available for run subset selection",
      context_label = "gene/transcript only"
    ))
  }

  base <- metadata
  if (is_scalar_informative(study)) {
    study_value <- study
    base <- base[get("study") == study_value]
  }
  exact <- filter_metadata_by_context(copy(base), ctx)
  case_runs <- condition_runs(exact, case_condition)
  control_runs <- condition_runs(exact, control_conditions)
  context_label <- if (!is.null(ctx) && "branch_context_label" %in% names(ctx) &&
                       is_scalar_informative(ctx$branch_context_label)) {
    scalar_text(ctx$branch_context_label)
  } else {
    paste(study, case_condition, sep = " / ")
  }
  note <- "Best branch context"

  if (!length(case_runs) && !length(control_runs) && nrow(base)) {
    case_runs <- condition_runs(base, case_condition)
    control_runs <- condition_runs(base, control_conditions)
    note <- "Study and condition subset; context fields did not match metadata"
  }
  if (!length(case_runs) && !length(control_runs)) {
    note <- "Gene/transcript link only; no matching case/control runs"
  }

  list(
    case_runs = case_runs,
    control_runs = control_runs,
    note = note,
    context_label = context_label,
    case_condition = case_condition,
    control_conditions = paste(control_conditions, collapse = "; "),
    study = study
  )
}

make_observatory_link_row <- function(row, metadata, gene_labels) {
  gene_label <- gene_labels[tx_id == row$tx_id[1], observatory_gene_label][1]
  if (!is_scalar_informative(gene_label)) gene_label <- row$gene_symbol[1]
  subset <- get_observatory_subset(row, metadata)

  selection_ids <- character()
  plot_selections <- list()
  data_table_selections <- list()
  labels <- list()
  if (length(subset$control_runs)) {
    selection_ids <- c(selection_ids, "1")
    plot_selections[["1"]] <- subset$control_runs
    data_table_selections[["1"]] <- subset$control_runs
    labels[["1"]] <- paste0("Control: ", subset$control_conditions)
  }
  if (length(subset$case_runs)) {
    case_id <- as.character(length(selection_ids) + 1L)
    selection_ids <- c(selection_ids, case_id)
    plot_selections[[case_id]] <- subset$case_runs
    data_table_selections[[case_id]] <- subset$case_runs
    labels[[case_id]] <- paste0("Case: ", subset$case_condition)
  }
  if (!length(selection_ids)) {
    selection_ids <- "1"
    plot_selections[["1"]] <- character()
    data_table_selections[["1"]] <- character()
    labels[["1"]] <- "No matched case/control subset"
  }
  if (!("3" %in% selection_ids)) {
    selection_ids <- unique(c(selection_ids, "3"))
    plot_selections[["3"]] <- character()
    data_table_selections[["3"]] <- character()
    labels[["3"]] <- "All merged"
  }

  browser <- list(
    gene = gene_label,
    tx = scalar_text(row$tx_id),
    frames_type = "columns",
    kmer = 1L,
    extendLeaders = 30L,
    extendTrailers = 30L,
    viewMode = FALSE,
    other_tx = FALSE,
    collapsed_introns = FALSE,
    collapsed_introns_width = 30L,
    genomic_region = "",
    zoom_range = "",
    customSequence = "",
    go = TRUE
  )
  state <- rdg_observatory_state_from_inputs(
    selected_experiment = "all_samples-Homo_sapiens",
    color_by = c("tissue", "cell_line"),
    view = "browser",
    browser = browser,
    selections = list(
      index = selection_ids,
      plot_selections = plot_selections,
      data_table_selections = data_table_selections,
      labels = labels,
      active_selection_id = "3"
    )
  )
  url <- paste0("https://ribocrypt.org/#Observatory?obs_state=",
                rdg_make_observatory_url_state_param(state))
  note <- paste0(
    subset$note,
    "; study=", subset$study %||% "",
    "; context=", subset$context_label %||% "",
    "; case n=", length(subset$case_runs),
    "; control n=", length(subset$control_runs)
  )
  data.table(
    observatory_url = url,
    observatory_link_note = note,
    observatory_gene_label = gene_label,
    observatory_case_n = length(subset$case_runs),
    observatory_control_n = length(subset$control_runs),
    observatory_context_label = subset$context_label %||% ""
  )
}

parse_partial_pooling_labels <- function(label_text) {
  label_text <- scalar_text(label_text)
  if (!is_scalar_informative(label_text)) return(data.table())
  rows <- unlist(strsplit(label_text, "\\s*;\\s*", perl = TRUE),
                 use.names = FALSE)
  rows <- rows[is_informative_text(rows)]
  out <- rbindlist(lapply(rows, function(label) {
    parts <- clean_text(strsplit(label, "\\s*\\|\\s*", perl = TRUE)[[1]])
    if (length(parts) < 5L) return(NULL)
    control <- sub("^vs\\s+", "", parts[5], ignore.case = TRUE)
    data.table(
      study = parts[1],
      CELL_LINE = parts[2],
      TISSUE = parts[3],
      case_condition = parts[4],
      control_conditions = control
    )
  }), fill = TRUE)
  if (!nrow(out)) return(out)
  unique(out)
}

filter_metadata_by_partial_label <- function(meta, label) {
  if (!nrow(meta) || !nrow(label)) return(meta)
  for (col in c("study", "CELL_LINE", "TISSUE")) {
    if (col %in% names(meta) && col %in% names(label) &&
        is_scalar_informative(label[[col]][1])) {
      value <- scalar_text(label[[col]][1])
      meta <- meta[get(col) == value]
    }
  }
  meta
}

partial_pooling_label_runs <- function(labels, metadata, condition_field) {
  if (!nrow(labels) || !nrow(metadata) ||
      !condition_field %in% names(labels)) {
    return(character())
  }
  unique(unlist(lapply(seq_len(nrow(labels)), function(i) {
    label <- labels[i]
    condition_values <- split_condition_values(label[[condition_field]][1])
    meta <- filter_metadata_by_partial_label(metadata, label)
    condition_runs(meta, condition_values)
  }), use.names = FALSE))
}

get_partial_pooling_observatory_subset <- function(row, metadata) {
  context_labels <- parse_partial_pooling_labels(row$context_top_labels)
  complement_labels <- parse_partial_pooling_labels(row$complement_top_labels)
  context_axis <- scalar_text(row$context_axis)
  context_value <- scalar_text(row$context_value)
  branch_class <- scalar_text(row$branch_class)
  context_label <- paste0(context_axis, "=", context_value,
                          "; branch=", branch_class)
  list(
    context_case_runs = partial_pooling_label_runs(
      context_labels, metadata, "case_condition"
    ),
    context_control_runs = partial_pooling_label_runs(
      context_labels, metadata, "control_conditions"
    ),
    complement_case_runs = partial_pooling_label_runs(
      complement_labels, metadata, "case_condition"
    ),
    complement_control_runs = partial_pooling_label_runs(
      complement_labels, metadata, "control_conditions"
    ),
    context_label = context_label
  )
}

make_partial_pooling_observatory_link_row <- function(row, metadata,
                                                     gene_labels) {
  gene_label <- gene_labels[tx_id == row$tx_id[1], observatory_gene_label][1]
  if (!is_scalar_informative(gene_label)) gene_label <- row$gene_symbol[1]
  subset <- get_partial_pooling_observatory_subset(row, metadata)

  selection_ids <- character()
  plot_selections <- list()
  data_table_selections <- list()
  labels <- list()

  add_selection <- function(label, runs) {
    runs <- unique(as.character(runs))
    runs <- runs[!is.na(runs) & nzchar(runs)]
    if (!length(runs)) return(invisible(NULL))
    selection_id <- as.character(length(selection_ids) + 1L)
    selection_ids <<- c(selection_ids, selection_id)
    plot_selections[[selection_id]] <<- runs
    data_table_selections[[selection_id]] <<- runs
    labels[[selection_id]] <<- label
    invisible(NULL)
  }

  add_selection(paste0("Context case: ", subset$context_label),
                subset$context_case_runs)
  add_selection(paste0("Context control: ", subset$context_label),
                subset$context_control_runs)
  add_selection("Complement case", subset$complement_case_runs)
  add_selection("Complement control", subset$complement_control_runs)

  if (!length(selection_ids)) {
    selection_ids <- "1"
    plot_selections[["1"]] <- character()
    data_table_selections[["1"]] <- character()
    labels[["1"]] <- "No matched partial-pooling subset"
  }
  all_id <- as.character(length(selection_ids) + 1L)
  selection_ids <- unique(c(selection_ids, all_id))
  plot_selections[[all_id]] <- character()
  data_table_selections[[all_id]] <- character()
  labels[[all_id]] <- "All merged"

  browser <- list(
    gene = gene_label,
    tx = scalar_text(row$tx_id),
    frames_type = "columns",
    kmer = 1L,
    extendLeaders = 30L,
    extendTrailers = 30L,
    viewMode = FALSE,
    other_tx = FALSE,
    collapsed_introns = FALSE,
    collapsed_introns_width = 30L,
    genomic_region = "",
    zoom_range = "",
    customSequence = "",
    go = TRUE
  )
  state <- rdg_observatory_state_from_inputs(
    selected_experiment = "all_samples-Homo_sapiens",
    color_by = c("tissue", "cell_line"),
    view = "browser",
    browser = browser,
    selections = list(
      index = selection_ids,
      plot_selections = plot_selections,
      data_table_selections = data_table_selections,
      labels = labels,
      active_selection_id = all_id
    )
  )
  note <- paste0(
    "Partial pooling: ", subset$context_label,
    "; context case n=", length(subset$context_case_runs),
    "; context control n=", length(subset$context_control_runs),
    "; complement case n=", length(subset$complement_case_runs),
    "; complement control n=", length(subset$complement_control_runs)
  )
  data.table(
    partial_pooling_observatory_url = paste0(
      "https://ribocrypt.org/#Observatory?obs_state=",
      rdg_make_observatory_url_state_param(state)
    ),
    partial_pooling_observatory_link_note = note,
    partial_pooling_context_case_n = length(subset$context_case_runs),
    partial_pooling_context_control_n = length(subset$context_control_runs),
    partial_pooling_complement_case_n = length(subset$complement_case_runs),
    partial_pooling_complement_control_n =
      length(subset$complement_control_runs)
  )
}

metadata_info <- load_rdg_metadata()
rdg_metadata <- metadata_info$data
gene_label_map <- make_gene_label_map(cards)
link_rows <- rbindlist(lapply(seq_len(nrow(cards)), function(i) {
  make_observatory_link_row(cards[i], rdg_metadata, gene_label_map)
}), fill = TRUE)
cards[, names(link_rows) := link_rows]
review_workbench <- make_review_workbench()
review_workbench_source_choices <- c(
  "All queues" = "all",
  stats::setNames(sort(unique(review_workbench$queue_source)),
                  sort(unique(review_workbench$queue_source)))
)
review_workbench_gene_choices <- safe_choices(review_workbench$gene_symbol)
review_workbench_display_columns <- c(
  "Links", "review_destination", "gene_symbol", "tx_id", "design_family",
  "queue_source", "queue_rank", "queue_score",
  "target_label_class", "label_gap_class", "target_label_deficit",
  "target_remaining_after_drafts", "label_acquisition_score",
  "recommended_action", "status_label", "inspection_source_row",
  "active_review_order", "draft_label_class", "draft_saved_at",
  "feedback_label_class", "feedback_review_date",
  "batch_id", "batch_role", "claim_stage", "readiness",
  "protocol_pooling_need", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "pooling_proof_tier",
  "review_question", "review_reason", "flags", "observatory_link_note"
)
review_manual_status_choices <- c(
  "Not set" = "",
  "Positive" = "positive",
  "Weak positive" = "weak_positive",
  "Negative" = "negative",
  "Confounded" = "confounded",
  "Measurement problem" = "measurement_problem",
  "Not reviewable" = "not_reviewable"
)
review_browser_grade_choices <- c(
  "Not set" = "",
  "Strong support" = "strong_positive",
  "Weak support" = "weak_positive",
  "No visible signal" = "negative",
  "Confounded" = "confounded",
  "Measurement problem" = "measurement_problem",
  "Unclear" = "unclear",
  "Not reviewable" = "not_reviewable"
)
review_direction_choices <- c(
  "Not set" = "",
  "Verified" = "verified",
  "Opposite / mismatch" = "opposite",
  "Unclear / mixed" = "unclear",
  "Not applicable" = "not_applicable"
)
review_translon_choices <- c(
  "Not set" = "",
  "Structure verified" = "verified",
  "Unclear" = "unclear",
  "Wrong isoform/transcript" = "wrong_isoform",
  "Measurement problem" = "measurement_problem",
  "Not reviewable" = "not_reviewable"
)
review_confounder_choices <- c(
  "Not set" = "",
  "No confounder seen" = "no_confounder_seen",
  "Possible confounder" = "possible_confounder",
  "Confounded" = "confounded",
  "Wrong context" = "wrong_context"
)
review_negative_control_choices <- c(
  "Not set" = "",
  "Pass / no signal" = "pass",
  "Fail / signal detected" = "fail",
  "Unclear" = "unclear",
  "Not applicable" = "not_applicable"
)
review_next_action_choices <- c(
  "Not set" = "",
  "Merge to browser notes" = "merge_to_browser_notes",
  "Needs Observatory inspection" = "needs_observatory_inspection",
  "Needs translon fix" = "needs_translon_fix",
  "Needs context audit" = "needs_context_audit",
  "Needs metadata audit" = "needs_metadata_audit",
  "Reject / deprioritize" = "reject",
  "Revisit later" = "revisit_later"
)
label_dashboard_targets <- sort(unique(clean_text(c(
  dt_col(validation_label_gap_queue, "target_label_class", character()),
  dt_col(validation_label_gap_queue, "label_gap_class", character()),
  dt_col(validation_label_balance, "label_class", character())
))))
label_dashboard_targets <- label_dashboard_targets[
  is_informative_text(label_dashboard_targets)
]
label_dashboard_target_choices <- c(
  "All label targets" = "all",
  stats::setNames(label_dashboard_targets, label_dashboard_targets)
)
review_workbench_target_choices <- label_dashboard_target_choices

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background: #f4f6f8;
        color: #1f2933;
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
      }
      .container-fluid { max-width: 1500px; }
      .app-header {
        padding: 18px 2px 12px 2px;
        border-bottom: 1px solid #dce2ea;
        margin-bottom: 14px;
      }
      .app-header h1 {
        margin: 0;
        font-size: 26px;
        line-height: 1.2;
        font-weight: 650;
      }
      .app-subtitle {
        color: #52606d;
        max-width: 960px;
        margin-top: 5px;
      }
      .filters-shell, .panel {
        background: #fff;
        border: 1px solid #d9dde3;
        border-radius: 8px;
        box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
      }
      .filters-shell { padding: 12px 12px 14px 12px; }
      .filters-header {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 10px;
        padding-bottom: 10px;
        margin-bottom: 2px;
        border-bottom: 1px solid #edf1f5;
      }
      .filters-header .filter-title { margin: 0; }
      .filter-section {
        border-bottom: 1px solid #edf1f5;
        padding-bottom: 12px;
        margin-bottom: 12px;
      }
      .filter-section:last-child { border-bottom: 0; margin-bottom: 0; }
      .filter-title {
        color: #243b53;
        font-size: 12px;
        font-weight: 700;
        letter-spacing: 0.04em;
        margin: 0 0 8px 0;
        text-transform: uppercase;
      }
      details {
        border-top: 1px solid #edf1f5;
        padding: 9px 0 0 0;
        margin-top: 9px;
      }
      details summary {
        cursor: pointer;
        color: #243b53;
        font-size: 12px;
        font-weight: 700;
        letter-spacing: 0.04em;
        text-transform: uppercase;
        margin-bottom: 8px;
      }
      .filter-summary {
        padding: 9px 0 2px 0;
      }
      .filter-summary-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 7px;
        margin-bottom: 8px;
      }
      .filter-stat {
        background: #f7fafc;
        border: 1px solid #e5eaf0;
        border-radius: 7px;
        min-height: 48px;
        padding: 7px 8px;
      }
      .filter-stat b {
        display: block;
        color: #243b53;
        font-size: 16px;
        line-height: 1.05;
      }
      .filter-stat span {
        color: #52606d;
        font-size: 11px;
      }
      .filter-chip-list {
        display: flex;
        flex-wrap: wrap;
        gap: 5px;
      }
      .filter-chip {
        max-width: 100%;
        border: 1px solid #cbd5e1;
        border-radius: 999px;
        color: #243b53;
        background: #fff;
        font-size: 11px;
        line-height: 1.2;
        padding: 4px 7px;
        overflow-wrap: anywhere;
      }
      .filter-chip.muted {
        color: #6b7280;
        background: #f8fafc;
      }
      .panel {
        padding: 13px 15px;
        margin-bottom: 13px;
      }
      .metrics-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(132px, 1fr));
        gap: 8px;
      }
      .metric {
        padding: 9px 10px;
        background: #f7fafc;
        border: 1px solid #e5eaf0;
        border-radius: 7px;
        min-height: 58px;
      }
      .metric b { display: block; font-size: 18px; line-height: 1.1; }
      .metric span { color: #52606d; font-size: 12px; }
      .detail-grid {
        display: grid;
        grid-template-columns: minmax(170px, 230px) minmax(0, 1fr);
        gap: 6px 13px;
        font-size: 13px;
      }
      .detail-grid .key { color: #596273; font-weight: 600; }
      .detail-grid div { overflow-wrap: anywhere; }
      .warning { color: #8a4b00; font-weight: 600; }
      .small-note { color: #596273; font-size: 12px; }
      .observatory-link {
        color: #0b5cab;
        font-weight: 700;
        text-decoration: none;
      }
      .observatory-link:hover { text-decoration: underline; }
      .file-link {
        color: #0b5cab;
        text-decoration: none;
        overflow-wrap: anywhere;
      }
      .file-link:hover { text-decoration: underline; }
      .missing-file {
        color: #8a4b00;
        overflow-wrap: anywhere;
      }
      .file-status {
        color: #8a4b00;
        font-size: 12px;
        font-weight: 650;
        margin-left: 6px;
      }
      .primary-link {
        display: inline-block;
        margin: 2px 10px 8px 0;
        padding: 7px 10px;
        border-radius: 6px;
        background: #0b5cab;
        color: #fff !important;
        font-weight: 650;
        text-decoration: none !important;
      }
      .link-note {
        color: #52606d;
        font-size: 12px;
        overflow-wrap: anywhere;
      }
      .doc-status {
        display: inline-block;
        margin-left: 9px;
        color: #52606d;
        font-size: 12px;
        font-weight: 650;
      }
      .doc-status.stale { color: #8a4b00; }
      .report-controls .form-group { margin-bottom: 6px; }
      .report-controls-grid {
        display: grid;
        grid-template-columns: minmax(180px, 280px) minmax(260px, 1fr);
        gap: 12px;
        align-items: end;
      }
      @media (max-width: 900px) {
        .report-controls-grid { grid-template-columns: minmax(0, 1fr); }
      }
      .report-frame-shell {
        padding: 0;
        overflow: hidden;
      }
      .report-open-row {
        padding: 10px 12px 2px 12px;
        border-bottom: 1px solid #edf1f5;
      }
      .report-frame {
        display: block;
        width: 100%;
        height: calc(100vh - 255px);
        min-height: 720px;
        border: 0;
        background: #fff;
      }
      .table-toolbar {
        display: grid;
        grid-template-columns: minmax(190px, 260px) minmax(0, 1fr);
        gap: 12px;
        align-items: end;
        margin-bottom: 10px;
      }
      .table-toolbar .form-group { margin-bottom: 0; }
      .review-workbench-toolbar {
        grid-template-columns: repeat(4, minmax(160px, 1fr));
      }
      .review-workbench-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        align-items: center;
        grid-column: span 3;
      }
      .review-workbench-note {
        grid-column: 1 / -1;
      }
      .workbench-overview {
        display: grid;
        grid-template-columns: minmax(240px, 0.9fr) minmax(320px, 1.4fr);
        gap: 13px;
        align-items: stretch;
      }
      .workbench-overview .panel {
        min-width: 0;
      }
      .workbench-metric-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 8px;
      }
      .workbench-selected-summary h4,
      .workbench-metrics h4 {
        margin-top: 0;
      }
      .workbench-selected-title {
        display: flex;
        flex-wrap: wrap;
        align-items: baseline;
        justify-content: space-between;
        gap: 8px;
      }
      .workbench-route-pill {
        display: inline-block;
        border: 1px solid #cbd5e1;
        border-radius: 999px;
        color: #243b53;
        background: #f8fafc;
        font-size: 12px;
        font-weight: 650;
        padding: 4px 8px;
      }
      .workbench-selected-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        align-items: center;
        margin: 10px 0 8px 0;
      }
      .workbench-detail-grid {
        display: grid;
        grid-template-columns: minmax(120px, 180px) minmax(0, 1fr);
        gap: 5px 12px;
        font-size: 13px;
      }
      .workbench-detail-grid .key {
        color: #596273;
        font-weight: 650;
      }
      .workbench-detail-grid div {
        overflow-wrap: anywhere;
      }
      .label-dashboard-controls {
        display: grid;
        grid-template-columns: minmax(180px, 260px) minmax(170px, 220px)
          auto minmax(0, 1fr);
        gap: 12px;
        align-items: end;
      }
      .label-dashboard-controls .form-group { margin-bottom: 0; }
      .label-dashboard-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
        gap: 13px;
        align-items: start;
      }
      .label-dashboard-grid .panel { min-width: 0; }
      .label-action-row {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        align-items: center;
        margin: 0 0 10px 0;
      }
      .label-file-links {
        display: flex;
        flex-wrap: wrap;
        gap: 8px 14px;
        align-items: center;
      }
      .label-file-links a {
        font-weight: 650;
      }
      .label-status-pill {
        display: inline-block;
        border: 1px solid #cbd5e1;
        border-radius: 999px;
        background: #fff;
        color: #243b53;
        font-size: 11px;
        line-height: 1.2;
        padding: 4px 8px;
        margin: 2px 5px 2px 0;
      }
      @media (max-width: 900px) {
        .table-toolbar { grid-template-columns: minmax(0, 1fr); }
        .workbench-overview { grid-template-columns: minmax(0, 1fr); }
        .workbench-metric-grid {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        .label-dashboard-controls,
        .label-dashboard-grid { grid-template-columns: minmax(0, 1fr); }
      }
      .inspection-source-row {
        display: grid;
        grid-template-columns: minmax(260px, 1fr) minmax(0, 1fr);
        gap: 12px;
        align-items: end;
      }
      .inspection-source-row .form-group { margin-bottom: 0; }
      @media (max-width: 900px) {
        .inspection-source-row { grid-template-columns: minmax(0, 1fr); }
      }
      .inspection-asset-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
        gap: 12px;
        margin-bottom: 13px;
      }
      .inspection-asset {
        background: #fff;
        border: 1px solid #d9dde3;
        border-radius: 8px;
        overflow: hidden;
      }
      .inspection-asset-header {
        padding: 9px 11px;
        border-bottom: 1px solid #edf1f5;
      }
      .inspection-asset-header b {
        display: block;
        color: #243b53;
      }
      .inspection-asset img {
        display: block;
        width: 100%;
        max-height: 560px;
        object-fit: contain;
        background: #fff;
      }
      .inspection-table-stack h4 {
        margin: 18px 0 8px 0;
      }
      .review-label-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(180px, 1fr));
        gap: 12px;
        align-items: end;
      }
      .review-label-grid .form-group { margin-bottom: 0; }
      .review-label-wide {
        grid-column: 1 / -1;
      }
      .review-label-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        align-items: center;
        margin-top: 12px;
      }
      .review-label-status {
        background: #f7fafc;
        border: 1px solid #e5eaf0;
        border-radius: 7px;
        padding: 9px 10px;
        margin-bottom: 12px;
      }
      @media (max-width: 1000px) {
        .review-label-grid { grid-template-columns: minmax(0, 1fr); }
        .review-label-wide { grid-column: auto; }
      }
      .selected-title {
        display: flex;
        align-items: baseline;
        justify-content: space-between;
        gap: 12px;
        flex-wrap: wrap;
      }
      .selected-title h4 { margin-top: 0; }
      .selected-summary-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(132px, 1fr));
        gap: 8px;
        margin: 8px 0 13px 0;
      }
      .selected-summary-item {
        background: #f7fafc;
        border: 1px solid #e5eaf0;
        border-radius: 7px;
        min-height: 58px;
        padding: 9px 10px;
      }
      .selected-summary-item b {
        display: block;
        color: #243b53;
        font-size: 15px;
        line-height: 1.12;
        overflow-wrap: anywhere;
      }
      .selected-summary-item span {
        color: #52606d;
        font-size: 12px;
      }
      .visual-selection-panel .selected-summary-grid {
        margin-bottom: 8px;
      }
	      .visual-selection-actions {
	        display: flex;
	        align-items: center;
	        gap: 8px;
	        flex-wrap: wrap;
	      }
	      .visual-selection-jumps {
	        margin: 8px 0 10px 0;
	      }
	      .visual-selection-jumps .btn {
	        margin-bottom: 4px;
	      }
	      .visual-toolbar-grid {
	        display: grid;
	        grid-template-columns: repeat(3, minmax(160px, 1fr));
	        gap: 12px;
	        align-items: end;
	      }
	      .visual-toolbar-grid .form-group { margin-bottom: 0; }
	      .visual-toolbar-note {
	        color: #52606d;
	        font-size: 12px;
	        margin-top: 8px;
	      }
	      @media (max-width: 900px) {
	        .visual-toolbar-grid { grid-template-columns: minmax(0, 1fr); }
	      }
	      .nav-tabs > li > a { color: #334e68; }
      .page-group > .nav-tabs {
        border-bottom-color: #d9dde3;
        margin-bottom: 8px;
      }
      .page-group > .nav-tabs > li > a {
        font-weight: 650;
      }
      .subtabs > .nav-pills {
        margin-bottom: 10px;
      }
      .subtabs > .nav-pills > li > a {
        border-radius: 6px;
        color: #334e68;
      }
      .subtabs > .nav-pills > li.active > a,
      .subtabs > .nav-pills > li.active > a:focus,
      .subtabs > .nav-pills > li.active > a:hover {
        background: #0b5cab;
      }
      .dataTables_wrapper { font-size: 12px; }
      .shiny-input-container { width: 100%; }
      .checkbox label, .radio label { font-weight: 400; }
    "))
  ),
  div(
    class = "app-header",
    h1("RDG Atlas Browser"),
    div(
      class = "app-subtitle",
      "Review dominant-state RDG candidates, model agreement, context evidence, ",
      "and direct RiboCrypt Observatory links for the best case/control subsets."
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(
        class = "filters-shell",
        div(
          class = "filters-header",
          div(class = "filter-title", "Filters"),
          actionButton("reset", "Reset", class = "btn btn-default btn-sm")
        ),
        div(class = "filter-summary", uiOutput("filter_summary")),
        tags$details(
          open = TRUE,
          tags$summary("Scope"),
          selectizeInput(
            "genes", "Genes (empty = all)",
            choices = safe_choices(cards$gene_symbol),
            selected = character(),
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          ),
          selectizeInput(
            "families", "Design families",
            choices = safe_choices(cards$design_family),
            selected = safe_choices(cards$design_family),
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          ),
          textInput("search", "Text search", value = "")
        ),
        tags$details(
          open = TRUE,
          tags$summary("Preset and Thresholds"),
          selectInput(
            "review_preset", "Preset",
            choices = c(
              "Review queue" = "review_queue",
              "All cards" = "all_cards",
              "Manually verified" = "validated",
              "Branch-usage shifts" = "branch_usage",
              "Context interactions" = "context_interaction",
              "Partial-pooling evidence" = "partial_pooling",
              "Pooling/protocol risk" = "pooling_hazard",
              "Model disagreements" = "model_disagreement"
            ),
            selected = "review_queue"
          ),
          sliderInput(
            "min_conf", "Minimum confidence",
            min = 0, max = 1, value = 0, step = 0.05
          ),
          sliderInput(
            "min_abs_pct", "Minimum absolute clean-CDS effect",
            min = 0, max = 250, value = 0, step = 5, post = "%"
          ),
          checkboxInput("review_only", "Review rows only", value = TRUE)
        ),
        tags$details(
          tags$summary("Evidence Classes"),
          checkboxGroupInput(
            "evidence", "Atlas evidence",
            choices = safe_choices(cards$atlas_evidence_class),
            selected = safe_choices(cards$atlas_evidence_class)
          ),
          checkboxGroupInput(
            "agreement_tiers", "Model agreement",
            choices = safe_choices(cards$model_agreement_tier),
            selected = safe_choices(cards$model_agreement_tier)
          )
        ),
        tags$details(
          tags$summary("Stability Classes"),
          checkboxGroupInput(
            "loo_classes", "Study omission",
            choices = safe_choices(cards$loo_consensus_class),
            selected = safe_choices(cards$loo_consensus_class)
          ),
          checkboxGroupInput(
            "transport_classes", "Grouped omission",
            choices = safe_choices(cards$transport_consensus_class),
            selected = safe_choices(cards$transport_consensus_class)
          )
        ),
        tags$details(
          tags$summary("Advanced Review Flags"),
          checkboxInput("manually_verified_only", "Manually verified only",
                        value = FALSE),
          checkboxInput("branch_usage_shift_only", "Branch-usage shift only",
                        value = FALSE),
          checkboxInput("hierarchical_interval_only",
                        "Hierarchical interval only", value = FALSE),
          checkboxInput("high_residual_only", "High hierarchical residual only",
                        value = FALSE),
          checkboxInput("context_improved_only",
                        "Context calibration improved only", value = FALSE),
          checkboxInput("context_model_improved_only",
                        "Context model improved only", value = FALSE),
          checkboxInput("context_interaction_only", "Context interaction only",
                        value = FALSE),
          checkboxInput("partial_pooling_only", "Partial-pooling evidence only",
                        value = FALSE),
          checkboxInput("pooling_hazard_only",
                        "Pooling/protocol risk only", value = FALSE),
          checkboxInput("model_disagreement_only", "Model disagreement only",
                        value = FALSE)
        )
      )
    ),
    mainPanel(
      width = 9,
      div(class = "panel", uiOutput("metrics")),
      div(
        class = "page-group",
        tabsetPanel(
          id = "page_group",
          tabPanel(
            "Review",
            br(),
            div(
              class = "subtabs",
              tabsetPanel(
                id = "review_pages",
                type = "pills",
                tabPanel(
                  "Cards", br(),
                  div(
                    class = "table-toolbar panel",
                    selectInput(
                      "card_table_view",
                      "Table view",
                      choices = c(
                        "Review summary" = "review",
                        "Model and stability" = "model",
                        "Context evidence" = "context",
                        "Full scan columns" = "full"
                      ),
                      selected = "review"
                    ),
                    div(class = "small-note", uiOutput("card_table_view_note"))
                  ),
                  DTOutput("cards_table")
                ),
                tabPanel(
                  "Label Dashboard", br(),
                  div(
                    class = "panel label-dashboard-controls",
                    selectInput(
                      "label_dashboard_target",
                      "Target label class",
                      choices = label_dashboard_target_choices,
                      selected = "all"
                    ),
                    checkboxInput(
                      "label_dashboard_open_only",
                      "Open rows only",
                      value = TRUE
                    ),
                    actionButton("start_label_dashboard_session",
                                 "Start Label Session",
                                 class = "btn btn-primary"),
                    actionButton("reload_review_label_drafts", "Reload Drafts",
                                 class = "btn btn-default"),
                    div(
                      div(class = "small-note", uiOutput("label_dashboard_note")),
                      uiOutput("label_dashboard_file_links")
                    )
                  ),
                  div(class = "panel", uiOutput("label_dashboard_metrics")),
                  div(
                    class = "label-dashboard-grid",
                    div(
                      class = "panel",
                      h4("Label Balance"),
                      DTOutput("label_balance_table")
                    ),
                    div(
                      class = "panel",
                      h4("Feedback Summary"),
                      DTOutput("label_feedback_summary_table")
                    )
                  ),
                  div(
                    class = "panel",
                    h4("Label Gap Queue"),
                    div(class = "small-note",
                        uiOutput("label_gap_selection_note")),
                    div(
                      class = "label-action-row",
                      actionButton("open_label_gap_review_labels",
                                   "Review Labels",
                                   class = "btn btn-primary btn-sm"),
                      actionButton("open_label_gap_inspection",
                                   "Inspection Assets",
                                   class = "btn btn-default btn-sm"),
                      actionButton("open_label_gap_workbench",
                                   "Review Workbench",
                                   class = "btn btn-default btn-sm")
                    ),
                    DTOutput("label_gap_table")
                  ),
                  div(
                    class = "label-dashboard-grid",
                    div(
                      class = "panel",
                      h4("Current Validation Labels"),
                      DTOutput("validation_label_set_table")
                    ),
                    div(
                      class = "panel",
                      h4("Saved Draft Labels"),
                      DTOutput("label_dashboard_drafts_table")
                    )
                  ),
                  div(
                    class = "label-dashboard-grid",
                    div(
                      class = "panel",
                      h4("Feedback Labels"),
                      DTOutput("review_feedback_labels_table")
                    ),
                    div(
                      class = "panel",
                      h4("Candidate Browser Notes"),
                      DTOutput("review_feedback_candidates_table")
                    )
                  )
                ),
                tabPanel(
                  "Review Workbench", br(),
                  div(
                    class = "table-toolbar panel review-workbench-toolbar",
                    selectInput(
                      "review_workbench_source",
                      "Queue",
                      choices = review_workbench_source_choices,
                      selected = "all"
                    ),
                    selectInput(
                      "review_workbench_target",
                      "Target label",
                      choices = review_workbench_target_choices,
                      selected = "all"
                    ),
                    selectizeInput(
                      "review_workbench_genes",
                      "Workbench genes",
                      choices = review_workbench_gene_choices,
                      selected = character(),
                      multiple = TRUE,
                      options = list(
                        placeholder = "Type gene, e.g. GABARAPL1"
                      )
                    ),
                    selectInput(
                      "review_workbench_scope",
                      "Scope",
                      choices = c(
                        "Current filters" = "filtered",
                        "Selected card" = "selected",
                        "All queued rows" = "all"
                      ),
                      selected = "filtered"
                    ),
                    checkboxInput(
                      "review_workbench_open_only",
                      "Open rows only",
                      value = TRUE
                    ),
                    checkboxInput(
                      "review_workbench_label_task_only",
                      "Label tasks only",
                      value = FALSE
                    ),
                    checkboxInput(
                      "review_workbench_inspection_only",
                      "Inspection packets only",
                      value = FALSE
                    ),
                    checkboxInput(
                      "review_workbench_protocol_only",
                      "Protocol audit only",
                      value = FALSE
                    ),
                    div(
                      class = "review-workbench-actions",
                      actionButton("start_review_workbench_session",
                                   "Start Session",
                                   class = "btn btn-primary"),
                      actionButton("start_selected_review_workbench_session",
                                   "Open Selected",
                                   class = "btn btn-default"),
                      actionButton("use_sidebar_genes_review_workbench",
                                   "Use Sidebar Genes",
                                   class = "btn btn-default"),
                      actionButton("use_selected_gene_review_workbench",
                                   "Use Selected Gene",
                                   class = "btn btn-default"),
                      actionButton("clear_review_workbench_genes",
                                   "Clear Genes",
                                   class = "btn btn-default")
                    ),
                    div(
                      class = "small-note review-workbench-note",
                      uiOutput("review_workbench_note")
                    )
                  ),
                  div(
                    class = "workbench-overview",
                    div(
                      class = "panel workbench-metrics",
                      h4("Current Workbench"),
                      uiOutput("review_workbench_metrics")
                    ),
                    div(
                      class = "panel workbench-selected-summary",
                      uiOutput("review_workbench_selected_summary")
                    )
                  ),
                  DTOutput("review_workbench_table")
                ),
                tabPanel(
                  "Review Labels", br(),
                  div(
                    class = "panel",
                    uiOutput("review_label_task_select_ui"),
                    uiOutput("review_label_status"),
                    uiOutput("review_label_inspection_actions"),
                    uiOutput("review_label_form")
                  ),
                  h4("Saved Draft Labels"),
                  DTOutput("review_label_drafts_table")
                ),
                tabPanel(
                  "Selected Card", br(),
                  div(
                    class = "table-toolbar panel",
                    selectInput(
                      "selected_card_view",
                      "Detail view",
                      choices = c(
                        "Evidence summary" = "summary",
                        "Model and stability" = "model",
                        "Context evidence" = "context",
                        "Full field dump" = "full"
                      ),
                      selected = "summary"
                    ),
                    div(class = "small-note", uiOutput("selected_card_view_note"))
                  ),
                  uiOutput("selected_card")
                ),
                tabPanel("Files", br(), uiOutput("file_panel"))
              )
            )
          ),
          tabPanel(
            "Visuals",
            br(),
            div(
	              class = "subtabs",
	              tabsetPanel(
	                id = "visual_pages",
	                type = "pills",
	                tabPanel("Effect Space", br(),
	                         div(
	                           class = "panel visual-toolbar",
	                           div(
	                             class = "visual-toolbar-grid",
	                             selectInput(
	                               "effect_y", "Y-axis",
	                               choices = effect_y_choices,
	                               selected = "atlas_confidence_score"
	                             ),
	                             selectInput(
	                               "effect_color", "Color",
	                               choices = effect_color_choices,
	                               selected = "atlas_evidence_class"
	                             ),
	                             selectInput(
	                               "effect_size", "Point size",
	                               choices = effect_size_choices,
	                               selected = "none"
	                             )
	                           ),
	                           div(class = "visual-toolbar-note",
	                               uiOutput("effect_space_note"))
	                         ),
	                         plotlyOutput("scatter", height = "620px")),
	                tabPanel("Heatmap", br(),
	                         plotlyOutput("heatmap", height = "720px")),
	                tabPanel("Branch Heatmap", br(),
	                         plotlyOutput("branch_heatmap", height = "720px"))
	              ),
	              uiOutput("visual_selected_card")
	            )
	          ),
          tabPanel(
            "Evidence",
            br(),
            div(
              class = "subtabs",
              tabsetPanel(
                id = "evidence_pages",
                type = "pills",
                tabPanel(
                  "Inspection Assets", br(),
                  div(
                    class = "panel inspection-source-row",
                    uiOutput("inspection_source_select_ui"),
                    div(class = "small-note", uiOutput("inspection_source_note"))
                  ),
                  uiOutput("inspection_assets"),
                  div(
                    class = "inspection-table-stack",
                    h4("Coverage Metrics"),
                    DTOutput("inspection_profile_metrics_table"),
                    h4("Feature Contrasts"),
                    DTOutput("inspection_feature_contrasts_table"),
                    h4("Replicate Summary"),
                    DTOutput("inspection_replicate_summary_table"),
                    h4("Sequence Features"),
                    DTOutput("inspection_sequence_features_table")
                  )
                ),
                tabPanel(
                  "Data Readiness", br(),
                  div(
                    class = "small-note",
                    "Mounted-data audit for the selected candidate and the ",
                    "global required FST pages. This reports whether the ",
                    "backing S-drive coverage, all-merged ORFik files, ",
                    "reference sequence, and static inspection packets are ",
                    "reachable; it is not a biological validation label."
                  ),
                  br(),
                  div(class = "panel", uiOutput("sdrive_readiness_metrics")),
                  div(
                    class = "panel",
                    h4("Selected Card Readiness"),
                    uiOutput("selected_sdrive_readiness")
                  ),
                  DTOutput("selected_sdrive_readiness_table"),
                  h4("Experiment And Reference Checks"),
                  DTOutput("sdrive_status_table"),
                  h4("Qualified Source Rows"),
                  DTOutput("sdrive_transcript_readiness_table"),
                  h4("Required FST Pages"),
                  DTOutput("sdrive_fst_page_readiness_table")
                ),
                tabPanel("Branch Contexts", br(),
                         DTOutput("branch_contexts_table")),
                tabPanel(
                  "Study Evidence", br(),
                  div(
                    class = "small-note",
                    "Per-study branch deltas behind the selected card. ",
                    "LOO support uses within-branch p-values as a robustness proxy, ",
                    "not a new independent validation dataset."
                  ),
                  br(),
                  DTOutput("study_evidence_table")
                ),
                tabPanel(
                  "Transport Evidence", br(),
                  div(
                    class = "small-note",
                    "Grouped omission removes all studies sharing a context, ",
                    "protocol, publication, or lab family. Direction retention ",
                    "tests whether the branch sign survives; support retention ",
                    "tests whether thresholded evidence survives."
                  ),
                  br(),
                  DTOutput("transport_evidence_table")
                ),
                tabPanel(
                  "Context Interactions", br(),
                  div(
                    class = "small-note",
                    "Branch effects in one tissue or cell-line context compared with ",
                    "the remaining viral-infection contexts. Interaction delta is ",
                    "context branch delta minus complement branch delta."
                  ),
                  br(),
                  DTOutput("context_interaction_table")
                ),
                tabPanel(
                  "Partial Pooling", br(),
                  div(
                    class = "small-note",
                    "Shrinkage-aware context rows. Pooled support means the context ",
                    "effect survives empirical-Bayes shrinkage across sibling contexts; ",
                    "single-context rows are reviewable but not pooled evidence."
                  ),
                  br(),
                  DTOutput("partial_pooling_table")
                )
              )
            )
          ),
          tabPanel("Reports", br(), report_browser_ui()),
          tabPanel("Tutorial", br(), uiOutput("tutorial_frame")),
          tabPanel("Scientific Note", br(), uiOutput("scientific_note_frame"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  selected_card_key <- reactiveVal(NULL)
  selected_inspection_source_row <- reactiveVal(NA_integer_)
  selected_review_task_order_hint <- reactiveVal(NA_integer_)
  review_label_drafts_state <- reactiveVal(read_review_label_drafts())
  cards_proxy <- dataTableProxy("cards_table", session = session)
  syncing_cards_table <- reactiveVal(FALSE)

  observeEvent(input$reset, {
    selected_card_key(NULL)
    selected_inspection_source_row(NA_integer_)
    selected_review_task_order_hint(NA_integer_)
    updateSelectInput(session, "review_preset", selected = "review_queue")
    updateSelectizeInput(session, "genes",
                         selected = character())
    updateSelectizeInput(session, "families",
                         selected = safe_choices(cards$design_family))
    updateCheckboxGroupInput(session, "evidence",
                             selected = safe_choices(cards$atlas_evidence_class))
    updateCheckboxGroupInput(session, "agreement_tiers",
                             selected = safe_choices(cards$model_agreement_tier))
    updateCheckboxGroupInput(session, "loo_classes",
                             selected = safe_choices(cards$loo_consensus_class))
    updateCheckboxGroupInput(
      session, "transport_classes",
      selected = safe_choices(cards$transport_consensus_class)
    )
    updateSliderInput(session, "min_conf", value = 0)
    updateSliderInput(session, "min_abs_pct", value = 0)
    updateTextInput(session, "search", value = "")
    updateCheckboxInput(session, "manually_verified_only", value = FALSE)
    updateCheckboxInput(session, "branch_usage_shift_only", value = FALSE)
    updateCheckboxInput(session, "hierarchical_interval_only", value = FALSE)
    updateCheckboxInput(session, "high_residual_only", value = FALSE)
    updateCheckboxInput(session, "context_improved_only", value = FALSE)
    updateCheckboxInput(session, "context_model_improved_only", value = FALSE)
    updateCheckboxInput(session, "context_interaction_only", value = FALSE)
    updateCheckboxInput(session, "partial_pooling_only", value = FALSE)
    updateCheckboxInput(session, "pooling_hazard_only", value = FALSE)
    updateCheckboxInput(session, "model_disagreement_only", value = FALSE)
    updateCheckboxInput(session, "review_only", value = TRUE)
  })

  observeEvent(input$review_preset, {
    preset <- input$review_preset %||% "review_queue"
    updateCheckboxInput(session, "manually_verified_only", value = FALSE)
    updateCheckboxInput(session, "branch_usage_shift_only", value = FALSE)
    updateCheckboxInput(session, "hierarchical_interval_only", value = FALSE)
    updateCheckboxInput(session, "high_residual_only", value = FALSE)
    updateCheckboxInput(session, "context_improved_only", value = FALSE)
    updateCheckboxInput(session, "context_model_improved_only", value = FALSE)
    updateCheckboxInput(session, "context_interaction_only", value = FALSE)
    updateCheckboxInput(session, "partial_pooling_only", value = FALSE)
    updateCheckboxInput(session, "pooling_hazard_only", value = FALSE)
    updateCheckboxInput(session, "model_disagreement_only", value = FALSE)
    updateSliderInput(session, "min_conf", value = 0)
    updateSliderInput(session, "min_abs_pct", value = 0)
    updateCheckboxInput(session, "review_only", value = preset != "all_cards")

    if (identical(preset, "validated")) {
      updateCheckboxInput(session, "review_only", value = FALSE)
      updateCheckboxInput(session, "manually_verified_only", value = TRUE)
    } else if (identical(preset, "branch_usage")) {
      updateCheckboxInput(session, "branch_usage_shift_only", value = TRUE)
    } else if (identical(preset, "context_interaction")) {
      updateCheckboxInput(session, "context_interaction_only", value = TRUE)
    } else if (identical(preset, "partial_pooling")) {
      updateCheckboxInput(session, "partial_pooling_only", value = TRUE)
    } else if (identical(preset, "pooling_hazard")) {
      updateCheckboxInput(session, "pooling_hazard_only", value = TRUE)
    } else if (identical(preset, "model_disagreement")) {
      updateCheckboxInput(session, "model_disagreement_only", value = TRUE)
    }
  }, ignoreInit = TRUE)

  filtered_cards <- reactive({
    dt <- copy(cards)
    if (length(input$genes)) dt <- dt[gene_symbol %chin% input$genes]
    if (length(input$families)) dt <- dt[design_family %chin% input$families]
    if (length(input$evidence)) dt <- dt[atlas_evidence_class %chin% input$evidence]
    if (length(input$agreement_tiers)) {
      dt <- dt[model_agreement_tier %chin% input$agreement_tiers]
    }
    if (length(input$loo_classes)) {
      dt <- dt[loo_consensus_class %chin% input$loo_classes]
    }
    if (length(input$transport_classes)) {
      dt <- dt[transport_consensus_class %chin% input$transport_classes]
    }
    dt <- dt[atlas_confidence_score >= input$min_conf]
    dt <- dt[abs_atlas_clean_cds_pct_change >= input$min_abs_pct]
    if (isTRUE(input$review_only)) {
      dt <- dt[atlas_evidence_class != "exploratory_low_support" |
                 atlas_confidence_score >= 0.45 |
                 (!is.na(sparse_onoff_score) & sparse_onoff_score >= 0.5)]
    }
    if (isTRUE(input$manually_verified_only)) {
      dt <- dt[is.finite(browser_validation_score) &
                 browser_validation_score > 0]
    }
    if (isTRUE(input$branch_usage_shift_only)) {
      keep <- grepl("branch_usage_shift",
                    dt$branch_usage_top_support_class,
                    fixed = TRUE)
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$hierarchical_interval_only)) {
      keep <- dt$hierarchical_evidence_class %chin% c(
        "hierarchical_calibrated_atlas_ready",
        "hierarchical_calibrated_replicated_interval_review",
        "hierarchical_calibrated_single_design_interval_review",
        "hierarchical_atlas_ready",
        "hierarchical_replicated_interval_review",
        "hierarchical_single_design_interval_review"
      )
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$high_residual_only)) {
      keep <- dt$residual_risk_class == "high_residual_heterogeneity"
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$context_improved_only)) {
      keep <- dt$context_calibration_status ==
        "context_reduces_residual_error"
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$context_model_improved_only)) {
      keep <- dt$context_model_status ==
        "context_model_reduces_residual_error"
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$context_interaction_only)) {
      keep <- dt$context_interaction_gene_class %chin% c(
        "context_direction_reversal",
        "context_specific_supported_shift",
        "context_magnitude_modulated"
      )
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$partial_pooling_only)) {
      keep <- dt$partial_pooling_gene_class %chin% c(
        "partial_pooled_context_supported",
        "partial_pooled_context_review",
        "single_context_raw_signal_not_pooled",
        "raw_context_signal_fragile_after_pooling"
      )
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$pooling_hazard_only)) {
      keep <- dt$pooling_hazard_class %chin% c(
        "high_pooling_confounding_hazard",
        "moderate_pooling_confounding_hazard"
      ) |
        dt$pooling_any_start_site_enriched_protocol == TRUE |
        dt$pooling_proof_tier %chin% c(
          "C_pooled_recurrence_or_transport",
          "D_context_specific_or_partial_pooled"
        )
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    if (isTRUE(input$model_disagreement_only)) {
      keep <- dt$model_disagreement_flag == TRUE
      keep[is.na(keep)] <- FALSE
      dt <- dt[keep]
    }
    query <- trimws(input$search)
    if (nzchar(query)) {
      hay <- paste(dt$gene_symbol, dt$design_family, dt$best_case_condition,
                   dt$top_shifted_feature_id, dt$warnings,
                   dt$browser_review_reason, dt$browser_validation_status,
                   dt$browser_validation_notes,
                   dt$branch_usage_top_feature_id,
                   dt$branch_usage_top_support_class,
                   dt$hierarchical_evidence_class,
                   dt$residual_risk_class,
                   dt$residual_worst_context_key,
                   dt$context_calibration_status,
                   dt$context_calibration_top_worse_context,
                   dt$context_model_status,
                   dt$context_model_worst_context,
                   dt$model_agreement_tier,
                   dt$model_support_signature,
                   dt$model_disagreement_flags,
                   dt$model_agreement_top_context,
                   dt$loo_consensus_class,
                   dt$loo_agreement_class,
                   dt$loo_support_robust_branch_names,
                   dt$loo_fragile_branch_names,
                   dt$loo_worst_omitted_studies,
                   dt$transport_consensus_class,
                   dt$transport_agreement_class,
                   dt$transport_evaluable_axis_names,
                   dt$transport_fragile_axis_names,
                   dt$transport_worst_omitted_groups,
                   dt$transport_worst_omitted_studies,
                   dt$context_interaction_gene_class,
                   dt$context_interaction_top_class,
                   dt$context_interaction_top_axis,
                   dt$context_interaction_top_value,
	                   dt$context_interaction_top_branch_class,
	                   dt$context_interaction_context_study_names,
	                   dt$context_interaction_complement_study_names,
	                   dt$context_interaction_summary,
	                   dt$partial_pooling_gene_class,
	                   dt$partial_pooling_top_class,
	                   dt$partial_pooling_top_axis,
	                   dt$partial_pooling_top_value,
	                   dt$partial_pooling_top_branch_class,
	                   dt$partial_pooling_top_context_summary,
	                   dt$pooling_proof_tier,
	                   dt$pooling_claim_language,
	                   dt$pooling_hazard_class,
	                   dt$pooling_hazard_sources,
	                   dt$pooling_top_protocol_classes,
	                   dt$hierarchical_prior_source, sep = " ")
      dt <- dt[grepl(query, hay, ignore.case = TRUE)]
    }
    setorder(dt, -atlas_browser_review_priority_score,
             -atlas_confidence_score, -abs_atlas_clean_cds_pct_change,
             gene_symbol, design_family)
    dt
  })

  output$filter_summary <- renderUI({
    dt <- filtered_cards()
    total_cards <- nrow(cards)
    atlas_fraction <- if (total_cards > 0) {
      paste0(round(100 * nrow(dt) / total_cards, 1), "%")
    } else {
      "NA"
    }
    chips <- list()
    add_chip <- function(label, muted = FALSE) {
      class <- if (isTRUE(muted)) "filter-chip muted" else "filter-chip"
      chips[[length(chips) + 1L]] <<- tags$span(class = class, label)
    }
    add_choice_chip <- function(input_id, all_choices, label) {
      selected <- as.character(input[[input_id]] %||% character())
      if (!setequal(selected, all_choices)) {
        add_chip(paste0(length(selected), "/", length(all_choices), " ", label))
      }
    }
    count_label <- function(n, singular, plural = paste0(singular, "s")) {
      paste(n, if (n == 1) singular else plural)
    }

    preset_labels <- c(
      review_queue = "Review queue",
      all_cards = "All cards",
      validated = "Manually verified",
      branch_usage = "Branch-usage shifts",
      context_interaction = "Context interactions",
      partial_pooling = "Partial pooling",
      pooling_hazard = "Pooling/protocol risk",
      model_disagreement = "Model disagreements"
    )
    preset <- input$review_preset %||% "review_queue"
    preset_label <- unname(preset_labels[preset])
    if (!is_scalar_informative(preset_label)) preset_label <- preset
    add_chip(paste0("Preset: ", preset_label), muted = identical(preset, "review_queue"))

    if (length(input$genes)) add_chip(count_label(length(input$genes), "gene"))
    query <- trimws(input$search %||% "")
    if (nzchar(query)) add_chip(paste0("Search: ", query))
    add_choice_chip("families", safe_choices(cards$design_family), "families")
    add_choice_chip("evidence", safe_choices(cards$atlas_evidence_class),
                    "atlas evidence")
    add_choice_chip("agreement_tiers", safe_choices(cards$model_agreement_tier),
                    "agreement tiers")
    add_choice_chip("loo_classes", safe_choices(cards$loo_consensus_class),
                    "LOO classes")
    add_choice_chip("transport_classes",
                    safe_choices(cards$transport_consensus_class),
                    "transport classes")
    if (!is.null(input$min_conf) && is.finite(input$min_conf) &&
        input$min_conf > 0) {
      add_chip(paste0("Confidence >= ", fmt_num(input$min_conf, 2)))
    }
    if (!is.null(input$min_abs_pct) && is.finite(input$min_abs_pct) &&
        input$min_abs_pct > 0) {
      add_chip(paste0("|CDS effect| >= ", fmt_num(input$min_abs_pct, 0), "%"))
    }
    if (isTRUE(input$review_only)) add_chip("Review rows")
    flag_labels <- c(
      manually_verified_only = "Manually verified",
      branch_usage_shift_only = "Branch-usage shift",
      hierarchical_interval_only = "Hierarchical interval",
      high_residual_only = "High residual",
      context_improved_only = "Context calibration improved",
      context_model_improved_only = "Context model improved",
      context_interaction_only = "Context interaction",
      partial_pooling_only = "Partial pooling",
      pooling_hazard_only = "Pooling/protocol risk",
      model_disagreement_only = "Model disagreement"
    )
    for (flag in names(flag_labels)) {
      if (isTRUE(input[[flag]])) add_chip(unname(flag_labels[flag]))
    }
    if (!length(chips)) add_chip("No active gates", muted = TRUE)

    div(
      div(
        class = "filter-summary-grid",
        div(class = "filter-stat", tags$b(nrow(dt)), span("cards")),
        div(class = "filter-stat",
            tags$b(uniqueN(dt$gene_symbol)), span("genes")),
        div(class = "filter-stat",
            tags$b(uniqueN(dt$design_family)), span("families")),
        div(class = "filter-stat", tags$b(atlas_fraction), span("of atlas"))
      ),
      div(class = "filter-chip-list", chips)
    )
  })

  output$metrics <- renderUI({
    dt <- filtered_cards()
    div(
      class = "metrics-grid",
      div(class = "metric", tags$b(nrow(dt)), span("cards")),
      div(class = "metric", tags$b(uniqueN(dt$gene_symbol)), span("genes")),
      div(class = "metric", tags$b(uniqueN(dt$design_family)), span("design families")),
      div(class = "metric",
          tags$b(sum(dt$atlas_evidence_class ==
                       "atlas_ready_replicated_interval", na.rm = TRUE)),
          span("atlas-ready")),
      div(class = "metric",
          tags$b(fmt_num(median(dt$atlas_confidence_score, na.rm = TRUE))),
          span("median confidence")),
      div(class = "metric",
          tags$b(sum(is.finite(dt$browser_validation_score) &
                       dt$browser_validation_score > 0, na.rm = TRUE)),
          span("manually verified")),
      div(class = "metric",
          tags$b(sum(dt$model_agreement_tier %chin% c(
            "A_browser_validated_consensus",
            "B_replicated_branch_consensus",
            "C_multi_model_consensus"
          ), na.rm = TRUE)),
          span("model consensus")),
      div(class = "metric",
          tags$b(sum(dt$model_agreement_tier ==
                       "D_model_disagreement_review", na.rm = TRUE)),
          span("model disagreements")),
      div(class = "metric",
          tags$b(sum(dt$loo_consensus_class ==
                       "consensus_loo_support_robust", na.rm = TRUE)),
          span("LOO support-robust")),
      div(class = "metric",
          tags$b(sum(dt$loo_consensus_class ==
                       "consensus_loo_fragile", na.rm = TRUE)),
          span("LOO-fragile consensus")),
      div(class = "metric",
          tags$b(sum(dt$transport_consensus_class ==
                       "consensus_transport_support_robust", na.rm = TRUE)),
          span("transport robust")),
      div(class = "metric",
          tags$b(sum(dt$transport_consensus_class ==
                       "consensus_transport_fragile", na.rm = TRUE)),
          span("transport fragile")),
      div(class = "metric",
          tags$b(sum(grepl("branch_usage_shift",
                           dt$branch_usage_top_support_class,
                           fixed = TRUE), na.rm = TRUE)),
          span("branch shifts")),
      div(class = "metric",
          tags$b(sum(dt$hierarchical_evidence_class %chin% c(
            "hierarchical_calibrated_atlas_ready",
            "hierarchical_calibrated_replicated_interval_review",
            "hierarchical_calibrated_single_design_interval_review",
            "hierarchical_atlas_ready",
            "hierarchical_replicated_interval_review",
            "hierarchical_single_design_interval_review"
          ), na.rm = TRUE)),
          span("hierarchical intervals")),
      div(class = "metric",
          tags$b(sum(dt$residual_risk_class ==
                       "high_residual_heterogeneity", na.rm = TRUE)),
          span("high residual")),
      div(class = "metric",
          tags$b(sum(dt$context_calibration_status ==
                       "context_reduces_residual_error", na.rm = TRUE)),
          span("context improved")),
      div(class = "metric",
          tags$b(sum(dt$context_model_status ==
                       "context_model_reduces_residual_error",
                     na.rm = TRUE)),
          span("context model improved")),
	      div(class = "metric",
	          tags$b(sum(dt$context_interaction_gene_class %chin% c(
	            "context_direction_reversal",
	            "context_specific_supported_shift",
	            "context_magnitude_modulated"
	          ), na.rm = TRUE)),
	          span("context branch shifts")),
	      div(class = "metric",
	          tags$b(sum(dt$partial_pooling_gene_class ==
	                       "partial_pooled_context_supported", na.rm = TRUE)),
	          span("pooled context")),
	      div(class = "metric",
	          tags$b(sum(dt$partial_pooling_gene_class ==
	                       "single_context_raw_signal_not_pooled",
	                     na.rm = TRUE)),
	          span("single context")),
	      div(class = "metric",
	          tags$b(sum(dt$pooling_proof_tier ==
	                       "A_exact_strict_project_matched_control",
	                     na.rm = TRUE)),
	          span("exact-control proof")),
	      div(class = "metric",
	          tags$b(sum(dt$pooling_hazard_class %chin% c(
	            "high_pooling_confounding_hazard",
	            "moderate_pooling_confounding_hazard"
	          ), na.rm = TRUE)),
	          span("pooling risk"))
	    )
  })

  empty_dt <- function(message) data.table(message = message)

  label_file_link <- function(path, label) {
    info <- make_selected_file_link(path)
    if (is.null(info) || !isTRUE(info$exists) ||
        !is_scalar_informative(info$href)) {
      return(tags$span(class = "missing-file", paste(label, "missing")))
    }
    tags$a(
      class = "file-link",
      href = info$href,
      target = "_blank",
      rel = "noopener noreferrer",
      label
    )
  }

  focus_review_context <- function(row) {
    if (is.null(row) || !nrow(row)) return(invisible(FALSE))
    if ("card_key" %in% names(row) && is_scalar_informative(row$card_key[1])) {
      selected_card_key(row$card_key[1])
    } else if (all(c("gene_symbol", "tx_id", "design_family") %in% names(row))) {
      selected_card_key(paste(row$gene_symbol[1], row$tx_id[1],
                              row$design_family[1], sep = "\r"))
    }
    active_order <- if ("active_review_order" %in% names(row)) {
      suppressWarnings(as.integer(row$active_review_order[1]))
    } else {
      NA_integer_
    }
    if (is.finite(active_order) && !is.na(active_order)) {
      selected_review_task_order_hint(active_order)
      updateSelectInput(session, "review_label_task",
                        selected = as.character(active_order))
    }
    source_row <- if ("inspection_source_row" %in% names(row)) {
      suppressWarnings(as.integer(row$inspection_source_row[1]))
    } else {
      NA_integer_
    }
    if (is.finite(source_row) && !is.na(source_row)) {
      selected_inspection_source_row(source_row)
      updateSelectInput(session, "inspection_source_row",
                        selected = as.character(source_row))
    }
    invisible(TRUE)
  }

  label_dashboard_gap_rows <- reactive({
    dt <- copy(validation_label_gap_queue)
    if (!nrow(dt)) return(dt)
    dt <- attach_review_task_lookup(dt)
    dt <- add_card_key(dt)
    drafts <- review_label_drafts_state()
    if (nrow(drafts)) {
      draft_status <- unique(drafts[, .(
        active_review_order = suppressWarnings(as.integer(global_review_order)),
        draft_label_class = review_label_class_preview,
        draft_saved_at,
        draft_next_action = next_action
      )], by = "active_review_order")
      dt <- merge(dt, draft_status, by = "active_review_order",
                  all.x = TRUE, sort = FALSE)
    } else {
      dt[, draft_label_class := ""]
      dt[, draft_saved_at := ""]
      dt[, draft_next_action := ""]
    }
    if (nrow(review_feedback_labels) &&
        "global_review_order" %in% names(review_feedback_labels)) {
      feedback_status <- unique(review_feedback_labels[, .(
        active_review_order = suppressWarnings(as.integer(global_review_order)),
        feedback_label_class = review_label_class,
        feedback_reviewer = reviewer,
        feedback_review_date = review_date
      )], by = "active_review_order")
      dt <- merge(dt, feedback_status, by = "active_review_order",
                  all.x = TRUE, sort = FALSE)
    } else {
      dt[, feedback_label_class := ""]
      dt[, feedback_reviewer := ""]
      dt[, feedback_review_date := ""]
    }
    for (column in c("draft_label_class", "draft_saved_at",
                     "draft_next_action", "feedback_label_class",
                     "feedback_reviewer", "feedback_review_date")) {
      if (!column %in% names(dt)) dt[, (column) := ""]
    }
    balance <- label_balance_status(drafts)
    if (nrow(balance) && "target_label_class" %in% names(dt)) {
      balance_keep <- intersect(c(
        "label_class", "current_count", "pending_draft_count",
        "label_deficit", "remaining_after_drafts", "acquisition_priority"
      ), names(balance))
      balance <- balance[, ..balance_keep]
      setnames(
        balance,
        old = setdiff(names(balance), "label_class"),
        new = paste0("target_", setdiff(names(balance), "label_class"))
      )
      setnames(balance, "label_class", "target_label_class")
      dt <- merge(dt, balance, by = "target_label_class",
                  all.x = TRUE, sort = FALSE)
    }
    target <- input$label_dashboard_target %||% "all"
    if (!identical(target, "all")) {
      dt <- dt[
        clean_text(target_label_class) == target |
          clean_text(label_gap_class) == target |
          clean_text(draft_label_class) == target |
          clean_text(feedback_label_class) == target
      ]
    }
    if (isTRUE(input$label_dashboard_open_only)) {
      already <- as.logical(dt_col(dt, "already_labeled", FALSE))
      already[is.na(already)] <- FALSE
      dt <- dt[
        !already &
          !is_informative_text(draft_label_class) &
          !is_informative_text(feedback_label_class)
      ]
    }
    setorder(dt, label_gap_rank, -label_acquisition_score)
    dt
  })

  selected_label_gap_row <- reactive({
    dt <- label_dashboard_gap_rows()
    if (!nrow(dt)) return(data.table())
    idx <- input$label_gap_table_rows_selected
    if (length(idx) != 1L || idx > nrow(dt)) idx <- 1L
    dt[idx]
  })

  observeEvent(input$label_gap_table_rows_selected, {
    focus_review_context(selected_label_gap_row())
  }, ignoreInit = TRUE)

  observeEvent(input$reload_review_label_drafts, {
    review_label_drafts_state(read_review_label_drafts())
    showNotification("Draft review labels reloaded.", type = "message")
  }, ignoreInit = TRUE)

  output$label_dashboard_note <- renderUI({
    dt <- label_dashboard_gap_rows()
    tagList(
      tags$b(nrow(dt)), " displayed gap rows from ",
      tags$b(nrow(validation_label_gap_queue)), " total active-review gaps. ",
      "Selecting a row focuses the selected card and active-review task."
    )
  })

  output$label_dashboard_file_links <- renderUI({
    div(
      class = "label-file-links",
      label_file_link(review_label_drafts_file, "Draft CSV"),
      label_file_link(file.path(
        analysis_dir, "dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_gap_queue.csv"
      ), "Gap queue"),
      label_file_link(file.path(
        analysis_dir, "dominant_rdg_validation_label_bridge",
        "dominant_rdg_validation_label_balance.csv"
      ), "Label balance"),
      label_file_link(file.path(
        analysis_dir, "dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_labels.csv"
      ), "Feedback labels"),
      label_file_link(file.path(
        analysis_dir, "dominant_rdg_review_feedback",
        "dominant_rdg_review_feedback_browser_validation_notes_candidate.csv"
      ), "Candidate notes")
    )
  })

  output$label_dashboard_metrics <- renderUI({
    drafts <- review_label_drafts_state()
    gaps <- label_dashboard_gap_rows()
    missing_classes <- if (nrow(validation_label_balance) &&
                           "label_deficit" %in% names(validation_label_balance)) {
      sum(suppressWarnings(as.numeric(validation_label_balance$label_deficit)) > 0,
          na.rm = TRUE)
    } else {
      NA_integer_
    }
    div(
      class = "metrics-grid",
      div(class = "metric",
          tags$b(nrow(validation_label_set)),
          span("current bridge labels")),
      div(class = "metric",
          tags$b(nrow(review_feedback_labels)),
          span("ingested feedback labels")),
      div(class = "metric",
          tags$b(sum(is_informative_text(
            dt_col(drafts, "review_label_class_preview", "")
          ))),
          span("saved draft labels")),
      div(class = "metric",
          tags$b(nrow(gaps)),
          span("displayed open gaps")),
      div(class = "metric",
          tags$b(missing_classes),
          span("label classes below target")),
      div(class = "metric",
          tags$b(nrow(review_feedback_candidates)),
          span("candidate browser notes"))
    )
  })

  output$label_balance_table <- renderDT({
    balance <- label_balance_status(review_label_drafts_state())
    if (!nrow(balance)) {
      return(datatable(empty_dt("No label-balance table found."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    keep <- intersect(c(
      "label_class", "acquisition_role", "target_minimum",
      "current_count", "pending_draft_count", "label_deficit",
      "remaining_after_drafts", "deficit_fraction", "acquisition_priority"
    ), names(balance))
    out <- round_numeric_columns(balance[, ..keep], 3)
    datatable(out, rownames = FALSE, filter = "top",
              class = "stripe hover compact",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$label_feedback_summary_table <- renderDT({
    out <- copy(review_feedback_summary)
    if (!nrow(out)) out <- empty_dt("No review-feedback summary found.")
    datatable(out, rownames = FALSE, class = "stripe hover compact",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$label_gap_selection_note <- renderUI({
    row <- selected_label_gap_row()
    if (!nrow(row)) return("No label-gap row selected.")
    pills <- tagList(
      tags$span(class = "label-status-pill",
                paste0("rank ", row$label_gap_rank[1])),
      tags$span(class = "label-status-pill",
                scalar_text(row$target_label_class)),
      tags$span(class = "label-status-pill",
                scalar_text(row$label_gap_class)),
      if (is_informative_text(row$draft_label_class[1])) {
        tags$span(class = "label-status-pill",
                  paste("draft", row$draft_label_class[1]))
      } else {
        NULL
      }
    )
    tagList(
      tags$b(paste(row$gene_symbol[1], row$design_family[1], sep = " / ")),
      " | active task #", row$active_review_order[1],
      tags$br(),
      pills,
      tags$br(),
      tags$span(row$review_reason[1])
    )
  })

  output$label_gap_table <- renderDT({
    dt <- copy(label_dashboard_gap_rows())
    if (!nrow(dt)) {
      return(datatable(empty_dt("No label-gap rows match the current filter."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    dt <- make_queue_observatory_links(dt)
    keep <- intersect(c(
      "Links", "label_gap_rank", "active_review_order",
      "draft_label_class", "draft_saved_at", "feedback_label_class",
      "gene_symbol", "tx_id", "design_family", "batch_id", "batch_role",
      "label_gap_class", "target_label_class", "label_acquisition_score",
      "target_label_deficit", "target_pending_draft_count",
      "target_remaining_after_drafts",
      "validation_label_class", "motif_prior_review_class",
      "motif_top_concordance_class", "atlas_evidence_class",
      "atlas_clean_cds_pct_change", "consensus_class", "consensus_score",
      "review_reason", "observatory_link_note"
    ), names(dt))
    out <- round_numeric_columns(dt[, ..keep], 4)
    datatable(
      out,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      escape = FALSE,
      class = "stripe hover compact",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })

  output$validation_label_set_table <- renderDT({
    out <- copy(validation_label_set)
    if (!nrow(out)) {
      return(datatable(empty_dt("No validation labels found."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    keep <- intersect(c(
      "label_source", "gene_symbol", "tx_id", "design_family",
      "label_class", "label_status", "label_score", "label_weight",
      "reviewer", "review_date", "review_notes", "global_review_order",
      "batch_id", "batch_role"
    ), names(out))
    out <- round_numeric_columns(out[, ..keep], 4)
    datatable(out, rownames = FALSE, filter = "top",
              class = "stripe hover compact",
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$label_dashboard_drafts_table <- renderDT({
    out <- copy(review_label_drafts_state())
    if (!nrow(out)) {
      return(datatable(empty_dt("No draft labels saved yet."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    keep <- intersect(c(
      "global_review_order", "review_label_class_preview", "gene_symbol",
      "tx_id", "design_family", "batch_id", "batch_role",
      "target_label_class", "label_gap_class", "protocol_pooling_need",
      "pooling_hazard_class",
      "manual_review_status", "browser_signal_grade",
      "metadata_confounder_status", "negative_control_result",
      "reviewer", "review_date", "next_action", "draft_saved_at",
      "reviewer_notes"
    ), names(out))
    datatable(out[, ..keep], rownames = FALSE, filter = "top",
              class = "stripe hover compact",
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$review_feedback_labels_table <- renderDT({
    out <- copy(review_feedback_labels)
    if (!nrow(out)) {
      return(datatable(empty_dt("No ingested feedback labels yet."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    keep <- intersect(c(
      "global_review_order", "review_label_class", "calibration_weight",
      "usable_for_training", "gene_symbol", "tx_id", "design_family",
      "batch_id", "batch_role", "target_label_class", "label_gap_class",
      "protocol_pooling_need", "pooling_hazard_class", "manual_review_status",
      "browser_signal_grade", "metadata_confounder_status",
      "negative_control_result", "reviewer", "review_date", "next_action",
      "reviewer_notes"
    ), names(out))
    out <- round_numeric_columns(out[, ..keep], 4)
    datatable(out, rownames = FALSE, filter = "top",
              class = "stripe hover compact",
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$review_feedback_candidates_table <- renderDT({
    out <- copy(review_feedback_candidates)
    if (!nrow(out)) {
      return(datatable(empty_dt("No candidate browser-validation notes yet."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    keep <- intersect(c(
      "gene_symbol", "tx_id", "design_family", "browser_validation_status",
      "browser_validation_label", "validation_score", "global_review_order",
      "batch_id", "review_label_class", "calibration_weight",
      "validated_by", "validation_date", "validation_notes"
    ), names(out))
    out <- round_numeric_columns(out[, ..keep], 4)
    datatable(out, rownames = FALSE, filter = "top",
              class = "stripe hover compact",
              options = list(pageLength = 8, scrollX = TRUE))
  })

  output$report_select_ui <- renderUI({
    report_category <- input$report_category %||% "dominant_rdgs"
    selectizeInput(
      "report_path",
      "Report",
      choices = report_choice_groups(report_category),
      selected = first_report_choice(report_category),
      options = list(placeholder = "Type to find a report")
    )
  })

  output$report_category_note <- renderUI({
    report_category <- input$report_category %||% "dominant_rdgs"
    report_category_note(report_category)
  })

  output$report_frame <- renderUI({
    render_report_frame(input$report_path)
  })

  output$tutorial_frame <- renderUI({
    render_static_html_frame(
      "inst/shiny/rdg_atlas_browser/docs/tutorial.html",
      source_rel_path = "md/dominant_rdg_atlas_browser_tutorial.md",
      missing_message = "Tutorial HTML is missing. Regenerate it from the tutorial Markdown.",
      open_label = "Open tutorial"
    )
  })

  output$scientific_note_frame <- renderUI({
    render_static_html_frame(
      "inst/shiny/rdg_atlas_browser/docs/scientific_note.html",
      source_rel_path = "md/dominant_cell_state_scientific_note.md",
      missing_message = "Scientific-note HTML is missing. Regenerate it from the scientific note Markdown.",
      open_label = "Open scientific note"
    )
  })

  card_table_view_note <- function(view) {
    switch(
      view,
      review = paste(
        "Compact index for picking cards. Use the detail tabs after selecting",
        "a row."
      ),
      model = paste(
        "Model agreement, study-omission, transport, hierarchy, and residual",
        "diagnostics for the same filtered cards."
      ),
      context = paste(
        "Branch-usage, context-interaction, partial-pooling, and context-model",
        "signals for the same filtered cards."
      ),
      full = paste(
        "Broader scan view with the main review columns. Selected Card keeps",
        "the exhaustive field dump."
      ),
      ""
    )
  }

  card_table_columns <- function(view) {
    base <- c(
      "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
      "atlas_confidence_score", "atlas_browser_review_priority_score",
      "atlas_clean_cds_pct_change",
      "atlas_clean_cds_pct_lower", "atlas_clean_cds_pct_upper",
      "strict_supported_designs", "strict_supported_studies",
      "top_shifted_feature_id", "top_shifted_feature_type",
      "browser_validation_status", "best_case_condition", "best_study"
    )
    review <- c(
      base,
      "model_agreement_tier", "loo_consensus_class",
      "transport_consensus_class", "context_interaction_gene_class",
      "partial_pooling_gene_class", "hierarchical_evidence_class",
      "pooling_proof_tier", "pooling_hazard_class",
      "residual_risk_class", "warnings"
    )
    model <- c(
      base,
      "model_agreement_tier", "model_agreement_class",
      "model_agreement_score", "model_agreement_review_priority",
      "model_support_signature", "model_disagreement_flags",
      "loo_consensus_class", "loo_stability_score",
      "transport_consensus_class", "transport_stability_score",
      "hierarchical_evidence_class", "hierarchical_priority_score",
      "hierarchical_clean_cds_pct_change",
      "hierarchical_prior_source", "hierarchical_evidence_weight",
      "residual_risk_class", "context_calibration_status",
      "context_model_status", "warnings"
    )
    context <- c(
      base,
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_hazard_class", "pooling_hazard_score",
      "pooling_hazard_sources", "pooling_top_protocol_classes",
      "branch_usage_top_feature_id", "branch_usage_top_feature_class",
      "branch_usage_top_support_class",
      "branch_usage_top_weighted_usage_delta",
      "context_interaction_gene_class",
      "context_interaction_top_class",
      "context_interaction_top_axis",
      "context_interaction_top_value",
      "context_interaction_top_branch_class",
      "context_interaction_delta",
      "context_interaction_delta_q",
      "partial_pooling_gene_class",
      "partial_pooling_top_class",
      "partial_pooling_top_axis",
      "partial_pooling_top_value",
      "partial_pooling_top_branch_class",
      "partial_pooling_top_delta",
      "partial_pooling_top_q",
      "context_calibration_status", "context_model_status",
      "warnings"
    )
    full <- c(
      base,
      "branch_usage_top_feature_id", "branch_usage_top_feature_class",
      "branch_usage_top_support_class",
      "branch_usage_top_weighted_usage_delta",
      "model_agreement_tier", "model_agreement_class",
      "model_agreement_score", "model_agreement_review_priority",
      "model_support_signature", "model_disagreement_flags",
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_hazard_class", "pooling_hazard_score",
      "pooling_hazard_sources", "pooling_top_protocol_classes",
      "loo_consensus_class", "loo_stability_score",
      "transport_consensus_class", "transport_stability_score",
      "context_interaction_gene_class",
      "context_interaction_top_class",
      "context_interaction_top_axis",
      "context_interaction_top_value",
      "context_interaction_top_branch_class",
	      "context_interaction_delta",
	      "context_interaction_delta_q",
	      "partial_pooling_gene_class",
	      "partial_pooling_top_class",
	      "partial_pooling_top_axis",
	      "partial_pooling_top_value",
	      "partial_pooling_top_branch_class",
	      "partial_pooling_top_delta",
	      "partial_pooling_top_q",
	      "hierarchical_evidence_class", "hierarchical_priority_score",
      "hierarchical_clean_cds_pct_change",
      "hierarchical_prior_source", "hierarchical_evidence_weight",
      "residual_risk_class", "context_calibration_status",
      "context_model_status",
      "warnings"
    )
    unique(switch(
      view,
      model = model,
      context = context,
      full = full,
      review
    ))
  }

  output$card_table_view_note <- renderUI({
    card_table_view_note(input$card_table_view %||% "review")
  })

  table_dt <- reactive({
    dt <- filtered_cards()
    keep <- card_table_columns(input$card_table_view %||% "review")
    keep <- intersect(keep, names(dt))
    out <- dt[, ..keep]
    link_labels <- as.character(htmltools::htmlEscape(dt$gene_symbol))
    link_urls <- as.character(htmltools::htmlEscape(dt$observatory_url,
                                                    attribute = TRUE))
    link_notes <- as.character(htmltools::htmlEscape(dt$observatory_link_note,
                                                     attribute = TRUE))
    out[, Links := sprintf(
      "<a class=\"observatory-link\" href=\"%s\" target=\"_blank\" rel=\"noopener noreferrer\" title=\"%s\">%s</a>",
      link_urls, link_notes, link_labels
    )]
    if ("gene_symbol" %in% names(out)) out[, gene_symbol := NULL]
    setcolorder(out, c("Links", setdiff(names(out), "Links")))
    pct_cols <- intersect(c("atlas_clean_cds_pct_change",
                            "atlas_clean_cds_pct_lower",
                            "atlas_clean_cds_pct_upper",
                            "hierarchical_clean_cds_pct_change",
                            "hierarchical_calibrated_clean_cds_pct_lower",
                            "hierarchical_calibrated_clean_cds_pct_upper"), names(out))
    for (col in pct_cols) set(out, j = col, value = round(out[[col]], 2))
    num_cols <- intersect(c("atlas_confidence_score",
                            "atlas_browser_review_priority_score",
                            "atlas_direction_fraction", "atlas_i2",
                            "sparse_onoff_score",
                            "branch_usage_top_score",
                            "branch_usage_top_weighted_usage_delta",
                            "branch_usage_top_q",
                            "grouped_branch_usage_top_weighted_usage_delta",
                            "count_aware_top_weighted_usage_delta",
                            "joint_top_branch_delta",
                            "hierarchical_joint_top_delta",
                            "hierarchical_joint_top_priority",
                            "hierarchical_joint_l1_effect",
                            "dm_top_delta",
                            "dm_top_priority",
                            "dm_l1_effect",
                            "model_agreement_score",
                            "model_agreement_review_priority",
                            "loo_stability_score",
                            "loo_review_priority",
                            "loo_min_direction_fraction",
                            "loo_min_support_retention",
                            "transport_stability_score",
                            "transport_review_priority",
                            "transport_min_direction_fraction",
                            "transport_min_support_retention",
                            "transport_min_remaining_studies",
                            "transport_coverage_fraction",
                            "transport_axis_coverage_fraction",
                            "context_interaction_top_score",
                            "context_interaction_review_priority",
	                            "context_interaction_context_delta",
	                            "context_interaction_complement_delta",
	                            "context_interaction_delta",
	                            "context_interaction_delta_q",
	                            "partial_pooling_review_priority",
	                            "partial_pooling_component",
	                            "partial_pooling_top_delta",
	                            "partial_pooling_top_q",
	                            "partial_pooling_top_weight",
	                            "partial_pooling_top_shrinkage_fraction",
	                            "pooling_proof_component",
	                            "pooling_hazard_score",
	                            "pooling_exact_context_fraction",
	                            "pooling_joint_strict_fraction",
	                            "hierarchical_priority_score",
                            "hierarchical_clean_cds_pct_change",
                            "hierarchical_evidence_weight",
                            "residual_heterogeneity_score",
                            "residual_rmse_log2fc",
                            "context_calibration_rmse_delta_log2fc",
                            "context_calibration_fraction_improved",
                            "context_model_rmse_delta_log2fc",
                            "context_model_fraction_improved"), names(out))
    for (col in num_cols) set(out, j = col, value = round(out[[col]], 3))
    out
  })

  output$cards_table <- renderDT({
    dt <- table_dt()
    priority_col <- match("atlas_browser_review_priority_score", names(dt)) - 1L
    if (!is.finite(priority_col)) priority_col <- 4L
    datatable(
      dt,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      escape = FALSE,
      class = "stripe hover compact",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(priority_col, "desc"))
      )
    )
  })

  observeEvent(input$cards_table_rows_selected, {
    if (isTRUE(syncing_cards_table())) return()
    dt <- filtered_cards()
    idx <- input$cards_table_rows_selected
    if (!nrow(dt) || length(idx) != 1 || idx > nrow(dt)) return()
    key <- dt$card_key[idx]
    if (!identical(key, selected_card_key())) selected_card_key(key)
  }, ignoreInit = TRUE)

  observeEvent(selected_card_key(), {
    key <- selected_card_key()
    dt <- filtered_cards()
    idx <- if (is_scalar_informative(key) && nrow(dt)) {
      which(dt$card_key == key)[1]
    } else {
      NA_integer_
    }
    syncing_cards_table(TRUE)
    on.exit(syncing_cards_table(FALSE), add = TRUE)
    if (is.finite(idx) && !is.na(idx)) {
      selectRows(cards_proxy, idx)
    } else {
      selectRows(cards_proxy, NULL)
    }
  }, ignoreInit = TRUE)

  review_workbench_filtered <- reactive({
    dt <- copy(review_workbench)
    if (!nrow(dt)) return(dt)
    drafts <- review_label_drafts_state()
    if (nrow(drafts)) {
      draft_status <- unique(drafts[, .(
        active_review_order = global_review_order,
        draft_label_class = review_label_class_preview,
        draft_saved_at
      )], by = "active_review_order")
      dt <- merge(dt, draft_status, by = "active_review_order",
                  all.x = TRUE, sort = FALSE)
    } else {
      dt[, draft_label_class := ""]
      dt[, draft_saved_at := ""]
    }
    if (nrow(review_feedback_labels) &&
        "global_review_order" %in% names(review_feedback_labels)) {
      feedback_status <- unique(review_feedback_labels[, .(
        active_review_order = suppressWarnings(as.integer(global_review_order)),
        feedback_label_class = review_label_class,
        feedback_review_date = review_date
      )], by = "active_review_order")
      dt <- merge(dt, feedback_status, by = "active_review_order",
                  all.x = TRUE, sort = FALSE)
    } else {
      dt[, feedback_label_class := ""]
      dt[, feedback_review_date := ""]
    }
    for (column in c("draft_label_class", "draft_saved_at",
                     "feedback_label_class", "feedback_review_date")) {
      if (!column %in% names(dt)) dt[, (column) := ""]
    }
    source <- input$review_workbench_source %||% "all"
    if (!identical(source, "all")) dt <- dt[queue_source == source]
    workbench_genes <- as.character(input$review_workbench_genes %||%
                                      character())
    workbench_genes <- workbench_genes[is_informative_text(workbench_genes)]
    if (length(workbench_genes)) {
      dt <- dt[gene_symbol %chin% workbench_genes]
    }
    balance <- label_balance_status(drafts)
    if (nrow(balance) && "target_label_class" %in% names(dt)) {
      balance_keep <- intersect(c(
        "label_class", "current_count", "pending_draft_count",
        "label_deficit", "remaining_after_drafts", "acquisition_priority"
      ), names(balance))
      balance <- balance[, ..balance_keep]
      setnames(
        balance,
        old = setdiff(names(balance), "label_class"),
        new = paste0("target_", setdiff(names(balance), "label_class"))
      )
      setnames(balance, "label_class", "target_label_class")
      dt <- merge(dt, balance, by = "target_label_class",
                  all.x = TRUE, sort = FALSE)
    }
    target <- input$review_workbench_target %||% "all"
    if (!identical(target, "all")) {
      for (column in c("target_label_class", "label_gap_class",
                       "draft_label_class", "feedback_label_class")) {
        if (!column %in% names(dt)) dt[, (column) := ""]
      }
      dt <- dt[
        clean_text(target_label_class) == target |
          clean_text(label_gap_class) == target |
          clean_text(draft_label_class) == target |
          clean_text(feedback_label_class) == target
      ]
    }
    if (isTRUE(input$review_workbench_open_only) && "is_open" %in% names(dt)) {
      dt <- dt[is.na(is_open) | is_open == TRUE]
      dt <- dt[
        !is_informative_text(draft_label_class) &
          !is_informative_text(feedback_label_class)
      ]
    }
    if (isTRUE(input$review_workbench_label_task_only)) {
      dt <- dt[workbench_finite_integer_flag(dt, "active_review_order")]
    }
    if (isTRUE(input$review_workbench_inspection_only)) {
      dt <- dt[workbench_finite_integer_flag(dt, "inspection_source_row")]
    }
    if (isTRUE(input$review_workbench_protocol_only)) {
      dt <- dt[review_workbench_protocol_audit_mask(dt)]
    }
    scope <- input$review_workbench_scope %||% "filtered"
    if (identical(scope, "filtered")) {
      dt <- dt[card_key %chin% filtered_cards()$card_key]
    } else if (identical(scope, "selected")) {
      row <- selected_row()
      if (is.null(row)) {
        dt <- dt[0]
      } else {
        dt <- dt[card_key == row$card_key[1]]
      }
    }
    dt[, review_destination := vapply(
      review_workbench_row_destination(dt),
      review_workbench_destination_label,
      character(1)
    )]
    setorder(dt, queue_source, queue_rank, -queue_score)
    dt
  })

  output$review_workbench_note <- renderUI({
    dt <- review_workbench_filtered()
    target <- input$review_workbench_target %||% "all"
    target_label <- if (!identical(target, "all")) {
      paste0(" Target: ", target, ". ")
    } else {
      ""
    }
    scope <- input$review_workbench_scope %||% "filtered"
    workbench_genes <- as.character(input$review_workbench_genes %||%
                                      character())
    workbench_genes <- unique(workbench_genes[
      is_informative_text(workbench_genes)
    ])
    sidebar_genes <- as.character(input$genes %||% character())
    sidebar_genes <- unique(sidebar_genes[is_informative_text(sidebar_genes)])
    short_gene_list <- function(x, max_items = 4L) {
      x <- unique(x[is_informative_text(x)])
      if (!length(x)) return("")
      label <- paste(utils::head(x, max_items), collapse = ", ")
      if (length(x) > max_items) {
        label <- paste0(label, " +", length(x) - max_items, " more")
      }
      label
    }
    scope_label <- switch(
      scope,
      filtered = "Current filters",
      selected = "Selected card",
      all = "All queued rows",
      scope
    )
    local_gates <- c(
      if (isTRUE(input$review_workbench_open_only)) "open rows",
      if (isTRUE(input$review_workbench_label_task_only)) "label tasks",
      if (isTRUE(input$review_workbench_inspection_only)) "inspection packets",
      if (isTRUE(input$review_workbench_protocol_only)) "protocol audits"
    )
    gate_label <- if (length(local_gates)) {
      paste0(" Gates: ", paste(local_gates, collapse = ", "), ". ")
    } else {
      ""
    }
    gene_scope_class <- ""
    gene_scope <- if (length(workbench_genes)) {
      gene_label <- short_gene_list(workbench_genes)
      if (identical(scope, "filtered") && length(sidebar_genes)) {
        paste0(
          "Gene scope: Workbench genes (", gene_label,
          ") intersect current sidebar filters, including sidebar Genes (",
          short_gene_list(sidebar_genes), ")."
        )
      } else if (identical(scope, "filtered")) {
        paste0(
          "Gene scope: Workbench genes (", gene_label,
          ") within current sidebar filters."
        )
      } else if (identical(scope, "selected")) {
        paste0(
          "Gene scope: Workbench genes (", gene_label,
          ") within the selected-card scope."
        )
      } else {
        paste0(
          "Gene scope: Workbench genes (", gene_label,
          ") across all queued rows."
        )
      }
    } else if (identical(scope, "filtered") && length(sidebar_genes)) {
      paste0("Gene scope: sidebar Genes (", short_gene_list(sidebar_genes), ").")
    } else if (identical(scope, "filtered")) {
      "Gene scope: current left-sidebar filters."
    } else if (identical(scope, "selected")) {
      "Gene scope: selected card only."
    } else if (length(sidebar_genes)) {
      gene_scope_class <- "warning"
      paste0(
        "Gene scope: all queued rows; sidebar Genes (",
        short_gene_list(sidebar_genes), ") are ignored."
      )
    } else {
      "Gene scope: all queued rows."
    }
    protocol_rows <- if ("protocol_pooling_need" %in% names(dt)) {
      sum(suppressWarnings(as.numeric(dt$protocol_pooling_need)) >= 0.55,
          na.rm = TRUE)
    } else {
      0L
    }
    open_label_rows <- nrow(review_label_open_rows(dt))
    tagList(
      tags$b(nrow(dt)), " queue rows; ",
      tags$b(uniqueN(dt$card_key)), " cards. ",
      target_label,
      gate_label,
      tags$b(open_label_rows), " open label tasks. ",
      tags$b(protocol_rows), " protocol/pooling-audit rows. ",
      "Scope: ", tags$b(scope_label), ". ",
      tags$span(class = gene_scope_class, gene_scope), " ",
      "Selecting a row updates the selected card and inspection context."
    )
  })

  output$review_workbench_table <- renderDT({
    dt <- copy(review_workbench_filtered())
    if (!nrow(dt)) {
      dt <- data.table(message = "No queue rows match the current scope.")
      return(datatable(dt, rownames = FALSE, options = list(dom = "t")))
    }
    dt <- make_queue_observatory_links(dt)
    keep <- intersect(review_workbench_display_columns, names(dt))
    out <- dt[, ..keep]
    out <- round_numeric_columns(out, 4)
    datatable(
      out,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      escape = FALSE,
      class = "stripe hover compact",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(match("queue_rank", names(out)) - 1L, "asc"))
      )
    )
  })

  observeEvent(input$review_workbench_table_rows_selected, {
    dt <- review_workbench_filtered()
    idx <- input$review_workbench_table_rows_selected
    if (!nrow(dt) || length(idx) != 1 || idx > nrow(dt)) return()
    focus_review_context(dt[idx])
  }, ignoreInit = TRUE)

  selected_review_workbench_row <- reactive({
    dt <- review_workbench_filtered()
    if (!nrow(dt)) return(data.table())
    idx <- input$review_workbench_table_rows_selected
    if (length(idx) != 1L || idx > nrow(dt)) return(data.table())
    dt[idx]
  })

  output$review_workbench_metrics <- renderUI({
    dt <- review_workbench_filtered()
    destination_counts <- review_workbench_destination_counts(dt)
    destination_count <- function(destination) {
      destination_name <- destination
      value <- destination_counts[destination == destination_name, rows]
      if (length(value)) value[1] else 0L
    }
    protocol_rows <- sum(review_workbench_protocol_audit_mask(dt), na.rm = TRUE)
    high_hazard_rows <- if ("pooling_hazard_class" %in% names(dt)) {
      sum(dt$pooling_hazard_class %chin% c(
        "high_pooling_confounding_hazard",
        "moderate_pooling_confounding_hazard"
      ), na.rm = TRUE)
    } else {
      0L
    }
    exact_proof_rows <- if ("pooling_proof_tier" %in% names(dt)) {
      sum(dt$pooling_proof_tier ==
            "A_exact_strict_project_matched_control", na.rm = TRUE)
    } else {
      0L
    }
    metric <- function(value, label) {
      div(class = "metric", tags$b(value), span(label))
    }
    div(
      class = "workbench-metric-grid",
      metric(nrow(dt), "queue rows"),
      metric(uniqueN(dt$card_key), "cards"),
      metric(nrow(review_label_open_rows(dt)), "open label tasks"),
      metric(destination_count("review_labels"), "open to labels"),
      metric(destination_count("inspection_assets"), "open to inspection"),
      metric(destination_count("selected_card"), "open to card"),
      metric(protocol_rows, "protocol audits"),
      metric(high_hazard_rows, "pooling risk"),
      metric(exact_proof_rows, "exact proof")
    )
  })

  output$review_workbench_selected_summary <- renderUI({
    row <- selected_review_workbench_row()
    if (!nrow(row)) {
      return(tagList(
        h4("Selected Row"),
        div(
          class = "small-note",
          "Select a Workbench row to see where it will open, why it was queued, ",
          "and which review or inspection context it carries."
        )
      ))
    }
    row_value <- function(column, default = "not set") {
      if (!column %in% names(row)) return(default)
      value <- scalar_text(row[[column]][1])
      if (is_scalar_informative(value)) value else default
    }
    row_number <- function(column, digits = 3, default = "not set") {
      if (!column %in% names(row)) return(default)
      value <- suppressWarnings(as.numeric(row[[column]][1]))
      if (!is.finite(value) || is.na(value)) return(default)
      fmt_num(value, digits)
    }
    destination <- review_workbench_selected_destination(row)
    destination_label <- review_workbench_destination_label(destination)
    active_task <- if (workbench_finite_integer_flag(row, "active_review_order")) {
      as.character(suppressWarnings(as.integer(row$active_review_order[1])))
    } else {
      "none"
    }
    inspection_source <- if (workbench_finite_integer_flag(row, "inspection_source_row")) {
      as.character(suppressWarnings(as.integer(row$inspection_source_row[1])))
    } else {
      "none"
    }
    observatory_url <- row_value("observatory_url", default = "")
    observatory_link <- if (is_scalar_informative(observatory_url)) {
      tags$a(
        class = "primary-link",
        href = observatory_url,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open Observatory"
      )
    }
    tagList(
      div(
        class = "workbench-selected-title",
        h4(paste(row_value("gene_symbol"), row_value("design_family"))),
        span(class = "workbench-route-pill", destination_label)
      ),
      div(
        class = "workbench-selected-actions",
        actionButton("open_selected_workbench_summary",
                     "Open Row",
                     class = "btn btn-primary btn-sm"),
        actionButton("use_selected_workbench_row_gene",
                     "Gene Session",
                     class = "btn btn-default btn-sm"),
        observatory_link
      ),
      div(
        class = "workbench-detail-grid",
        div(class = "key", "Transcript"), div(row_value("tx_id")),
        div(class = "key", "Queue"), div(row_value("queue_source")),
        div(class = "key", "Queue rank"), div(row_number("queue_rank", 1)),
        div(class = "key", "Active task"), div(active_task),
        div(class = "key", "Inspection source"), div(inspection_source),
        div(class = "key", "Target label"), div(row_value("target_label_class")),
        div(class = "key", "Label gap"), div(row_value("label_gap_class")),
        div(class = "key", "Recommended action"), div(row_value("recommended_action")),
        div(class = "key", "Pooling proof"), div(row_value("pooling_proof_tier")),
        div(class = "key", "Pooling hazard"), div(row_value("pooling_hazard_class")),
        div(class = "key", "Protocol need"), div(row_number("protocol_pooling_need")),
        div(class = "key", "Review reason"), div(row_value("review_reason")),
        div(class = "key", "Link note"), div(row_value("observatory_link_note"))
      )
    )
  })

  selected_review_task_rows <- reactive({
    row <- selected_row()
    if (is.null(row) || !nrow(review_batch_template)) return(data.table())
    out <- review_batch_template[
      gene_symbol == row$gene_symbol[1] &
        tx_id == row$tx_id[1] &
        design_family == row$design_family[1]
    ]
    if (!nrow(out)) return(out)
    setorder(out, global_review_order)
    out
  })

  preferred_review_task_order <- reactive({
    rows <- selected_review_task_rows()
    if (!nrow(rows)) return(NA_integer_)
    active_order <- selected_review_task_order_hint()
    if (length(active_order) == 1L && is.finite(active_order) &&
        active_order %in% rows$global_review_order) {
      return(as.integer(active_order))
    }
    rows$global_review_order[1]
  })

  selected_review_task_order <- reactive({
    rows <- selected_review_task_rows()
    if (!nrow(rows)) return(NA_integer_)
    input_order <- suppressWarnings(as.integer(input$review_label_task))
    if (length(input_order) == 1L && is.finite(input_order) &&
        input_order %in% rows$global_review_order) {
      return(input_order)
    }
    preferred_review_task_order()
  })

  selected_review_task <- reactive({
    rows <- selected_review_task_rows()
    task_order <- selected_review_task_order()
    if (!nrow(rows) || !is.finite(task_order)) return(data.table())
    rows[global_review_order == task_order][1]
  })

  selected_review_task_context <- reactive({
    task <- selected_review_task()
    if (!nrow(task) || !"global_review_order" %in% names(task)) {
      return(data.table())
    }
    task_order <- suppressWarnings(as.integer(task$global_review_order[1]))
    out <- review_workbench[
      queue_source == "Review batch" &
        active_review_order == task_order
    ]
    if (!nrow(out)) return(data.table())
    out[1]
  })

  selected_review_draft <- reactive({
    task_order <- selected_review_task_order()
    drafts <- review_label_drafts_state()
    if (!nrow(drafts) || !is.finite(task_order)) return(data.table())
    drafts[global_review_order == task_order][1]
  })

  draft_value <- function(draft, task, field) {
    if (nrow(draft) && field %in% names(draft) &&
        is_scalar_informative(draft[[field]][1])) {
      return(scalar_text(draft[[field]][1]))
    }
    if (nrow(task) && field %in% names(task) &&
        is_scalar_informative(task[[field]][1])) {
      return(scalar_text(task[[field]][1]))
    }
    ""
  }

  output$review_label_task_select_ui <- renderUI({
    rows <- selected_review_task_rows()
    validate(need(nrow(rows) > 0,
                  "No active-review task exists for the selected card."))
    labels <- vapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i]
      paste0(
        "#", row$global_review_order,
        " | ", row$batch_id,
        " | ", row$batch_role,
        " | ", row$action_label
      )
    }, character(1))
    div(
      class = "inspection-source-row",
      selectInput(
        "review_label_task",
        "Active review task",
        choices = stats::setNames(as.character(rows$global_review_order),
                                  labels),
        selected = as.character(selected_review_task_order())
      ),
      div(
        class = "small-note",
        "Draft labels are saved separately and consumed by the feedback ingest; ",
        "the linked review template is not edited by the app."
      )
    )
  })

  output$review_label_status <- renderUI({
    task <- selected_review_task()
    validate(need(nrow(task) > 0,
                  "Select a card with an active-review task."))
    context <- selected_review_task_context()
    draft <- selected_review_draft()
    label <- if (nrow(draft)) {
      draft_value(draft, task, "review_label_class_preview")
    } else {
      ""
    }
    target <- if (nrow(context)) scalar_text(context$target_label_class) else ""
    gap <- if (nrow(context)) scalar_text(context$label_gap_class) else ""
    hazard <- if (nrow(context)) scalar_text(context$pooling_hazard_class) else ""
    pooling_need <- if (nrow(context) &&
                        "protocol_pooling_need" %in% names(context) &&
                        is.finite(suppressWarnings(
                          as.numeric(context$protocol_pooling_need[1])
                        ))) {
      fmt_num(as.numeric(context$protocol_pooling_need[1]), 2)
    } else {
      ""
    }
    div(
      class = "review-label-status",
      tags$b(paste0(task$gene_symbol[1], " / ", task$design_family[1])),
      " | task #", task$global_review_order[1],
      " | ", task$batch_id[1],
      " | ", task$batch_role[1],
      if (is_scalar_informative(label)) {
        tagList(" | draft label: ", tags$b(label))
      } else {
        " | no draft label saved"
      },
      tags$br(),
      if (is_scalar_informative(target) || is_scalar_informative(gap) ||
          is_scalar_informative(hazard)) {
        tagList(
          tags$span(class = "label-status-pill",
                    paste("target", target)),
          tags$span(class = "label-status-pill",
                    paste("gap", gap)),
          if (is_scalar_informative(pooling_need)) {
            tags$span(class = "label-status-pill",
                      paste("protocol/pooling", pooling_need))
          },
          if (is_scalar_informative(hazard)) {
            tags$span(class = "label-status-pill", hazard)
          },
          tags$br()
        )
      } else {
        NULL
      },
      tags$span(class = "small-note", task$review_question[1]),
      if (nrow(draft) && is_scalar_informative(draft$draft_saved_at[1])) {
        tagList(tags$br(), tags$span(
          class = "small-note",
          paste("Saved:", draft$draft_saved_at[1])
        ))
      } else {
        NULL
      }
    )
  })

  review_label_link_value <- function(field) {
    context <- selected_review_task_context()
    task <- selected_review_task()
    row <- selected_row()
    for (source in list(context, task, row)) {
      if (!is.null(source) && nrow(source) && field %in% names(source) &&
          is_scalar_informative(source[[field]][1])) {
        return(scalar_text(source[[field]][1]))
      }
    }
    ""
  }

  output$review_label_inspection_actions <- renderUI({
    task <- selected_review_task()
    validate(need(nrow(task) > 0,
                  "Select a card with an active-review task."))
    observatory_url <- review_label_link_value("observatory_url")
    observatory_note <- review_label_link_value("observatory_link_note")
    observatory_link <- if (is_scalar_informative(observatory_url)) {
      tags$a(
        class = "primary-link",
        href = observatory_url,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open in Observatory"
      )
    } else {
      tags$span(class = "missing-file", "No Observatory link")
    }
    div(
      class = "review-label-wide",
      div(
        class = "label-action-row",
        observatory_link,
        actionButton("show_inspection_from_label", "Inspection Assets",
                     class = "btn btn-default btn-sm"),
        actionButton("show_selected_card_from_label", "Selected Card",
                     class = "btn btn-default btn-sm"),
        actionButton("show_workbench_from_label", "Review Workbench",
                     class = "btn btn-default btn-sm"),
        actionButton("show_files_from_label", "Files",
                     class = "btn btn-default btn-sm")
      ),
      if (is_scalar_informative(observatory_note)) {
        div(class = "link-note", observatory_note)
      } else {
        NULL
      }
    )
  })

  output$review_label_form <- renderUI({
    task <- selected_review_task()
    validate(need(nrow(task) > 0,
                  "Select a card with an active-review task."))
    draft <- selected_review_draft()
    div(
      class = "review-label-grid",
      div(
        class = "review-label-wide label-action-row",
        actionButton("quick_label_positive", "Positive",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_weak_positive", "Weak Positive",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_negative", "Negative",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_confounded", "Confounded",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_measurement", "Measurement Problem",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_not_reviewable", "Not Reviewable",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_negative_control_pass",
                     "Negative-Control Pass",
                     class = "btn btn-default btn-sm"),
        actionButton("quick_label_negative_control_fail",
                     "Negative-Control Fail",
                     class = "btn btn-default btn-sm")
      ),
      selectInput(
        "draft_manual_review_status",
        "Manual review status",
        choices = review_manual_status_choices,
        selected = draft_value(draft, task, "manual_review_status")
      ),
      selectInput(
        "draft_browser_signal_grade",
        "Browser signal grade",
        choices = review_browser_grade_choices,
        selected = draft_value(draft, task, "browser_signal_grade")
      ),
      selectInput(
        "draft_negative_control_result",
        "Negative-control result",
        choices = review_negative_control_choices,
        selected = draft_value(draft, task, "negative_control_result")
      ),
      selectInput(
        "draft_branch_direction_verified",
        "Branch direction",
        choices = review_direction_choices,
        selected = draft_value(draft, task, "branch_direction_verified")
      ),
      selectInput(
        "draft_clean_cds_direction_verified",
        "Clean-CDS direction",
        choices = review_direction_choices,
        selected = draft_value(draft, task, "clean_cds_direction_verified")
      ),
      selectInput(
        "draft_context_match_verified",
        "Context match",
        choices = review_direction_choices,
        selected = draft_value(draft, task, "context_match_verified")
      ),
      selectInput(
        "draft_translon_structure_status",
        "Translon / isoform status",
        choices = review_translon_choices,
        selected = draft_value(draft, task, "translon_structure_status")
      ),
      selectInput(
        "draft_metadata_confounder_status",
        "Metadata confounder",
        choices = review_confounder_choices,
        selected = draft_value(draft, task, "metadata_confounder_status")
      ),
      selectInput(
        "draft_next_action",
        "Next action",
        choices = review_next_action_choices,
        selected = draft_value(draft, task, "next_action")
      ),
      textInput(
        "draft_reviewer",
        "Reviewer",
        value = {
          value <- draft_value(draft, task, "reviewer")
          if (is_scalar_informative(value)) value else Sys.info()[["user"]]
        }
      ),
      textInput(
        "draft_review_date",
        "Review date",
        value = {
          value <- draft_value(draft, task, "review_date")
          if (is_scalar_informative(value)) value else as.character(Sys.Date())
        }
      ),
      div(
        class = "small-note",
        uiOutput("review_label_preview")
      ),
      div(
        class = "review-label-wide",
        textAreaInput(
          "draft_reviewer_notes",
          "Reviewer notes",
          value = draft_value(draft, task, "reviewer_notes"),
          rows = 4,
          placeholder = paste(
            "Record what was seen in Inspection Assets / Observatory,",
            "including coverage breadth, replicate consistency, RDG-flow match,",
            "and why the row should or should not train the next model."
          )
        )
      ),
      div(
        class = "review-label-wide review-label-actions",
        actionButton("save_review_label_draft", "Save Draft Label",
                     class = "btn btn-primary"),
        actionButton("save_review_label_draft_next", "Save & Next",
                     class = "btn btn-default"),
        actionButton("delete_review_label_draft", "Delete Draft",
                     class = "btn btn-default"),
        uiOutput("review_label_draft_file_link")
      )
    )
  })

  output$review_label_preview <- renderUI({
    task <- selected_review_task()
    if (!nrow(task)) return("")
    label <- infer_review_label_class(
      input$draft_manual_review_status %||% "",
      input$draft_browser_signal_grade %||% "",
      input$draft_branch_direction_verified %||% "",
      input$draft_clean_cds_direction_verified %||% "",
      input$draft_context_match_verified %||% "",
      input$draft_translon_structure_status %||% "",
      input$draft_metadata_confounder_status %||% "",
      input$draft_negative_control_result %||% "",
      task$batch_role[1]
    )
    if (!is_scalar_informative(label)) label <- "no label yet"
    tagList("Preview label: ", tags$b(label))
  })

  output$review_label_draft_file_link <- renderUI({
    file_info <- make_selected_file_link(review_label_drafts_file)
    if (is.null(file_info) || !isTRUE(file_info$exists) ||
        !is_scalar_informative(file_info$href)) {
      return(tags$span(class = "small-note", "No draft file saved yet."))
    }
    tags$a(
      class = "file-link",
      href = file_info$href,
      target = "_blank",
      rel = "noopener noreferrer",
      "Open draft CSV"
    )
  })

  apply_quick_label <- function(manual = "", browser = "", negative = "",
                                branch = "", clean_cds = "", context = "",
                                translon = "", confounder = "",
                                next_action = "") {
    updateSelectInput(session, "draft_manual_review_status",
                      selected = manual)
    updateSelectInput(session, "draft_browser_signal_grade",
                      selected = browser)
    updateSelectInput(session, "draft_negative_control_result",
                      selected = negative)
    updateSelectInput(session, "draft_branch_direction_verified",
                      selected = branch)
    updateSelectInput(session, "draft_clean_cds_direction_verified",
                      selected = clean_cds)
    updateSelectInput(session, "draft_context_match_verified",
                      selected = context)
    updateSelectInput(session, "draft_translon_structure_status",
                      selected = translon)
    updateSelectInput(session, "draft_metadata_confounder_status",
                      selected = confounder)
    updateSelectInput(session, "draft_next_action",
                      selected = next_action)
  }

  observeEvent(input$quick_label_positive, {
    apply_quick_label(
      manual = "positive",
      browser = "strong_positive",
      branch = "verified",
      clean_cds = "verified",
      context = "verified",
      translon = "verified",
      confounder = "no_confounder_seen",
      negative = "not_applicable",
      next_action = "merge_to_browser_notes"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_weak_positive, {
    apply_quick_label(
      manual = "weak_positive",
      browser = "weak_positive",
      branch = "unclear",
      clean_cds = "verified",
      context = "verified",
      translon = "unclear",
      confounder = "no_confounder_seen",
      negative = "not_applicable",
      next_action = "needs_observatory_inspection"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_negative, {
    apply_quick_label(
      manual = "negative",
      browser = "negative",
      branch = "opposite",
      clean_cds = "opposite",
      context = "verified",
      translon = "verified",
      confounder = "no_confounder_seen",
      negative = "not_applicable",
      next_action = "reject"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_confounded, {
    apply_quick_label(
      manual = "confounded",
      browser = "confounded",
      branch = "unclear",
      clean_cds = "unclear",
      context = "opposite",
      translon = "unclear",
      confounder = "confounded",
      negative = "not_applicable",
      next_action = "needs_context_audit"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_measurement, {
    apply_quick_label(
      manual = "measurement_problem",
      browser = "measurement_problem",
      branch = "unclear",
      clean_cds = "unclear",
      context = "unclear",
      translon = "measurement_problem",
      confounder = "",
      negative = "not_applicable",
      next_action = "needs_translon_fix"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_not_reviewable, {
    apply_quick_label(
      manual = "not_reviewable",
      browser = "not_reviewable",
      branch = "not_applicable",
      clean_cds = "not_applicable",
      context = "not_applicable",
      translon = "not_reviewable",
      confounder = "",
      negative = "not_applicable",
      next_action = "revisit_later"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_negative_control_pass, {
    apply_quick_label(
      manual = "negative",
      browser = "negative",
      branch = "not_applicable",
      clean_cds = "not_applicable",
      context = "verified",
      translon = "verified",
      confounder = "no_confounder_seen",
      negative = "pass",
      next_action = "merge_to_browser_notes"
    )
  }, ignoreInit = TRUE)

  observeEvent(input$quick_label_negative_control_fail, {
    apply_quick_label(
      manual = "positive",
      browser = "strong_positive",
      branch = "verified",
      clean_cds = "verified",
      context = "verified",
      translon = "verified",
      confounder = "no_confounder_seen",
      negative = "fail",
      next_action = "needs_observatory_inspection"
    )
  }, ignoreInit = TRUE)

  save_current_review_label_draft <- function() {
    task <- selected_review_task()
    validate(need(nrow(task) > 0,
                  "No active-review task selected."))
    context <- selected_review_task_context()
    context_value <- function(field) {
      if (nrow(context) && field %in% names(context)) {
        return(scalar_text(context[[field]]))
      }
      if (field %in% names(task)) {
        return(scalar_text(task[[field]]))
      }
      ""
    }
    label <- infer_review_label_class(
      input$draft_manual_review_status %||% "",
      input$draft_browser_signal_grade %||% "",
      input$draft_branch_direction_verified %||% "",
      input$draft_clean_cds_direction_verified %||% "",
      input$draft_context_match_verified %||% "",
      input$draft_translon_structure_status %||% "",
      input$draft_metadata_confounder_status %||% "",
      input$draft_negative_control_result %||% "",
      task$batch_role[1]
    )
    draft <- data.table(
      global_review_order = as.integer(task$global_review_order[1]),
      batch_id = scalar_text(task$batch_id),
      batch_name = scalar_text(task$batch_name),
      batch_review_order = scalar_text(task$batch_review_order),
      batch_role = scalar_text(task$batch_role),
      gene_symbol = scalar_text(task$gene_symbol),
      tx_id = scalar_text(task$tx_id),
      design_family = scalar_text(task$design_family),
      review_packet_id = scalar_text(task$review_packet_id),
      action_label = scalar_text(task$action_label),
      claim_stage = scalar_text(task$claim_stage),
      review_question = scalar_text(task$review_question),
      review_reason = scalar_text(task$review_reason),
      observatory_link_note = scalar_text(task$observatory_link_note),
      target_label_class = context_value("target_label_class"),
      label_gap_class = context_value("label_gap_class"),
      label_acquisition_score = context_value("label_acquisition_score"),
      protocol_pooling_need = context_value("protocol_pooling_need"),
      pooling_proof_tier = context_value("pooling_proof_tier"),
      pooling_hazard_class = context_value("pooling_hazard_class"),
      pooling_hazard_score = context_value("pooling_hazard_score"),
      pooling_hazard_sources = context_value("pooling_hazard_sources"),
      manual_review_status = scalar_text(input$draft_manual_review_status %||% ""),
      browser_signal_grade = scalar_text(input$draft_browser_signal_grade %||% ""),
      branch_direction_verified =
        scalar_text(input$draft_branch_direction_verified %||% ""),
      clean_cds_direction_verified =
        scalar_text(input$draft_clean_cds_direction_verified %||% ""),
      context_match_verified =
        scalar_text(input$draft_context_match_verified %||% ""),
      translon_structure_status =
        scalar_text(input$draft_translon_structure_status %||% ""),
      metadata_confounder_status =
        scalar_text(input$draft_metadata_confounder_status %||% ""),
      negative_control_result =
        scalar_text(input$draft_negative_control_result %||% ""),
      reviewer_notes = scalar_text(input$draft_reviewer_notes %||% ""),
      reviewer = scalar_text(input$draft_reviewer %||% ""),
      review_date = scalar_text(input$draft_review_date %||% ""),
      next_action = scalar_text(input$draft_next_action %||% ""),
      review_label_class_preview = label,
      draft_saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    saved <- write_review_label_draft(draft)
    review_label_drafts_state(saved)
    list(
      drafts = saved,
      draft = draft,
      task_order = as.integer(task$global_review_order[1]),
      label = label
    )
  }

  observeEvent(input$save_review_label_draft, {
    save_current_review_label_draft()
    showNotification("Draft review label saved.", type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$save_review_label_draft_next, {
    candidate_rows <- review_workbench_filtered()
    saved <- save_current_review_label_draft()
    next_row <- choose_next_review_label_row(
      candidate_rows,
      current_order = saved$task_order
    )
    if (nrow(next_row)) {
      focus_review_context(next_row)
      jump_to_review_page("Review Labels")
      showNotification("Draft saved; advanced to next open label task.",
                       type = "message")
    } else {
      showNotification("Draft saved; no next open label task in this Workbench context.",
                       type = "message")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$delete_review_label_draft, {
    task_order <- selected_review_task_order()
    validate(need(is.finite(task_order), "No draft task selected."))
    review_label_drafts_state(delete_review_label_draft(task_order))
    showNotification("Draft review label deleted.", type = "warning")
  }, ignoreInit = TRUE)

  output$review_label_drafts_table <- renderDT({
    drafts <- copy(review_label_drafts_state())
    if (!nrow(drafts)) {
      drafts <- data.table(message = "No draft labels saved yet.")
      return(datatable(drafts, rownames = FALSE, options = list(dom = "t")))
    }
    row <- selected_row()
    if (!is.null(row)) {
      selected_drafts <- drafts[
        gene_symbol == row$gene_symbol[1] &
          tx_id == row$tx_id[1] &
          design_family == row$design_family[1]
      ]
      if (nrow(selected_drafts)) drafts <- selected_drafts
    }
    keep <- intersect(c(
      "global_review_order", "batch_id", "batch_role", "gene_symbol",
      "tx_id", "design_family", "target_label_class", "label_gap_class",
      "protocol_pooling_need", "pooling_hazard_class",
      "manual_review_status",
      "browser_signal_grade", "branch_direction_verified",
      "clean_cds_direction_verified", "context_match_verified",
      "translon_structure_status", "metadata_confounder_status",
      "negative_control_result", "review_label_class_preview",
      "reviewer", "review_date", "next_action", "reviewer_notes",
      "draft_saved_at"
    ), names(drafts))
    drafts <- drafts[, ..keep]
    datatable(
      drafts,
      rownames = FALSE,
      filter = "top",
      class = "stripe hover compact",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  observe_plot_selection <- function(source) {
    observeEvent(plotly::event_data("plotly_click", source = source,
                                    priority = "event"), {
      click <- plotly::event_data("plotly_click", source = source,
                                  priority = "event")
      key <- click[["key"]] %||% click[["customdata"]]
      key <- as.character(unlist(key, use.names = FALSE))
      key <- key[is_informative_text(key)]
      if (!length(key) && source %chin% c("effect_heatmap", "branch_heatmap")) {
        gene <- as.character(unlist(click[["y"]] %||% "", use.names = FALSE))[1]
        family <- as.character(unlist(click[["x"]] %||% "", use.names = FALSE))[1]
        dt <- filtered_cards()[gene_symbol == gene & design_family == family]
        if (identical(source, "branch_heatmap")) {
          dt <- dt[is.finite(joint_top_l1_shift)]
          if (nrow(dt)) {
            dt[, branch_heatmap_value := 100 * joint_top_l1_shift]
            setorder(dt, -branch_heatmap_value,
                     -atlas_browser_review_priority_score)
          }
        } else if (nrow(dt)) {
          setorder(dt, -atlas_browser_review_priority_score,
                   -atlas_confidence_score)
        }
        if (nrow(dt)) key <- dt$card_key[1]
      }
      if (length(key)) selected_card_key(key[1])
    }, ignoreInit = TRUE)
  }
  observe_plot_selection("effect_space")
  observe_plot_selection("effect_heatmap")
  observe_plot_selection("branch_heatmap")

  selected_row <- reactive({
    dt <- filtered_cards()
    if (!nrow(dt)) return(NULL)
    key <- selected_card_key()
    if (is_scalar_informative(key)) {
      hit <- which(dt$card_key == key)[1]
      if (is.finite(hit) && !is.na(hit)) return(dt[hit])
    }
    idx <- input$cards_table_rows_selected
    if (length(idx) != 1 || idx > nrow(dt)) idx <- 1L
    dt[idx]
  })

  visual_selected_key <- function() {
    key <- selected_card_key()
    if (is_scalar_informative(key)) return(key)
    row <- selected_row()
    if (is.null(row)) return("")
    row$card_key[1]
  }

  selected_inspection_rows <- reactive({
    inspection_rows_for_card(selected_row())
  })

  selected_inspection_source <- reactive({
    rows <- selected_inspection_rows()
    if (!nrow(rows)) return(NA_integer_)
    input_source <- suppressWarnings(as.integer(input$inspection_source_row))
    if (length(input_source) == 1L && is.finite(input_source) &&
        input_source %in% rows$source_row) {
      return(input_source)
    }
    preferred <- suppressWarnings(as.integer(selected_inspection_source_row()))
    if (length(preferred) == 1L && is.finite(preferred) &&
        preferred %in% rows$source_row) {
      return(preferred)
    }
    rows$source_row[1]
  })

  selected_inspection_manifest_row <- reactive({
    rows <- selected_inspection_rows()
    source_id <- selected_inspection_source()
    if (!nrow(rows) || !is.finite(source_id)) return(data.table())
    out <- rows[get("source_row") == source_id]
    if (!nrow(out)) return(data.table())
    out[1]
  })

  inspection_source_label <- function(row) {
    paste0(
      "#", row$source_row,
      " | ", row$study,
      " | ", row$tissue,
      " | ", row$cell_line,
      " | ", row$case_condition,
      " vs ", row$control_conditions
    )
  }

  output$inspection_source_select_ui <- renderUI({
    rows <- selected_inspection_rows()
    validate(need(nrow(rows) > 0,
                  "No qualified inspection assets for the selected card."))
    labels <- vapply(seq_len(nrow(rows)), function(i) {
      inspection_source_label(rows[i])
    }, character(1))
    choices <- stats::setNames(as.character(rows$source_row), labels)
    selectInput(
      "inspection_source_row",
      "Inspection source row",
      choices = choices,
      selected = as.character(selected_inspection_source())
    )
  })

  observeEvent(input$inspection_source_row, {
    source_row <- suppressWarnings(as.integer(input$inspection_source_row))
    if (length(source_row) == 1L && is.finite(source_row)) {
      selected_inspection_source_row(source_row)
    }
  }, ignoreInit = TRUE)

  output$inspection_source_note <- renderUI({
    row <- selected_inspection_manifest_row()
    if (!nrow(row)) return("No qualified inspection source selected.")
    tagList(
      tags$b(row$inspection_readiness[1]), "; ",
      "branch min count ", tags$b(row$branch_min_case_control_counts[1]), "; ",
      "display min count ", tags$b(row$display_min_exact_counts[1]), "; ",
      "top feature min count ", tags$b(row$top_feature_min_exact_counts[1]), "; ",
      "mRNA sequence ", tags$b(if (isTRUE(row$mrna_sequence_available[1])) "available" else "missing"),
      if (is_informative_text(row$flags[1])) tagList(tags$br(), row$flags[1]) else NULL
    )
  })

  inspection_source_table <- function(dt) {
    source_id <- selected_inspection_source()
    if (!nrow(dt) || !is.finite(source_id)) return(data.table())
    dt[get("source_row") == source_id]
  }

  output$inspection_assets <- renderUI({
    row <- selected_inspection_manifest_row()
    validate(need(nrow(row) > 0,
                  "No qualified inspection assets for the selected card."))
    index_link <- make_selected_file_link(
      "dominant_rdg_qualified_inspection/qualified_review_index.html"
    )
    index_node <- if (!is.null(index_link) && isTRUE(index_link$exists)) {
      tags$a(
        class = "primary-link",
        href = index_link$href,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open inspection index"
      )
    } else {
      NULL
    }
    asset_fields <- c(
      profile_png = "Coverage profile",
      feature_allocation_png = "Feature allocation",
      rdg_flow_png = "RDG flow annotation",
      replicate_diagnostics_png = "Replicate diagnostics"
    )
    asset_cards <- lapply(names(asset_fields), function(field) {
      file_info <- make_selected_file_link(row[[field]][1])
      if (is.null(file_info)) return(NULL)
      body <- if (isTRUE(file_info$exists) &&
                  is_scalar_informative(file_info$href)) {
        tags$a(
          href = file_info$href,
          target = "_blank",
          rel = "noopener noreferrer",
          tags$img(src = file_info$href, loading = "lazy",
                   alt = unname(asset_fields[field]))
        )
      } else {
        div(class = "missing-file", file_info$abs_path)
      }
      div(
        class = "inspection-asset",
        div(
          class = "inspection-asset-header",
          tags$b(unname(asset_fields[field])),
          div(class = "link-note", file_info$abs_path)
        ),
        body
      )
    })
    tagList(
      div(
        class = "panel",
        div(
          class = "selected-title",
          h4(paste(row$gene_symbol, row$design_family,
                   paste0("source row ", row$source_row), sep = " / ")),
          index_node
        ),
        div(class = "link-note", row$expected_browser_pattern[1]),
        div(
          class = "selected-summary-grid",
          div(class = "selected-summary-item",
              tags$b(row$inspection_readiness[1]),
              span("readiness")),
          div(class = "selected-summary-item",
              tags$b(row$top_branch_class[1]),
              span("top branch")),
          div(class = "selected-summary-item",
              tags$b(fmt_num(row$top_branch_delta[1], 3)),
              span("branch delta")),
          div(class = "selected-summary-item",
              tags$b(row$top_feature_id[1]),
              span("top feature")),
          div(class = "selected-summary-item",
              tags$b(row$context_count_gate[1]),
              span("context count gate")),
          div(class = "selected-summary-item",
              tags$b(row$replicate_gate[1]),
              span("replicate gate"))
        )
      ),
      div(class = "inspection-asset-grid", asset_cards)
    )
  })

  output$inspection_profile_metrics_table <- renderDT({
    out <- inspection_source_table(qualified_profile_metrics)
    keep <- intersect(c(
      "group_id", "group_label", "group_role", "group_scope", "run_count",
      "display_width", "display_total_counts", "max_position_count",
      "max_position_fraction", "top3_position_fraction", "nonzero_positions"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No profile metrics for the selected source row.")
    }
    out <- round_numeric_columns(out, 4)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$inspection_feature_contrasts_table <- renderDT({
    out <- inspection_source_table(qualified_feature_contrasts)
    if (nrow(out) && "contrast_id" %in% names(out)) {
      exact <- out[contrast_id == "exact_case_vs_exact_control"]
      if (nrow(exact)) out <- exact
    }
    if (nrow(out)) {
      out[, abs_allocation_delta := abs(allocation_prop_delta)]
      out[, abs_prop_delta := abs(prop_delta)]
      setorder(out, -abs_allocation_delta, -abs_prop_delta)
    }
    keep <- intersect(c(
      "contrast_id", "feature_id", "feature_type", "category",
      "numerator_raw", "denominator_raw", "raw_log2_ratio", "prop_delta",
      "allocation_prop_delta", "numerator_max_run_fraction",
      "denominator_max_run_fraction"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No feature contrasts for the selected source row.")
    }
    out <- round_numeric_columns(out, 4)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 12, scrollX = TRUE))
  })

  output$inspection_replicate_summary_table <- renderDT({
    out <- inspection_source_table(qualified_replicate_summary)
    keep <- intersect(c(
      "group_id", "group_label", "group_role", "exact_run_count",
      "exact_measured_run_count", "exact_missing_run_count",
      "exact_signal_run_count", "exact_signal_run_fraction",
      "exact_total_top_feature_counts", "exact_median_top_feature_counts",
      "exact_max_top_feature_counts", "exact_max_run_fraction"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No replicate summary for the selected source row.")
    }
    out <- round_numeric_columns(out, 4)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$inspection_sequence_features_table <- renderDT({
    out <- inspection_source_table(qualified_sequence_features)
    if (nrow(out)) {
      setorder(out, tx_start, tx_end, feature_id)
    }
    keep <- intersect(c(
      "feature_id", "feature_type", "category", "tx_start", "tx_end",
      "feature_bases", "cds_overlap_bases", "start_codon", "stop_codon",
      "start_codon_class", "frame_relative_to_cds", "peptide_length_aa",
      "kozak_strength", "coding_potential_proxy", "start_context_sequence",
      "mrna_sequence_available", "mrna_sequence_length",
      "display_sequence_available", "display_sequence_length"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No sequence features for the selected source row.")
    }
    out <- round_numeric_columns(out, 4)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 12, scrollX = TRUE))
  })

  sdrive_metric_value <- function(metric_name, default = "NA") {
    if (!nrow(sdrive_readiness_summary_metrics) ||
        !"metric" %in% names(sdrive_readiness_summary_metrics)) {
      return(default)
    }
    hit <- sdrive_readiness_summary_metrics[metric == metric_name]
    if (!nrow(hit) || !"value" %in% names(hit)) return(default)
    scalar_text(hit$value[1], default)
  }

  selected_sdrive_readiness_dt <- reactive({
    row <- selected_row()
    key_cols <- c("gene_symbol", "tx_id", "design_family")
    if (!nrow(sdrive_transcript_readiness) || is.null(row) ||
        !all(key_cols %in% names(sdrive_transcript_readiness))) {
      return(data.table())
    }
    out <- sdrive_transcript_readiness[
      gene_symbol == row$gene_symbol[1] &
        tx_id == row$tx_id[1] &
        design_family == row$design_family[1]
    ]
    if (!nrow(out)) return(out)
    source_id <- selected_inspection_source()
    if (is.finite(source_id) && "source_row" %in% names(out) &&
        source_id %in% out$source_row) {
      setorder(out, source_row)
      return(out[source_row == source_id])
    }
    if ("source_row" %in% names(out)) setorder(out, source_row)
    out
  })

  output$sdrive_readiness_metrics <- renderUI({
    metric_box <- function(metric, label) {
      div(class = "metric", tags$b(sdrive_metric_value(metric)), span(label))
    }
    summary_link <- make_selected_file_link(
      "dominant_rdg_sdrive_readiness/sdrive_readiness_summary.md"
    )
    summary_node <- if (!is.null(summary_link) &&
                        isTRUE(summary_link$exists)) {
      tags$a(
        class = "primary-link",
        href = summary_link$href,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open audit summary"
      )
    } else {
      span(class = "small-note", "Audit summary has not been generated yet.")
    }
    tagList(
      div(
        class = "selected-title",
        h4("S-Drive Backing Data"),
        summary_node
      ),
      div(
        class = "metrics-grid",
        metric_box("sdrive_mount_ready", "S drive mounted"),
        metric_box("experiment_loaded", "ORFik experiment"),
        metric_box("required_fst_pages_present", "FST pages present"),
        metric_box("required_fst_pages_missing", "FST pages missing"),
        metric_box("sdrive_ready_source_rows", "ready source rows"),
        metric_box("coverage_reviewable_source_rows", "coverage-reviewable rows"),
        metric_box("spike_risk_source_rows", "spike-risk rows"),
        metric_box("not_claimable_static_source_rows", "not-claimable rows")
      )
    )
  })

  output$selected_sdrive_readiness <- renderUI({
    out <- selected_sdrive_readiness_dt()
    validate(need(nrow(out) > 0,
                  "No S-drive readiness row for the selected card."))
    row <- out[1]
    value <- function(field, default = "NA") {
      if (!field %in% names(row)) return(default)
      x <- row[[field]][1]
      if (is.logical(x)) return(if (isTRUE(x)) "TRUE" else "FALSE")
      if (is.numeric(x)) return(fmt_num(x, 3))
      if (is.na(x) || !nzchar(as.character(x))) return(default)
      as.character(x)
    }
    tagList(
      div(
        class = "selected-summary-grid",
        div(class = "selected-summary-item",
            tags$b(value("sdrive_readiness_class")),
            span("data readiness")),
        div(class = "selected-summary-item",
            tags$b(value("review_support_class")),
            span("coverage support")),
        div(class = "selected-summary-item",
            tags$b(paste0(value("required_fst_pages_present"), "/",
                          value("required_fst_page_count"))),
            span("FST pages")),
        div(class = "selected-summary-item",
            tags$b(value("all_merged_core_ready")),
            span("all-merged files")),
        div(class = "selected-summary-item",
            tags$b(value("reference_ready")),
            span("reference")),
        div(class = "selected-summary-item",
            tags$b(value("mrna_sequence_available")),
            span("mRNA sequence")),
        div(class = "selected-summary-item",
            tags$b(value("top_feature_min_exact_counts")),
            span("top feature count")),
        div(class = "selected-summary-item",
            tags$b(value("top_feature_max_run_fraction")),
            span("top feature run max"))
      ),
      if (is_informative_text(value("flags", ""))) {
        div(class = "small-note", tags$b("Flags: "), value("flags", ""))
      } else {
        NULL
      }
    )
  })

  output$selected_sdrive_readiness_table <- renderDT({
    out <- selected_sdrive_readiness_dt()
    keep <- intersect(c(
      "source_row", "gene_symbol", "tx_id", "design_family", "study",
      "tissue", "cell_line", "case_condition", "control_conditions",
      "inspection_readiness", "review_support_class",
      "sdrive_readiness_class", "static_inspection_status",
      "fst_pages_ready", "required_fst_page_count",
      "required_fst_pages_present", "required_fst_pages_missing",
      "mrna_sequence_available", "display_sequence_available",
      "reference_ready", "all_merged_core_ready", "rdg_feature_count",
      "leader_uorf_features", "overlapping_uorf_features",
      "internal_orf_features", "nte_extension_features",
      "branch_min_case_control_counts", "display_min_exact_counts",
      "top_feature_min_exact_counts", "max_exact_position_fraction",
      "top_feature_max_run_fraction", "top_feature_in_rdg_annotation",
      "has_internal_or_nte_architecture", "flags",
      "missing_required_fst_page_basenames"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No S-drive readiness row for the selected card.")
    }
    out <- round_numeric_columns(out, 4)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 5, scrollX = TRUE))
  })

  output$sdrive_status_table <- renderDT({
    out <- copy(sdrive_experiment_status)
    keep <- intersect(c(
      "check", "status", "ok", "value", "path", "detail",
      "size_bytes", "mtime", "generated_at"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No S-drive experiment audit has been generated.")
    }
    out <- round_numeric_columns(out, 2)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })

  output$sdrive_transcript_readiness_table <- renderDT({
    out <- copy(sdrive_transcript_readiness)
    keep <- intersect(c(
      "source_row", "gene_symbol", "tx_id", "design_family",
      "inspection_readiness", "review_support_class",
      "sdrive_readiness_class", "study", "tissue", "cell_line",
      "case_condition", "control_conditions", "required_fst_page_count",
      "required_fst_pages_missing", "mrna_sequence_available",
      "display_sequence_available", "branch_min_case_control_counts",
      "display_min_exact_counts", "top_feature_min_exact_counts",
      "top_feature_max_run_fraction", "flags"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No transcript readiness audit has been generated.")
    }
    out <- round_numeric_columns(out, 4)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })

  output$sdrive_fst_page_readiness_table <- renderDT({
    out <- copy(sdrive_fst_page_readiness)
    keep <- intersect(c(
      "fst_basename", "direction", "chr", "part", "n_genes", "genes",
      "n_transcripts", "readiness_status", "local_file_exists_now",
      "file_size_mb", "file_mtime", "local_path"
    ), names(out))
    out <- if (nrow(out)) out[, ..keep] else {
      data.table(message = "No required FST page audit has been generated.")
    }
    out <- round_numeric_columns(out, 2)
    datatable(out, rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })

  add_heatmap_selection_trace <- function(plot, key_mat) {
    row <- selected_row()
    if (is.null(row) || !length(key_mat)) return(plot)
    x_value <- row$design_family[1]
    y_value <- row$gene_symbol[1]
    if (!(x_value %chin% colnames(key_mat)) || !(y_value %chin% rownames(key_mat))) {
      return(plot)
    }
    add_trace(
      plot,
      x = x_value,
      y = y_value,
      type = "scatter",
      mode = "markers",
      inherit = FALSE,
      showlegend = FALSE,
      hoverinfo = "skip",
      marker = list(
        symbol = "square-open",
        size = 19,
        color = "#111827",
        line = list(color = "#111827", width = 3)
      )
    )
  }

  output$visual_selected_card <- renderUI({
    row <- selected_row()
    validate(need(!is.null(row), "No card selected."))
    row_text <- function(field, digits = 3) {
      if (!field %in% names(row)) return("NA")
      value <- row[[field]][1]
      if (is.numeric(value)) value <- fmt_num(value, digits)
      if (is.na(value) || !nzchar(as.character(value))) return("NA")
      as.character(value)
    }
    row_pct <- function(field) {
      if (!field %in% names(row)) return("NA")
      value <- row[[field]][1]
      if (!is.numeric(value) || is.na(value)) return(row_text(field))
      fmt_pct(value)
    }
    link <- if ("observatory_url" %in% names(row) &&
                is_scalar_informative(row$observatory_url)) {
      tags$a(
        class = "primary-link",
        href = row$observatory_url[1],
        target = "_blank",
        rel = "noopener noreferrer",
        "Open in Observatory"
      )
    } else {
      NULL
    }
    div(
      class = "panel visual-selection-panel",
      div(
        class = "selected-title",
        h4(paste(row$gene_symbol, row$design_family, sep = " / ")),
        div(
          class = "visual-selection-actions",
          actionButton("show_selected_card_from_visual", "Selected Card",
                       class = "btn btn-default btn-sm"),
          link
        )
	      ),
	      div(class = "link-note", row$observatory_link_note[1]),
	      div(
	        class = "visual-selection-actions visual-selection-jumps",
	        actionButton("show_cards_from_visual", "Cards",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_workbench_from_visual", "Review Workbench",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_files_from_visual", "Files",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_inspection_from_visual", "Inspection Assets",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_branch_contexts_from_visual", "Branch Contexts",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_study_evidence_from_visual", "Study Evidence",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_transport_evidence_from_visual",
	                     "Transport Evidence",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_context_interactions_from_visual",
	                     "Context Interactions",
	                     class = "btn btn-default btn-sm"),
	        actionButton("show_partial_pooling_from_visual", "Partial Pooling",
	                     class = "btn btn-default btn-sm")
	      ),
	      div(
	        class = "selected-summary-grid",
        div(class = "selected-summary-item",
            tags$b(row_text("atlas_confidence_score", 2)),
            span("confidence")),
        div(class = "selected-summary-item",
            tags$b(row_pct("atlas_clean_cds_pct_change")),
            span("clean-CDS effect")),
        div(class = "selected-summary-item",
            tags$b(row_text("model_agreement_tier")),
            span("model agreement")),
        div(class = "selected-summary-item",
            tags$b(row_text("transport_consensus_class")),
            span("grouped omission")),
        div(class = "selected-summary-item",
            tags$b(row_text("joint_top_branch_class")),
            span("top branch")),
        div(class = "selected-summary-item",
            tags$b(row_text("context_interaction_gene_class")),
            span("context interaction"))
      )
    )
  })

  jump_to_review_page <- function(page) {
    updateTabsetPanel(session, "page_group", selected = "Review")
    updateTabsetPanel(session, "review_pages", selected = page)
  }
  jump_to_evidence_page <- function(page) {
    updateTabsetPanel(session, "page_group", selected = "Evidence")
    updateTabsetPanel(session, "evidence_pages", selected = page)
  }

  start_review_label_session <- function(rows,
                                         empty_message = "No open label task in this context.") {
    next_row <- choose_next_review_label_row(rows)
    if (!nrow(next_row)) {
      showNotification(empty_message, type = "warning")
      return(invisible(FALSE))
    }
    focus_review_context(next_row)
    jump_to_review_page("Review Labels")
    showNotification(
      paste0(
        "Started label session at task #",
        next_row$active_review_order[1], "."
      ),
      type = "message"
    )
    invisible(TRUE)
  }

  open_selected_workbench_row <- function(row = selected_review_workbench_row()) {
    if (!nrow(row)) {
      showNotification("Select a Workbench row first.", type = "warning")
      return(invisible(FALSE))
    }
    destination <- review_workbench_selected_destination(row)
    focus_review_context(row)
    active_order <- if ("active_review_order" %in% names(row)) {
      suppressWarnings(as.integer(row$active_review_order[1]))
    } else {
      NA_integer_
    }
    if (identical(destination, "review_labels")) {
      jump_to_review_page("Review Labels")
      showNotification(
        paste0("Opened selected review task #", active_order, "."),
        type = "message"
      )
    } else if (identical(destination, "inspection_assets")) {
      jump_to_evidence_page("Inspection Assets")
      showNotification(
        "Selected row has no active label task; opened Inspection Assets.",
        type = "message"
      )
    } else {
      jump_to_review_page("Selected Card")
      showNotification(
        "Selected row has no active label task or inspection packet; opened Selected Card.",
        type = "message"
      )
    }
    invisible(TRUE)
  }

  set_workbench_gene_session_from_row <- function(row = selected_review_workbench_row()) {
    if (!nrow(row) || !is_scalar_informative(row$gene_symbol[1])) {
      showNotification("Select a Workbench row with a gene first.",
                       type = "warning")
      return(invisible(FALSE))
    }
    gene <- row$gene_symbol[1]
    updateSelectizeInput(session, "review_workbench_genes",
                         selected = gene)
    updateSelectInput(session, "review_workbench_scope",
                      selected = "all")
    showNotification(
      paste0("Workbench gene set to ", gene, "; scope set to all queued rows."),
      type = "message"
    )
    invisible(TRUE)
  }

  observeEvent(input$show_cards_from_visual, {
    jump_to_review_page("Cards")
  })
  observeEvent(input$show_workbench_from_visual, {
    jump_to_review_page("Review Workbench")
  })
  observeEvent(input$show_selected_card_from_label, {
    jump_to_review_page("Selected Card")
  }, ignoreInit = TRUE)
  observeEvent(input$show_workbench_from_label, {
    jump_to_review_page("Review Workbench")
  }, ignoreInit = TRUE)
  observeEvent(input$show_files_from_label, {
    jump_to_review_page("Files")
  }, ignoreInit = TRUE)
  observeEvent(input$show_inspection_from_label, {
    jump_to_evidence_page("Inspection Assets")
  }, ignoreInit = TRUE)
  observeEvent(input$use_sidebar_genes_review_workbench, {
    genes <- as.character(input$genes %||% character())
    genes <- unique(genes[is_informative_text(genes)])
    available <- genes[genes %chin% review_workbench_gene_choices]
    if (!length(available)) {
      showNotification(
        "No sidebar Genes are available in the Workbench queues.",
        type = "warning"
      )
      return(invisible(FALSE))
    }
    updateSelectizeInput(session, "review_workbench_genes",
                         selected = available)
    updateSelectInput(session, "review_workbench_scope",
                      selected = "all")
    showNotification(
      "Workbench genes set from the sidebar; scope set to all queued rows.",
      type = "message"
    )
  }, ignoreInit = TRUE)
  observeEvent(input$use_selected_gene_review_workbench, {
    row <- selected_row()
    if (is.null(row) || !nrow(row) ||
        !is_scalar_informative(row$gene_symbol[1])) {
      showNotification("No selected card gene is available.", type = "warning")
      return(invisible(FALSE))
    }
    gene <- row$gene_symbol[1]
    if (!gene %chin% review_workbench_gene_choices) {
      showNotification(
        paste0("Selected gene ", gene, " has no Workbench queue rows."),
        type = "warning"
      )
      return(invisible(FALSE))
    }
    updateSelectizeInput(session, "review_workbench_genes",
                         selected = gene)
    updateSelectInput(session, "review_workbench_scope",
                      selected = "all")
    showNotification(
      paste0("Workbench gene set to ", gene, "; scope set to all queued rows."),
      type = "message"
    )
  }, ignoreInit = TRUE)
  observeEvent(input$clear_review_workbench_genes, {
    updateSelectizeInput(session, "review_workbench_genes",
                         selected = character())
    showNotification("Workbench gene filter cleared.", type = "message")
  }, ignoreInit = TRUE)
  observeEvent(input$start_label_dashboard_session, {
    target <- input$label_dashboard_target %||% "all"
    updateSelectInput(session, "review_workbench_source",
                      selected = "Label gap")
    updateSelectInput(session, "review_workbench_target",
                      selected = target)
    updateSelectInput(session, "review_workbench_scope",
                      selected = "all")
    updateCheckboxInput(session, "review_workbench_open_only",
                        value = TRUE)
    start_review_label_session(
      label_dashboard_gap_rows(),
      empty_message = "No open label-gap task matches the Label Dashboard filter."
    )
  }, ignoreInit = TRUE)
  observeEvent(input$start_review_workbench_session, {
    start_review_label_session(
      review_workbench_filtered(),
      empty_message = paste(
        "No open label task matches the current Workbench filters.",
        "Select a row and use Open Selected for inspection-only rows."
      )
    )
  }, ignoreInit = TRUE)
  observeEvent(input$start_selected_review_workbench_session, {
    open_selected_workbench_row()
  }, ignoreInit = TRUE)
  observeEvent(input$open_selected_workbench_summary, {
    open_selected_workbench_row()
  }, ignoreInit = TRUE)
  observeEvent(input$use_selected_workbench_row_gene, {
    set_workbench_gene_session_from_row()
  }, ignoreInit = TRUE)
  observeEvent(input$open_label_gap_review_labels, {
    focus_review_context(selected_label_gap_row())
    jump_to_review_page("Review Labels")
  })
  observeEvent(input$open_label_gap_inspection, {
    focus_review_context(selected_label_gap_row())
    jump_to_evidence_page("Inspection Assets")
  })
  observeEvent(input$open_label_gap_workbench, {
    focus_review_context(selected_label_gap_row())
    updateSelectInput(session, "review_workbench_source",
                      selected = "Label gap")
    updateSelectInput(session, "review_workbench_scope",
                      selected = "selected")
    updateCheckboxInput(session, "review_workbench_open_only",
                        value = FALSE)
    jump_to_review_page("Review Workbench")
  })
  observeEvent(input$show_selected_card_from_visual, {
    jump_to_review_page("Selected Card")
  })
  observeEvent(input$show_files_from_visual, {
    jump_to_review_page("Files")
  })
  observeEvent(input$show_inspection_from_visual, {
    jump_to_evidence_page("Inspection Assets")
  })
  observeEvent(input$show_branch_contexts_from_visual, {
    jump_to_evidence_page("Branch Contexts")
  })
  observeEvent(input$show_study_evidence_from_visual, {
    jump_to_evidence_page("Study Evidence")
  })
  observeEvent(input$show_transport_evidence_from_visual, {
    jump_to_evidence_page("Transport Evidence")
  })
  observeEvent(input$show_context_interactions_from_visual, {
    jump_to_evidence_page("Context Interactions")
  })
  observeEvent(input$show_partial_pooling_from_visual, {
    jump_to_evidence_page("Partial Pooling")
  })
  observeEvent(input$show_workbench_from_selected, {
    jump_to_review_page("Review Workbench")
  })
  observeEvent(input$show_inspection_from_selected, {
    jump_to_evidence_page("Inspection Assets")
  })

  output$heatmap <- renderPlotly({
    dt <- filtered_cards()
    validate(need(nrow(dt) > 0, "No cards match the current filters."))
    dt <- dt[order(-atlas_browser_review_priority_score)][seq_len(min(.N, 160))]
    gene_order <- dt[, .(score = max(atlas_confidence_score, na.rm = TRUE)),
                     by = gene_symbol][order(score), gene_symbol]
    family_order <- dt[, .(score = max(atlas_confidence_score, na.rm = TRUE)),
                       by = design_family][order(-score), design_family]
    cell_dt <- dt[order(-atlas_browser_review_priority_score,
                        -atlas_confidence_score), .(
      heatmap_value = mean(atlas_clean_cds_pct_change, na.rm = TRUE),
      heatmap_key = card_key[1],
      heatmap_hover = paste0(
        "Gene: ", gene_symbol[1],
        "<br>Design: ", design_family[1],
        "<br>Mean clean-CDS effect: ",
        fmt_pct(mean(atlas_clean_cds_pct_change, na.rm = TRUE)),
        "<br>Representative transcript: ", tx_id[1],
        "<br>Representative confidence: ",
        fmt_num(atlas_confidence_score[1]),
        "<br>Representative evidence: ", atlas_evidence_class[1],
        "<br>Cards in cell: ", .N,
        "<br>Top feature: ", top_shifted_feature_id[1],
        "<br>Warnings:<br>", fmt_warning_hover(warnings[1])
      )
    ), by = .(gene_symbol, design_family)]
    mat <- dcast(cell_dt, gene_symbol ~ design_family,
                 value.var = "heatmap_value")
    keys <- dcast(cell_dt, gene_symbol ~ design_family,
                  value.var = "heatmap_key")
    hover <- dcast(cell_dt, gene_symbol ~ design_family,
                   value.var = "heatmap_hover")
    z <- as.matrix(mat[, ..family_order])
    rownames(z) <- mat$gene_symbol
    key_mat <- as.matrix(keys[, ..family_order])
    rownames(key_mat) <- keys$gene_symbol
    txt <- as.matrix(hover[, ..family_order])
    rownames(txt) <- hover$gene_symbol
    keep_rows <- intersect(gene_order, rownames(z))
    z <- z[keep_rows, , drop = FALSE]
    key_mat <- key_mat[keep_rows, , drop = FALSE]
    txt <- txt[keep_rows, , drop = FALSE]
    z <- pmax(pmin(z, 300), -100)
    p <- plot_ly(
      source = "effect_heatmap",
      x = colnames(z),
      y = rownames(z),
      z = z,
      customdata = key_mat,
      text = txt,
      type = "heatmap",
      zmin = -100,
      zmax = 300,
      colorscale = list(
        list(0, "#2b6cb0"),
        list(0.25, "white"),
        list(1, "#8b0000")
      ),
      hovertemplate = paste("%{text}", "<extra></extra>")
    )
    p <- add_heatmap_selection_trace(p, key_mat)
    p %>%
      layout(
        xaxis = list(title = "", tickangle = 35),
        yaxis = list(title = ""),
        margin = list(l = 90, b = 120)
      )
  })

  output$branch_heatmap <- renderPlotly({
    dt <- filtered_cards()
    validate(need(nrow(dt) > 0, "No cards match the current filters."))
    validate(need("joint_top_l1_shift" %in% names(dt),
                  "No joint branch-allocation columns found in atlas cards."))
    dt <- dt[is.finite(joint_top_l1_shift)]
    validate(need(nrow(dt) > 0,
                  "No joint branch-allocation rows match the current filters."))
    dt[, branch_heatmap_value := 100 * joint_top_l1_shift]
    dt[, branch_hover := paste0(
      "Gene: ", gene_symbol,
      "<br>Design: ", design_family,
      "<br>Joint branch shift: ", round(branch_heatmap_value, 1), " pp",
      "<br>Top branch: ", joint_top_branch_class,
      "<br>Top branch delta: ", fmt_pp(joint_top_branch_delta),
      "<br>Joint review: ", joint_top_review_class,
      "<br>Joint context: ", joint_top_context,
      "<br>Clean-CDS effect: ", fmt_pct(atlas_clean_cds_pct_change),
      "<br>Warnings:<br>", fmt_warning_hover(warnings)
    )]
    dt <- dt[order(-atlas_browser_review_priority_score,
                   -branch_heatmap_value)][seq_len(min(.N, 160))]
    cell_dt <- dt[order(-branch_heatmap_value,
                        -atlas_browser_review_priority_score), .(
      branch_heatmap_value = safe_max(branch_heatmap_value),
      branch_heatmap_key = card_key[1],
      branch_hover = paste0(branch_hover[1], "<br>Cards in cell: ", .N)
    ), by = .(gene_symbol, design_family)]
    gene_order <- dt[, .(score = max(branch_heatmap_value, na.rm = TRUE)),
                     by = gene_symbol][order(score), gene_symbol]
    family_order <- dt[, .(score = max(branch_heatmap_value, na.rm = TRUE)),
                       by = design_family][order(-score), design_family]
    mat <- dcast(cell_dt, gene_symbol ~ design_family,
                 value.var = "branch_heatmap_value")
    keys <- dcast(cell_dt, gene_symbol ~ design_family,
                  value.var = "branch_heatmap_key")
    hover <- dcast(cell_dt, gene_symbol ~ design_family,
                   value.var = "branch_hover")
    z <- as.matrix(mat[, ..family_order])
    rownames(z) <- mat$gene_symbol
    key_mat <- as.matrix(keys[, ..family_order])
    rownames(key_mat) <- keys$gene_symbol
    txt <- as.matrix(hover[, ..family_order])
    rownames(txt) <- hover$gene_symbol
    keep_rows <- intersect(gene_order, rownames(z))
    z <- z[keep_rows, , drop = FALSE]
    key_mat <- key_mat[keep_rows, , drop = FALSE]
    txt <- txt[keep_rows, , drop = FALSE]
    z <- pmin(z, 100)
    p <- plot_ly(
      source = "branch_heatmap",
      x = colnames(z),
      y = rownames(z),
      z = z,
      customdata = key_mat,
      text = txt,
      type = "heatmap",
      zmin = 0,
      zmax = 100,
      colorscale = list(
        list(0, "white"),
        list(0.5, "#f4a582"),
        list(1, "#7f0000")
      ),
      hovertemplate = paste("%{text}", "<extra></extra>")
    )
    p <- add_heatmap_selection_trace(p, key_mat)
    p %>%
      layout(
        xaxis = list(title = "", tickangle = 35),
        yaxis = list(title = ""),
        margin = list(l = 90, b = 120)
      )
  })

  output$effect_space_note <- renderUI({
    dt <- filtered_cards()
    y_col <- safe_choice_value(input$effect_y, effect_y_choices,
                               "atlas_confidence_score")
    color_col <- safe_choice_value(input$effect_color, effect_color_choices,
                                   "atlas_evidence_class")
    size_col <- safe_choice_value(input$effect_size, effect_size_choices,
                                  "none")
    finite_n <- if (nrow(dt)) {
      sum(is.finite(dt$atlas_clean_cds_pct_change) & is.finite(dt[[y_col]]))
    } else {
      0L
    }
    tagList(
      tags$b(finite_n), " / ", nrow(dt), " finite x/y cards; ",
      "y = ", choice_label(effect_y_choices, y_col), "; ",
      "color = ", choice_label(effect_color_choices, color_col), "; ",
      "size = ", choice_label(effect_size_choices, size_col)
    )
  })

  output$scatter <- renderPlotly({
    dt <- filtered_cards()
    validate(need(nrow(dt) > 0, "No cards match the current filters."))
    y_col <- safe_choice_value(input$effect_y, effect_y_choices,
                               "atlas_confidence_score")
    color_col <- safe_choice_value(input$effect_color, effect_color_choices,
                                   "atlas_evidence_class")
    size_col <- safe_choice_value(input$effect_size, effect_size_choices,
                                  "none")
    y_label <- choice_label(effect_y_choices, y_col)
    color_label <- choice_label(effect_color_choices, color_col)
    size_label <- choice_label(effect_size_choices, size_col)

    dt[, plot_x := atlas_clean_cds_pct_change]
    dt[, plot_y := as.numeric(get(y_col))]
    dt[, plot_color := clean_text(get(color_col))]
    dt[!is_informative_text(plot_color), plot_color := "not annotated"]
    if (identical(size_col, "none")) {
      dt[, plot_size := 8]
    } else {
      size_value <- as.numeric(dt[[size_col]])
      if (grepl("delta", size_col, fixed = TRUE)) size_value <- abs(size_value)
      dt[, plot_size := scale_marker_size(size_value)]
    }
    dt <- dt[is.finite(plot_x) & is.finite(plot_y)]
    validate(need(nrow(dt) > 0,
                  "No cards have finite x/y values for the selected view."))
    dt[, plot_hover := paste0(
      gene_symbol, " / ", design_family,
      "<br>CDS effect: ", fmt_pct(atlas_clean_cds_pct_change),
      "<br>", y_label, ": ", fmt_num(plot_y, 3),
      "<br>", color_label, ": ", plot_color,
      "<br>Point size: ", size_label,
      "<br>Confidence: ", fmt_num(atlas_confidence_score),
      "<br>Review priority: ", fmt_num(atlas_browser_review_priority_score),
      "<br>Browser validation: ", browser_validation_status,
      "<br>Branch usage: ", branch_usage_top_feature_id,
      " / ", branch_usage_top_support_class,
      "<br>Grouped branch: ", grouped_branch_usage_top_feature_id,
      " / ", grouped_branch_usage_top_support_class,
      " / ", grouped_branch_usage_top_context,
      "<br>Count-aware branch: ", count_aware_top_feature_id,
      " / ", count_aware_top_support_class,
      " / ", count_aware_top_context,
      "<br>Joint branch: ", joint_top_branch_class,
      " / ", joint_top_review_class,
      " / delta ", fmt_num(joint_top_branch_delta, 3),
      "<br>Joint context: ", joint_top_context,
      "<br>Hier. joint: ", hierarchical_joint_top_branch_class,
      " / ", hierarchical_joint_top_support_class,
      " / delta ", fmt_num(hierarchical_joint_top_delta, 3),
      "<br>Hier. joint context: ", hierarchical_joint_top_context,
      "<br>DM branch: ", dm_top_branch_class,
      " / ", dm_top_support_class,
      " / delta ", fmt_num(dm_top_delta, 3),
      "<br>DM context: ", dm_top_context,
      "<br>Model agreement: ", model_agreement_tier,
      " / ", fmt_num(model_agreement_score, 3),
      "<br>Model support: ", model_support_signature,
      "<br>Model disagreement: ", model_disagreement_flags,
      "<br>LOO stability: ", loo_consensus_class,
      " / ", fmt_num(loo_stability_score, 3),
      "<br>LOO minimum direction/support retention: ",
      fmt_num(loo_min_direction_fraction, 3), " / ",
      fmt_num(loo_min_support_retention, 3),
      "<br>Worst omitted studies: ", loo_worst_omitted_studies,
      "<br>Grouped transport: ", transport_consensus_class,
      " / ", fmt_num(transport_stability_score, 3),
      "<br>Transport minimum direction/support retention: ",
      fmt_num(transport_min_direction_fraction, 3), " / ",
      fmt_num(transport_min_support_retention, 3),
      "<br>Transport branch/axis coverage: ",
      fmt_num(transport_coverage_fraction, 3), " / ",
      fmt_num(transport_axis_coverage_fraction, 3),
      "<br>Worst omitted groups: ", transport_worst_omitted_groups,
      "<br>Context interaction: ", context_interaction_gene_class,
      " / ", context_interaction_top_axis, "=",
      context_interaction_top_value,
      " / ", context_interaction_top_branch_class,
      " / delta ", fmt_num(context_interaction_delta, 3),
      "<br>Partial pooling: ", partial_pooling_gene_class,
      " / ", partial_pooling_top_axis, "=",
      partial_pooling_top_value,
      " / ", partial_pooling_top_branch_class,
      " / pooled delta ", fmt_num(partial_pooling_top_delta, 3),
      " / q ", fmt_num(partial_pooling_top_q, 3),
      "<br>Hierarchical: ", hierarchical_evidence_class,
      " / ", fmt_num(hierarchical_priority_score),
      "<br>Context model: ", context_model_status,
      "<br>Top feature: ", top_shifted_feature_id,
      "<br>Warnings:<br>", fmt_warning_hover(warnings)
    )]
    selected_key <- selected_card_key()
    if (!is_scalar_informative(selected_key)) {
      row <- selected_row()
      if (!is.null(row)) selected_key <- row$card_key[1]
    }
    selected_dt <- if (is_scalar_informative(selected_key)) {
      dt[card_key == selected_key]
    } else {
      dt[0]
    }
    if (nrow(selected_dt)) {
      selected_dt[, plot_selected_size := pmin(plot_size + 8, 30)]
    }
    y_axis <- list(title = y_label, autorange = TRUE, zeroline = TRUE,
                   automargin = TRUE)
    nonnegative_y <- y_col %chin% c(
      "atlas_confidence_score", "atlas_browser_review_priority_score",
      "strict_supported_studies", "strict_supported_designs", "atlas_i2",
      "sparse_onoff_score", "joint_top_l1_shift", "model_agreement_score",
      "loo_stability_score", "transport_stability_score",
      "partial_pooling_top_q", "hierarchical_priority_score",
      "residual_heterogeneity_score"
    )
    if (isTRUE(nonnegative_y)) y_axis$rangemode <- "tozero"

    p <- plot_ly(
      dt,
      source = "effect_space",
      x = ~plot_x,
      y = ~plot_y,
      color = ~plot_color,
      colors = effect_palette,
      key = ~card_key,
      customdata = ~card_key,
      type = "scatter",
      mode = "markers",
      text = ~plot_hover,
      hoverinfo = "text",
      marker = list(
        size = dt$plot_size,
        opacity = 0.72,
        line = list(color = "rgba(31,41,55,0.25)", width = 0.5)
      )
    )
    if (nrow(selected_dt)) {
      p <- add_trace(
        p,
        data = selected_dt,
        x = ~plot_x,
        y = ~plot_y,
        key = ~card_key,
        customdata = ~card_key,
        type = "scatter",
        mode = "markers",
        text = ~plot_hover,
        hoverinfo = "text",
        inherit = FALSE,
        showlegend = FALSE,
        marker = list(
          size = selected_dt$plot_selected_size,
          symbol = "circle-open",
          color = "#111827",
          line = list(color = "#111827", width = 3)
        )
      )
    }
    p %>%
      layout(
        xaxis = list(title = "Clean-CDS buffering effect (%)",
                     zeroline = TRUE, automargin = TRUE),
        yaxis = y_axis,
        legend = list(orientation = "h", y = -0.25),
        margin = list(l = 80, r = 20, b = 120, t = 15)
      )
  })

  branch_contexts_dt <- reactive({
    row <- selected_row()
    if (!nrow(branch_contexts) || is.null(row)) return(data.table())
    out <- branch_contexts[
      gene_symbol == row$gene_symbol[1] &
        design_family == row$design_family[1]
    ]
    if ("tx_id" %in% names(out) && "tx_id" %in% names(row)) {
      out <- out[tx_id == row$tx_id[1]]
    }
    if (!nrow(out)) return(out)
    keep <- c(
      "study", "AUTHOR", "case_condition", "TISSUE", "CELL_LINE", "GENE",
      "branch_context_rank", "branch_context_label",
      "Cancer_type", "control_conditions", "n_case_runs_merged",
      "n_control_runs_merged", "case_branch_total_counts",
      "control_branch_total_counts", "joint_review_class",
      "joint_review_priority", "joint_top_branch_class",
      "joint_top_branch_delta", "joint_top_branch_q",
      "joint_l1_allocation_shift", "joint_delta_clean_CDS",
      "joint_delta_leader_uORF", "joint_delta_overlapping_uORF",
      "joint_case_prop_clean_CDS", "joint_control_prop_clean_CDS",
      "joint_case_prop_leader_uORF", "joint_control_prop_leader_uORF",
      "joint_case_prop_overlapping_uORF",
      "joint_control_prop_overlapping_uORF", "joint_evaluable",
      "joint_strict", "match_scope"
    )
    keep <- intersect(keep, names(out))
    out <- out[, ..keep]
    num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    for (col in num_cols) set(out, j = col, value = round(out[[col]], 4))
    if ("joint_review_priority" %in% names(out)) {
      setorder(out, -joint_review_priority)
    }
    out
  })

  output$branch_contexts_table <- renderDT({
    out <- branch_contexts_dt()
    if (!nrow(out)) {
      out <- data.table(message = "No joint branch contexts for the selected card.")
    }
    datatable(
      out,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })

  study_evidence_dt <- reactive({
    row <- selected_row()
    if (!nrow(study_evidence) || is.null(row)) return(data.table())
    out <- study_evidence[
      gene_symbol == row$gene_symbol[1] &
        design_family == row$design_family[1]
    ]
    if ("tx_id" %in% names(out) && "tx_id" %in% names(row)) {
      out <- out[tx_id == row$tx_id[1]]
    }
    if (!nrow(out)) return(out)
    keep <- c(
      "study", "study_evaluable_branches", "study_both_model_branches",
      "study_same_direction_branches", "study_opposite_direction_branches",
      "study_top_branch_class", "study_top_evidence_class",
      "study_top_mean_delta", "study_top_hierarchical_delta",
      "study_top_dm_delta", "study_top_priority", "study_top_context",
      "study_branch_summary"
    )
    keep <- intersect(keep, names(out))
    out <- out[, ..keep]
    num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    for (col in num_cols) set(out, j = col, value = round(out[[col]], 4))
    if ("study_top_priority" %in% names(out)) {
      setorder(out, -study_top_priority, study)
    }
    out
  })

  output$study_evidence_table <- renderDT({
    out <- study_evidence_dt()
    if (!nrow(out)) {
      out <- data.table(
        message = "No study-level model evidence for the selected card."
      )
    }
    datatable(
      out,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  transport_evidence_dt <- reactive({
    row <- selected_row()
    if (!nrow(transport_evidence) || is.null(row)) return(data.table())
    out <- transport_evidence[
      gene_symbol == row$gene_symbol[1] &
        design_family == row$design_family[1]
    ]
    if ("tx_id" %in% names(out) && "tx_id" %in% names(row)) {
      out <- out[tx_id == row$tx_id[1]]
    }
    if (!nrow(out)) return(out)
    keep <- c(
      "transport_axis", "transport_axis_class",
      "transport_axis_evaluable_branches",
      "transport_axis_support_robust_branches",
      "transport_axis_direction_fragile_branches",
      "transport_axis_min_direction_fraction",
      "transport_axis_min_support_retention",
      "transport_axis_min_remaining_studies",
      "transport_axis_worst_omitted_groups"
    )
    keep <- intersect(keep, names(out))
    out <- out[, ..keep]
    num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    for (col in num_cols) set(out, j = col, value = round(out[[col]], 4))
    setorder(out, transport_axis)
    out
  })

  output$transport_evidence_table <- renderDT({
    out <- transport_evidence_dt()
    if (!nrow(out)) {
      out <- data.table(
        message = "No grouped-omission transport evidence for the selected card."
      )
    }
    datatable(
      out,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  context_interaction_dt <- reactive({
    row <- selected_row()
    if (!nrow(context_interactions) || is.null(row)) return(data.table())
    out <- context_interactions[
      gene_symbol == row$gene_symbol[1] &
        design_family == row$design_family[1]
    ]
    if ("tx_id" %in% names(out) && "tx_id" %in% names(row)) {
      out <- out[tx_id == row$tx_id[1]]
    }
    if (!nrow(out)) return(out)
    keep <- c(
      "context_axis", "context_value", "branch_class",
      "context_interaction_class", "context_review_priority",
      "context_interaction_score", "context_studies",
      "complement_studies", "context_delta", "context_delta_lower",
      "context_delta_upper", "complement_delta",
      "complement_delta_lower", "complement_delta_upper",
      "interaction_delta", "interaction_delta_q",
      "context_support_class", "complement_support_class",
      "context_study_names", "complement_study_names",
      "context_top_labels", "complement_top_labels"
    )
    keep <- intersect(keep, names(out))
    out <- out[, ..keep]
    num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    for (col in num_cols) set(out, j = col, value = round(out[[col]], 4))
    if ("context_review_priority" %in% names(out)) {
      setorder(out, -context_review_priority,
               context_axis, context_value, branch_class)
    }
    out
  })

  output$context_interaction_table <- renderDT({
    out <- context_interaction_dt()
    if (!nrow(out)) {
      out <- data.table(
        message = "No context-interaction evidence for the selected card."
      )
    }
    datatable(
      out,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  partial_pooling_dt <- reactive({
    row <- selected_row()
    if (!nrow(partial_pooling) || is.null(row)) return(data.table())
    out <- partial_pooling[
      gene_symbol == row$gene_symbol[1] &
        design_family == row$design_family[1]
    ]
    if ("tx_id" %in% names(out) && "tx_id" %in% names(row)) {
      out <- out[tx_id == row$tx_id[1]]
    }
    if (!nrow(out)) return(out)
    link_rows <- rbindlist(lapply(seq_len(nrow(out)), function(i) {
      make_partial_pooling_observatory_link_row(
        out[i], rdg_metadata, gene_label_map
      )
    }), fill = TRUE)
    out[, names(link_rows) := link_rows]
    link_labels <- as.character(htmltools::htmlEscape(paste(
      out$gene_symbol, out$branch_class, out$context_axis,
      out$context_value, sep = " / "
    )))
    link_urls <- as.character(htmltools::htmlEscape(
      out$partial_pooling_observatory_url, attribute = TRUE
    ))
    link_notes <- as.character(htmltools::htmlEscape(
      out$partial_pooling_observatory_link_note, attribute = TRUE
    ))
    out[, Links := sprintf(
      "<a class=\"observatory-link\" href=\"%s\" target=\"_blank\" rel=\"noopener noreferrer\" title=\"%s\">%s</a>",
      link_urls, link_notes, link_labels
    )]
    keep <- c(
      "Links",
      "partial_pooling_context_case_n", "partial_pooling_context_control_n",
      "partial_pooling_complement_case_n",
      "partial_pooling_complement_control_n",
      "context_axis", "context_value", "branch_class",
      "partial_pooled_context_class", "partial_pooling_survival_class",
      "partial_pooling_review_priority", "interaction_delta",
      "interaction_delta_se", "interaction_delta_q",
      "pooled_interaction_delta", "pooled_interaction_lower",
      "pooled_interaction_upper", "pooled_interaction_q",
      "pooling_group_contexts", "partial_pooling_weight",
      "shrinkage_fraction", "pooled_same_direction_as_raw",
      "raw_context_supported", "context_studies", "complement_studies",
      "context_interaction_class"
    )
    keep <- intersect(keep, names(out))
    out <- out[, ..keep]
    num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    for (col in num_cols) set(out, j = col, value = round(out[[col]], 4))
    if ("partial_pooling_review_priority" %in% names(out)) {
      setorder(out, -partial_pooling_review_priority,
               context_axis, context_value, branch_class)
    }
    out
  })

  output$partial_pooling_table <- renderDT({
    out <- partial_pooling_dt()
    if (!nrow(out)) {
      out <- data.table(
        message = "No partial-pooling context evidence for the selected card."
      )
    }
    datatable(
      out,
      rownames = FALSE,
      filter = "top",
      escape = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  output$selected_card_view_note <- renderUI({
    switch(
      input$selected_card_view %||% "summary",
      model = "Model agreement, study omission, grouped omission, hierarchy, and residual fields.",
      context = "Branch allocation, context-interaction, partial-pooling, and context-model fields.",
      full = "Complete selected-card field dump.",
      "Compact first-pass review fields for the selected candidate."
    )
  })

  output$selected_card <- renderUI({
    row <- selected_row()
    validate(need(!is.null(row), "No card selected."))
    full_fields <- c(
      "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
      "observatory_gene_label", "observatory_context_label",
      "observatory_case_n", "observatory_control_n",
      "observatory_link_note",
      "atlas_confidence_score", "atlas_browser_review_priority_score",
      "atlas_clean_cds_pct_change",
      "atlas_clean_cds_pct_lower", "atlas_clean_cds_pct_upper",
      "atlas_direction_fraction", "atlas_i2", "strict_supported_designs",
      "strict_supported_studies", "top_shifted_feature_id",
      "top_shifted_feature_type", "feature_id_uORF_ouORF",
      "weighted_relative_use_delta_uORF_ouORF", "feature_id_iORF",
      "weighted_relative_use_delta_iORF", "sparse_feature_id",
      "sparse_onoff_pattern", "sparse_onoff_score",
      "browser_validation_status", "browser_validation_label",
      "browser_validation_notes", "best_case_condition",
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_proof_component", "pooling_hazard_class",
      "pooling_hazard_score", "pooling_hazard_sources",
      "pooling_context_rows", "pooling_exact_context_rows",
      "pooling_relaxed_context_rows", "pooling_joint_strict_rows",
      "pooling_joint_evaluable_rows", "pooling_exact_context_fraction",
      "pooling_joint_strict_fraction", "pooling_joint_evaluable_fraction",
      "pooling_context_wt_like_control_rows",
      "pooling_context_strong_signal_rows",
      "pooling_top_protocol_classes", "pooling_top_inhibitors",
      "pooling_top_tissues", "pooling_top_cell_lines",
      "pooling_top_case_conditions", "pooling_top_controls",
      "pooling_any_start_site_enriched_protocol",
      "pooling_any_relaxed_context",
      "branch_usage_top_feature_id", "branch_usage_top_feature_class",
      "branch_usage_top_support_class", "branch_usage_top_score",
      "branch_usage_top_weighted_usage_delta", "branch_usage_top_log_or",
      "branch_usage_top_q", "branch_usage_top_i2",
                "branch_usage_top_direction_fraction",
                "branch_usage_top_context",
                "grouped_branch_usage_top_feature_id",
                "grouped_branch_usage_top_feature_class",
                "grouped_branch_usage_top_support_class",
                "grouped_branch_usage_top_weighted_usage_delta",
                "grouped_branch_usage_top_context",
                "count_aware_top_feature_id",
                "count_aware_top_feature_class",
                "count_aware_top_support_class",
                "count_aware_top_weighted_usage_delta",
                "count_aware_top_context",
                "joint_top_review_class", "joint_top_branch_class",
                "joint_top_branch_delta", "joint_top_context",
                "hierarchical_joint_top_support_class",
                "hierarchical_joint_top_priority",
                "hierarchical_joint_top_branch_class",
                "hierarchical_joint_top_delta",
                "hierarchical_joint_top_delta_lower",
                "hierarchical_joint_top_delta_upper",
                "hierarchical_joint_top_delta_q",
                "hierarchical_joint_l1_effect",
                "hierarchical_joint_supported_branch_count",
                "hierarchical_joint_top_studies",
                "hierarchical_joint_top_context_rows",
                "hierarchical_joint_top_direction_fraction",
                "hierarchical_joint_top_i2",
                "hierarchical_joint_top_context",
                "dm_top_support_class", "dm_top_priority",
                "dm_top_branch_class", "dm_top_delta",
                "dm_top_delta_lower", "dm_top_delta_upper",
                "dm_top_delta_q", "dm_l1_effect",
                "dm_supported_branch_count", "dm_top_studies",
                "dm_top_context_rows", "dm_top_direction_fraction",
                "dm_top_i2", "dm_top_prior_source", "dm_top_context",
                "model_agreement_tier", "model_agreement_class",
                "model_agreement_score", "model_agreement_review_priority",
                "model_support_count", "replicated_model_count",
                "model_support_signature", "model_disagreement_flag",
                "model_disagreement_flags", "branch_model_consensus",
                "branch_model_direction_disagreement",
                "different_branch_support",
                "allocation_output_disagreement",
                "high_branch_heterogeneity",
                "common_supported_branches",
                "agreeing_common_branches",
                "disagreeing_common_branches",
                "common_supported_branch_names",
                "agreeing_common_branch_names",
                "disagreeing_common_branch_names",
                "clean_cds_output_allocation_alignment",
                "top_model_disagreement_branch",
                "top_model_disagreement_abs_difference",
                "model_agreement_top_context",
                "loo_consensus_class", "loo_agreement_class",
                "loo_stability_score", "loo_review_priority",
                "loo_relevant_supported_branches",
                "loo_support_robust_branches",
                "loo_direction_robust_support_fragile_branches",
                "loo_direction_fragile_branches",
                "loo_not_evaluable_branches",
                "loo_min_direction_fraction",
                "loo_min_support_retention", "loo_min_model_studies",
                "loo_support_robust_branch_names",
                "loo_fragile_branch_names", "loo_worst_omitted_studies",
                "transport_consensus_class", "transport_agreement_class",
                "transport_stability_score", "transport_review_priority",
                "transport_relevant_supported_branches",
                "transport_evaluable_supported_branches",
                "transport_support_robust_branches",
                "transport_direction_robust_support_fragile_branches",
                "transport_direction_fragile_branches",
                "transport_not_evaluable_branches",
                "transport_min_direction_fraction",
                "transport_min_support_retention",
                "transport_min_remaining_studies",
                "transport_coverage_fraction",
                "transport_axis_coverage_fraction",
                "transport_evaluable_axis_names",
                "transport_support_robust_axis_names",
                "transport_fragile_axis_names",
                "transport_worst_omitted_groups",
                "transport_worst_omitted_studies",
                "context_interaction_gene_class",
                "context_interaction_top_class",
                "context_interaction_top_score",
                "context_interaction_review_priority",
                "context_interaction_top_axis",
                "context_interaction_top_value",
                "context_interaction_top_branch_class",
                "context_interaction_context_delta",
                "context_interaction_context_delta_lower",
                "context_interaction_context_delta_upper",
                "context_interaction_complement_delta",
                "context_interaction_complement_delta_lower",
                "context_interaction_complement_delta_upper",
                "context_interaction_delta",
                "context_interaction_delta_q",
	                "context_interaction_context_studies",
	                "context_interaction_complement_studies",
	                "context_interaction_context_study_names",
	                "context_interaction_complement_study_names",
	                "context_interaction_summary",
	                "partial_pooling_gene_class",
	                "partial_pooling_top_class",
	                "partial_pooling_top_survival_class",
	                "partial_pooling_review_priority",
	                "partial_pooling_component",
	                "partial_pooling_top_axis",
	                "partial_pooling_top_value",
	                "partial_pooling_top_branch_class",
	                "partial_pooling_top_delta",
	                "partial_pooling_top_lower",
	                "partial_pooling_top_upper",
	                "partial_pooling_top_q",
	                "partial_pooling_top_raw_delta",
	                "partial_pooling_top_weight",
	                "partial_pooling_top_shrinkage_fraction",
	                "partial_pooling_rows",
	                "partial_pooling_pooled_supported_rows",
	                "partial_pooling_pooled_review_rows",
	                "partial_pooling_raw_supported_rows",
	                "partial_pooling_raw_survived_rows",
	                "partial_pooling_raw_shrunk_rows",
	                "partial_pooling_single_context_raw_signal_rows",
	                "partial_pooling_pooled_emergent_rows",
	                "partial_pooling_top_context_summary",
	                "partial_pooling_context_values",
	                "hierarchical_evidence_class", "hierarchical_priority_score",
                "hierarchical_clean_cds_pct_change",
                "hierarchical_clean_cds_pct_lower",
                "hierarchical_clean_cds_pct_upper",
                "hierarchical_calibrated_clean_cds_pct_lower",
                "hierarchical_calibrated_clean_cds_pct_upper",
                "hierarchical_calibrated_ci_excludes_zero",
                "hierarchical_conformal_q95_log2fc",
                "hierarchical_calibration_scope",
                "hierarchical_prior_source", "hierarchical_evidence_weight",
                "hierarchical_loo_mae_log2fc",
                "hierarchical_loo_interval_coverage",
                "hierarchical_loo_calibrated_interval_coverage",
                "residual_risk_class", "residual_heterogeneity_score",
                "residual_predictions", "residual_mae_log2fc",
                "residual_rmse_log2fc", "residual_fraction_over_design_q90",
                "residual_direction_accuracy",
                "residual_worst_abs_error_log2fc",
                "residual_worst_context_key",
                "residual_top_studies", "residual_top_cell_lines",
                "residual_top_tissues", "residual_top_conditions",
                "context_calibration_status",
                "context_calibration_predictions",
                "context_calibration_peer_corrected_fraction",
                "context_calibration_original_rmse_log2fc",
                "context_calibration_context_rmse_log2fc",
                "context_calibration_rmse_delta_log2fc",
                "context_calibration_original_q95_log2fc",
                "context_calibration_context_q95_log2fc",
                "context_calibration_q95_delta_log2fc",
                "context_calibration_fraction_improved",
                "context_calibration_top_offset_level",
                "context_calibration_top_worse_context",
                "context_calibration_top_improved_contexts",
                "context_model_status",
                "context_model_predictions",
                "context_model_original_rmse_log2fc",
                "context_model_rmse_log2fc",
                "context_model_rmse_delta_log2fc",
                "context_model_original_q95_log2fc",
                "context_model_q95_log2fc",
                "context_model_q95_delta_log2fc",
                "context_model_fraction_improved",
                "context_model_mean_abs_offset_log2fc",
                "context_model_worst_context",
                "context_model_top_improved_contexts",
                "best_control_conditions", "best_study", "browser_review_reason",
      "browser_review_question", "warnings", "nearest_designs"
    )
    summary_fields <- c(
      "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
      "browser_validation_status", "browser_validation_label",
      "atlas_confidence_score", "atlas_browser_review_priority_score",
      "atlas_clean_cds_pct_change",
      "atlas_clean_cds_pct_lower", "atlas_clean_cds_pct_upper",
      "strict_supported_designs", "strict_supported_studies",
      "best_study", "best_case_condition", "best_control_conditions",
      "top_shifted_feature_id", "top_shifted_feature_type",
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_hazard_class", "pooling_hazard_score",
      "pooling_hazard_sources", "pooling_top_protocol_classes",
      "branch_usage_top_feature_id", "branch_usage_top_feature_class",
      "branch_usage_top_support_class",
      "branch_usage_top_weighted_usage_delta",
      "joint_top_review_class", "joint_top_branch_class",
      "joint_top_branch_delta", "joint_top_context",
      "model_agreement_tier", "model_agreement_class",
      "model_agreement_score", "model_support_signature",
      "loo_consensus_class", "transport_consensus_class",
      "context_interaction_gene_class", "context_interaction_top_class",
      "partial_pooling_gene_class", "partial_pooling_top_class",
      "hierarchical_evidence_class", "residual_risk_class",
      "context_calibration_status", "context_model_status",
      "browser_review_reason", "browser_review_question",
      "warnings", "nearest_designs"
    )
    model_fields <- c(
      "gene_symbol", "tx_id", "design_family",
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_proof_component", "pooling_hazard_class",
      "pooling_hazard_score", "pooling_hazard_sources",
      "pooling_exact_context_fraction",
      "pooling_joint_strict_fraction",
      "pooling_top_protocol_classes",
      "model_agreement_tier", "model_agreement_class",
      "model_agreement_score", "model_agreement_review_priority",
      "model_support_count", "replicated_model_count",
      "model_support_signature", "model_disagreement_flag",
      "model_disagreement_flags", "branch_model_consensus",
      "branch_model_direction_disagreement",
      "different_branch_support", "allocation_output_disagreement",
      "high_branch_heterogeneity", "common_supported_branch_names",
      "agreeing_common_branch_names", "disagreeing_common_branch_names",
      "clean_cds_output_allocation_alignment",
      "top_model_disagreement_branch",
      "top_model_disagreement_abs_difference",
      "model_agreement_top_context",
      "loo_consensus_class", "loo_agreement_class",
      "loo_stability_score", "loo_review_priority",
      "loo_min_direction_fraction", "loo_min_support_retention",
      "loo_min_model_studies", "loo_support_robust_branch_names",
      "loo_fragile_branch_names", "loo_worst_omitted_studies",
      "transport_consensus_class", "transport_agreement_class",
      "transport_stability_score", "transport_review_priority",
      "transport_min_direction_fraction",
      "transport_min_support_retention",
      "transport_min_remaining_studies",
      "transport_coverage_fraction",
      "transport_axis_coverage_fraction",
      "transport_evaluable_axis_names",
      "transport_support_robust_axis_names",
      "transport_fragile_axis_names",
      "transport_worst_omitted_groups",
      "transport_worst_omitted_studies",
      "hierarchical_evidence_class", "hierarchical_priority_score",
      "hierarchical_clean_cds_pct_change",
      "hierarchical_calibrated_clean_cds_pct_lower",
      "hierarchical_calibrated_clean_cds_pct_upper",
      "hierarchical_prior_source", "hierarchical_evidence_weight",
      "residual_risk_class", "residual_heterogeneity_score",
      "residual_rmse_log2fc", "residual_direction_accuracy",
      "residual_worst_abs_error_log2fc",
      "residual_worst_context_key", "residual_top_studies",
      "residual_top_cell_lines", "residual_top_tissues",
      "residual_top_conditions", "warnings"
    )
    context_fields <- c(
      "gene_symbol", "tx_id", "design_family",
      "observatory_gene_label", "observatory_context_label",
      "observatory_case_n", "observatory_control_n",
      "observatory_link_note", "best_study", "best_case_condition",
      "best_control_conditions",
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_proof_component", "pooling_hazard_class",
      "pooling_hazard_score", "pooling_hazard_sources",
      "pooling_exact_context_fraction",
      "pooling_joint_strict_fraction",
      "pooling_top_protocol_classes",
      "branch_usage_top_feature_id", "branch_usage_top_feature_class",
      "branch_usage_top_support_class", "branch_usage_top_score",
      "branch_usage_top_weighted_usage_delta",
      "branch_usage_top_direction_fraction",
      "branch_usage_top_context",
      "grouped_branch_usage_top_feature_id",
      "grouped_branch_usage_top_feature_class",
      "grouped_branch_usage_top_support_class",
      "grouped_branch_usage_top_weighted_usage_delta",
      "grouped_branch_usage_top_context",
      "count_aware_top_feature_id",
      "count_aware_top_feature_class",
      "count_aware_top_support_class",
      "count_aware_top_weighted_usage_delta",
      "count_aware_top_context",
      "joint_top_review_class", "joint_top_branch_class",
      "joint_top_branch_delta", "joint_top_context",
      "hierarchical_joint_top_support_class",
      "hierarchical_joint_top_priority",
      "hierarchical_joint_top_branch_class",
      "hierarchical_joint_top_delta",
      "hierarchical_joint_top_delta_lower",
      "hierarchical_joint_top_delta_upper",
      "hierarchical_joint_top_delta_q",
      "hierarchical_joint_top_context",
      "dm_top_support_class", "dm_top_priority",
      "dm_top_branch_class", "dm_top_delta",
      "dm_top_delta_lower", "dm_top_delta_upper", "dm_top_delta_q",
      "dm_top_prior_source", "dm_top_context",
      "context_interaction_gene_class",
      "context_interaction_top_class",
      "context_interaction_top_score",
      "context_interaction_review_priority",
      "context_interaction_top_axis",
      "context_interaction_top_value",
      "context_interaction_top_branch_class",
      "context_interaction_context_delta",
      "context_interaction_complement_delta",
      "context_interaction_delta",
      "context_interaction_delta_q",
      "context_interaction_context_study_names",
      "context_interaction_complement_study_names",
      "context_interaction_summary",
      "partial_pooling_gene_class",
      "partial_pooling_top_class",
      "partial_pooling_top_survival_class",
      "partial_pooling_review_priority",
      "partial_pooling_top_axis",
      "partial_pooling_top_value",
      "partial_pooling_top_branch_class",
      "partial_pooling_top_delta",
      "partial_pooling_top_lower",
      "partial_pooling_top_upper",
      "partial_pooling_top_q",
      "partial_pooling_top_raw_delta",
      "partial_pooling_top_weight",
      "partial_pooling_top_shrinkage_fraction",
      "partial_pooling_top_context_summary",
      "partial_pooling_context_values",
      "context_calibration_status",
      "context_calibration_rmse_delta_log2fc",
      "context_calibration_fraction_improved",
      "context_calibration_top_worse_context",
      "context_calibration_top_improved_contexts",
      "context_model_status", "context_model_rmse_delta_log2fc",
      "context_model_fraction_improved",
      "context_model_worst_context",
      "context_model_top_improved_contexts", "warnings"
    )
    selected_view <- input$selected_card_view %||% "summary"
    fields <- switch(
      selected_view,
      model = model_fields,
      context = context_fields,
      full = full_fields,
      summary_fields
    )
    fields <- intersect(fields, names(row))
    link <- if ("observatory_url" %in% names(row) &&
                is_scalar_informative(row$observatory_url)) {
      tags$a(
        class = "primary-link",
        href = row$observatory_url[1],
        target = "_blank",
        rel = "noopener noreferrer",
        "Open in Observatory"
      )
    } else {
      NULL
    }
    row_text <- function(field, digits = 3) {
      if (!field %in% names(row)) return("NA")
      value <- row[[field]][1]
      if (is.numeric(value)) value <- fmt_num(value, digits)
      if (is.na(value) || !nzchar(as.character(value))) return("NA")
      as.character(value)
    }
    row_pct <- function(field) {
      if (!field %in% names(row)) return("NA")
      value <- row[[field]][1]
      if (!is.numeric(value) || is.na(value)) return(row_text(field))
      fmt_pct(value)
    }
    selected_summary <- div(
      class = "selected-summary-grid",
      div(class = "selected-summary-item",
          tags$b(row_text("atlas_confidence_score", 2)),
          span("confidence")),
      div(class = "selected-summary-item",
          tags$b(row_pct("atlas_clean_cds_pct_change")),
          span("clean-CDS effect")),
      div(class = "selected-summary-item",
          tags$b(row_text("model_agreement_tier")),
          span("model agreement")),
      div(class = "selected-summary-item",
          tags$b(row_text("loo_consensus_class")),
          span("study omission")),
      div(class = "selected-summary-item",
          tags$b(row_text("transport_consensus_class")),
          span("grouped omission")),
      div(class = "selected-summary-item",
          tags$b(row_text("joint_top_branch_class")),
          span("top branch")),
      div(class = "selected-summary-item",
          tags$b(row_text("context_interaction_gene_class")),
          span("context interaction")),
      div(class = "selected-summary-item",
          tags$b(row_text("partial_pooling_gene_class")),
          span("partial pooling")),
      div(class = "selected-summary-item",
          tags$b(row_text("pooling_proof_tier")),
          span("pooling proof")),
      div(class = "selected-summary-item",
          tags$b(row_text("pooling_hazard_class")),
          span("pooling hazard"))
    )
    div(
      class = "panel",
      div(
        class = "selected-title",
        h4(paste(row$gene_symbol, row$design_family, sep = " / ")),
        div(
          class = "visual-selection-actions",
          actionButton("show_workbench_from_selected", "Review Workbench",
                       class = "btn btn-default btn-sm"),
          actionButton("show_inspection_from_selected", "Inspection Assets",
                       class = "btn btn-default btn-sm"),
          link
        )
      ),
      div(class = "link-note", row$observatory_link_note[1]),
      br(),
      selected_summary,
      div(class = "detail-grid",
          lapply(fields, function(field) {
            value <- row[[field]][1]
            if (is.numeric(value)) value <- fmt_num(value, 3)
            if (is.na(value) || !nzchar(as.character(value))) value <- "NA"
            tagList(
              div(class = "key", field),
              div(class = if (field == "warnings") "warning" else NULL,
                  as.character(value))
            )
          }))
    )
  })

  output$file_panel <- renderUI({
    row <- selected_row()
    validate(need(!is.null(row), "No card selected."))
    file_fields <- intersect(c(
      "rdg_png_file", "state_rdg_png_file", "relative_usage_matrix_file",
      "matched_control_table_file"
    ), names(row))
    links <- lapply(file_fields, function(field) {
      rel <- row[[field]][1]
      if (is.na(rel) || !nzchar(rel)) return(NULL)
      file_info <- make_selected_file_link(rel)
      if (is.null(file_info)) return(NULL)
      file_node <- if (isTRUE(file_info$exists) &&
                       is_scalar_informative(file_info$href)) {
        tags$a(
          class = "file-link",
          href = file_info$href,
          target = "_blank",
          rel = "noopener noreferrer",
          title = file_info$abs_path,
          file_info$abs_path
        )
      } else {
        tagList(
          tags$span(class = "missing-file", file_info$abs_path),
          tags$span(class = "file-status", "missing")
        )
      }
      tags$li(tags$code(field), ": ", file_node)
    })
    div(
      class = "panel",
      h4("Selected-card file pointers"),
      tags$ul(links),
      div(class = "small-note",
          "Existing files open through the local Shiny session. ",
          "Missing files are marked.")
    )
  })
}

shinyApp(ui, server, options = list("launch.browser" = ifelse(interactive(), TRUE, FALSE)))
