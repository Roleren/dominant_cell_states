#!/usr/bin/env Rscript

# Rank merged-replicate RDG rows where clean-CDS allocation changes relative
# to uORF allocation, then annotate what is still missing before interpretation.

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
output_dir <- file.path(analysis_dir, "dominant_rdg_merged_replicate_translon_shift")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

clip01 <- function(x) pmin(1, pmax(0, x))

safe_fread <- function(path) {
  if (!file.exists(path)) return(data.table())
  data.table::fread(path, showProgress = FALSE, nThread = 1)
}

read_required <- function(path, label) {
  dt <- safe_fread(path)
  if (nrow(dt) == 0) {
    stop("Missing or empty ", label, ": ", path, call. = FALSE)
  }
  dt
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
  data.table::fifelse(
    y %chin% c("true", "t", "1", "yes", "y"),
    TRUE,
    data.table::fifelse(
      y %chin% c("false", "f", "0", "no", "n"),
      FALSE,
      default
    )
  )
}

has_text <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(trimws(x))
}

text_has <- function(x, pattern) {
  grepl(pattern, as.character(x), ignore.case = TRUE, perl = TRUE)
}

safe_max <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  max(x)
}

safe_min <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  min(x)
}

collapse_unique <- function(x, max_items = 8) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  if (!length(x)) return("")
  paste(utils::head(x, max_items), collapse = "; ")
}

merge_new_columns <- function(x, y, by) {
  if (nrow(y) == 0) return(x)
  add <- setdiff(names(y), names(x))
  if (!length(add)) return(x)
  merge(x, y[, c(by, add), with = FALSE], by = by, all.x = TRUE,
        sort = FALSE)
}

q_score <- function(q) {
  q <- suppressWarnings(as.numeric(q))
  q[!is.finite(q) | q <= 0] <- 1
  clip01(-log10(pmax(q, 1e-300)) / 6)
}

compact_text <- function(x, max_chars = 240) {
  x <- gsub("\\s+", " ", as.character(x))
  x[is.na(x)] <- ""
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 3), "...")
  x
}

joint_file <- file.path(
  analysis_dir,
  "dominant_rdg_joint_branch_allocation",
  "dominant_rdg_joint_branch_allocation_rows.csv"
)
atlas_file <- file.path(
  analysis_dir,
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_cards.csv"
)
evidence_file <- file.path(
  analysis_dir,
  "dominant_rdg_evidence_gaps",
  "dominant_rdg_evidence_gap_priorities.csv"
)
sdrive_file <- file.path(
  analysis_dir,
  "dominant_rdg_sdrive_readiness",
  "sdrive_transcript_readiness.csv"
)
inspection_file <- file.path(
  analysis_dir,
  "dominant_rdg_qualified_inspection",
  "qualified_inspection_manifest.csv"
)
disagreement_file <- file.path(
  analysis_dir,
  "dominant_rdg_allocation_output_disagreements",
  "dominant_rdg_allocation_output_disagreement_audit.csv"
)

rows_out <- file.path(output_dir, "merged_replicate_translon_shift_rows.csv")
gene_out <- file.path(output_dir, "merged_replicate_translon_shift_gene_design_summary.csv")
queue_out <- file.path(output_dir, "merged_replicate_translon_shift_review_queue.csv")
candidate_out <- file.path(output_dir, "merged_replicate_translon_shift_candidate_queue.csv")
pooling_out <- file.path(output_dir, "merged_replicate_translon_shift_pooling_summary.csv")
metrics_out <- file.path(output_dir, "merged_replicate_translon_shift_summary_metrics.csv")
report_out <- file.path(output_dir, "merged_replicate_translon_shift_report.md")

joint <- read_required(joint_file, "joint branch-allocation rows")
keys3 <- c("gene_symbol", "tx_id", "design_family")
required <- c(
  keys3, "study", "case_condition", "control_conditions", "match_scope",
  "n_case_runs_merged", "n_control_runs_merged",
  "case_branch_total_counts", "control_branch_total_counts",
  "joint_case_count_clean_CDS", "joint_control_count_clean_CDS",
  "joint_case_count_leader_uORF", "joint_control_count_leader_uORF",
  "joint_case_count_overlapping_uORF",
  "joint_control_count_overlapping_uORF",
  "joint_delta_clean_CDS", "joint_delta_leader_uORF",
  "joint_delta_overlapping_uORF", "joint_l1_allocation_shift",
  "joint_chisq_q", "joint_min_runs_merged", "joint_passes_total_gate",
  "joint_passes_class_gate_both", "joint_evaluable", "joint_strict",
  "joint_review_priority", "joint_review_class"
)
missing_required <- setdiff(required, names(joint))
if (length(missing_required)) {
  stop("Joint branch-allocation rows are missing columns: ",
       paste(missing_required, collapse = ", "), call. = FALSE)
}

dt <- copy(joint)

numeric_columns <- intersect(c(
  "n_case_runs_merged", "n_control_runs_merged",
  "case_branch_total_counts", "control_branch_total_counts",
  "joint_case_count_clean_CDS", "joint_control_count_clean_CDS",
  "joint_case_count_leader_uORF", "joint_control_count_leader_uORF",
  "joint_case_count_overlapping_uORF",
  "joint_control_count_overlapping_uORF",
  "joint_delta_clean_CDS", "joint_delta_leader_uORF",
  "joint_delta_overlapping_uORF", "joint_l1_allocation_shift",
  "joint_chisq_q", "joint_min_runs_merged", "joint_review_priority"
), names(dt))
for (column in numeric_columns) {
  set(dt, j = column, value = num_col(dt, column))
}

dt[, clean_cds_abs_delta := abs(joint_delta_clean_CDS)]
dt[, clean_cds_min_count := pmin(joint_case_count_clean_CDS,
                                 joint_control_count_clean_CDS)]
dt[, leader_uorf_abs_delta := abs(joint_delta_leader_uORF)]
dt[, overlapping_uorf_abs_delta := abs(joint_delta_overlapping_uORF)]
dt[, use_overlapping_uorf := overlapping_uorf_abs_delta > leader_uorf_abs_delta]
dt[, dominant_uorf_branch_class := data.table::fifelse(
  use_overlapping_uorf, "overlapping_uORF", "leader_uORF"
)]
dt[, dominant_uorf_delta := data.table::fifelse(
  use_overlapping_uorf,
  joint_delta_overlapping_uORF,
  joint_delta_leader_uORF
)]
dt[, dominant_uorf_abs_delta := abs(dominant_uorf_delta)]
dt[, dominant_uorf_case_count := data.table::fifelse(
  use_overlapping_uorf,
  joint_case_count_overlapping_uORF,
  joint_case_count_leader_uORF
)]
dt[, dominant_uorf_control_count := data.table::fifelse(
  use_overlapping_uorf,
  joint_control_count_overlapping_uORF,
  joint_control_count_leader_uORF
)]
dt[, dominant_uorf_min_count := pmin(dominant_uorf_case_count,
                                     dominant_uorf_control_count)]
dt[, clean_uorf_delta_gap := joint_delta_clean_CDS - dominant_uorf_delta]
dt[, clean_uorf_contrast_abs := abs(clean_uorf_delta_gap)]
dt[, clean_uorf_opposition := sign(joint_delta_clean_CDS) !=
     sign(dominant_uorf_delta) &
     clean_cds_abs_delta >= 0.02 &
     dominant_uorf_abs_delta >= 0.02]

dt[, clean_uorf_shift_class := data.table::fcase(
  clean_uorf_opposition & joint_delta_clean_CDS > 0 &
    dominant_uorf_delta < 0,
  "clean_CDS_up_uORF_down",
  clean_uorf_opposition & joint_delta_clean_CDS < 0 &
    dominant_uorf_delta > 0,
  "clean_CDS_down_uORF_up",
  clean_cds_abs_delta >= 0.05 & dominant_uorf_abs_delta < 0.02,
  "clean_CDS_shift_without_uORF_opposition",
  dominant_uorf_abs_delta >= 0.05 & clean_cds_abs_delta < 0.02,
  "uORF_shift_without_clean_CDS_opposition",
  default = "weak_or_no_clean_uORF_shift"
)]

dt[, merged_replicate_support_class := data.table::fcase(
  joint_min_runs_merged >= 4, "strong_merged_replicate_support",
  joint_min_runs_merged >= 2, "merged_replicate_support",
  default = "single_run_or_unmerged"
)]
dt[, pooling_scope_class := data.table::fcase(
  match_scope == "exact_context" & joint_min_runs_merged >= 2,
  "same_condition_merged_replicates",
  match_scope == "exact_context",
  "same_condition_single_or_unmerged",
  text_has(match_scope, "relaxed|study|cell|tissue|context"),
  "relaxed_context_pooling",
  default = "design_family_or_unclassified_pooling"
)]
dt[, exact_merged_contrast := match_scope == "exact_context" &
     joint_min_runs_merged >= 2]
dt[, clean_uorf_two_sided_count_support :=
     clean_cds_min_count >= 30 & dominant_uorf_min_count >= 30]

atlas <- safe_fread(atlas_file)
if (nrow(atlas) > 0 && all(keys3 %in% names(atlas))) {
  atlas_keep <- intersect(c(
    keys3,
    "atlas_confidence_score", "atlas_browser_review_priority_score",
    "atlas_evidence_class", "atlas_clean_cds_pct_change",
    "abs_atlas_clean_cds_pct_change", "weighted_clean_cds_pct_change",
    "browser_validation_label", "browser_validation_status",
    "browser_validation_score", "browser_signal_grade",
    "transcript_qc_class", "cds_exon_qc_status",
    "n_measured_features", "n_measured_uorf_ouorf",
    "n_measured_internal_orf", "n_measured_nte_ntt",
    "model_agreement_tier", "loo_consensus_class",
    "transport_consensus_class", "context_interaction_gene_class",
    "partial_pooling_gene_class",
    "pooling_proof_tier", "pooling_claim_language",
    "pooling_proof_component", "pooling_hazard_class",
    "pooling_hazard_score", "pooling_hazard_sources",
    "pooling_top_protocol_classes", "pooling_exact_context_fraction",
    "pooling_joint_strict_fraction",
    "pooling_any_start_site_enriched_protocol",
    "pooling_any_relaxed_context",
    "warnings", "warning_items"
  ), names(atlas))
  atlas_small <- unique(atlas[, atlas_keep, with = FALSE], by = keys3)
  dt <- merge_new_columns(dt, atlas_small, by = keys3)
}

evidence <- safe_fread(evidence_file)
if (nrow(evidence) > 0 && all(keys3 %in% names(evidence))) {
  evidence_keep <- intersect(c(
    keys3, "acquisition_rank", "action_label", "claim_stage",
    "acquisition_priority_score", "counterfactual_readiness_score",
    "evidence_gap_score", "biological_value_score",
    "browser_need", "measurement_need", "replication_need",
    "robustness_need", "context_need", "metadata_confounder_need",
    "protocol_pooling_need", "therapeutic_need",
    "negative_control_score", "review_reason", "review_packet_id",
    "selected_for_primary_review", "review_priority_tier",
    "validation_priority_score"
  ), names(evidence))
  evidence_small <- unique(evidence[, evidence_keep, with = FALSE], by = keys3)
  dt <- merge_new_columns(dt, evidence_small, by = keys3)
}

sdrive <- safe_fread(sdrive_file)
if (nrow(sdrive) > 0 && all(keys3 %in% names(sdrive))) {
  sdrive[, sdrive_ready_flag := chr_col(.SD, "sdrive_readiness_class") ==
           "sdrive_ready"]
  sdrive[, static_coverage_reviewable_flag :=
           chr_col(.SD, "review_support_class") == "coverage_reviewable"]
  sdrive[, fst_pages_ready_flag := bool_col(.SD, "fst_pages_ready")]
  sdrive_summary <- sdrive[, .(
    sdrive_ready_rows = sum(sdrive_ready_flag, na.rm = TRUE),
    static_coverage_reviewable_rows = sum(
      static_coverage_reviewable_flag, na.rm = TRUE
    ),
    fst_pages_ready_rows = sum(fst_pages_ready_flag, na.rm = TRUE),
    best_branch_min_case_control_counts = safe_max(
      num_col(.SD, "branch_min_case_control_counts")
    ),
    best_display_min_exact_counts = safe_max(
      num_col(.SD, "display_min_exact_counts")
    ),
    best_top_feature_min_exact_counts = safe_max(
      num_col(.SD, "top_feature_min_exact_counts")
    ),
    worst_top_feature_max_run_fraction = safe_max(
      num_col(.SD, "top_feature_max_run_fraction")
    ),
    mrna_sequence_ready_rows = sum(
      bool_col(.SD, "mrna_sequence_available"), na.rm = TRUE
    ),
    top_feature_in_rdg_annotation_rows = sum(
      bool_col(.SD, "top_feature_in_rdg_annotation"), na.rm = TRUE
    ),
    sdrive_flags = collapse_unique(chr_col(.SD, "flags"))
  ), by = keys3]
  dt <- merge_new_columns(dt, sdrive_summary, by = keys3)
}

inspection <- safe_fread(inspection_file)
if (nrow(inspection) > 0 && all(keys3 %in% names(inspection))) {
  inspection_summary <- inspection[, .(
    qualified_inspection_rows = .N,
    qualified_profile_ready_rows = sum(
      chr_col(.SD, "profile_status") == "profile_ready", na.rm = TRUE
    ),
    qualified_static_ready_rows = sum(
      chr_col(.SD, "inspection_readiness") == "ready", na.rm = TRUE
    ),
    qualified_mrna_sequence_rows = sum(
      bool_col(.SD, "mrna_sequence_available"), na.rm = TRUE
    ),
    qualified_rdg_annotation_rows = sum(
      bool_col(.SD, "top_feature_in_rdg_annotation"), na.rm = TRUE
    ),
    qualified_flags = collapse_unique(chr_col(.SD, "flags"))
  ), by = keys3]
  dt <- merge_new_columns(dt, inspection_summary, by = keys3)
}

disagreement <- safe_fread(disagreement_file)
if (nrow(disagreement) > 0 && all(keys3 %in% names(disagreement))) {
  disagreement_keep <- intersect(c(
    keys3,
    "allocation_output_disagreement_flag",
    "clean_cds_output_allocation_alignment",
    "clean_cds_output_pct",
    "clean_cds_output_pct_lower",
    "clean_cds_output_pct_upper",
    "clean_cds_output_log2fc",
    "clean_cds_allocation_delta",
    "leader_uorf_allocation_delta",
    "overlapping_uorf_allocation_delta",
    "upstream_allocation_delta",
    "output_allocation_direction_pair",
    "total_load_compensation_proxy",
    "upstream_shift_proxy",
    "disagreement_audit_class",
    "disagreement_review_priority"
  ), names(disagreement))
  disagreement_small <- unique(
    disagreement[, disagreement_keep, with = FALSE],
    by = keys3
  )
  dt <- merge_new_columns(dt, disagreement_small, by = keys3)
}

fill_zero_cols <- c(
  "pooling_hazard_score", "protocol_pooling_need", "browser_need",
  "measurement_need", "replication_need", "robustness_need",
  "context_need", "metadata_confounder_need",
  "counterfactual_readiness_score", "acquisition_priority_score",
  "sdrive_ready_rows", "static_coverage_reviewable_rows",
  "fst_pages_ready_rows", "best_branch_min_case_control_counts",
  "best_display_min_exact_counts", "best_top_feature_min_exact_counts",
  "worst_top_feature_max_run_fraction", "mrna_sequence_ready_rows",
  "top_feature_in_rdg_annotation_rows", "qualified_inspection_rows",
  "qualified_profile_ready_rows", "qualified_static_ready_rows",
  "qualified_mrna_sequence_rows", "qualified_rdg_annotation_rows",
  "clean_cds_output_pct", "clean_cds_output_pct_lower",
  "clean_cds_output_pct_upper", "clean_cds_output_log2fc",
  "clean_cds_allocation_delta", "leader_uorf_allocation_delta",
  "overlapping_uorf_allocation_delta", "upstream_allocation_delta",
  "total_load_compensation_proxy", "upstream_shift_proxy",
  "disagreement_review_priority"
)
for (column in setdiff(fill_zero_cols, names(dt))) {
  dt[, (column) := 0]
}
for (column in intersect(fill_zero_cols, names(dt))) {
  set(dt, j = column, value = num_col(dt, column))
}

fill_text_cols <- c(
  "pooling_proof_tier", "pooling_claim_language", "pooling_hazard_class",
  "pooling_hazard_sources", "pooling_top_protocol_classes",
  "browser_validation_label", "browser_validation_status",
  "browser_signal_grade", "transcript_qc_class", "cds_exon_qc_status",
  "model_agreement_tier", "loo_consensus_class",
  "transport_consensus_class", "context_interaction_gene_class",
  "partial_pooling_gene_class", "action_label", "claim_stage",
  "review_reason", "review_priority_tier", "sdrive_flags",
  "qualified_flags", "clean_cds_output_allocation_alignment",
  "output_allocation_direction_pair", "disagreement_audit_class"
)
for (column in setdiff(fill_text_cols, names(dt))) {
  dt[, (column) := ""]
}
for (column in intersect(fill_text_cols, names(dt))) {
  set(dt, j = column, value = chr_col(dt, column))
}

for (column in c("pooling_any_start_site_enriched_protocol",
                 "pooling_any_relaxed_context")) {
  if (!column %in% names(dt)) dt[, (column) := FALSE]
  if (column %in% names(dt)) set(dt, j = column, value = bool_col(dt, column))
}
if (!"allocation_output_disagreement_flag" %in% names(dt)) {
  dt[, allocation_output_disagreement_flag := FALSE]
}
dt[, allocation_output_disagreement_flag :=
     bool_col(.SD, "allocation_output_disagreement_flag")]

dt[, missing_merged_replicates := joint_min_runs_merged < 2]
dt[, missing_branch_total_counts := !bool_col(.SD, "joint_passes_total_gate")]
dt[, missing_two_sided_uorf_counts :=
     !bool_col(.SD, "joint_passes_class_gate_both") |
     !clean_uorf_two_sided_count_support]
dt[, missing_clean_uorf_opposition := !clean_uorf_opposition]
dt[, missing_exact_matched_context := match_scope != "exact_context"]
dt[, missing_browser_label := !has_text(chr_col(.SD, "browser_validation_label"))]
dt[, missing_transport_support := !text_has(
  chr_col(.SD, "transport_consensus_class"),
  "support-robust|direction-robust"
)]
dt[, missing_loo_support := !text_has(
  chr_col(.SD, "loo_consensus_class"),
  "support-robust|direction-robust"
)]
dt[, missing_partial_pooling_support := text_has(
  chr_col(.SD, "partial_pooling_gene_class"),
  "not|single|fragile|non"
) | !has_text(chr_col(.SD, "partial_pooling_gene_class"))]
dt[, protocol_pooling_audit_needed :=
     num_col(.SD, "pooling_hazard_score") >= 0.55 |
     num_col(.SD, "protocol_pooling_need") >= 0.55 |
     bool_col(.SD, "pooling_any_start_site_enriched_protocol") |
     bool_col(.SD, "pooling_any_relaxed_context")]
dt[, missing_static_coverage_support :=
     num_col(.SD, "static_coverage_reviewable_rows") <= 0 &
     num_col(.SD, "qualified_static_ready_rows") <= 0]
dt[, single_run_dominance_static_asset :=
     num_col(.SD, "worst_top_feature_max_run_fraction") > 0.65]
dt[, missing_sdrive_backing_data := num_col(.SD, "sdrive_ready_rows") <= 0]
dt[, missing_mrna_sequence := num_col(.SD, "mrna_sequence_ready_rows") <= 0 &
     num_col(.SD, "qualified_mrna_sequence_rows") <= 0]
dt[, missing_rdg_annotation_track :=
     num_col(.SD, "top_feature_in_rdg_annotation_rows") <= 0 &
     num_col(.SD, "qualified_rdg_annotation_rows") <= 0]

flag_columns <- c(
  "missing_merged_replicates", "missing_branch_total_counts",
  "missing_two_sided_uorf_counts", "missing_clean_uorf_opposition",
  "missing_exact_matched_context", "missing_browser_label",
  "missing_transport_support", "missing_loo_support",
  "missing_partial_pooling_support", "protocol_pooling_audit_needed",
  "missing_static_coverage_support", "single_run_dominance_static_asset",
  "missing_sdrive_backing_data", "missing_mrna_sequence",
  "missing_rdg_annotation_track"
)
flag_df <- as.data.frame(dt[, flag_columns, with = FALSE])
dt[, missing_evidence_flags := apply(flag_df, 1, function(x) {
  hits <- names(flag_df)[which(as.logical(x) & !is.na(x))]
  if (!length(hits)) "none" else paste(hits, collapse = ";")
})]

dt[, opposition_score_component := data.table::fifelse(
  clean_uorf_opposition,
  clip01(clean_uorf_contrast_abs / 0.35),
  0.35 * clip01(clean_uorf_contrast_abs / 0.35)
)]
dt[, l1_score_component := clip01(joint_l1_allocation_shift / 0.55)]
dt[, q_score_component := q_score(joint_chisq_q)]
dt[, merged_score_component := clip01((joint_min_runs_merged - 1) / 3)]
dt[, count_score_component := clip01(
  log1p(pmin(clean_cds_min_count, dominant_uorf_min_count)) / log1p(300)
)]
dt[, exact_context_component := data.table::fifelse(
  match_scope == "exact_context", 1, 0.35
)]
dt[, static_review_component := data.table::fifelse(
  num_col(.SD, "static_coverage_reviewable_rows") > 0 |
    num_col(.SD, "qualified_static_ready_rows") > 0,
  1,
  data.table::fifelse(num_col(.SD, "sdrive_ready_rows") > 0, 0.55, 0)
)]
dt[, model_support_component := pmax(
  data.table::fifelse(text_has(chr_col(.SD, "transport_consensus_class"),
                               "support-robust"), 1,
                      data.table::fifelse(text_has(
                        chr_col(.SD, "transport_consensus_class"),
                        "direction-robust"
                      ), 0.65, 0)),
  data.table::fifelse(text_has(chr_col(.SD, "loo_consensus_class"),
                               "support-robust"), 0.8,
                      data.table::fifelse(text_has(
                        chr_col(.SD, "loo_consensus_class"),
                        "direction-robust"
                      ), 0.5, 0))
)]
dt[, design_focus_component := data.table::fifelse(
  text_has(design_family, "viral|infection|starvation|ISR|stress"),
  1,
  0.65
)]
dt[, hazard_penalty :=
     0.10 * clip01(num_col(.SD, "pooling_hazard_score")) +
     0.08 * as.numeric(protocol_pooling_audit_needed) +
     0.06 * as.numeric(missing_exact_matched_context) +
     0.05 * as.numeric(missing_two_sided_uorf_counts) +
     0.08 * as.numeric(missing_static_coverage_support) +
     0.05 * as.numeric(missing_sdrive_backing_data) +
     0.03 * as.numeric(missing_mrna_sequence) +
     0.03 * as.numeric(missing_rdg_annotation_track) +
     0.04 * as.numeric(single_run_dominance_static_asset)]

dt[, discovery_priority_score := clip01(
  0.23 * opposition_score_component +
    0.16 * l1_score_component +
    0.11 * q_score_component +
    0.11 * merged_score_component +
    0.12 * count_score_component +
    0.08 * exact_context_component +
    0.08 * static_review_component +
    0.07 * model_support_component +
    0.04 * design_focus_component -
    hazard_penalty
)]

inspection_audit_available <- any(
  num_col(dt, "sdrive_ready_rows") > 0 |
    num_col(dt, "qualified_inspection_rows") > 0,
  na.rm = TRUE
)
dt[, inspection_possible := num_col(.SD, "sdrive_ready_rows") > 0 |
     num_col(.SD, "static_coverage_reviewable_rows") > 0 |
     num_col(.SD, "qualified_inspection_rows") > 0 |
     num_col(.SD, "qualified_static_ready_rows") > 0]
if (!inspection_audit_available) {
  dt[, inspection_possible := TRUE]
}

dt[, claim_readiness_class := data.table::fcase(
  exact_merged_contrast & clean_uorf_opposition &
    clean_uorf_two_sided_count_support &
    !missing_static_coverage_support & !protocol_pooling_audit_needed &
    !single_run_dominance_static_asset & !missing_browser_label,
  "claim_candidate_exact_merged_review_ready",
  protocol_pooling_audit_needed,
  "pooled_or_protocol_confounded",
  missing_merged_replicates | missing_branch_total_counts |
    missing_two_sided_uorf_counts,
  "measurement_limited",
  missing_static_coverage_support | missing_sdrive_backing_data |
    missing_mrna_sequence | missing_rdg_annotation_track,
  "inspection_asset_limited",
  exact_merged_contrast & clean_uorf_opposition &
    clean_uorf_two_sided_count_support &
    missing_browser_label,
  "review_ready_missing_browser_label",
  clean_uorf_opposition,
  "biological_review_lead",
  default = "low_contrast_or_background"
)]

setorder(dt, -discovery_priority_score, -clean_uorf_opposition,
         -clean_uorf_contrast_abs, gene_symbol, design_family, study)
dt[, merged_translon_shift_rank := seq_len(.N)]

queue <- dt[
  discovery_priority_score >= 0.40 &
    joint_min_runs_merged >= 2 &
    inspection_possible &
    (clean_uorf_opposition | clean_uorf_contrast_abs >= 0.10)
]
setorder(queue, -discovery_priority_score, -clean_uorf_opposition,
         -clean_uorf_contrast_abs, gene_symbol, design_family, study)
queue[, review_queue_rank := seq_len(.N)]
queue[, review_focus := data.table::fcase(
  clean_uorf_shift_class == "clean_CDS_up_uORF_down",
  "CDS gain with uORF branch loss",
  clean_uorf_shift_class == "clean_CDS_down_uORF_up",
  "CDS loss with uORF branch gain",
  clean_cds_abs_delta >= dominant_uorf_abs_delta,
  "CDS shift dominates uORF contrast",
  default = "uORF shift dominates CDS contrast"
)]
queue[, review_instruction := compact_text(paste(
  "Inspect merged replicate contrast",
  paste0(case_condition, " vs ", control_conditions),
  "in", study,
  "for", review_focus,
  "using exact coverage, RDG-flow annotation, and replicate diagnostics.",
  "Missing:", missing_evidence_flags
))]

dt[, in_review_queue := merged_translon_shift_rank %in%
     queue$merged_translon_shift_rank]

candidate_keys <- c(keys3, "clean_uorf_shift_class")
candidate_queue <- data.table()
if (nrow(queue) > 0) {
  candidate_queue <- queue[, {
    top <- .SD[which.max(discovery_priority_score)]
    .(
      candidate_rows = .N,
      candidate_studies = uniqueN(study),
      candidate_case_conditions = uniqueN(case_condition),
      candidate_control_conditions = uniqueN(control_conditions),
      exact_context_rows = sum(match_scope == "exact_context", na.rm = TRUE),
      exact_merged_rows = sum(exact_merged_contrast, na.rm = TRUE),
      strict_rows = sum(bool_col(.SD, "joint_strict"), na.rm = TRUE),
      evaluable_rows = sum(bool_col(.SD, "joint_evaluable"), na.rm = TRUE),
      max_discovery_priority_score = max(discovery_priority_score,
                                         na.rm = TRUE),
      median_discovery_priority_score = stats::median(
        discovery_priority_score, na.rm = TRUE
      ),
      max_clean_uorf_contrast_abs = max(clean_uorf_contrast_abs,
                                        na.rm = TRUE),
      median_clean_uorf_contrast_abs = stats::median(
        clean_uorf_contrast_abs, na.rm = TRUE
      ),
      median_clean_cds_delta = stats::median(joint_delta_clean_CDS,
                                             na.rm = TRUE),
      median_dominant_uorf_delta = stats::median(dominant_uorf_delta,
                                                 na.rm = TRUE),
      max_joint_min_runs_merged = max(joint_min_runs_merged, na.rm = TRUE),
      min_joint_min_runs_merged = min(joint_min_runs_merged, na.rm = TRUE),
      max_clean_cds_min_count = max(clean_cds_min_count, na.rm = TRUE),
      max_dominant_uorf_min_count = max(dominant_uorf_min_count,
                                        na.rm = TRUE),
      protocol_or_pooling_audit_rows = sum(
        protocol_pooling_audit_needed, na.rm = TRUE
      ),
      measurement_limited_rows = sum(
        claim_readiness_class == "measurement_limited", na.rm = TRUE
      ),
      inspection_asset_limited_rows = sum(
        claim_readiness_class == "inspection_asset_limited", na.rm = TRUE
      ),
      browser_label_missing_rows = sum(missing_browser_label, na.rm = TRUE),
      static_coverage_missing_rows = sum(
        missing_static_coverage_support, na.rm = TRUE
      ),
      single_run_dominance_rows = sum(
        single_run_dominance_static_asset, na.rm = TRUE
      ),
      claim_candidate_rows = sum(
        claim_readiness_class ==
          "claim_candidate_exact_merged_review_ready",
        na.rm = TRUE
      ),
      allocation_output_disagreement_rows = sum(
        bool_col(.SD, "allocation_output_disagreement_flag"), na.rm = TRUE
      ),
      max_disagreement_review_priority = max(
        num_col(.SD, "disagreement_review_priority"), na.rm = TRUE
      ),
      clean_cds_output_pct = top$clean_cds_output_pct[1],
      clean_cds_output_allocation_alignment =
        top$clean_cds_output_allocation_alignment[1],
      output_allocation_direction_pair =
        top$output_allocation_direction_pair[1],
      disagreement_audit_class = top$disagreement_audit_class[1],
      top_review_queue_rank = top$review_queue_rank[1],
      top_study = top$study[1],
      top_case_condition = top$case_condition[1],
      top_control_conditions = top$control_conditions[1],
      top_match_scope = top$match_scope[1],
      top_clean_cds_delta = top$joint_delta_clean_CDS[1],
      top_dominant_uorf_branch_class = top$dominant_uorf_branch_class[1],
      top_dominant_uorf_delta = top$dominant_uorf_delta[1],
      top_claim_readiness_class = top$claim_readiness_class[1],
      top_missing_evidence_flags = top$missing_evidence_flags[1],
      studies = collapse_unique(study, 10),
      case_conditions = collapse_unique(case_condition, 10),
      claim_readiness_classes = collapse_unique(claim_readiness_class, 6),
      missing_evidence_flags = collapse_unique(missing_evidence_flags, 8)
    )
  }, by = candidate_keys]

  direction_summary <- queue[, .(
    queue_rows_for_gene_design = .N,
    direction_classes_for_gene_design = uniqueN(clean_uorf_shift_class),
    direction_class_list = collapse_unique(clean_uorf_shift_class, 6)
  ), by = keys3]
  candidate_queue <- merge(
    candidate_queue,
    direction_summary,
    by = keys3,
    all.x = TRUE,
    sort = FALSE
  )
  candidate_queue[, opposite_direction_rows :=
                    queue_rows_for_gene_design - candidate_rows]
  candidate_queue[, recurrence_component := clip01(
    0.45 * clip01((pmin(candidate_rows, 5) - 1) / 4) +
      0.35 * clip01((pmin(candidate_studies, 4) - 1) / 3) +
      0.20 * clip01((pmin(candidate_case_conditions, 4) - 1) / 3)
  )]
  candidate_queue[, exact_scope_component := data.table::fifelse(
    exact_context_rows > 0, exact_context_rows / pmax(candidate_rows, 1), 0
  )]
  candidate_queue[, output_context_component := data.table::fifelse(
    allocation_output_disagreement_rows > 0 |
      has_text(clean_cds_output_allocation_alignment),
    1, 0
  )]
  candidate_queue[, static_component := data.table::fifelse(
    static_coverage_missing_rows == 0, 1,
    data.table::fifelse(inspection_asset_limited_rows < candidate_rows,
                        0.55, 0)
  )]
  candidate_queue[, candidate_priority_score := clip01(
    0.52 * max_discovery_priority_score +
      0.16 * recurrence_component +
      0.10 * exact_scope_component +
      0.07 * output_context_component +
      0.10 * static_component +
      0.05 * clip01(max_disagreement_review_priority) -
      0.07 * (protocol_or_pooling_audit_rows / pmax(candidate_rows, 1)) -
      0.06 * (measurement_limited_rows / pmax(candidate_rows, 1)) -
      0.05 * (inspection_asset_limited_rows / pmax(candidate_rows, 1)) -
      0.05 * as.numeric(opposite_direction_rows > 0)
  )]
  candidate_queue[, recurrence_class := data.table::fcase(
    candidate_studies >= 2 & exact_context_rows >= 2,
    "multi_study_exact_recurrence",
    candidate_case_conditions >= 2 & exact_context_rows >= 2,
    "multi_condition_exact_recurrence",
    exact_context_rows >= 1,
    "single_exact_context",
    default = "pooled_or_relaxed_only"
  )]
  candidate_queue[, pooling_benefit_class := data.table::fcase(
    candidate_studies >= 2,
    "recurrence_across_studies",
    candidate_case_conditions >= 2,
    "condition_series_within_study",
    max_joint_min_runs_merged >= 4,
    "strong_replicate_depth",
    default = "single_context_depth_gain"
  )]
  candidate_queue[, pooling_loss_class := data.table::fcase(
    protocol_or_pooling_audit_rows > 0,
    "protocol_or_start_site_bias_risk",
    opposite_direction_rows > 0,
    "context_direction_conflict",
    measurement_limited_rows > 0,
    "measurement_count_limit",
    inspection_asset_limited_rows > 0,
    "coverage_asset_limit",
    browser_label_missing_rows > 0,
    "browser_label_missing",
    default = "no_major_pooling_loss_flag"
  )]
  candidate_queue[, candidate_interpretation_class := data.table::fcase(
    claim_candidate_rows > 0,
    "claim_candidate_after_browser_review",
    protocol_or_pooling_audit_rows > 0,
    "protocol_pooling_audit_lead",
    opposite_direction_rows > 0,
    "context_direction_conflict_lead",
    measurement_limited_rows > 0,
    "measurement_limited_candidate",
    inspection_asset_limited_rows > 0,
    "coverage_asset_limited_candidate",
    browser_label_missing_rows > 0,
    "browser_validation_candidate",
    default = "biological_review_candidate"
  )]
  candidate_queue[, next_review_action := data.table::fcase(
    protocol_or_pooling_audit_rows > 0,
    "audit_protocol_pooling_and_start_site_bias",
    opposite_direction_rows > 0,
    "split_contexts_before_claim",
    measurement_limited_rows > 0,
    "rebuild_or_review_translon_measurement",
    inspection_asset_limited_rows > 0 | static_coverage_missing_rows > 0,
    "generate_or_open_static_coverage_assets",
    browser_label_missing_rows > 0,
    "browser_validate_exact_contrast",
    default = "write_claim_candidate_review"
  )]
  setorder(candidate_queue, -candidate_priority_score,
           -candidate_rows, -candidate_studies, gene_symbol, design_family)
  candidate_queue[, candidate_queue_rank := seq_len(.N)]
  setcolorder(candidate_queue, c(
    "candidate_queue_rank", keys3, "clean_uorf_shift_class",
    setdiff(names(candidate_queue),
            c("candidate_queue_rank", keys3, "clean_uorf_shift_class"))
  ))
}

pooling_summary <- dt[, .(
  rows = .N,
  review_queue_rows = sum(in_review_queue, na.rm = TRUE),
  genes = uniqueN(gene_symbol),
  gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
  clean_uorf_opposition_rows = sum(clean_uorf_opposition, na.rm = TRUE),
  exact_merged_rows = sum(exact_merged_contrast, na.rm = TRUE),
  strict_rows = sum(bool_col(.SD, "joint_strict"), na.rm = TRUE),
  median_clean_uorf_contrast_abs = stats::median(
    clean_uorf_contrast_abs, na.rm = TRUE
  ),
  max_discovery_priority_score = max(discovery_priority_score, na.rm = TRUE),
  protocol_or_pooling_audit_rows = sum(
    protocol_pooling_audit_needed, na.rm = TRUE
  ),
  measurement_limited_rows = sum(
    claim_readiness_class == "measurement_limited", na.rm = TRUE
  ),
  inspection_asset_limited_rows = sum(
    claim_readiness_class == "inspection_asset_limited", na.rm = TRUE
  )
), by = .(design_family, pooling_scope_class, clean_uorf_shift_class,
          claim_readiness_class)]
setorder(pooling_summary, design_family, -review_queue_rows,
         -max_discovery_priority_score)

gene_summary <- dt[, {
  top <- .SD[which.max(discovery_priority_score)]
  .(
    merged_shift_rows = .N,
    review_queue_rows = sum(discovery_priority_score >= 0.40 &
                              joint_min_runs_merged >= 2 &
                              inspection_possible &
                              (clean_uorf_opposition |
                                 clean_uorf_contrast_abs >= 0.10),
                            na.rm = TRUE),
    exact_merged_opposition_rows = sum(exact_merged_contrast &
                                         clean_uorf_opposition, na.rm = TRUE),
    claim_candidate_rows = sum(
      claim_readiness_class ==
        "claim_candidate_exact_merged_review_ready",
      na.rm = TRUE
    ),
    protocol_or_pooling_audit_rows = sum(
      protocol_pooling_audit_needed, na.rm = TRUE
    ),
    max_discovery_priority_score = max(discovery_priority_score, na.rm = TRUE),
    top_study = top$study[1],
    top_case_condition = top$case_condition[1],
    top_control_conditions = top$control_conditions[1],
    top_match_scope = top$match_scope[1],
    top_shift_class = top$clean_uorf_shift_class[1],
    top_claim_readiness_class = top$claim_readiness_class[1],
    top_clean_cds_delta = top$joint_delta_clean_CDS[1],
    top_dominant_uorf_branch_class = top$dominant_uorf_branch_class[1],
    top_dominant_uorf_delta = top$dominant_uorf_delta[1],
    top_missing_evidence_flags = top$missing_evidence_flags[1],
    collapsed_candidate_rows = if (nrow(candidate_queue)) {
      candidate_queue[
        gene_symbol == .BY$gene_symbol &
          tx_id == .BY$tx_id &
          design_family == .BY$design_family,
        .N
      ]
    } else {
      0L
    }
  )
}, by = keys3]
setorder(gene_summary, -max_discovery_priority_score, gene_symbol,
         design_family)
gene_summary[, merged_shift_gene_design_rank := seq_len(.N)]

metrics <- data.table(
  metric = c(
    "joint_rows",
    "exact_merged_rows",
    "clean_uorf_opposition_rows",
    "exact_merged_opposition_rows",
    "review_queue_rows",
    "candidate_queue_rows",
    "candidate_gene_designs",
    "multi_study_exact_recurrence_candidates",
    "context_direction_conflict_candidates",
    "protocol_pooling_audit_candidates",
    "claim_candidate_exact_merged_review_ready_rows",
    "review_ready_missing_browser_label_rows",
    "pooled_or_protocol_confounded_rows",
    "measurement_limited_rows",
    "inspection_asset_limited_rows",
    "clean_CDS_up_uORF_down_rows",
    "clean_CDS_down_uORF_up_rows",
    "viral_infection_review_queue_rows",
    "nutrient_starvation_review_queue_rows",
    "tumor_cancer_context_review_queue_rows",
    "static_coverage_reviewable_rows",
    "sdrive_ready_rows"
  ),
  value = as.character(c(
    nrow(dt),
    sum(dt$exact_merged_contrast, na.rm = TRUE),
    sum(dt$clean_uorf_opposition, na.rm = TRUE),
    sum(dt$exact_merged_contrast & dt$clean_uorf_opposition, na.rm = TRUE),
    nrow(queue),
    nrow(candidate_queue),
    uniqueN(paste(candidate_queue$gene_symbol,
                  candidate_queue$tx_id,
                  candidate_queue$design_family)),
    sum(candidate_queue$recurrence_class ==
          "multi_study_exact_recurrence", na.rm = TRUE),
    sum(candidate_queue$opposite_direction_rows > 0, na.rm = TRUE),
    sum(candidate_queue$protocol_or_pooling_audit_rows > 0, na.rm = TRUE),
    sum(dt$claim_readiness_class ==
          "claim_candidate_exact_merged_review_ready", na.rm = TRUE),
    sum(dt$claim_readiness_class ==
          "review_ready_missing_browser_label", na.rm = TRUE),
    sum(dt$claim_readiness_class ==
          "pooled_or_protocol_confounded", na.rm = TRUE),
    sum(dt$claim_readiness_class == "measurement_limited", na.rm = TRUE),
    sum(dt$claim_readiness_class == "inspection_asset_limited", na.rm = TRUE),
    sum(dt$clean_uorf_shift_class == "clean_CDS_up_uORF_down", na.rm = TRUE),
    sum(dt$clean_uorf_shift_class == "clean_CDS_down_uORF_up", na.rm = TRUE),
    sum(queue$design_family == "viral_infection", na.rm = TRUE),
    sum(queue$design_family == "nutrient_starvation", na.rm = TRUE),
    sum(queue$design_family == "tumor_cancer_context", na.rm = TRUE),
    sum(num_col(dt, "static_coverage_reviewable_rows") > 0 |
          num_col(dt, "qualified_static_ready_rows") > 0, na.rm = TRUE),
    sum(num_col(dt, "sdrive_ready_rows") > 0, na.rm = TRUE)
  ))
)

data.table::fwrite(dt, rows_out)
data.table::fwrite(gene_summary, gene_out)
data.table::fwrite(queue, queue_out)
data.table::fwrite(candidate_queue, candidate_out)
data.table::fwrite(pooling_summary, pooling_out)
data.table::fwrite(metrics, metrics_out)

plot_dt <- dt[clean_uorf_opposition | discovery_priority_score >= 0.40]
if (nrow(plot_dt) > 0) {
  if (nrow(plot_dt) > 5000) plot_dt <- plot_dt[seq_len(5000)]
  scatter <- ggplot(
    plot_dt,
    aes(x = joint_delta_clean_CDS, y = dominant_uorf_delta,
        color = claim_readiness_class, size = discovery_priority_score)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey65") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey65") +
    geom_abline(slope = -1, intercept = 0, linewidth = 0.25,
                linetype = "dashed", color = "grey55") +
    geom_point(alpha = 0.65) +
    scale_size_continuous(range = c(0.8, 3.2), limits = c(0, 1)) +
    labs(
      title = "Merged replicate clean-CDS versus strongest uORF shift",
      x = "Clean-CDS allocation delta",
      y = "Strongest uORF allocation delta",
      color = "Claim gate",
      size = "Priority"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(figure_dir, "merged_replicate_translon_shift_clean_uorf_scatter.png"),
    scatter, width = 8.5, height = 6.2, dpi = 160
  )
}

top_queue <- queue[seq_len(min(.N, 30))]
if (nrow(top_queue) > 0) {
  top_queue[, plot_label := paste(
    paste0("#", review_queue_rank),
    gene_symbol,
    design_family,
    study,
    sep = " | "
  )]
  top_queue[, plot_label := factor(plot_label, levels = rev(plot_label))]
  bar <- ggplot(
    top_queue,
    aes(x = plot_label, y = discovery_priority_score,
        fill = clean_uorf_shift_class)
  ) +
    geom_col(width = 0.74) +
    coord_flip() +
    labs(
      title = "Top merged-replicate translon-shift review rows",
      x = NULL,
      y = "Discovery priority",
      fill = "Shift class"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(figure_dir, "merged_replicate_translon_shift_priority_top.png"),
    bar, width = 9, height = 7.2, dpi = 160
  )
}

top_candidate <- candidate_queue[seq_len(min(.N, 30))]
if (nrow(top_candidate) > 0) {
  top_candidate[, plot_label := paste(
    paste0("#", candidate_queue_rank),
    gene_symbol,
    design_family,
    clean_uorf_shift_class,
    sep = " | "
  )]
  top_candidate[, plot_label := factor(plot_label, levels = rev(plot_label))]
  candidate_plot <- ggplot(
    top_candidate,
    aes(x = plot_label, y = candidate_priority_score,
        fill = candidate_interpretation_class)
  ) +
    geom_col(width = 0.74) +
    coord_flip() +
    labs(
      title = "Collapsed merged-translon candidate queue",
      x = NULL,
      y = "Candidate priority",
      fill = "Interpretation"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(figure_dir, "merged_replicate_translon_shift_candidate_top.png"),
    candidate_plot, width = 10, height = 7.4, dpi = 160
  )
}

top_lines <- character()
if (nrow(queue) > 0) {
  top_rows <- queue[seq_len(min(.N, 20))]
  top_lines <- paste0(
    seq_len(nrow(top_rows)), ". `", top_rows$gene_symbol, " | ",
    top_rows$design_family, " | ", top_rows$study, "`: ",
    top_rows$case_condition, " vs ", top_rows$control_conditions,
    "; clean-CDS delta ",
    sprintf("%.3f", top_rows$joint_delta_clean_CDS),
    ", ", top_rows$dominant_uorf_branch_class, " delta ",
    sprintf("%.3f", top_rows$dominant_uorf_delta),
    "; ", top_rows$claim_readiness_class,
    "; missing `", top_rows$missing_evidence_flags, "`."
  )
}

candidate_lines <- character()
if (nrow(candidate_queue) > 0) {
  top_candidates <- candidate_queue[seq_len(min(.N, 20))]
  candidate_lines <- paste0(
    seq_len(nrow(top_candidates)), ". `",
    top_candidates$gene_symbol, " | ",
    top_candidates$design_family, " | ",
    top_candidates$clean_uorf_shift_class, "`: score ",
    sprintf("%.3f", top_candidates$candidate_priority_score),
    "; ", top_candidates$candidate_rows, " row(s), ",
    top_candidates$candidate_studies, " study/studies, ",
    top_candidates$candidate_case_conditions, " case condition(s); ",
    top_candidates$candidate_interpretation_class,
    "; next `", top_candidates$next_review_action, "`."
  )
}

metric_value <- function(name) metrics[metric == name, value][1]
report <- c(
  "# Merged-Replicate Translon Shift Review",
  "",
  paste0("Input joint rows: ", metric_value("joint_rows"), "."),
  paste0("Exact same-condition merged rows: ",
         metric_value("exact_merged_rows"), "."),
  paste0("Clean-CDS/uORF opposition rows: ",
         metric_value("clean_uorf_opposition_rows"), "."),
  paste0("Exact merged opposition rows: ",
         metric_value("exact_merged_opposition_rows"), "."),
  paste0("Review queue rows: ", metric_value("review_queue_rows"), "."),
  paste0("Collapsed candidate rows: ",
         metric_value("candidate_queue_rows"), "."),
  "",
  "## Interpretation",
  "",
  paste(
    "This table is a claim gate for the current discovery goal: cases where",
    "merged replicates show changed clean-CDS allocation relative to uORF",
    "allocation. The strongest class is an exact matched condition contrast",
    "with at least two runs on both sides, two-sided clean-CDS/uORF counts,",
    "coverage assets, and no immediate protocol/pooling hazard. Rows outside",
    "that class remain useful review leads, but their missing-evidence flags",
    "should be resolved before biological claims."
  ),
  "",
  "## Top Collapsed Candidates",
  "",
  if (length(candidate_lines)) {
    candidate_lines
  } else {
    "No collapsed candidates passed the review queue gate."
  },
  "",
  "## Top Review Rows",
  "",
  if (length(top_lines)) top_lines else "No rows passed the review queue gate.",
  "",
  "## Output Files",
  "",
  paste0("- `", basename(rows_out), "`: all joint rows with clean-CDS/uORF ",
         "contrast, claim-readiness, and missing-evidence flags."),
  paste0("- `", basename(queue_out), "`: prioritized manual/Codex review queue."),
  paste0("- `", basename(candidate_out), "`: collapsed candidate queue with ",
         "recurrence, pooling tradeoff, and next-review-action fields."),
  paste0("- `", basename(pooling_out), "`: design/pooling-scope summary."),
  paste0("- `", basename(gene_out), "`: gene/transcript/design summary."),
  paste0("- `", basename(metrics_out), "`: compact run metrics."),
  "- `figures/merged_replicate_translon_shift_clean_uorf_scatter.png`: clean-CDS versus strongest uORF delta.",
  "- `figures/merged_replicate_translon_shift_priority_top.png`: top queue rows.",
  "- `figures/merged_replicate_translon_shift_candidate_top.png`: top collapsed candidates."
)
writeLines(report, report_out)

message("Wrote merged-replicate translon shift outputs to: ", output_dir)
