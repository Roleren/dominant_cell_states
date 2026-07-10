#!/usr/bin/env Rscript

# Build browser-review dossiers from the RDG consensus atlas.
# This step is deliberately not another discovery model. It turns the current
# consensus queue into concrete validation packets: what to open, what pattern
# to expect, what context to check first, and what warnings should be considered
# before treating a row as biology.

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
output_dir <- file.path(analysis_dir, "dominant_rdg_validation_dossiers")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

consensus_file <- file.path(
  analysis_dir,
  "dominant_rdg_consensus",
  "dominant_rdg_consensus_candidates.csv"
)
review_file <- file.path(
  analysis_dir,
  "dominant_rdg_consensus",
  "dominant_rdg_consensus_review_queue.csv"
)
evidence_file <- file.path(
  analysis_dir,
  "dominant_rdg_consensus",
  "dominant_rdg_consensus_evidence_matrix.csv"
)
branch_context_file <- file.path(
  analysis_dir,
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_branch_contexts.csv"
)
perturb_file <- file.path(
  analysis_dir,
  "dominant_rdg_therapeutic_perturbation_predictions",
  "dominant_rdg_therapeutic_perturbation_gene_predictions.csv"
)

candidate_out <- file.path(
  output_dir,
  "dominant_rdg_validation_dossier_candidates.csv"
)
context_out <- file.path(
  output_dir,
  "dominant_rdg_validation_dossier_branch_contexts.csv"
)
template_out <- file.path(
  output_dir,
  "dominant_rdg_validation_dossier_manual_review_template.csv"
)
gene_out <- file.path(
  output_dir,
  "dominant_rdg_validation_dossier_gene_summary.csv"
)
metrics_out <- file.path(
  output_dir,
  "dominant_rdg_validation_dossier_summary_metrics.csv"
)
report_out <- file.path(
  output_dir,
  "dominant_rdg_validation_dossier_report.md"
)

if (!file.exists(consensus_file)) {
  stop("Missing consensus candidate file: ", consensus_file, call. = FALSE)
}
if (!file.exists(review_file)) {
  stop("Missing consensus review queue file: ", review_file, call. = FALSE)
}
if (!file.exists(evidence_file)) {
  stop("Missing consensus evidence matrix file: ", evidence_file, call. = FALSE)
}
if (!file.exists(branch_context_file)) {
  stop("Missing atlas branch-context file: ", branch_context_file,
       call. = FALSE)
}

read_optional <- function(path) {
  if (!file.exists(path)) return(data.table())
  fread(path, showProgress = FALSE)
}

num_col <- function(dt, col, default = NA_real_) {
  if (!col %chin% names(dt)) return(rep(default, nrow(dt)))
  x <- suppressWarnings(as.numeric(dt[[col]]))
  x[!is.finite(x)] <- default
  x
}

chr_col <- function(dt, col, default = "") {
  if (!col %chin% names(dt)) return(rep(default, nrow(dt)))
  x <- as.character(dt[[col]])
  x[is.na(x)] <- default
  x
}

bounded01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- 0
  pmin(pmax(x, 0), 1)
}

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

short_label <- function(x, width = 44) {
  vapply(as.character(x), function(s) {
    s <- if (is.na(s)) "" else s
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1))
}

collapse_unique <- function(x, n = 4) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  if (length(x) > n) {
    paste0(paste(x[seq_len(n)], collapse = "; "), "; +", length(x) - n,
           " more")
  } else {
    paste(x, collapse = "; ")
  }
}

signed_pct <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), sprintf("%+.1f%%", x), "unknown")
}

signed_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), sprintf("%+.3f", x), "unknown")
}

first_nonempty <- function(...) {
  values <- list(...)
  n <- length(values[[1]])
  out <- rep("", n)
  for (x in values) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    take <- !nzchar(out) & nzchar(x)
    out[take] <- x[take]
  }
  out
}

save_plot <- function(plot, name, width, height) {
  png_file <- file.path(figure_dir, paste0(name, ".png"))
  pdf_file <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220,
         bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  if (exists("dominant_save_ggplotly", mode = "function") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      plot,
      file.path(figure_dir, paste0(name, ".html")),
      title = name
    )
  }
}

keys <- c("gene_symbol", "tx_id", "design_family")

candidate <- fread(consensus_file, showProgress = FALSE)
review <- fread(review_file, showProgress = FALSE)
evidence <- fread(evidence_file, showProgress = FALSE)
branch_context <- fread(branch_context_file, showProgress = FALSE)
perturb <- read_optional(perturb_file)

required_candidate <- c(
  "consensus_rank", keys, "consensus_class", "consensus_score",
  "atlas_clean_cds_pct_change", "atlas_direction", "branch_component",
  "model_component", "robustness_component", "browser_component",
  "context_component", "pooling_component", "pooling_hazard_component",
  "aux_context_component", "therapeutic_component",
  "fragility_component", "warning_penalty", "consensus_review_question"
)
missing_candidate <- setdiff(required_candidate, names(candidate))
if (length(missing_candidate)) {
  stop("Consensus candidates are missing columns: ",
       paste(missing_candidate, collapse = ", "), call. = FALSE)
}

candidate[, is_viral_design := grepl("viral", design_family,
                                     ignore.case = TRUE)]
candidate[, is_isr_design := grepl("isr|stress", design_family,
                                   ignore.case = TRUE)]
candidate[, has_browser_support := num_col(.SD, "browser_component", 0) > 0]
candidate[, has_branch_support := num_col(.SD, "branch_component", 0) >= 0.5]
candidate[, has_robustness_support :=
            num_col(.SD, "robustness_component", 0) >= 0.5]
candidate[, has_context_warning :=
            num_col(.SD, "context_component", 0) >= 0.5 |
            num_col(.SD, "aux_context_component", 0) >= 0.5]
candidate[, has_exact_pooling_proof :=
            chr_col(.SD, "pooling_proof_tier") ==
            "A_exact_strict_project_matched_control"]
candidate[, has_pooling_protocol_hazard :=
            num_col(.SD, "pooling_hazard_component", 0) >= 0.40]
candidate[, has_start_site_protocol :=
            grepl("start_site_protocol_present|start_site_enriched",
                  paste(chr_col(.SD, "pooling_hazard_sources"),
                        chr_col(.SD, "pooling_top_protocol_classes"),
                        chr_col(.SD, "warnings")),
                  ignore.case = TRUE)]
candidate[, has_therapeutic_link :=
            num_col(.SD, "therapeutic_component", 0) >= 0.5]

candidate[, selected_for_primary_review :=
            consensus_rank <= 100 |
            consensus_class %chin% c(
              "A_browser_validated_consensus",
              "B_replicated_robust_consensus",
              "C_multi_model_branch_consensus",
              "D_context_specific_review",
              "E_buffering_compensation_review"
            ) |
            has_browser_support |
            (is_viral_design & consensus_score >= 0.35 & has_branch_support) |
            (has_pooling_protocol_hazard & consensus_score >= 0.35) |
            (has_therapeutic_link & consensus_score >= 0.38)]

candidate[, review_priority_tier := fcase(
  consensus_class == "A_browser_validated_consensus",
  "anchor_positive_control",
  consensus_class == "B_replicated_robust_consensus",
  "immediate_browser_validation",
  consensus_class == "C_multi_model_branch_consensus",
  "branch_mechanism_validation",
  consensus_class == "D_context_specific_review",
  "context_specific_validation",
  consensus_class == "E_buffering_compensation_review",
  "buffering_compensation_review",
  consensus_class == "F_therapeutic_translation_review" &
    selected_for_primary_review,
  "therapeutic_followup_review",
  has_pooling_protocol_hazard & selected_for_primary_review,
  "pooling_protocol_confounder_review",
  consensus_class == "I_fragile_or_confounded_review",
  "fragility_or_confounder_review",
  consensus_class == "H_low_support_background",
  "defer_low_support",
  default = "secondary_atlas_review"
)]

candidate[, review_priority_weight := fcase(
  review_priority_tier == "anchor_positive_control", 1.00,
  review_priority_tier == "immediate_browser_validation", 0.95,
  review_priority_tier == "branch_mechanism_validation", 0.86,
  review_priority_tier == "context_specific_validation", 0.80,
  review_priority_tier == "buffering_compensation_review", 0.74,
  review_priority_tier == "therapeutic_followup_review", 0.68,
  review_priority_tier == "pooling_protocol_confounder_review", 0.62,
  review_priority_tier == "secondary_atlas_review", 0.55,
  review_priority_tier == "fragility_or_confounder_review", 0.36,
  default = 0.20
)]
candidate[, validation_priority_score := bounded01(
  0.72 * consensus_score +
    0.12 * review_priority_weight +
    0.06 * as.numeric(has_browser_support) +
    0.05 * as.numeric(is_viral_design) +
    0.04 * as.numeric(has_exact_pooling_proof) +
    0.05 * as.numeric(has_therapeutic_link) -
    0.04 * num_col(.SD, "pooling_hazard_component", 0) -
    0.08 * num_col(.SD, "warning_penalty", 0)
)]

candidate[, main_branch_class := first_nonempty(
  chr_col(.SD, "joint_top_branch_class"),
  chr_col(.SD, "hierarchical_joint_top_branch_class"),
  chr_col(.SD, "dm_top_branch_class"),
  chr_col(.SD, "top_shifted_feature_type")
)]
candidate[, main_branch_delta := fifelse(
  nzchar(chr_col(.SD, "joint_top_branch_class")) &
    is.finite(num_col(.SD, "joint_top_branch_delta")),
  num_col(.SD, "joint_top_branch_delta"),
  fifelse(
    nzchar(chr_col(.SD, "hierarchical_joint_top_branch_class")) &
      is.finite(num_col(.SD, "hierarchical_joint_top_delta")),
    num_col(.SD, "hierarchical_joint_top_delta"),
    num_col(.SD, "dm_top_delta")
  )
)]
candidate[!is.finite(main_branch_delta), main_branch_delta := NA_real_]

candidate[, expected_cds_pattern := paste0(
  "Clean CDS ", ifelse(num_col(.SD, "atlas_clean_cds_pct_change", 0) >= 0,
                       "higher", "lower"),
  " in case by ", signed_pct(num_col(.SD, "atlas_clean_cds_pct_change")),
  " (", chr_col(.SD, "atlas_direction", "unknown_direction"), ")."
)]
candidate[, expected_branch_pattern := paste0(
  "Top branch: ",
  ifelse(nzchar(main_branch_class), main_branch_class, "unknown"),
  " delta ", signed_num(main_branch_delta),
  ifelse(nzchar(chr_col(.SD, "top_shifted_feature_id")),
         paste0("; feature ", chr_col(.SD, "top_shifted_feature_id"),
                " / ", chr_col(.SD, "top_shifted_feature_type")),
         ""),
  "."
)]
candidate[, expected_browser_pattern := paste(
  expected_cds_pattern,
  expected_branch_pattern,
  ifelse(nzchar(chr_col(.SD, "pooling_claim_language")),
         paste0("Claim tier: ", chr_col(.SD, "pooling_claim_language"),
                "; hazard: ",
                chr_col(.SD, "pooling_hazard_class", "unknown"), "."),
         ""),
  ifelse(nzchar(chr_col(.SD, "best_case_condition")),
         paste0("Open ", chr_col(.SD, "best_case_condition"),
                " versus ", chr_col(.SD, "best_control_conditions"),
                " in ", chr_col(.SD, "best_study"), "."),
         "")
)]

candidate[, warning_items := gsub(";\\s*", " | ", chr_col(.SD, "warnings"))]
candidate[, consensus_warning_items :=
            gsub(";\\s*", " | ", chr_col(.SD, "consensus_warning_summary"))]
candidate[, review_packet_id := paste0(
  "RDG", sprintf("%04d", consensus_rank), "_",
  safe_id(gene_symbol), "_", safe_id(design_family)
)]
candidate[, primary_review_rank := NA_integer_]
setorder(candidate, -selected_for_primary_review,
         -validation_priority_score, consensus_rank)
candidate[selected_for_primary_review == TRUE,
          primary_review_rank := seq_len(.N)]
setorder(candidate, consensus_rank)

packet_cols <- c(
  "review_packet_id", "primary_review_rank", "selected_for_primary_review",
  "review_priority_tier", "validation_priority_score",
  "consensus_rank", "gene_symbol", "tx_id", "design_family",
  "consensus_class", "consensus_score", "evidence_tags",
  "atlas_evidence_class", "atlas_clean_cds_pct_change",
  "atlas_clean_cds_pct_lower", "atlas_clean_cds_pct_upper",
  "atlas_direction", "n_studies", "strict_supported_studies",
  "total_case_samples", "total_control_samples", "best_case_condition",
  "best_control_conditions", "best_study", "main_branch_class",
  "main_branch_delta", "top_shifted_feature_id", "top_shifted_feature_type",
  "top_shifted_feature_family", "model_agreement_tier",
  "loo_consensus_class", "transport_consensus_class",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_proof_component", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "pooling_top_protocol_classes", "pooling_exact_context_fraction",
  "pooling_joint_strict_fraction",
  "context_interaction_gene_class", "context_interaction_summary",
  "partial_pooling_gene_class", "browser_validation_label",
  "browser_validation_context", "aux_top_state", "aux_top_overlap_class",
  "aux_audit_warning", "therapeutic_gene_label", "perturbation_label",
  "perturbation_support_class", "half_normalization_cds_pct_change",
  "validation_assay", "expected_browser_pattern",
  "consensus_review_question", "consensus_warning_items", "warning_items",
  "atlas_component", "branch_component", "model_component",
  "robustness_component", "browser_component", "output_effect_component",
  "context_component", "pooling_component", "pooling_hazard_component",
  "aux_context_component", "therapeutic_component",
  "fragility_component", "warning_penalty"
)
packet_cols <- intersect(packet_cols, names(candidate))
candidate_packet <- copy(candidate[, ..packet_cols])
numeric_packet <- names(candidate_packet)[
  vapply(candidate_packet, is.numeric, logical(1))
]
for (col in numeric_packet) {
  set(candidate_packet, j = col,
      value = round(candidate_packet[[col]], if (col %chin% c(
        "consensus_rank", "primary_review_rank", "n_studies",
        "strict_supported_studies", "total_case_samples",
        "total_control_samples"
      )) 0 else 4))
}

branch_required <- c(keys, "branch_context_rank", "branch_context_label")
missing_branch <- setdiff(branch_required, names(branch_context))
if (length(missing_branch)) {
  stop("Branch-context table is missing columns: ",
       paste(missing_branch, collapse = ", "), call. = FALSE)
}

primary_keys <- unique(candidate[
  selected_for_primary_review == TRUE,
  c(keys, "review_packet_id", "primary_review_rank",
    "validation_priority_score"),
  with = FALSE
])

context_packet <- merge(
  branch_context,
  primary_keys,
  by = keys,
  all = FALSE,
  sort = FALSE
)
context_packet <- context_packet[branch_context_rank <= 3]
context_packet[, context_count_gate := fcase(
  num_col(.SD, "case_branch_total_counts", 0) >= 30 &
    num_col(.SD, "control_branch_total_counts", 0) >= 30,
  "passes_30_count_gate_both_sides",
  num_col(.SD, "case_branch_total_counts", 0) >= 30 |
    num_col(.SD, "control_branch_total_counts", 0) >= 30,
  "passes_30_count_gate_one_side",
  default = "low_branch_counts"
)]
context_packet[, replicate_gate := fcase(
  num_col(.SD, "n_case_runs_merged", 0) >= 2 &
    num_col(.SD, "n_control_runs_merged", 0) >= 2,
  "merged_replicates_both_sides",
  num_col(.SD, "n_case_runs_merged", 0) >= 2 |
    num_col(.SD, "n_control_runs_merged", 0) >= 2,
  "merged_replicates_one_side",
  default = "single_or_unmerged_context"
)]
context_packet[, expected_context_pattern := paste0(
  branch_context_label,
  " | top branch ", chr_col(.SD, "joint_top_branch_class", "unknown"),
  " delta ", signed_num(num_col(.SD, "joint_top_branch_delta")),
  " | clean CDS delta ", signed_num(num_col(.SD, "joint_delta_clean_CDS")),
  " | counts case/control ",
  round(num_col(.SD, "case_branch_total_counts", 0)), "/",
  round(num_col(.SD, "control_branch_total_counts", 0))
)]
context_cols <- c(
  "review_packet_id", "primary_review_rank", keys,
  "branch_context_rank", "branch_context_label", "study", "AUTHOR",
  "case_condition", "control_conditions", "TISSUE", "CELL_LINE", "GENE",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
  "match_scope", "n_case_runs_merged", "n_control_runs_merged",
  "case_branch_total_counts", "control_branch_total_counts",
  "context_count_gate", "replicate_gate", "joint_review_class",
  "joint_review_priority", "joint_top_branch_class",
  "joint_top_branch_delta", "joint_top_branch_q",
  "joint_l1_allocation_shift", "joint_delta_clean_CDS",
  "joint_delta_leader_uORF", "joint_delta_overlapping_uORF",
  "joint_case_prop_clean_CDS", "joint_control_prop_clean_CDS",
  "joint_case_prop_leader_uORF", "joint_control_prop_leader_uORF",
  "joint_case_prop_overlapping_uORF",
  "joint_control_prop_overlapping_uORF",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_hazard_class", "pooling_hazard_score",
  "pooling_hazard_sources", "pooling_top_protocol_classes",
  "context_interaction_gene_class", "context_interaction_summary",
  "partial_pooling_gene_class", "partial_pooling_top_class",
  "warnings", "expected_context_pattern"
)
context_cols <- intersect(context_cols, names(context_packet))
context_packet <- context_packet[, ..context_cols]
for (col in names(context_packet)[vapply(context_packet, is.numeric,
                                         logical(1))]) {
  set(context_packet, j = col, value = round(context_packet[[col]], 4))
}
setorder(context_packet, primary_review_rank, branch_context_rank)

top_context <- context_packet[
  branch_context_rank == 1,
  .(review_packet_id, branch_context_label, expected_context_pattern,
    context_count_gate, replicate_gate)
]
template <- merge(
  candidate[
    selected_for_primary_review == TRUE,
    .(
      review_packet_id, primary_review_rank, consensus_rank,
      gene_symbol, tx_id, design_family, review_priority_tier,
      validation_priority_score, consensus_class, consensus_score,
      pooling_proof_tier, pooling_claim_language, pooling_hazard_class,
      pooling_hazard_score, pooling_hazard_sources,
      expected_browser_pattern, consensus_review_question,
      consensus_warning_items, warning_items
    )
  ],
  top_context,
  by = "review_packet_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(template, primary_review_rank)
template[, validation_status := ""]
template[, browser_signal_grade := ""]
template[, branch_direction_verified := ""]
template[, clean_cds_direction_verified := ""]
template[, context_match_verified := ""]
template[, likely_failure_mode := ""]
template[, reviewer_notes := ""]
template[, reviewer := ""]
template[, review_date := ""]
template[, next_action := ""]
template_cols <- c(
  "review_packet_id", "primary_review_rank", "consensus_rank",
  "gene_symbol", "tx_id", "design_family", "review_priority_tier",
  "validation_priority_score", "consensus_class", "consensus_score",
  "pooling_proof_tier", "pooling_claim_language", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "branch_context_label", "context_count_gate", "replicate_gate",
  "expected_browser_pattern", "expected_context_pattern",
  "consensus_review_question", "consensus_warning_items", "warning_items",
  "validation_status", "browser_signal_grade",
  "branch_direction_verified", "clean_cds_direction_verified",
  "context_match_verified", "likely_failure_mode", "reviewer_notes",
  "reviewer", "review_date", "next_action"
)
template <- template[, ..template_cols]
for (col in names(template)[vapply(template, is.numeric, logical(1))]) {
  set(template, j = col, value = round(template[[col]], 4))
}

if (nrow(perturb) && all(c(keys, "therapeutic_hypothesis") %chin%
                         names(perturb))) {
  perturb_summary <- perturb[
    ,
    .(
      perturbation_hypotheses = collapse_unique(therapeutic_hypothesis, 3),
      perturbation_labels = collapse_unique(hypothesis_label, 3),
      max_half_normalization_abs_pct = {
        values <- abs(num_col(.SD, "half_normalization_cds_pct_change"))
        if (all(!is.finite(values))) NA_real_ else max(values, na.rm = TRUE)
      }
    ),
    by = gene_symbol
  ]
  perturb_summary[!is.finite(max_half_normalization_abs_pct),
                  max_half_normalization_abs_pct := NA_real_]
} else {
  perturb_summary <- data.table(
    gene_symbol = character(),
    perturbation_hypotheses = character(),
    perturbation_labels = character(),
    max_half_normalization_abs_pct = numeric()
  )
}

gene_summary <- candidate[
  ,
  .(
    best_consensus_rank = min(consensus_rank),
    best_primary_review_rank = {
      rank_values <- as.numeric(primary_review_rank)
      if (all(is.na(rank_values))) NA_real_ else min(rank_values, na.rm = TRUE)
    },
    best_consensus_score = max(consensus_score, na.rm = TRUE),
    best_validation_priority_score =
      max(validation_priority_score, na.rm = TRUE),
    n_designs = .N,
    n_primary_review_designs = sum(selected_for_primary_review, na.rm = TRUE),
    n_viral_designs = sum(is_viral_design, na.rm = TRUE),
    n_exact_pooling_proof_designs =
      sum(has_exact_pooling_proof, na.rm = TRUE),
    n_pooling_protocol_hazard_designs =
      sum(has_pooling_protocol_hazard, na.rm = TRUE),
    has_browser_support = any(has_browser_support, na.rm = TRUE),
    has_immediate_candidate = any(
      review_priority_tier %chin% c(
        "anchor_positive_control",
        "immediate_browser_validation",
        "branch_mechanism_validation",
        "context_specific_validation"
      ),
      na.rm = TRUE
    ),
    best_design_family = design_family[which.max(validation_priority_score)],
    best_consensus_class = consensus_class[which.max(validation_priority_score)],
    best_review_priority_tier =
      review_priority_tier[which.max(validation_priority_score)],
    best_expected_pattern =
      expected_browser_pattern[which.max(validation_priority_score)]
  ),
  by = gene_symbol
]
gene_summary[!is.finite(best_primary_review_rank),
             best_primary_review_rank := NA_real_]
gene_summary <- merge(gene_summary, perturb_summary, by = "gene_symbol",
                      all.x = TRUE, sort = FALSE)
gene_summary[, gene_review_class := fcase(
  has_browser_support, "browser_supported_gene",
  has_immediate_candidate & n_viral_designs > 0, "viral_candidate_gene",
  has_immediate_candidate, "immediate_review_gene",
  n_primary_review_designs > 0, "secondary_review_gene",
  default = "background_gene"
)]
for (col in names(gene_summary)[vapply(gene_summary, is.numeric, logical(1))]) {
  set(gene_summary, j = col, value = round(gene_summary[[col]], 4))
}
setorder(gene_summary, best_consensus_rank)

tier_counts <- candidate_packet[
  ,
  .(
    n_rows = .N,
    n_primary = sum(selected_for_primary_review, na.rm = TRUE)
  ),
  by = review_priority_tier
]
metrics <- rbindlist(list(
  data.table(metric = "candidate_rows", value = nrow(candidate_packet)),
  data.table(metric = "primary_review_rows",
             value = sum(candidate$selected_for_primary_review, na.rm = TRUE)),
  data.table(metric = "manual_template_rows", value = nrow(template)),
  data.table(metric = "branch_context_rows", value = nrow(context_packet)),
  data.table(metric = "genes_in_dossier", value = uniqueN(candidate$gene_symbol)),
  data.table(metric = "primary_review_genes",
             value = uniqueN(candidate[
               selected_for_primary_review == TRUE,
               gene_symbol
             ])),
  data.table(metric = "viral_primary_review_rows",
             value = sum(candidate$selected_for_primary_review &
                           candidate$is_viral_design, na.rm = TRUE)),
  data.table(metric = "browser_supported_primary_rows",
             value = sum(candidate$selected_for_primary_review &
                           candidate$has_browser_support, na.rm = TRUE)),
  data.table(metric = "pooling_exact_primary_rows",
             value = sum(candidate$selected_for_primary_review &
                           candidate$has_exact_pooling_proof, na.rm = TRUE)),
  data.table(metric = "pooling_protocol_hazard_primary_rows",
             value = sum(candidate$selected_for_primary_review &
                           candidate$has_pooling_protocol_hazard,
                         na.rm = TRUE)),
  data.table(metric = "contexts_pass_30_count_gate_both_sides",
             value = sum(context_packet$context_count_gate ==
                           "passes_30_count_gate_both_sides",
                         na.rm = TRUE)),
  data.table(metric = "contexts_merged_replicates_both_sides",
             value = sum(context_packet$replicate_gate ==
                           "merged_replicates_both_sides",
                         na.rm = TRUE)),
  tier_counts[, .(
    metric = paste0("tier_", safe_id(review_priority_tier)),
    value = n_rows
  )]
), fill = TRUE)
metrics[, value := as.numeric(value)]

fwrite(candidate_packet, candidate_out)
fwrite(context_packet, context_out)
fwrite(template, template_out)
fwrite(gene_summary, gene_out)
fwrite(metrics, metrics_out)

primary_plot_dt <- candidate[
  selected_for_primary_review == TRUE
][order(primary_review_rank)][seq_len(min(.N, 60))]
primary_plot_dt[, plot_label := factor(
  paste0(primary_review_rank, ". ", gene_symbol, " / ", design_family),
  levels = rev(paste0(primary_review_rank, ". ", gene_symbol, " / ",
                      design_family))
)]
  tier_palette <- c(
  anchor_positive_control = "#D95F02",
  immediate_browser_validation = "#2B6CB0",
  branch_mechanism_validation = "#1B9E77",
  context_specific_validation = "#7570B3",
  buffering_compensation_review = "#E7298A",
  therapeutic_followup_review = "#66A61E",
  pooling_protocol_confounder_review = "#A23E48",
  secondary_atlas_review = "#A6761D",
  fragility_or_confounder_review = "#666666",
  defer_low_support = "#BDBDBD"
)
present_tiers <- unique(as.character(primary_plot_dt$review_priority_tier))
palette_values <- tier_palette[present_tiers]
missing_tiers <- present_tiers[is.na(palette_values)]
if (length(missing_tiers)) {
  fallback <- grDevices::hcl.colors(length(missing_tiers), "Dark 3")
  names(fallback) <- missing_tiers
  palette_values <- c(palette_values[!is.na(palette_values)], fallback)
}
p_priority <- ggplot(primary_plot_dt,
                     aes(validation_priority_score, plot_label,
                         fill = review_priority_tier)) +
  geom_col(width = 0.78, color = "grey20", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%.2f", validation_priority_score)),
            hjust = -0.1, size = 2.6) +
  scale_x_continuous(limits = c(0, max(primary_plot_dt$validation_priority_score,
                                       na.rm = TRUE) + 0.08),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "RDG Validation Dossier Priority",
    subtitle = "Rows selected from the consensus atlas for direct browser review",
    x = "Validation priority score",
    y = NULL,
    fill = "Review tier"
  ) +
  scale_fill_manual(
    values = palette_values,
    labels = function(x) gsub("_", " ", x, fixed = TRUE)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title.position = "plot"
  )
save_plot(p_priority, "rdg_validation_dossier_priority", 10, 12.3)

workload <- candidate[
  ,
  .N,
  by = .(selected_for_primary_review, review_priority_tier)
]
workload[, selected_for_primary_review := ifelse(selected_for_primary_review,
                                                "Primary review",
                                                "Not primary")]
p_workload <- ggplot(workload,
                     aes(review_priority_tier, N,
                         fill = selected_for_primary_review)) +
  geom_col(position = "stack", width = 0.72, color = "grey25",
           linewidth = 0.15) +
  coord_flip() +
  labs(
    title = "Validation Workload by Review Tier",
    subtitle = "Primary rows are intended for immediate browser inspection",
    x = NULL,
    y = "Rows",
    fill = "Scope"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot")
save_plot(p_workload, "rdg_validation_dossier_workload", 8, 5.6)

evidence_primary <- evidence[
  candidate[selected_for_primary_review == TRUE],
  on = keys,
  nomatch = 0
]
evidence_primary <- evidence_primary[
  consensus_rank <= max(primary_plot_dt$consensus_rank, na.rm = TRUE) |
    review_packet_id %chin% primary_plot_dt$review_packet_id
]
evidence_primary <- evidence_primary[
  review_packet_id %chin% primary_plot_dt$review_packet_id
]
rank_col <- if ("i.primary_review_rank" %chin% names(evidence_primary)) {
  "i.primary_review_rank"
} else {
  "primary_review_rank"
}
evidence_primary[, candidate_label := factor(
  paste0(get(rank_col), ". ", gene_symbol, " / ", design_family),
  levels = rev(paste0(primary_plot_dt$primary_review_rank, ". ",
                      primary_plot_dt$gene_symbol, " / ",
                      primary_plot_dt$design_family))
)]
evidence_primary[, component_label := factor(
  component_label,
  levels = unique(component_label)
)]
p_evidence <- ggplot(evidence_primary,
                     aes(component_label, candidate_label,
                         fill = component_score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2B6CB0", mid = "white",
                       high = "#991B1B", midpoint = 0.5,
                       limits = c(0, 1), oob = scales::squish) +
  labs(
    title = "Evidence Components for Validation Packets",
    subtitle = "Blue = weak component, white = intermediate, red = strong",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.title.position = "plot"
  )
save_plot(p_evidence, "rdg_validation_dossier_evidence_heatmap", 12, 12)

context_plot <- context_packet[
  review_packet_id %chin% primary_plot_dt$review_packet_id
]
if (nrow(context_plot)) {
  context_plot[, candidate_label := factor(
    paste0(primary_review_rank, ". ", gene_symbol, " / ", design_family),
    levels = rev(paste0(primary_plot_dt$primary_review_rank, ". ",
                        primary_plot_dt$gene_symbol, " / ",
                        primary_plot_dt$design_family))
  )]
  context_plot[, context_axis := paste0("context ", branch_context_rank)]
  p_context <- ggplot(context_plot,
                      aes(context_axis, candidate_label,
                          fill = joint_top_branch_delta)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = joint_top_branch_class), size = 2.4) +
    scale_fill_gradient2(low = "#2B6CB0", mid = "white",
                         high = "#991B1B", midpoint = 0,
                         na.value = "grey90") +
    labs(
      title = "Top Matched Branch Contexts",
      subtitle = "Fill is top-branch delta in case minus control",
      x = NULL,
      y = NULL,
      fill = "Delta"
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot")
  save_plot(p_context, "rdg_validation_dossier_branch_contexts", 8.5, 12)
}

md_table <- function(dt, cols, n = 15) {
  if (!nrow(dt)) return(character())
  dt <- copy(dt[seq_len(min(.N, n)), ..cols])
  for (col in names(dt)) dt[[col]] <- gsub("\\|", "/", as.character(dt[[col]]))
  header <- paste0("| ", paste(names(dt), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(dt)), collapse = " | "), " |")
  rows <- apply(dt, 1, function(x) {
    paste0("| ", paste(x, collapse = " | "), " |")
  })
  c(header, sep, rows)
}

top_report <- candidate_packet[
  selected_for_primary_review == TRUE
][order(primary_review_rank)][
  ,
  .(
    rank = primary_review_rank,
    gene = gene_symbol,
    design = design_family,
    tier = review_priority_tier,
    score = sprintf("%.2f", validation_priority_score),
    consensus = sprintf("%.2f", consensus_score),
    expected = expected_browser_pattern
  )
]
top_report[, expected := short_label(expected, 70)]

report_lines <- c(
  "# RDG Validation Dossiers",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Purpose",
  "",
  paste(
    "This folder converts the RDG consensus atlas into concrete browser",
    "review packets. The score is a practical review priority, not a new",
    "statistical significance value. A row is high priority when multiple",
    "layers agree, when the effect is robust across omission tests, when",
    "the context is biologically important, or when it is useful as a",
    "perturbation/positive-control candidate. It also carries the pooling",
    "proof tier and pooling/protocol hazard class so reviewers know whether",
    "a row is an exact local claim candidate or a confounded review lead."
  ),
  "",
  "## Summary Metrics",
  "",
  md_table(metrics, c("metric", "value"), n = nrow(metrics)),
  "",
  "## Highest Priority Browser Packets",
  "",
  md_table(top_report, names(top_report), n = 25),
  "",
  "## Review Protocol",
  "",
  "1. Open the first context listed for each packet.",
  "2. Check whether clean-CDS direction matches the expected pattern.",
  "3. Check whether the listed branch shift is visible in the RDG/coverage view.",
  "4. Check the pooling proof tier and pooling hazard class before wording the claim.",
  "5. Treat auxiliary-state and pooling/protocol warnings as possible context drivers, not as disqualifiers.",
  "6. Mark the manual template as positive, weak positive, negative, confounded, or not reviewable.",
  "7. Feed validated positives and negatives back into motif calibration.",
  "",
  "## Practical Interpretation",
  "",
  paste(
    "DDIT3 remains the anchor-positive control because browser support and",
    "multiple model layers agree. PPP1R15A and ATF4 remain the strongest",
    "viral branch candidates. IFIH1 remains a lung/context-specific viral",
    "hypothesis and should be reviewed by comparing lung and non-lung",
    "infected contexts directly."
  )
)
writeLines(report_lines, report_out, useBytes = TRUE)

cat("RDG validation dossiers:\n")
cat("  candidate rows: ", nrow(candidate_packet), "\n", sep = "")
cat("  primary review rows: ",
    sum(candidate$selected_for_primary_review, na.rm = TRUE), "\n", sep = "")
cat("  branch context rows: ", nrow(context_packet), "\n", sep = "")
cat("  genes: ", uniqueN(candidate$gene_symbol), "\n", sep = "")
cat("  outputs: ", output_dir, "\n", sep = "")

stopifnot(nrow(candidate_packet) == nrow(candidate))
stopifnot(nrow(template) == sum(candidate$selected_for_primary_review,
                                na.rm = TRUE))
stopifnot(all(candidate_packet$validation_priority_score >= 0 &
                candidate_packet$validation_priority_score <= 1))
stopifnot(all(template$primary_review_rank == seq_len(nrow(template))))
