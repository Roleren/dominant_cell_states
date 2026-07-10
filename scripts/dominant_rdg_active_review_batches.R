#!/usr/bin/env Rscript

# Build concrete active-review batches from the RDG evidence-gap queue.
# This is a human-facing layer: small balanced batches, not another discovery
# model. Each row gets one review job and blank columns for manual feedback.

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
output_dir <- file.path(analysis_dir, "dominant_rdg_review_batches")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

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

make_row_key <- function(dt) {
  paste(chr_col(dt, "gene_symbol"),
        chr_col(dt, "tx_id"),
        chr_col(dt, "design_family"),
        sep = "||")
}

compact_text <- function(x, max_chars = 220) {
  x <- gsub("\\s+", " ", as.character(x))
  x[is.na(x)] <- ""
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 3), "...")
  x
}

read_required <- function(path, label) {
  dt <- safe_fread(path)
  if (nrow(dt) == 0) {
    stop("Missing or empty ", label, ": ", path, call. = FALSE)
  }
  dt
}

priority_file <- file.path(
  analysis_dir,
  "dominant_rdg_evidence_gaps",
  "dominant_rdg_evidence_gap_priorities.csv"
)
negative_file <- file.path(
  analysis_dir,
  "dominant_rdg_evidence_gaps",
  "dominant_rdg_evidence_gap_negative_controls.csv"
)
candidates_file <- file.path(
  analysis_dir,
  "dominant_rdg_validation_dossiers",
  "dominant_rdg_validation_dossier_candidates.csv"
)
contexts_file <- file.path(
  analysis_dir,
  "dominant_rdg_validation_dossiers",
  "dominant_rdg_validation_dossier_branch_contexts.csv"
)

batches_out <- file.path(output_dir, "dominant_rdg_review_batches.csv")
template_out <- file.path(output_dir, "dominant_rdg_review_batch_template.csv")
summary_out <- file.path(output_dir, "dominant_rdg_review_batch_summary.csv")
gene_balance_out <- file.path(output_dir, "dominant_rdg_review_batch_gene_balance.csv")
metrics_out <- file.path(output_dir, "dominant_rdg_review_batch_metrics.csv")
report_out <- file.path(output_dir, "dominant_rdg_review_batch_report.md")

priority <- read_required(priority_file, "evidence-gap priorities")
negative <- read_required(negative_file, "evidence-gap negative controls")
candidates <- read_required(candidates_file, "validation dossier candidates")
contexts <- safe_fread(contexts_file)

keys3 <- c("gene_symbol", "tx_id", "design_family")
if (!all(keys3 %in% names(priority))) {
  stop("Priority table lacks required columns: ", paste(keys3, collapse = ", "),
       call. = FALSE)
}

priority[, row_key := make_row_key(.SD)]
negative[, row_key := make_row_key(.SD)]

candidate_keep <- intersect(
  c("review_packet_id", keys3, "primary_review_rank",
    "expected_browser_pattern", "consensus_review_question",
    "context_interaction_summary", "aux_audit_warning",
    "half_normalization_cds_pct_change", "validation_assay",
    "evidence_tags", "consensus_warning_items",
    "pooling_proof_tier", "pooling_claim_language",
    "pooling_hazard_class", "pooling_hazard_score",
    "pooling_hazard_sources", "pooling_top_protocol_classes"),
  names(candidates)
)
candidate_info <- unique(candidates[, candidate_keep, with = FALSE])
if ("review_packet_id" %in% names(priority) &&
    "review_packet_id" %in% names(candidate_info)) {
  add_cols <- setdiff(names(candidate_info), names(priority))
  if (length(add_cols)) {
    priority <- merge(
      priority,
      candidate_info[, c("review_packet_id", add_cols), with = FALSE],
      by = "review_packet_id",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

context_info <- data.table()
if (nrow(contexts) > 0 && "review_packet_id" %in% names(contexts)) {
  setorder(contexts, review_packet_id, branch_context_rank)
  context_keep <- intersect(
    c("review_packet_id", "branch_context_label", "study", "AUTHOR",
      "case_condition", "control_conditions", "TISSUE", "CELL_LINE", "GENE",
      "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
      "match_scope", "n_case_runs_merged", "n_control_runs_merged",
      "case_branch_total_counts", "control_branch_total_counts",
      "context_count_gate", "replicate_gate", "joint_review_class",
      "joint_review_priority", "joint_top_branch_class",
      "joint_top_branch_delta", "joint_l1_allocation_shift",
      "joint_delta_clean_CDS", "joint_delta_leader_uORF",
      "joint_delta_overlapping_uORF", "expected_context_pattern",
      "pooling_proof_tier", "pooling_claim_language",
      "pooling_hazard_class", "pooling_hazard_score",
      "pooling_hazard_sources", "pooling_top_protocol_classes"),
    names(contexts)
  )
  context_info <- contexts[, context_keep, with = FALSE][
    , .SD[1], by = review_packet_id
  ]
}

take_balanced <- function(pool, n, used = character(), max_per_gene = 2,
                          max_per_design = 8) {
  if (nrow(pool) == 0 || n <= 0) return(pool[0])
  pool <- copy(pool)
  if (!"row_key" %in% names(pool)) pool[, row_key := make_row_key(.SD)]
  pool <- pool[!row_key %chin% used]
  if (nrow(pool) == 0) return(pool[0])
  setorder(pool, -acquisition_priority_score,
           -counterfactual_readiness_score, gene_symbol, design_family)
  selected <- integer()
  gene_counts <- list()
  design_counts <- list()
  get_count <- function(x, nm) {
    val <- x[[nm]]
    if (is.null(val)) 0L else val
  }
  for (i in seq_len(nrow(pool))) {
    gene <- as.character(pool$gene_symbol[i])
    design <- as.character(pool$design_family[i])
    if (get_count(gene_counts, gene) >= max_per_gene) next
    if (get_count(design_counts, design) >= max_per_design) next
    selected <- c(selected, i)
    gene_counts[[gene]] <- get_count(gene_counts, gene) + 1L
    design_counts[[design]] <- get_count(design_counts, design) + 1L
    if (length(selected) >= n) break
  }
  pool[selected]
}

batches <- list()
used_keys <- character()

add_batch <- function(pool, batch_id, batch_name, batch_goal, batch_role,
                      n, max_per_gene = 2, max_per_design = 8,
                      reuse_rows = FALSE) {
  selected <- take_balanced(
    pool = pool,
    n = n,
    used = if (reuse_rows) character() else used_keys,
    max_per_gene = max_per_gene,
    max_per_design = max_per_design
  )
  if (nrow(selected) == 0) return(invisible(NULL))
  selected[, `:=`(
    batch_id = batch_id,
    batch_name = batch_name,
    batch_goal = batch_goal,
    batch_role = batch_role,
    batch_review_order = seq_len(.N)
  )]
  batches[[length(batches) + 1L]] <<- selected
  if (!reuse_rows) used_keys <<- unique(c(used_keys, selected$row_key))
  invisible(NULL)
}

viral_or_stress <- grepl(
  "viral|interferon|isr|stress|nutrient|damage",
  chr_col(priority, "design_family"),
  ignore.case = TRUE
)

add_batch(
  priority[action_label == "browser_validate_now" & viral_or_stress],
  "B01",
  "Browser validation: viral/stress anchors",
  paste(
    "Validate the strongest acquisition rows in the post-viral/ISR/stress",
    "space before adding more model complexity."
  ),
  "candidate_positive",
  n = 22,
  max_per_gene = 2,
  max_per_design = 14
)

add_batch(
  priority[
    claim_stage == "near_counterfactual_candidate" |
      action_label == "perturbation_followup" |
      (counterfactual_readiness_score >= 65 & acquisition_priority_score >= 45)
  ],
  "B02",
  "Counterfactual readiness and perturbation follow-up",
  paste(
    "Rows closest to a protein-output prediction; review whether they are",
    "ready for perturbation framing or still need context qualification."
  ),
  "counterfactual_review",
  n = 18,
  max_per_gene = 3,
  max_per_design = 10,
  reuse_rows = TRUE
)

add_batch(
  priority[
    action_label == "protocol_pooling_audit" |
      protocol_pooling_need >= 0.55 |
      num_col(priority, "pooling_hazard_score") >= 0.40 |
      grepl("pooling/protocol|protocol", review_reason, ignore.case = TRUE)
  ],
  "B03",
  "Pooling and protocol confounder audit",
  paste(
    "Rows where inhibitor, fraction, relaxed pooling, start-site enrichment,",
    "or mixed protocol context could explain the RDG signal."
  ),
  "protocol_pooling_audit",
  n = 20,
  max_per_gene = 2,
  max_per_design = 8
)

add_batch(
  priority[action_label == "metadata_context_audit"],
  "B04",
  "Context and metadata confounder audit",
  paste(
    "High-value rows where state, tissue, cell line, tumor status, or immune",
    "composition could explain the RDG signal."
  ),
  "context_audit",
  n = 28,
  max_per_gene = 2,
  max_per_design = 8
)

measurement_pool <- priority[
  action_label == "add_or_verify_translon" |
    measurement_need >= 0.58 |
    grepl("translon|measurement|hidden", review_reason, ignore.case = TRUE)
]
add_batch(
  measurement_pool,
  "B05",
  "Translon and measurement review",
  paste(
    "Rows where wrong/missing uORF, NTE/NTT, iORF, isoform, or low-count",
    "measurement could dominate interpretation."
  ),
  "measurement_review",
  n = 24,
  max_per_gene = 3,
  max_per_design = 10
)

add_batch(
  priority[action_label == "add_context_replication"],
  "B06",
  "Replication and transport follow-up",
  paste(
    "Rows that are interesting enough to track but not yet robust across",
    "study/context omission."
  ),
  "replication_followup",
  n = 28,
  max_per_gene = 2,
  max_per_design = 8
)

negative[, `:=`(
  action_label = "negative_control_review",
  claim_stage = "falsification_control",
  review_reason = "low-effect row with adequate measurement support; use as a falsification control"
)]
add_batch(
  negative,
  "B07",
  "Negative controls and falsification",
  paste(
    "Low-effect but measurable rows matched across design families; these",
    "should not show strong browser branch changes."
  ),
  "negative_control",
  n = 36,
  max_per_gene = 2,
  max_per_design = 5,
  reuse_rows = TRUE
)

if (!length(batches)) {
  stop("No review-batch rows selected from evidence-gap priorities.",
       call. = FALSE)
}

batch_dt <- rbindlist(batches, fill = TRUE)
setorder(batch_dt, batch_id, batch_review_order)
batch_dt[, global_review_order := seq_len(.N)]

if (nrow(context_info) > 0 && "review_packet_id" %in% names(batch_dt)) {
  add_cols <- setdiff(names(context_info), names(batch_dt))
  batch_dt <- merge(
    batch_dt,
    context_info[, c("review_packet_id", add_cols), with = FALSE],
    by = "review_packet_id",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(batch_dt, batch_id, batch_review_order)
}

batch_dt[, expected_browser_pattern := compact_text(
  chr_col(.SD, "expected_browser_pattern"), 260
)]
batch_dt[, expected_context_pattern := compact_text(
  chr_col(.SD, "expected_context_pattern"), 260
)]
batch_dt[, review_question := compact_text(
  fifelse(has_text(chr_col(.SD, "consensus_review_question")),
          chr_col(.SD, "consensus_review_question"),
          chr_col(.SD, "review_reason")),
  260
)]
batch_dt[, review_instruction := fifelse(
  batch_role == "negative_control",
  "Check that the branch/CDS pattern is weak or absent; flag if a strong signal appears.",
  fifelse(
    batch_role == "measurement_review",
    "Inspect isoform and translon structure before interpreting branch direction.",
    fifelse(
      batch_role == "protocol_pooling_audit",
      "Check exact matched controls, inhibitor/fraction, relaxed pooling, and whether the visible branch signal survives local same-protocol comparisons.",
      fifelse(
      batch_role == "context_audit",
      "Inspect whether the signal follows condition or a confounded context/state.",
      fifelse(
        batch_role == "replication_followup",
        "Look for the same direction across comparable studies or contexts.",
        "Validate expected clean-CDS and branch direction in RiboCrypt."
      )
      )
    )
  )
)]

output_cols <- intersect(c(
  "global_review_order", "batch_id", "batch_name", "batch_review_order",
  "batch_role", "action_label", "claim_stage", "gene_symbol", "tx_id",
  "design_family", "acquisition_rank", "acquisition_priority_score",
  "counterfactual_readiness_score", "evidence_gap_score",
  "biological_value_score", "negative_control_score", "review_reason",
  "protocol_pooling_need", "pooling_proof_tier", "pooling_claim_language",
  "pooling_hazard_class", "pooling_hazard_score",
  "pooling_hazard_sources", "pooling_top_protocol_classes",
  "pooling_exact_context_fraction", "pooling_joint_strict_fraction",
  "pooling_any_start_site_enriched_protocol",
  "pooling_any_relaxed_context", "review_instruction",
  "review_packet_id", "primary_review_rank",
  "review_priority_tier", "consensus_rank", "consensus_class",
  "consensus_score", "atlas_evidence_class", "atlas_clean_cds_pct_change",
  "main_branch_class", "main_branch_delta", "top_shifted_feature_id",
  "top_shifted_feature_type", "top_shifted_feature_family",
  "model_agreement_tier", "loo_consensus_class",
  "transport_consensus_class", "context_interaction_gene_class",
  "partial_pooling_gene_class", "aux_top_state", "aux_top_overlap_class",
  "therapeutic_gene_label", "perturbation_label",
  "perturbation_support_class", "half_normalization_cds_pct_change",
  "n_studies", "strict_supported_studies", "total_case_samples",
  "total_control_samples", "best_case_condition", "best_control_conditions",
  "best_study", "branch_context_label", "study", "AUTHOR",
  "case_condition", "control_conditions", "TISSUE", "CELL_LINE", "GENE",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
  "match_scope", "n_case_runs_merged", "n_control_runs_merged",
  "case_branch_total_counts", "control_branch_total_counts",
  "context_count_gate", "replicate_gate", "joint_review_class",
  "joint_review_priority", "joint_top_branch_class",
  "joint_top_branch_delta", "joint_l1_allocation_shift",
  "joint_delta_clean_CDS", "joint_delta_leader_uORF",
  "joint_delta_overlapping_uORF", "transcript_qc_class",
  "cds_exon_qc_status", "n_measured_features",
  "n_signal_supported_features", "n_low_count_features",
  "n_manual_features", "n_overlapping_features", "n_internal_features",
  "expected_browser_pattern", "expected_context_pattern",
  "review_question", "warning_item_count", "warning_items", "warnings",
  "batch_goal"
), names(batch_dt))

batch_out <- batch_dt[, output_cols, with = FALSE]
fwrite(batch_out, batches_out)

template <- copy(batch_out)
template[, `:=`(
  manual_review_status = "",
  browser_signal_grade = "",
  branch_direction_verified = "",
  clean_cds_direction_verified = "",
  context_match_verified = "",
  translon_structure_status = "",
  metadata_confounder_status = "",
  negative_control_result = "",
  reviewer_notes = "",
  reviewer = "",
  review_date = "",
  next_action = ""
)]
template_cols <- c(
  "global_review_order", "batch_id", "batch_name", "batch_review_order",
  "batch_role", "gene_symbol", "tx_id", "design_family",
  "review_instruction", "expected_browser_pattern",
  "expected_context_pattern", "review_question",
  "manual_review_status", "browser_signal_grade",
  "branch_direction_verified", "clean_cds_direction_verified",
  "context_match_verified", "translon_structure_status",
  "metadata_confounder_status", "negative_control_result",
  "reviewer_notes", "reviewer", "review_date", "next_action",
  setdiff(names(template), c(
    "global_review_order", "batch_id", "batch_name", "batch_review_order",
    "batch_role", "gene_symbol", "tx_id", "design_family",
    "review_instruction", "expected_browser_pattern",
    "expected_context_pattern", "review_question",
    "manual_review_status", "browser_signal_grade",
    "branch_direction_verified", "clean_cds_direction_verified",
    "context_match_verified", "translon_structure_status",
    "metadata_confounder_status", "negative_control_result",
    "reviewer_notes", "reviewer", "review_date", "next_action"
  ))
)
setcolorder(template, intersect(template_cols, names(template)))
fwrite(template, template_out)

batch_summary <- batch_out[, .(
  n_rows = .N,
  n_genes = uniqueN(gene_symbol),
  n_design_families = uniqueN(design_family),
  median_acquisition_priority = round(
    median(num_col(.SD, "acquisition_priority_score"), na.rm = TRUE), 2
  ),
  median_counterfactual_readiness = round(
    median(num_col(.SD, "counterfactual_readiness_score"), na.rm = TRUE), 2
  ),
  top_gene_designs = paste(
    paste(head(gene_symbol, 8), head(design_family, 8), sep = " | "),
    collapse = "; "
  ),
  batch_goal = batch_goal[1]
), by = .(batch_id, batch_name, batch_role)]
setorder(batch_summary, batch_id)
fwrite(batch_summary, summary_out)

gene_balance <- batch_out[, .(
  n_rows = .N,
  batches = paste(unique(batch_id), collapse = "; "),
  roles = paste(unique(batch_role), collapse = "; "),
  max_acquisition_priority = round(
    max(num_col(.SD, "acquisition_priority_score"), na.rm = TRUE), 2
  ),
  max_counterfactual_readiness = round(
    max(num_col(.SD, "counterfactual_readiness_score"), na.rm = TRUE), 2
  ),
  design_families = paste(unique(design_family), collapse = "; ")
), by = gene_symbol]
setorder(gene_balance, -max_acquisition_priority, gene_symbol)
fwrite(gene_balance, gene_balance_out)

metrics <- data.table(
  metric = c(
    "review_batch_rows", "review_batch_genes", "review_batches",
    "candidate_positive_rows", "counterfactual_review_rows",
    "protocol_pooling_audit_rows", "context_audit_rows",
    "measurement_review_rows",
    "replication_followup_rows", "negative_control_rows",
    "median_acquisition_priority",
    "median_counterfactual_readiness"
  ),
  value = c(
    nrow(batch_out),
    uniqueN(batch_out$gene_symbol),
    uniqueN(batch_out$batch_id),
    sum(batch_out$batch_role == "candidate_positive"),
    sum(batch_out$batch_role == "counterfactual_review"),
    sum(batch_out$batch_role == "protocol_pooling_audit"),
    sum(batch_out$batch_role == "context_audit"),
    sum(batch_out$batch_role == "measurement_review"),
    sum(batch_out$batch_role == "replication_followup"),
    sum(batch_out$batch_role == "negative_control"),
    round(median(num_col(batch_out, "acquisition_priority_score"),
                 na.rm = TRUE), 2),
    round(median(num_col(batch_out, "counterfactual_readiness_score"),
                 na.rm = TRUE), 2)
  )
)
fwrite(metrics, metrics_out)

count_plot <- ggplot(
  batch_out,
  aes(x = batch_id, fill = batch_role)
) +
  geom_bar(width = 0.74) +
  labs(
    title = "RDG Active-Review Batches",
    subtitle = "Each batch is a different acquisition job, not a combined discovery score.",
    x = "Batch",
    y = "Rows",
    fill = "Review role"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(figure_dir, "review_batch_action_counts.png"),
       count_plot, width = 8, height = 4.8, dpi = 160)
ggsave(file.path(figure_dir, "review_batch_action_counts.pdf"),
       count_plot, width = 8, height = 4.8)

plot_dt <- copy(batch_out)
plot_dt[, plot_label := fifelse(
  global_review_order <= 28,
  paste(gene_symbol, design_family, sep = "\n"),
  ""
)]
batch_shape_ids <- sort(unique(as.character(plot_dt$batch_id)))
batch_shape_values <- c(16, 17, 15, 3, 7, 8, 4, 18, 0, 1)
batch_shape_values <- batch_shape_values[seq_along(batch_shape_ids)]
names(batch_shape_values) <- batch_shape_ids
priority_plot <- ggplot(
  plot_dt,
  aes(
    x = counterfactual_readiness_score,
    y = acquisition_priority_score,
    color = batch_role,
    shape = batch_id,
    label = plot_label
  )
) +
  geom_point(alpha = 0.82, size = 2.7) +
  geom_text(check_overlap = TRUE, size = 2.3, vjust = -0.7,
            show.legend = FALSE) +
  scale_shape_manual(values = batch_shape_values, drop = FALSE) +
  labs(
    title = "Review-Batch Effect Space",
    subtitle = "Readiness is current predictive usability; priority is value of acquiring the missing evidence.",
    x = "Counterfactual readiness",
    y = "Acquisition priority",
    color = "Review role",
    shape = "Batch"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(figure_dir, "review_batch_priority_readiness.png"),
       priority_plot, width = 10, height = 7, dpi = 160)
ggsave(file.path(figure_dir, "review_batch_priority_readiness.pdf"),
       priority_plot, width = 10, height = 7)
if (exists("dominant_save_ggplotly") && requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    priority_plot,
    file.path(figure_dir, "review_batch_priority_readiness.html"),
    title = "review_batch_priority_readiness"
  )
}

top_genes <- head(gene_balance$gene_symbol, 45)
matrix_dt <- batch_out[gene_symbol %chin% top_genes, .(
  max_priority = max(num_col(.SD, "acquisition_priority_score"), na.rm = TRUE),
  role_label = paste(unique(batch_role), collapse = "\n"),
  n_rows = .N
), by = .(gene_symbol, batch_id)]
matrix_dt[, gene_symbol := factor(gene_symbol, levels = rev(top_genes))]
matrix_plot <- ggplot(
  matrix_dt,
  aes(x = batch_id, y = gene_symbol, fill = max_priority)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = n_rows), size = 2.2) +
  scale_fill_gradient(low = "#f2f2f2", high = "#8b0000",
                      name = "Max priority") +
  labs(
    title = "Gene Coverage Across Review Batches",
    subtitle = "Cell labels are number of selected rows for that gene and batch.",
    x = "Batch",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
ggsave(file.path(figure_dir, "review_batch_gene_matrix.png"),
       matrix_plot, width = 8.5, height = 9.5, dpi = 160)
ggsave(file.path(figure_dir, "review_batch_gene_matrix.pdf"),
       matrix_plot, width = 8.5, height = 9.5)

top_rows <- batch_out[seq_len(min(.N, 20))]
report_lines <- c(
  "# RDG Active-Review Batches",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Purpose",
  "",
  paste(
    "This layer converts the evidence-gap queue into a small, balanced",
    "active-learning review plan. It is meant to guide browser validation,",
    "negative-control review, protocol/pooling audits, context audits,",
    "translon fixes, and replication follow-up."
  ),
  "",
  "## Batch Summary",
  "",
  paste0(
    "- ", batch_summary$batch_id, " / ", batch_summary$batch_name,
    ": ", batch_summary$n_rows, " rows, ", batch_summary$n_genes,
    " genes. Goal: ", batch_summary$batch_goal
  ),
  "",
  "## First Rows",
  "",
  paste0(
    seq_len(nrow(top_rows)), ". ",
    top_rows$batch_id, " | ", top_rows$gene_symbol, " | ",
    top_rows$design_family, " | role=", top_rows$batch_role,
    " | priority=", top_rows$acquisition_priority_score,
    " | readiness=", top_rows$counterfactual_readiness_score,
    " | instruction=", top_rows$review_instruction
  ),
  "",
  "## How To Use",
  "",
  paste(
    "Use the template CSV for manual review. Fill the review columns and",
    "keep the original batch/order columns unchanged, so later scripts can",
    "learn from positive, weak, negative, confounded, and not-reviewable",
    "labels."
  )
)
writeLines(report_lines, report_out, useBytes = TRUE)

message("Wrote active-review batches: ", batches_out)
message("Wrote active-review template: ", template_out)
