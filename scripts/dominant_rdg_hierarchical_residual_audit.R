#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) return(here)
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) && dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
output_dir <- file.path(analysis_dir, "dominant_rdg_hierarchical_residual_audit")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

read_dt <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing required input: ", path, call. = FALSE)
    return(data.table())
  }
  fread(path, showProgress = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

collapse_unique <- function(x, n = 6L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

first_nonempty <- function(x, default = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else default
}

mean_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

median_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

max_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

rmse_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else sqrt(mean(x^2))
}

quantile_safe <- function(x, prob = 0.95, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE,
                             type = 8, na.rm = TRUE))
}

bool_col <- function(dt, col, default = FALSE) {
  if (!col %in% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[col]]
  if (is.logical(x)) return(fifelse(is.na(x), default, x))
  if (is.numeric(x)) return(fifelse(is.na(x), default, x != 0))
  y <- tolower(as.character(x))
  fifelse(is.na(y), default, y %chin% c("true", "t", "1", "yes"))
}

normalize_text <- function(x) {
  x <- tolower(as.character(x))
  x[is.na(x)] <- ""
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "unknown"
  x
}

message("Hierarchical residual heterogeneity audit:")
message("  1. Attach leave-one-study-out residuals to matched-design context")
message("  2. Summarize residual risk by gene, study, cell line, tissue, and condition")
message("  3. Export atlas-ready residual-risk fields")

cds_file <- file.path(
  analysis_dir, "dominant_next_model",
  "dominant_next_model_matched_control_cds_buffering.csv"
)
loo_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical",
  "dominant_rdg_hierarchical_leave_one_study_calibration.csv"
)
gene_design_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical",
  "dominant_rdg_hierarchical_gene_design_effects.csv"
)
interval_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical",
  "dominant_rdg_hierarchical_interval_calibration.csv"
)

cds <- read_dt(cds_file, required = TRUE)
loo <- read_dt(loo_file, required = TRUE)
gene_design <- read_dt(gene_design_file, required = TRUE)
interval_calibration <- read_dt(interval_file)

if (!"design_family" %in% names(cds)) cds[, design_family := perturbation_class]
cds[is.na(design_family) | !nzchar(design_family), design_family := "unknown"]
cds <- cds[feature_id == "clean_CDS" | feature_type == "clean_CDS"]

cds[, effect_log2fc := fifelse(
  is.finite(safe_num(shrunk_matched_fpkm_log2_fc)),
  safe_num(shrunk_matched_fpkm_log2_fc),
  safe_num(matched_fpkm_log2_fc_pseudocount1)
)]
cds[, strict_support := bool_col(.SD, "matched_feature_passes_sample_gate") &
      bool_col(.SD, "matched_feature_passes_count_gate_both") &
      safe_num(n_case_samples) >= 2 &
      safe_num(n_control_samples) >= 2 &
      safe_num(min_feature_raw_counts) >= 30]
cds[, supported_support := strict_support |
      (bool_col(.SD, "matched_feature_passes_count_gate_any") &
         safe_num(min_matched_samples) >= 2 &
         safe_num(max_feature_raw_counts) >= 30)]
cds[, condition_token := normalize_text(case_condition)]
cds[, context_key := paste(study, CELL_LINE, TISSUE, condition_token, sep = " | ")]

context_by_study <- cds[, .(
  matched_rows_in_context = .N,
  strict_rows_in_context = sum(strict_support, na.rm = TRUE),
  supported_rows_in_context = sum(supported_support, na.rm = TRUE),
  mean_row_effect_log2fc = mean_safe(effect_log2fc),
  median_row_effect_log2fc = median_safe(effect_log2fc),
  top_abs_row_effect_log2fc = max_safe(abs(effect_log2fc)),
  cell_lines = collapse_unique(CELL_LINE, 6L),
  tissues = collapse_unique(TISSUE, 6L),
  authors = collapse_unique(AUTHOR, 4L),
  case_conditions = collapse_unique(case_condition, 8L),
  control_conditions = collapse_unique(control_conditions, 8L),
  condition_tokens = collapse_unique(condition_token, 8L),
  matched_designs = collapse_unique(matched_design_signature, 5L),
  total_case_samples = sum(safe_num(n_case_samples), na.rm = TRUE),
  total_control_samples = sum(safe_num(n_control_samples), na.rm = TRUE),
  total_case_counts = sum(safe_num(case_sum_raw_counts), na.rm = TRUE),
  total_control_counts = sum(safe_num(control_sum_raw_counts), na.rm = TRUE)
), by = .(gene_symbol, tx_id, design_family, heldout_study = study)]

loo_context <- merge(
  loo,
  context_by_study,
  by = c("gene_symbol", "tx_id", "design_family", "heldout_study"),
  all.x = TRUE
)

if (!"hierarchical_conformal_q95_log2fc" %in% names(loo_context)) {
  loo_context[, hierarchical_conformal_q95_log2fc := NA_real_]
}
if (nrow(interval_calibration)) {
  interval_small <- interval_calibration[, .(
    design_family,
    residual_audit_q90_log2fc = hierarchical_conformal_q90_log2fc,
    residual_audit_q95_log2fc = hierarchical_conformal_q95_log2fc,
    residual_audit_calibration_scope = hierarchical_calibration_scope
  )]
  loo_context <- merge(loo_context, interval_small,
                       by = "design_family", all.x = TRUE)
}
if (!"residual_audit_q95_log2fc" %in% names(loo_context)) {
  loo_context[, residual_audit_q95_log2fc := quantile_safe(prediction_abs_error_log2fc, 0.95)]
  loo_context[, residual_audit_q90_log2fc := quantile_safe(prediction_abs_error_log2fc, 0.90)]
  loo_context[, residual_audit_calibration_scope := "global_from_loo"]
}
loo_context[!is.finite(residual_audit_q95_log2fc),
            residual_audit_q95_log2fc := quantile_safe(prediction_abs_error_log2fc, 0.95)]
loo_context[!is.finite(residual_audit_q90_log2fc),
            residual_audit_q90_log2fc := quantile_safe(prediction_abs_error_log2fc, 0.90)]

loo_context[, residual_direction := fifelse(
  prediction_error_log2fc > 0, "heldout_higher_than_predicted",
  fifelse(prediction_error_log2fc < 0, "heldout_lower_than_predicted",
          "no_error")
)]
loo_context[, residual_abs_over_q90 := prediction_abs_error_log2fc >= residual_audit_q90_log2fc]
loo_context[, residual_abs_over_q95 := prediction_abs_error_log2fc >= residual_audit_q95_log2fc]
loo_context[, residual_signed_bias_class := fifelse(
  mean_safe(prediction_error_log2fc) > 0, "under_predicted_heldout",
  fifelse(mean_safe(prediction_error_log2fc) < 0, "over_predicted_heldout",
          "balanced")
), by = .(gene_symbol, tx_id, design_family)]

summarise_residuals <- function(dt, by_cols) {
  dt[, {
    worst_idx <- which.max(prediction_abs_error_log2fc)
    if (!length(worst_idx) || is.na(worst_idx)) worst_idx <- 1L
    data.table(
      residual_predictions = .N,
      residual_mae_log2fc = mean(prediction_abs_error_log2fc, na.rm = TRUE),
      residual_rmse_log2fc = rmse_safe(prediction_error_log2fc),
      residual_median_abs_error_log2fc = median_safe(prediction_abs_error_log2fc),
      residual_q90_abs_error_log2fc = quantile_safe(prediction_abs_error_log2fc, 0.90),
      residual_q95_abs_error_log2fc = quantile_safe(prediction_abs_error_log2fc, 0.95),
      residual_mean_signed_error_log2fc = mean(prediction_error_log2fc, na.rm = TRUE),
      residual_fraction_over_design_q90 = mean(residual_abs_over_q90, na.rm = TRUE),
      residual_fraction_over_design_q95 = mean(residual_abs_over_q95, na.rm = TRUE),
      residual_direction_accuracy = mean(sign(heldout_log2fc) == sign(predicted_log2fc), na.rm = TRUE),
      residual_interval_coverage = mean(prediction_interval_contains_heldout, na.rm = TRUE),
      residual_calibrated_interval_coverage = if ("calibrated_prediction_interval_contains_heldout" %in% names(.SD)) {
        mean(calibrated_prediction_interval_contains_heldout, na.rm = TRUE)
      } else NA_real_,
      residual_worst_abs_error_log2fc = prediction_abs_error_log2fc[worst_idx],
      residual_worst_signed_error_log2fc = prediction_error_log2fc[worst_idx],
      residual_worst_study = first_nonempty(heldout_study[worst_idx]),
      residual_worst_cell_line = first_nonempty(cell_lines[worst_idx]),
      residual_worst_tissue = first_nonempty(tissues[worst_idx]),
      residual_worst_case_condition = first_nonempty(case_conditions[worst_idx]),
      residual_worst_control_condition = first_nonempty(control_conditions[worst_idx]),
      residual_worst_context_key = paste(
        first_nonempty(heldout_study[worst_idx]),
        first_nonempty(cell_lines[worst_idx]),
        first_nonempty(tissues[worst_idx]),
        first_nonempty(case_conditions[worst_idx]),
        sep = " | "
      ),
      residual_top_studies = collapse_unique(heldout_study[order(-prediction_abs_error_log2fc)], 5L),
      residual_top_cell_lines = collapse_unique(cell_lines[order(-prediction_abs_error_log2fc)], 5L),
      residual_top_tissues = collapse_unique(tissues[order(-prediction_abs_error_log2fc)], 5L),
      residual_top_conditions = collapse_unique(case_conditions[order(-prediction_abs_error_log2fc)], 5L)
    )
  }, by = by_cols]
}

gene_design_residual <- summarise_residuals(
  loo_context,
  c("gene_symbol", "tx_id", "design_family")
)

# Merge hierarchy-level fields that help interpret residual risk.
hier_keep <- intersect(names(gene_design), c(
  "gene_symbol", "tx_id", "design_family",
  "hierarchical_evidence_class", "hierarchical_uncalibrated_evidence_class",
  "hierarchical_priority_score", "hierarchical_clean_cds_pct_change",
  "hierarchical_calibrated_clean_cds_pct_lower",
  "hierarchical_calibrated_clean_cds_pct_upper",
  "hierarchical_conformal_q95_log2fc", "hierarchical_calibration_scope",
  "hierarchical_loo_mae_log2fc", "hierarchical_loo_rmse_log2fc",
  "hierarchical_loo_calibrated_interval_coverage",
  "hierarchical_loo_direction_accuracy", "n_studies", "n_strict_studies",
  "n_rows", "n_strict_rows"
))
gene_design_residual <- merge(
  gene_design_residual,
  gene_design[, ..hier_keep],
  by = c("gene_symbol", "tx_id", "design_family"),
  all.x = TRUE
)

gene_design_residual[, residual_relative_to_family_q95 :=
                       residual_rmse_log2fc / pmax(hierarchical_conformal_q95_log2fc, 1e-6)]
gene_design_residual[, residual_heterogeneity_score := pmin(
  1,
  0.35 * pmin(residual_mae_log2fc / 0.5, 1) +
    0.25 * pmin(residual_rmse_log2fc / 0.8, 1) +
    0.20 * pmin(residual_fraction_over_design_q90 / 0.25, 1) +
    0.10 * as.numeric(is.finite(residual_direction_accuracy) &
                        residual_direction_accuracy < 0.55) +
    0.10 * as.numeric(is.finite(residual_interval_coverage) &
                        residual_interval_coverage < 0.60)
)]
gene_design_residual[, residual_risk_class := fifelse(
  residual_predictions < 2,
  "insufficient_loo",
  fifelse(residual_heterogeneity_score >= 0.70,
          "high_residual_heterogeneity",
          fifelse(residual_heterogeneity_score >= 0.45,
                  "moderate_residual_heterogeneity",
                  "low_residual_heterogeneity"))
)]
setorder(gene_design_residual, -residual_heterogeneity_score,
         -residual_rmse_log2fc, gene_symbol, design_family)

context_dims <- list(
  study = c("design_family", "heldout_study"),
  cell_line = c("design_family", "cell_lines"),
  tissue = c("design_family", "tissues"),
  condition = c("design_family", "condition_tokens")
)
context_summaries <- rbindlist(lapply(names(context_dims), function(dim_name) {
  by_cols <- context_dims[[dim_name]]
  out <- summarise_residuals(loo_context, by_cols)
  out[, residual_context_dimension := dim_name]
  if (dim_name == "study") setnames(out, "heldout_study", "residual_context_value")
  if (dim_name == "cell_line") setnames(out, "cell_lines", "residual_context_value")
  if (dim_name == "tissue") setnames(out, "tissues", "residual_context_value")
  if (dim_name == "condition") setnames(out, "condition_tokens", "residual_context_value")
  out
}), fill = TRUE)
context_summaries[, residual_context_score := pmin(
  1,
  0.45 * pmin(residual_mae_log2fc / 0.5, 1) +
    0.30 * pmin(residual_rmse_log2fc / 0.8, 1) +
    0.25 * pmin(residual_fraction_over_design_q90 / 0.25, 1)
)]
setorder(context_summaries, -residual_context_score,
         -residual_predictions, design_family, residual_context_dimension)

viral_context <- context_summaries[design_family == "viral_infection"]
viral_gene_design <- gene_design_residual[design_family == "viral_infection"]

loo_context_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_residual_audit_loo_context.csv"
)
gene_design_file_out <- file.path(
  output_dir, "dominant_rdg_hierarchical_residual_audit_gene_design.csv"
)
context_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_residual_audit_context_summary.csv"
)
viral_context_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_residual_audit_viral_contexts.csv"
)
viral_gene_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_residual_audit_viral_gene_design.csv"
)
summary_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_residual_audit_summary_metrics.csv"
)

fwrite(loo_context, loo_context_file)
fwrite(gene_design_residual, gene_design_file_out)
fwrite(context_summaries, context_file)
fwrite(viral_context, viral_context_file)
fwrite(viral_gene_design, viral_gene_file)

summary_metrics <- data.table(
  metric = c(
    "loo_context_rows",
    "gene_design_residual_rows",
    "high_residual_gene_design_rows",
    "moderate_residual_gene_design_rows",
    "viral_gene_design_rows",
    "viral_high_residual_gene_design_rows",
    "viral_median_mae_log2fc",
    "viral_median_rmse_log2fc",
    "top_context_rows",
    "viral_context_rows"
  ),
  value = c(
    nrow(loo_context),
    nrow(gene_design_residual),
    sum(gene_design_residual$residual_risk_class ==
          "high_residual_heterogeneity", na.rm = TRUE),
    sum(gene_design_residual$residual_risk_class ==
          "moderate_residual_heterogeneity", na.rm = TRUE),
    nrow(viral_gene_design),
    sum(viral_gene_design$residual_risk_class ==
          "high_residual_heterogeneity", na.rm = TRUE),
    median_safe(viral_gene_design$residual_mae_log2fc),
    median_safe(viral_gene_design$residual_rmse_log2fc),
    nrow(context_summaries),
    nrow(viral_context)
  )
)
fwrite(summary_metrics, summary_file)

# Plots
if (nrow(gene_design_residual)) {
  p1 <- ggplot(
    gene_design_residual,
    aes(x = residual_mae_log2fc, y = residual_rmse_log2fc,
        color = design_family,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Design: ", design_family,
          "<br>Risk: ", residual_risk_class,
          "<br>MAE: ", round(residual_mae_log2fc, 3),
          "<br>RMSE: ", round(residual_rmse_log2fc, 3),
          "<br>Worst: ", residual_worst_context_key
        ))
  ) +
    geom_point(alpha = 0.70, size = 2) +
    geom_vline(xintercept = 0.5, color = "grey65", linetype = "dashed",
               linewidth = 0.25) +
    geom_hline(yintercept = 0.8, color = "grey65", linetype = "dashed",
               linewidth = 0.25) +
    labs(
      title = "Hierarchical Residual Heterogeneity by Gene and Design Family",
      subtitle = "Leave-one-study-out residual error after hierarchical shrinkage",
      x = "LOO MAE (log2FC)", y = "LOO RMSE (log2FC)"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_residual_gene_design_scatter.png"),
         p1, width = 9.5, height = 6.2, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_residual_gene_design_scatter.pdf"),
         p1, width = 9.5, height = 6.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p1,
      file.path(figure_dir, "hierarchical_residual_gene_design_scatter.html"),
      title = "Hierarchical residual gene-design scatter"
    )
  }
}

if (nrow(viral_context)) {
  viral_plot <- viral_context[
    residual_predictions >= 3
  ][order(-residual_context_score)][seq_len(min(.N, 30))]
  if (nrow(viral_plot)) {
    viral_plot[, label := paste(residual_context_dimension,
                                residual_context_value, sep = ": ")]
    viral_plot[, label := factor(label, levels = rev(label))]
    p2 <- ggplot(
      viral_plot,
      aes(x = residual_context_score, y = label,
          fill = residual_mae_log2fc,
          text = paste0(
            "Context: ", label,
            "<br>Predictions: ", residual_predictions,
            "<br>MAE: ", round(residual_mae_log2fc, 3),
            "<br>RMSE: ", round(residual_rmse_log2fc, 3),
            "<br>Worst gene/study: ", residual_worst_context_key
          ))
    ) +
      geom_col(color = "grey25", linewidth = 0.15) +
      scale_fill_gradient(low = "#e8f1ff", high = "#8b0000",
                          name = "MAE") +
      labs(
        title = "Viral-Infection Hierarchical Residual Contexts",
        subtitle = "Contexts where broad viral-infection predictions fail most",
        x = "Residual context score", y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            panel.grid.major.y = element_blank())
    ggsave(file.path(figure_dir, "hierarchical_residual_viral_contexts.png"),
           p2, width = 10.5, height = 7.2, dpi = 180)
    ggsave(file.path(figure_dir, "hierarchical_residual_viral_contexts.pdf"),
           p2, width = 10.5, height = 7.2)
    if (exists("dominant_save_ggplotly") &&
        requireNamespace("plotly", quietly = TRUE)) {
      dominant_save_ggplotly(
        p2,
        file.path(figure_dir, "hierarchical_residual_viral_contexts.html"),
        title = "Hierarchical residual viral contexts"
      )
    }
  }
}

message("Saved hierarchical residual audit outputs in: ", output_dir)
print(summary_metrics)
