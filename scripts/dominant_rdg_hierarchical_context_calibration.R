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
output_dir <- file.path(analysis_dir, "dominant_rdg_hierarchical_context_calibration")
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

rmse_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else sqrt(mean(x^2))
}

quantile_safe <- function(x, prob = 0.95) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, type = 8,
                             na.rm = TRUE, names = FALSE))
}

collapse_unique <- function(x, n = 6L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

metric_row <- function(dt, scope) {
  data.table(
    scope = scope,
    n_predictions = nrow(dt),
    original_mae_log2fc = mean_safe(dt$prediction_abs_error_log2fc),
    context_mae_log2fc = mean_safe(dt$context_abs_error_log2fc),
    original_rmse_log2fc = rmse_safe(dt$prediction_error_log2fc),
    context_rmse_log2fc = rmse_safe(dt$context_prediction_error_log2fc),
    original_q90_abs_error_log2fc =
      quantile_safe(dt$prediction_abs_error_log2fc, 0.90),
    context_q90_abs_error_log2fc =
      quantile_safe(dt$context_abs_error_log2fc, 0.90),
    original_q95_abs_error_log2fc =
      quantile_safe(dt$prediction_abs_error_log2fc, 0.95),
    context_q95_abs_error_log2fc =
      quantile_safe(dt$context_abs_error_log2fc, 0.95),
    original_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$predicted_log2fc),
           na.rm = TRUE),
    context_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$context_predicted_log2fc),
           na.rm = TRUE),
    peer_corrected_fraction = mean(dt$context_peer_count >= 5, na.rm = TRUE)
  )
}

add_peer_offset <- function(dt, by_cols, prefix) {
  cols <- c(by_cols, "gene_symbol")
  group <- dt[, .(
    group_n = .N,
    group_sum_error = sum(prediction_error_log2fc, na.rm = TRUE),
    group_mean_error = mean(prediction_error_log2fc, na.rm = TRUE),
    group_median_abs_error = median_safe(prediction_abs_error_log2fc)
  ), by = by_cols]
  gene_group <- dt[, .(
    same_gene_n = .N,
    same_gene_sum_error = sum(prediction_error_log2fc, na.rm = TRUE)
  ), by = cols]
  out <- merge(dt[, c(by_cols, "gene_symbol"), with = FALSE],
               group, by = by_cols, all.x = TRUE)
  out <- merge(out, gene_group, by = cols, all.x = TRUE)
  out[, peer_n := group_n - same_gene_n]
  out[, peer_offset := (group_sum_error - same_gene_sum_error) /
        pmax(peer_n, 1)]
  out[peer_n <= 0 | !is.finite(peer_offset), peer_offset := NA_real_]
  out[, peer_shrinkage := fifelse(
    is.finite(peer_offset),
    pmin(1, peer_n / (peer_n + 10)),
    NA_real_
  )]
  out[, peer_offset_shrunk := peer_offset * peer_shrinkage]
  keep <- c(cols, "peer_n", "peer_offset", "peer_shrinkage",
            "peer_offset_shrunk", "group_n", "group_mean_error",
            "group_median_abs_error")
  out <- out[, ..keep]
  setnames(
    out,
    c("peer_n", "peer_offset", "peer_shrinkage", "peer_offset_shrunk",
      "group_n", "group_mean_error", "group_median_abs_error"),
    paste0(prefix, c("_peer_n", "_offset", "_shrinkage", "_offset_shrunk",
                     "_group_n", "_mean_error", "_median_abs_error"))
  )
  unique(out)
}

message("Hierarchical context calibration:")
message("  1. Use peer genes to estimate study/cell/tissue/condition residual offsets")
message("  2. Apply shrinkage-corrected context offsets to held-out predictions")
message("  3. Measure whether broad design-family residuals shrink")

loo_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical_residual_audit",
  "dominant_rdg_hierarchical_residual_audit_loo_context.csv"
)
gene_design_residual_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical_residual_audit",
  "dominant_rdg_hierarchical_residual_audit_gene_design.csv"
)

loo <- read_dt(loo_file, required = TRUE)
gene_design_residual <- read_dt(gene_design_residual_file, required = TRUE)

required <- c(
  "gene_symbol", "tx_id", "design_family", "heldout_study",
  "heldout_log2fc", "predicted_log2fc", "prediction_error_log2fc",
  "prediction_abs_error_log2fc", "cell_lines", "tissues", "condition_tokens"
)
missing <- setdiff(required, names(loo))
if (length(missing)) {
  stop("LOO context table is missing required columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

for (column in c("cell_lines", "tissues", "condition_tokens", "heldout_study")) {
  loo[is.na(get(column)) | !nzchar(get(column)), (column) := "unknown"]
}

loo[, context_exact := paste(
  design_family, heldout_study, cell_lines, tissues, condition_tokens,
  sep = " | "
)]
loo[, context_cell_tissue_condition := paste(
  design_family, cell_lines, tissues, condition_tokens,
  sep = " | "
)]
loo[, context_cell_tissue := paste(design_family, cell_lines, tissues,
                                   sep = " | ")]
loo[, context_condition := paste(design_family, condition_tokens, sep = " | ")]

offset_specs <- list(
  exact = c("design_family", "heldout_study", "cell_lines", "tissues",
            "condition_tokens"),
  study = c("design_family", "heldout_study"),
  cell_tissue_condition = c("design_family", "cell_lines", "tissues",
                            "condition_tokens"),
  cell_tissue = c("design_family", "cell_lines", "tissues"),
  condition = c("design_family", "condition_tokens"),
  design_family = c("design_family")
)

loo_offsets <- copy(loo)
offset_tables <- list()
for (name in names(offset_specs)) {
  prefix <- paste0("context_", name)
  offset <- add_peer_offset(loo, offset_specs[[name]], prefix)
  offset_tables[[name]] <- copy(offset)
  merge_by <- c(offset_specs[[name]], "gene_symbol")
  loo_offsets <- merge(loo_offsets, offset, by = merge_by, all.x = TRUE)
}

candidate_levels <- c("exact", "study", "cell_tissue_condition",
                      "cell_tissue", "condition", "design_family")
loo_offsets[, `:=`(
  context_offset_level = "none",
  context_peer_count = 0L,
  context_offset_log2fc = 0,
  context_offset_shrinkage = 0
)]
for (level in candidate_levels) {
  peer_col <- paste0("context_", level, "_peer_n")
  offset_col <- paste0("context_", level, "_offset_shrunk")
  shrink_col <- paste0("context_", level, "_shrinkage")
  take <- loo_offsets$context_offset_level == "none" &
    is.finite(loo_offsets[[offset_col]]) &
    loo_offsets[[peer_col]] >= 5
  if (any(take)) {
    loo_offsets[take, `:=`(
      context_offset_level = level,
      context_peer_count = get(peer_col),
      context_offset_log2fc = get(offset_col),
      context_offset_shrinkage = get(shrink_col)
    )]
  }
}

loo_offsets[, context_predicted_log2fc :=
              predicted_log2fc + context_offset_log2fc]
loo_offsets[, context_prediction_error_log2fc :=
              heldout_log2fc - context_predicted_log2fc]
loo_offsets[, context_abs_error_log2fc :=
              abs(context_prediction_error_log2fc)]
loo_offsets[, context_error_delta_log2fc :=
              prediction_abs_error_log2fc - context_abs_error_log2fc]
loo_offsets[, context_improved :=
              context_abs_error_log2fc < prediction_abs_error_log2fc]
loo_offsets[, context_improvement_fraction := fifelse(
  is.finite(prediction_abs_error_log2fc) & prediction_abs_error_log2fc > 0,
  context_error_delta_log2fc / prediction_abs_error_log2fc,
  NA_real_
)]
loo_offsets[, context_correction_class := fifelse(
  context_peer_count < 5,
  "no_peer_context",
  fifelse(context_error_delta_log2fc > 0.05,
          "context_improved",
          fifelse(context_error_delta_log2fc < -0.05,
                  "context_worse",
                  "context_neutral"))
)]

summary_rows <- list(metric_row(loo_offsets, "all"))
summary_rows <- c(
  summary_rows,
  lapply(
    split(loo_offsets, loo_offsets$design_family),
    function(x) metric_row(x, paste0("design_family:", x$design_family[[1]]))
  )
)
summary <- rbindlist(summary_rows, fill = TRUE)
summary[, `:=`(
  mae_delta_log2fc = original_mae_log2fc - context_mae_log2fc,
  rmse_delta_log2fc = original_rmse_log2fc - context_rmse_log2fc,
  q95_delta_log2fc = original_q95_abs_error_log2fc -
    context_q95_abs_error_log2fc
)]

gene_context <- loo_offsets[, .(
  context_calibration_predictions = .N,
  context_calibration_peer_corrected_fraction =
    mean(context_peer_count >= 5, na.rm = TRUE),
  context_calibration_original_mae_log2fc =
    mean_safe(prediction_abs_error_log2fc),
  context_calibration_context_mae_log2fc =
    mean_safe(context_abs_error_log2fc),
  context_calibration_original_rmse_log2fc =
    rmse_safe(prediction_error_log2fc),
  context_calibration_context_rmse_log2fc =
    rmse_safe(context_prediction_error_log2fc),
  context_calibration_original_q95_log2fc =
    quantile_safe(prediction_abs_error_log2fc, 0.95),
  context_calibration_context_q95_log2fc =
    quantile_safe(context_abs_error_log2fc, 0.95),
  context_calibration_mean_delta_log2fc =
    mean_safe(context_error_delta_log2fc),
  context_calibration_median_delta_log2fc =
    median_safe(context_error_delta_log2fc),
  context_calibration_fraction_improved =
    mean(context_improved, na.rm = TRUE),
  context_calibration_top_offset_level =
    names(sort(table(context_offset_level), decreasing = TRUE))[1],
  context_calibration_top_worse_context =
    paste(
      heldout_study[which.max(context_abs_error_log2fc)],
      cell_lines[which.max(context_abs_error_log2fc)],
      tissues[which.max(context_abs_error_log2fc)],
      condition_tokens[which.max(context_abs_error_log2fc)],
      sep = " | "
    ),
  context_calibration_top_worse_abs_error_log2fc =
    max(context_abs_error_log2fc, na.rm = TRUE),
  context_calibration_top_improved_contexts =
    collapse_unique(
      paste(heldout_study, cell_lines, tissues, condition_tokens, sep = " | ")
        [order(-context_error_delta_log2fc)],
      5L
    )
), by = .(gene_symbol, tx_id, design_family)]

gene_context[, context_calibration_mae_delta_log2fc :=
               context_calibration_original_mae_log2fc -
               context_calibration_context_mae_log2fc]
gene_context[, context_calibration_rmse_delta_log2fc :=
               context_calibration_original_rmse_log2fc -
               context_calibration_context_rmse_log2fc]
gene_context[, context_calibration_q95_delta_log2fc :=
               context_calibration_original_q95_log2fc -
               context_calibration_context_q95_log2fc]
gene_context[, context_calibration_status := fifelse(
  context_calibration_peer_corrected_fraction < 0.5,
  "low_peer_context_support",
  fifelse(context_calibration_rmse_delta_log2fc > 0.05,
          "context_reduces_residual_error",
          fifelse(context_calibration_rmse_delta_log2fc < -0.05,
                  "context_increases_residual_error",
                  "context_neutral"))
)]

if (nrow(gene_design_residual)) {
  keep <- intersect(names(gene_design_residual), c(
    "gene_symbol", "tx_id", "design_family", "residual_risk_class",
    "residual_heterogeneity_score", "residual_mae_log2fc",
    "residual_rmse_log2fc", "residual_worst_context_key"
  ))
  gene_context <- merge(
    gene_context,
    gene_design_residual[, ..keep],
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE
  )
}

context_offset_summary <- rbindlist(lapply(names(offset_tables), function(level) {
  dt <- copy(offset_tables[[level]])
  peer_col <- paste0("context_", level, "_peer_n")
  offset_col <- paste0("context_", level, "_offset")
  shrink_col <- paste0("context_", level, "_shrinkage")
  data.table(
    context_offset_level = level,
    rows = nrow(dt),
    usable_rows = sum(dt[[peer_col]] >= 5 & is.finite(dt[[offset_col]]),
                      na.rm = TRUE),
    median_peer_count = median_safe(dt[[peer_col]]),
    median_abs_offset_log2fc = median_safe(abs(dt[[offset_col]])),
    q95_abs_offset_log2fc = quantile_safe(abs(dt[[offset_col]]), 0.95),
    median_shrinkage = median_safe(dt[[shrink_col]])
  )
}), fill = TRUE)

loo_file_out <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_calibration_loo.csv"
)
gene_file_out <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_calibration_gene_design.csv"
)
summary_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_calibration_summary.csv"
)
offset_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_calibration_offset_summary.csv"
)
viral_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_calibration_viral_gene_design.csv"
)

fwrite(loo_offsets, loo_file_out)
fwrite(gene_context, gene_file_out)
fwrite(summary, summary_file)
fwrite(context_offset_summary, offset_file)
fwrite(gene_context[design_family == "viral_infection"], viral_file)

if (nrow(summary)) {
  plot_dt <- summary[grepl("^design_family:", scope)]
  plot_dt[, design_family := sub("^design_family:", "", scope)]
  plot_dt <- melt(
    plot_dt,
    id.vars = "design_family",
    measure.vars = c("original_rmse_log2fc", "context_rmse_log2fc"),
    variable.name = "metric",
    value.name = "rmse_log2fc"
  )
  plot_dt[, metric := fifelse(
    metric == "original_rmse_log2fc",
    "hierarchical", "context-corrected"
  )]
  p <- ggplot(
    plot_dt,
    aes(x = reorder(design_family, rmse_log2fc, FUN = max),
        y = rmse_log2fc, fill = metric,
        text = paste0(
          "Design: ", design_family,
          "<br>Model: ", metric,
          "<br>RMSE: ", round(rmse_log2fc, 3)
        ))
  ) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c("hierarchical" = "#6b7280",
                                 "context-corrected" = "#2563eb")) +
    labs(
      title = "Context Calibration of Hierarchical RDG Effects",
      subtitle = "Peer-gene context offsets reduce residual error if metadata context is informative",
      x = NULL, y = "Leave-one-study-out RMSE (log2FC)",
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_context_calibration_rmse.png"),
         p, width = 9.5, height = 6.8, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_context_calibration_rmse.pdf"),
         p, width = 9.5, height = 6.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "hierarchical_context_calibration_rmse.html"),
      title = "Hierarchical context calibration RMSE"
    )
  }
}

viral_gene <- gene_context[design_family == "viral_infection"]
if (nrow(viral_gene)) {
  viral_plot <- viral_gene[
    order(-abs(context_calibration_rmse_delta_log2fc))
  ][seq_len(min(.N, 30))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p2 <- ggplot(
    viral_plot,
    aes(x = context_calibration_rmse_delta_log2fc, y = gene_symbol,
        fill = context_calibration_status,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Status: ", context_calibration_status,
          "<br>RMSE delta: ",
          round(context_calibration_rmse_delta_log2fc, 3),
          "<br>Original RMSE: ",
          round(context_calibration_original_rmse_log2fc, 3),
          "<br>Context RMSE: ",
          round(context_calibration_context_rmse_log2fc, 3),
          "<br>Worst after correction: ",
          context_calibration_top_worse_context
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(
      "context_reduces_residual_error" = "#0f766e",
      "context_increases_residual_error" = "#b91c1c",
      "context_neutral" = "#6b7280",
      "low_peer_context_support" = "#9ca3af"
    )) +
    labs(
      title = "Viral-Infection Context Calibration by Gene",
      subtitle = "Positive values mean same-context peer genes improve prediction",
      x = "RMSE reduction (log2FC)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_context_calibration_viral_genes.png"),
         p2, width = 9.5, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_context_calibration_viral_genes.pdf"),
         p2, width = 9.5, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir, "hierarchical_context_calibration_viral_genes.html"),
      title = "Viral gene context calibration"
    )
  }
}

message("Saved hierarchical context-calibration outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
