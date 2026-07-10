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
output_dir <- file.path(analysis_dir, "dominant_rdg_hierarchical_context_model")
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

clean_factor <- function(x, prefix = "unknown") {
  x <- tolower(as.character(x))
  x[is.na(x) | !nzchar(x)] <- prefix
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- prefix
  factor(x)
}

metric_row <- function(dt, scope) {
  data.table(
    scope = scope,
    n_predictions = nrow(dt),
    original_mae_log2fc = mean_safe(dt$prediction_abs_error_log2fc),
    context_model_mae_log2fc = mean_safe(dt$context_model_abs_error_log2fc),
    original_rmse_log2fc = rmse_safe(dt$prediction_error_log2fc),
    context_model_rmse_log2fc = rmse_safe(dt$context_model_error_log2fc),
    original_q90_abs_error_log2fc =
      quantile_safe(dt$prediction_abs_error_log2fc, 0.90),
    context_model_q90_abs_error_log2fc =
      quantile_safe(dt$context_model_abs_error_log2fc, 0.90),
    original_q95_abs_error_log2fc =
      quantile_safe(dt$prediction_abs_error_log2fc, 0.95),
    context_model_q95_abs_error_log2fc =
      quantile_safe(dt$context_model_abs_error_log2fc, 0.95),
    original_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$predicted_log2fc),
           na.rm = TRUE),
    context_model_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$context_model_predicted_log2fc),
           na.rm = TRUE)
  )
}

message("Hierarchical context-aware residual model:")
message("  1. Train residual model with gene, context, and gene/context terms")
message("  2. Evaluate by leave-one-study-out without using held-out study rows")
message("  3. Compare against original hierarchy and peer-offset context calibration")

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required for the context-aware residual model.",
       call. = FALSE)
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required for sparse model matrices.",
       call. = FALSE)
}

loo_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical_residual_audit",
  "dominant_rdg_hierarchical_residual_audit_loo_context.csv"
)
residual_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical_residual_audit",
  "dominant_rdg_hierarchical_residual_audit_gene_design.csv"
)
context_calibration_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical_context_calibration",
  "dominant_rdg_hierarchical_context_calibration_gene_design.csv"
)

loo <- read_dt(loo_file, required = TRUE)
gene_residual <- read_dt(residual_file)
context_calibration <- read_dt(context_calibration_file)

required <- c(
  "gene_symbol", "tx_id", "design_family", "heldout_study",
  "heldout_log2fc", "predicted_log2fc", "prediction_error_log2fc",
  "prediction_abs_error_log2fc", "cell_lines", "tissues",
  "condition_tokens"
)
missing <- setdiff(required, names(loo))
if (length(missing)) {
  stop("LOO context table is missing required columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

model_dt <- copy(loo)
for (column in c("gene_symbol", "design_family", "heldout_study",
                 "cell_lines", "tissues", "condition_tokens")) {
  model_dt[is.na(get(column)) | !nzchar(get(column)), (column) := "unknown"]
}
model_dt[, `:=`(
  f_gene = clean_factor(gene_symbol, "gene_unknown"),
  f_design = clean_factor(design_family, "design_unknown"),
  f_study = clean_factor(heldout_study, "study_unknown"),
  f_cell = clean_factor(cell_lines, "cell_unknown"),
  f_tissue = clean_factor(tissues, "tissue_unknown"),
  f_condition = clean_factor(condition_tokens, "condition_unknown"),
  predicted_log2fc_num = safe_num(predicted_log2fc),
  abs_predicted_log2fc = abs(safe_num(predicted_log2fc)),
  log_matched_rows = log1p(safe_num(matched_rows_in_context)),
  log_strict_rows = log1p(safe_num(strict_rows_in_context)),
  log_case_counts = log1p(safe_num(total_case_counts)),
  log_control_counts = log1p(safe_num(total_control_counts))
)]
for (column in c("predicted_log2fc_num", "abs_predicted_log2fc",
                 "log_matched_rows", "log_strict_rows", "log_case_counts",
                 "log_control_counts")) {
  model_dt[!is.finite(get(column)), (column) := 0]
}

keep <- is.finite(model_dt$heldout_log2fc) &
  is.finite(model_dt$predicted_log2fc) &
  is.finite(model_dt$prediction_error_log2fc)
model_dt <- model_dt[keep]
if (!nrow(model_dt)) stop("No modelable LOO rows.", call. = FALSE)

residual_cap <- 3
model_dt[, residual_target := pmax(pmin(prediction_error_log2fc,
                                        residual_cap), -residual_cap)]

model_formula <- stats::as.formula(
  "~ 0 + f_gene + f_design + f_cell + f_tissue + f_condition +
   f_gene:f_design +
   f_design:f_cell + f_design:f_tissue + f_design:f_condition +
   f_gene:f_cell + f_gene:f_tissue +
   predicted_log2fc_num + abs_predicted_log2fc +
   log_matched_rows + log_strict_rows + log_case_counts + log_control_counts"
)

x_all <- Matrix::sparse.model.matrix(model_formula, data = model_dt)
y_all <- model_dt$residual_target

study_groups <- sort(unique(as.character(model_dt$heldout_study)))
fold_id <- match(as.character(model_dt$heldout_study), study_groups)

set.seed(1)
cvfit <- glmnet::cv.glmnet(
  x_all,
  y_all,
  family = "gaussian",
  alpha = 0,
  foldid = fold_id,
  standardize = TRUE,
  intercept = TRUE,
  nlambda = 80
)

predict_blocked_residuals <- function(lambda_value) {
  offset <- rep(NA_real_, nrow(model_dt))
  train_rows <- rep(NA_integer_, nrow(model_dt))
  for (study_id in study_groups) {
    test_idx <- which(as.character(model_dt$heldout_study) == study_id)
    train_idx <- which(as.character(model_dt$heldout_study) != study_id)
    if (length(test_idx) == 0 || length(train_idx) < 100) next
    fit <- glmnet::glmnet(
      x_all[train_idx, , drop = FALSE],
      y_all[train_idx],
      family = "gaussian",
      alpha = 0,
      lambda = lambda_value,
      standardize = TRUE,
      intercept = TRUE
    )
    offset[test_idx] <- as.numeric(stats::predict(
      fit,
      newx = x_all[test_idx, , drop = FALSE],
      s = lambda_value
    ))
    train_rows[test_idx] <- length(train_idx)
  }
  offset[!is.finite(offset)] <- 0
  data.table(
    context_model_residual_offset_log2fc = offset,
    context_model_train_rows = train_rows
  )
}

lambda_candidates <- data.table(
  lambda_name = c("lambda_min", "lambda_1se"),
  lambda_value = c(cvfit$lambda.min, cvfit$lambda.1se)
)
lambda_candidates <- lambda_candidates[is.finite(lambda_value)]
lambda_candidates <- unique(lambda_candidates, by = "lambda_value")
if (!nrow(lambda_candidates)) {
  stop("No finite glmnet lambda candidates were available.", call. = FALSE)
}

candidate_predictions <- vector("list", nrow(lambda_candidates))
lambda_selection <- rbindlist(lapply(seq_len(nrow(lambda_candidates)),
                                     function(i) {
  pred_dt <- predict_blocked_residuals(lambda_candidates$lambda_value[[i]])
  candidate_predictions[[i]] <<- pred_dt
  predicted <- model_dt$predicted_log2fc +
    pred_dt$context_model_residual_offset_log2fc
  err <- model_dt$heldout_log2fc - predicted
  data.table(
    lambda_name = lambda_candidates$lambda_name[[i]],
    lambda_value = lambda_candidates$lambda_value[[i]],
    mae_log2fc = mean_safe(abs(err)),
    rmse_log2fc = rmse_safe(err),
    q95_abs_error_log2fc = quantile_safe(abs(err), 0.95),
    direction_accuracy = mean(sign(model_dt$heldout_log2fc) ==
                                sign(predicted), na.rm = TRUE),
    mean_abs_offset_log2fc =
      mean_safe(abs(pred_dt$context_model_residual_offset_log2fc))
  )
}), fill = TRUE)
setorder(lambda_selection, rmse_log2fc, mae_log2fc, q95_abs_error_log2fc)
lambda_use <- lambda_selection$lambda_value[[1]]
lambda_label <- lambda_selection$lambda_name[[1]]
chosen_idx <- match(lambda_use, lambda_candidates$lambda_value)
chosen_pred <- candidate_predictions[[chosen_idx]]

model_dt[, `:=`(
  context_model_residual_offset_log2fc =
    chosen_pred$context_model_residual_offset_log2fc,
  context_model_train_rows = chosen_pred$context_model_train_rows,
  context_model_lambda = lambda_use,
  context_model_lambda_name = lambda_label
)]

if (!nrow(model_dt[is.finite(context_model_train_rows)])) {
  warning("No blocked context-model predictions had explicit train rows.")
}

model_dt[, context_model_predicted_log2fc :=
           predicted_log2fc + context_model_residual_offset_log2fc]
model_dt[, context_model_error_log2fc :=
           heldout_log2fc - context_model_predicted_log2fc]
model_dt[, context_model_abs_error_log2fc :=
           abs(context_model_error_log2fc)]
model_dt[, context_model_error_delta_log2fc :=
           prediction_abs_error_log2fc - context_model_abs_error_log2fc]
model_dt[, context_model_improved :=
           context_model_abs_error_log2fc < prediction_abs_error_log2fc]
model_dt[, context_model_improvement_fraction := fifelse(
  is.finite(prediction_abs_error_log2fc) & prediction_abs_error_log2fc > 0,
  context_model_error_delta_log2fc / prediction_abs_error_log2fc,
  NA_real_
)]
model_dt[, context_model_correction_class := fifelse(
  context_model_error_delta_log2fc > 0.05,
  "context_model_improved",
  fifelse(context_model_error_delta_log2fc < -0.05,
          "context_model_worse",
          "context_model_neutral")
)]

summary <- rbindlist(c(
  list(metric_row(model_dt, "all")),
  lapply(split(model_dt, model_dt$design_family), function(x) {
    metric_row(x, paste0("design_family:", x$design_family[[1]]))
  })
), fill = TRUE)
summary[, `:=`(
  mae_delta_log2fc = original_mae_log2fc - context_model_mae_log2fc,
  rmse_delta_log2fc = original_rmse_log2fc - context_model_rmse_log2fc,
  q95_delta_log2fc =
    original_q95_abs_error_log2fc - context_model_q95_abs_error_log2fc,
  direction_accuracy_delta =
    context_model_direction_accuracy - original_direction_accuracy
)]

gene_context <- model_dt[, {
  worst_idx <- which.max(context_model_abs_error_log2fc)
  if (!length(worst_idx) || is.na(worst_idx)) worst_idx <- 1L
  data.table(
    context_model_predictions = .N,
    context_model_original_mae_log2fc =
      mean_safe(prediction_abs_error_log2fc),
    context_model_mae_log2fc =
      mean_safe(context_model_abs_error_log2fc),
    context_model_original_rmse_log2fc =
      rmse_safe(prediction_error_log2fc),
    context_model_rmse_log2fc =
      rmse_safe(context_model_error_log2fc),
    context_model_original_q95_log2fc =
      quantile_safe(prediction_abs_error_log2fc, 0.95),
    context_model_q95_log2fc =
      quantile_safe(context_model_abs_error_log2fc, 0.95),
    context_model_mean_delta_log2fc =
      mean_safe(context_model_error_delta_log2fc),
    context_model_median_delta_log2fc =
      median_safe(context_model_error_delta_log2fc),
    context_model_fraction_improved =
      mean(context_model_improved, na.rm = TRUE),
    context_model_mean_abs_offset_log2fc =
      mean_safe(abs(context_model_residual_offset_log2fc)),
    context_model_q95_abs_offset_log2fc =
      quantile_safe(abs(context_model_residual_offset_log2fc), 0.95),
    context_model_worst_context =
      paste(heldout_study[worst_idx], cell_lines[worst_idx],
            tissues[worst_idx], condition_tokens[worst_idx], sep = " | "),
    context_model_worst_abs_error_log2fc =
      context_model_abs_error_log2fc[worst_idx],
    context_model_top_improved_contexts = collapse_unique(
      paste(heldout_study, cell_lines, tissues, condition_tokens, sep = " | ")
        [order(-context_model_error_delta_log2fc)],
      5L
    )
  )
}, by = .(gene_symbol, tx_id, design_family)]

gene_context[, `:=`(
  context_model_mae_delta_log2fc =
    context_model_original_mae_log2fc - context_model_mae_log2fc,
  context_model_rmse_delta_log2fc =
    context_model_original_rmse_log2fc - context_model_rmse_log2fc,
  context_model_q95_delta_log2fc =
    context_model_original_q95_log2fc - context_model_q95_log2fc
)]
gene_context[, context_model_status := fifelse(
  context_model_rmse_delta_log2fc > 0.05 &
    context_model_mae_delta_log2fc >= -0.02,
  "context_model_reduces_residual_error",
  fifelse(context_model_rmse_delta_log2fc < -0.05 |
            context_model_mae_delta_log2fc < -0.05,
          "context_model_increases_residual_error",
          "context_model_neutral")
)]

if (nrow(gene_residual)) {
  keep_cols <- intersect(names(gene_residual), c(
    "gene_symbol", "tx_id", "design_family", "residual_risk_class",
    "residual_heterogeneity_score", "residual_worst_context_key"
  ))
  gene_context <- merge(
    gene_context,
    gene_residual[, ..keep_cols],
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE
  )
}
if (nrow(context_calibration)) {
  keep_cols <- intersect(names(context_calibration), c(
    "gene_symbol", "tx_id", "design_family",
    "context_calibration_status",
    "context_calibration_rmse_delta_log2fc",
    "context_calibration_fraction_improved"
  ))
  gene_context <- merge(
    gene_context,
    context_calibration[, ..keep_cols],
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE
  )
}

final_fit <- glmnet::glmnet(
  x_all,
  y_all,
  family = "gaussian",
  alpha = 0,
  lambda = lambda_use,
  standardize = TRUE,
  intercept = TRUE
)
coef_mat <- as.matrix(stats::coef(final_fit, s = lambda_use))
coef_dt <- data.table(
  feature = rownames(coef_mat),
  coefficient = as.numeric(coef_mat[, 1])
)
coef_dt <- coef_dt[feature != "(Intercept)"]
coef_dt[, abs_coefficient := abs(coefficient)]
setorder(coef_dt, -abs_coefficient)
coef_dt <- coef_dt[seq_len(min(.N, 500))]

summary_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_summary.csv"
)
loo_file_out <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_loo.csv"
)
gene_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_gene_design.csv"
)
viral_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_viral_gene_design.csv"
)
coef_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_coefficients.csv"
)
lambda_selection_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_lambda_selection.csv"
)
metadata_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_context_model_metadata.csv"
)

fwrite(summary, summary_file)
fwrite(model_dt, loo_file_out)
fwrite(gene_context, gene_file)
fwrite(gene_context[design_family == "viral_infection"], viral_file)
fwrite(coef_dt, coef_file)
fwrite(lambda_selection, lambda_selection_file)
fwrite(data.table(
  model = "glmnet_ridge_context_residual_model",
  lambda_min = cvfit$lambda.min,
  lambda_1se = cvfit$lambda.1se,
  lambda_used = lambda_use,
  lambda_used_name = lambda_label,
  n_rows = nrow(model_dt),
  n_features = ncol(x_all),
  residual_cap_log2fc = residual_cap,
  blocked_cv_groups = length(study_groups)
), metadata_file)

if (nrow(summary)) {
  plot_dt <- summary[grepl("^design_family:", scope)]
  plot_dt[, design_family := sub("^design_family:", "", scope)]
  plot_dt <- melt(
    plot_dt,
    id.vars = "design_family",
    measure.vars = c("original_rmse_log2fc", "context_model_rmse_log2fc"),
    variable.name = "metric",
    value.name = "rmse_log2fc"
  )
  plot_dt[, metric := fifelse(
    metric == "original_rmse_log2fc",
    "hierarchical", "context-aware residual model"
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
    scale_fill_manual(values = c(
      "hierarchical" = "#6b7280",
      "context-aware residual model" = "#047857"
    )) +
    labs(
      title = "Context-Aware Hierarchical Residual Model",
      subtitle = "Leave-one-study-out residual correction trained on other studies",
      x = NULL, y = "RMSE (log2FC)", fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_context_model_rmse.png"),
         p, width = 9.5, height = 6.8, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_context_model_rmse.pdf"),
         p, width = 9.5, height = 6.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "hierarchical_context_model_rmse.html"),
      title = "Context-aware hierarchical residual model"
    )
  }
}

viral_gene <- gene_context[design_family == "viral_infection"]
if (nrow(viral_gene)) {
  viral_plot <- viral_gene[
    order(-abs(context_model_rmse_delta_log2fc))
  ][seq_len(min(.N, 30))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p2 <- ggplot(
    viral_plot,
    aes(x = context_model_rmse_delta_log2fc, y = gene_symbol,
        fill = context_model_status,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Status: ", context_model_status,
          "<br>RMSE delta: ", round(context_model_rmse_delta_log2fc, 3),
          "<br>MAE delta: ", round(context_model_mae_delta_log2fc, 3),
          "<br>Original RMSE: ",
          round(context_model_original_rmse_log2fc, 3),
          "<br>Context model RMSE: ",
          round(context_model_rmse_log2fc, 3),
          "<br>Worst context: ", context_model_worst_context
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(
      "context_model_reduces_residual_error" = "#0f766e",
      "context_model_increases_residual_error" = "#b91c1c",
      "context_model_neutral" = "#6b7280"
    )) +
    labs(
      title = "Viral-Infection Context-Aware Residual Model by Gene",
      subtitle = "Positive values mean the context-aware model improves held-out RMSE",
      x = "RMSE reduction (log2FC)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_context_model_viral_genes.png"),
         p2, width = 9.5, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_context_model_viral_genes.pdf"),
         p2, width = 9.5, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir, "hierarchical_context_model_viral_genes.html"),
      title = "Viral context-aware residual model"
    )
  }
}

message("Saved context-aware hierarchical residual model outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
