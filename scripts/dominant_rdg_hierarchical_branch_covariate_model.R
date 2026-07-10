#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        dir.exists(file.path(candidate, "scripts"))) {
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
output_dir <- file.path(analysis_dir,
                        "dominant_rdg_hierarchical_branch_covariates")
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

max_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) 0 else max(x)
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
    metadata_model_mae_log2fc =
      mean_safe(dt$metadata_model_abs_error_log2fc),
    branch_model_mae_log2fc =
      mean_safe(dt$branch_model_abs_error_log2fc),
    original_rmse_log2fc = rmse_safe(dt$prediction_error_log2fc),
    metadata_model_rmse_log2fc =
      rmse_safe(dt$metadata_model_error_log2fc),
    branch_model_rmse_log2fc =
      rmse_safe(dt$branch_model_error_log2fc),
    original_q90_abs_error_log2fc =
      quantile_safe(dt$prediction_abs_error_log2fc, 0.90),
    metadata_model_q90_abs_error_log2fc =
      quantile_safe(dt$metadata_model_abs_error_log2fc, 0.90),
    branch_model_q90_abs_error_log2fc =
      quantile_safe(dt$branch_model_abs_error_log2fc, 0.90),
    original_q95_abs_error_log2fc =
      quantile_safe(dt$prediction_abs_error_log2fc, 0.95),
    metadata_model_q95_abs_error_log2fc =
      quantile_safe(dt$metadata_model_abs_error_log2fc, 0.95),
    branch_model_q95_abs_error_log2fc =
      quantile_safe(dt$branch_model_abs_error_log2fc, 0.95),
    original_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$predicted_log2fc),
           na.rm = TRUE),
    metadata_model_direction_accuracy =
      mean(sign(dt$heldout_log2fc) ==
             sign(dt$metadata_model_predicted_log2fc), na.rm = TRUE),
    branch_model_direction_accuracy =
      mean(sign(dt$heldout_log2fc) ==
             sign(dt$branch_model_predicted_log2fc), na.rm = TRUE)
  )
}

branch_class_columns <- c(
  "train_delta", "train_abs_delta", "train_eval_rows", "train_strict_rows",
  "train_sig_rows", "train_studies"
)

build_branch_covariates <- function(model_dt, branch_dt) {
  classes <- data.table(
    branch_feature_class = c("clean_CDS", "leader_uORF",
                             "overlapping_uORF")
  )
  key_cols <- c("gene_symbol", "tx_id", "design_family",
                "heldout_study")
  keys <- unique(model_dt[, ..key_cols])
  if (!nrow(branch_dt)) {
    out <- copy(keys)
    for (prefix in branch_class_columns) {
      for (class_name in classes$branch_feature_class) {
        out[, (paste(prefix, class_name, sep = "_")) := 0]
      }
    }
    out[, `:=`(
      branch_covariate_any_training_rows = FALSE,
      branch_covariate_any_count_aware_support = FALSE,
      branch_covariate_max_abs_delta = 0,
      branch_covariate_total_strict_rows = 0,
      branch_covariate_total_evaluable_rows = 0
    )]
    return(out)
  }

  dt <- copy(branch_dt)
  dt <- dt[
    branch_feature_class %chin% classes$branch_feature_class &
      count_aware_evaluable == TRUE
  ]
  if (!nrow(dt)) {
    return(build_branch_covariates(model_dt, data.table()))
  }

  dt[, `:=`(
    delta = safe_num(count_aware_usage_delta),
    weight = pmax(1, safe_num(count_aware_row_weight)),
    q = safe_num(count_aware_usage_delta_q)
  )]
  dt[!is.finite(delta), delta := 0]
  dt[!is.finite(weight) | weight <= 0, weight := 1]
  dt[, sig_row := is.finite(q) & q <= 0.10 & abs(delta) >= 0.04]
  dt[, strict_row := count_aware_strict == TRUE]

  by_class <- c("gene_symbol", "tx_id", "design_family",
                "branch_feature_class")
  total <- dt[, .(
    total_weight = sum(weight),
    total_weighted_delta = sum(weight * delta),
    total_weighted_abs_delta = sum(weight * abs(delta)),
    total_eval_rows = .N,
    total_strict_rows = sum(strict_row, na.rm = TRUE),
    total_sig_rows = sum(sig_row, na.rm = TRUE),
    total_studies = uniqueN(study)
  ), by = by_class]
  heldout <- dt[, .(
    heldout_weight = sum(weight),
    heldout_weighted_delta = sum(weight * delta),
    heldout_weighted_abs_delta = sum(weight * abs(delta)),
    heldout_eval_rows = .N,
    heldout_strict_rows = sum(strict_row, na.rm = TRUE),
    heldout_sig_rows = sum(sig_row, na.rm = TRUE)
  ), by = c(by_class, "heldout_study" = "study")]

  keys[, tmp_join_key := 1L]
  classes[, tmp_join_key := 1L]
  grid <- merge(keys, classes, by = "tmp_join_key",
                allow.cartesian = TRUE)
  grid[, tmp_join_key := NULL]
  classes[, tmp_join_key := NULL]

  grid <- merge(grid, total, by = by_class, all.x = TRUE, sort = FALSE)
  grid <- merge(grid, heldout, by = c(by_class, "heldout_study"),
                all.x = TRUE, sort = FALSE)

  for (column in setdiff(names(grid), c(by_class, "heldout_study"))) {
    if (is.numeric(grid[[column]]) || is.integer(grid[[column]])) {
      grid[is.na(get(column)), (column) := 0]
    }
  }

  grid[, `:=`(
    train_weight = pmax(total_weight - heldout_weight, 0),
    train_weighted_delta =
      total_weighted_delta - heldout_weighted_delta,
    train_weighted_abs_delta =
      total_weighted_abs_delta - heldout_weighted_abs_delta,
    train_eval_rows = pmax(total_eval_rows - heldout_eval_rows, 0),
    train_strict_rows = pmax(total_strict_rows - heldout_strict_rows, 0),
    train_sig_rows = pmax(total_sig_rows - heldout_sig_rows, 0),
    train_studies =
      pmax(total_studies - as.integer(heldout_eval_rows > 0), 0)
  )]
  grid[, train_delta := fifelse(
    train_weight > 0,
    train_weighted_delta / train_weight,
    0
  )]
  grid[, train_abs_delta := fifelse(
    train_weight > 0,
    train_weighted_abs_delta / train_weight,
    0
  )]

  wide <- dcast(
    grid,
    gene_symbol + tx_id + design_family + heldout_study ~
      branch_feature_class,
    value.var = branch_class_columns,
    fill = 0
  )
  value_cols <- setdiff(names(wide), key_cols)
  for (column in value_cols) {
    wide[!is.finite(get(column)), (column) := 0]
  }

  wide[, branch_covariate_any_training_rows :=
         rowSums(.SD) > 0,
       .SDcols = grep("^train_eval_rows_", names(wide), value = TRUE)]
  wide[, branch_covariate_any_count_aware_support :=
         rowSums(.SD) > 0,
       .SDcols = grep("^train_sig_rows_", names(wide), value = TRUE)]
  wide[, branch_covariate_max_abs_delta :=
         do.call(pmax, c(.SD, na.rm = TRUE)),
       .SDcols = grep("^train_abs_delta_", names(wide), value = TRUE)]
  wide[, branch_covariate_total_strict_rows :=
         rowSums(.SD),
       .SDcols = grep("^train_strict_rows_", names(wide), value = TRUE)]
  wide[, branch_covariate_total_evaluable_rows :=
         rowSums(.SD),
       .SDcols = grep("^train_eval_rows_", names(wide), value = TRUE)]
  wide[]
}

fit_blocked_residual_model <- function(model_dt, formula, label) {
  x_all <- Matrix::sparse.model.matrix(formula, data = model_dt)
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

  predict_for_lambda <- function(lambda_value) {
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
    data.table(residual_offset_log2fc = offset, train_rows = train_rows)
  }

  lambdas <- data.table(
    model = label,
    lambda_name = c("lambda_min", "lambda_1se"),
    lambda_value = c(cvfit$lambda.min, cvfit$lambda.1se)
  )
  lambdas <- unique(lambdas[is.finite(lambda_value)], by = "lambda_value")
  if (!nrow(lambdas)) stop("No finite lambdas for model: ", label,
                           call. = FALSE)

  pred_list <- vector("list", nrow(lambdas))
  lambda_selection <- rbindlist(lapply(seq_len(nrow(lambdas)), function(i) {
    pred <- predict_for_lambda(lambdas$lambda_value[[i]])
    pred_list[[i]] <<- pred
    predicted <- model_dt$predicted_log2fc + pred$residual_offset_log2fc
    err <- model_dt$heldout_log2fc - predicted
    data.table(
      model = label,
      lambda_name = lambdas$lambda_name[[i]],
      lambda_value = lambdas$lambda_value[[i]],
      mae_log2fc = mean_safe(abs(err)),
      rmse_log2fc = rmse_safe(err),
      q95_abs_error_log2fc = quantile_safe(abs(err), 0.95),
      direction_accuracy = mean(sign(model_dt$heldout_log2fc) ==
                                  sign(predicted), na.rm = TRUE),
      mean_abs_offset_log2fc = mean_safe(abs(pred$residual_offset_log2fc))
    )
  }), fill = TRUE)
  setorder(lambda_selection, rmse_log2fc, mae_log2fc, q95_abs_error_log2fc)
  chosen <- lambda_selection[1]
  pred <- pred_list[[match(chosen$lambda_value, lambdas$lambda_value)]]

  final_fit <- glmnet::glmnet(
    x_all,
    y_all,
    family = "gaussian",
    alpha = 0,
    lambda = chosen$lambda_value,
    standardize = TRUE,
    intercept = TRUE
  )
  coef_mat <- as.matrix(stats::coef(final_fit, s = chosen$lambda_value))
  coef_dt <- data.table(
    model = label,
    feature = rownames(coef_mat),
    coefficient = as.numeric(coef_mat[, 1])
  )
  coef_dt <- coef_dt[feature != "(Intercept)"]
  coef_dt[, abs_coefficient := abs(coefficient)]
  setorder(coef_dt, -abs_coefficient)

  list(
    prediction = pred,
    lambda_selection = lambda_selection,
    coefficient = coef_dt,
    metadata = data.table(
      model = label,
      lambda_used = chosen$lambda_value,
      lambda_used_name = chosen$lambda_name,
      n_rows = nrow(model_dt),
      n_features = ncol(x_all),
      blocked_cv_groups = length(study_groups)
    )
  )
}

message("Hierarchical branch-covariate residual model:")
message("  1. Build leakage-controlled branch covariates excluding held-out study")
message("  2. Fit metadata-only and metadata+branch residual models")
message("  3. Compare held-out CDS-buffering prediction error")

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required.", call. = FALSE)
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required.", call. = FALSE)
}

loo_file <- file.path(
  analysis_dir, "dominant_rdg_hierarchical_residual_audit",
  "dominant_rdg_hierarchical_residual_audit_loo_context.csv"
)
branch_file <- file.path(
  analysis_dir, "dominant_rdg_branch_allocation",
  "dominant_rdg_branch_allocation_posterior_contrasts.csv"
)

loo <- read_dt(loo_file, required = TRUE)
branch_dt <- read_dt(branch_file, required = TRUE)

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
numeric_base <- c("predicted_log2fc_num", "abs_predicted_log2fc",
                  "log_matched_rows", "log_strict_rows",
                  "log_case_counts", "log_control_counts")
for (column in numeric_base) {
  model_dt[!is.finite(get(column)), (column) := 0]
}

keep <- is.finite(model_dt$heldout_log2fc) &
  is.finite(model_dt$predicted_log2fc) &
  is.finite(model_dt$prediction_error_log2fc)
model_dt <- model_dt[keep]
if (!nrow(model_dt)) stop("No modelable LOO rows.", call. = FALSE)

branch_cov <- build_branch_covariates(model_dt, branch_dt)
cov_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariates_training_covariates.csv"
)
fwrite(branch_cov, cov_file)

model_dt <- merge(
  model_dt,
  branch_cov,
  by = c("gene_symbol", "tx_id", "design_family", "heldout_study"),
  all.x = TRUE,
  sort = FALSE
)
branch_covariate_cols <- setdiff(names(branch_cov),
                                 c("gene_symbol", "tx_id",
                                   "design_family", "heldout_study"))
for (column in branch_covariate_cols) {
  if (is.logical(model_dt[[column]])) {
    model_dt[is.na(get(column)), (column) := FALSE]
    model_dt[, (column) := as.numeric(get(column))]
  } else {
    model_dt[is.na(get(column)) | !is.finite(safe_num(get(column))),
             (column) := 0]
  }
}

residual_cap <- 3
model_dt[, residual_target := pmax(pmin(prediction_error_log2fc,
                                        residual_cap), -residual_cap)]

base_terms <- paste(
  "0 + f_gene + f_design + f_cell + f_tissue + f_condition",
  "f_gene:f_design",
  "f_design:f_cell + f_design:f_tissue + f_design:f_condition",
  "f_gene:f_cell + f_gene:f_tissue",
  paste(numeric_base, collapse = " + "),
  sep = " + "
)
base_formula <- stats::as.formula(paste("~", base_terms))
branch_formula <- stats::as.formula(paste(
  "~", base_terms, "+", paste(branch_covariate_cols, collapse = " + ")
))

metadata_fit <- fit_blocked_residual_model(model_dt, base_formula,
                                           "metadata_context")
branch_fit <- fit_blocked_residual_model(model_dt, branch_formula,
                                         "metadata_context_plus_branch")

model_dt[, `:=`(
  metadata_model_residual_offset_log2fc =
    metadata_fit$prediction$residual_offset_log2fc,
  metadata_model_train_rows = metadata_fit$prediction$train_rows,
  branch_model_residual_offset_log2fc =
    branch_fit$prediction$residual_offset_log2fc,
  branch_model_train_rows = branch_fit$prediction$train_rows
)]
model_dt[, metadata_model_predicted_log2fc :=
           predicted_log2fc + metadata_model_residual_offset_log2fc]
model_dt[, branch_model_predicted_log2fc :=
           predicted_log2fc + branch_model_residual_offset_log2fc]
model_dt[, `:=`(
  metadata_model_error_log2fc =
    heldout_log2fc - metadata_model_predicted_log2fc,
  branch_model_error_log2fc =
    heldout_log2fc - branch_model_predicted_log2fc
)]
model_dt[, `:=`(
  metadata_model_abs_error_log2fc = abs(metadata_model_error_log2fc),
  branch_model_abs_error_log2fc = abs(branch_model_error_log2fc)
)]
model_dt[, `:=`(
  metadata_model_error_delta_log2fc =
    prediction_abs_error_log2fc - metadata_model_abs_error_log2fc,
  branch_model_error_delta_log2fc =
    prediction_abs_error_log2fc - branch_model_abs_error_log2fc,
  branch_vs_metadata_abs_error_delta_log2fc =
    metadata_model_abs_error_log2fc - branch_model_abs_error_log2fc
)]
model_dt[, branch_model_improved_vs_metadata :=
           branch_model_abs_error_log2fc < metadata_model_abs_error_log2fc]
model_dt[, branch_model_correction_class := fifelse(
  branch_vs_metadata_abs_error_delta_log2fc > 0.03,
  "branch_covariates_improve",
  fifelse(branch_vs_metadata_abs_error_delta_log2fc < -0.03,
          "branch_covariates_worse",
          "branch_covariates_neutral")
)]

summary <- rbindlist(c(
  list(metric_row(model_dt, "all")),
  lapply(split(model_dt, model_dt$design_family), function(x) {
    metric_row(x, paste0("design_family:", x$design_family[[1]]))
  })
), fill = TRUE)
summary[, `:=`(
  metadata_vs_original_mae_delta =
    original_mae_log2fc - metadata_model_mae_log2fc,
  branch_vs_original_mae_delta =
    original_mae_log2fc - branch_model_mae_log2fc,
  branch_vs_metadata_mae_delta =
    metadata_model_mae_log2fc - branch_model_mae_log2fc,
  metadata_vs_original_rmse_delta =
    original_rmse_log2fc - metadata_model_rmse_log2fc,
  branch_vs_original_rmse_delta =
    original_rmse_log2fc - branch_model_rmse_log2fc,
  branch_vs_metadata_rmse_delta =
    metadata_model_rmse_log2fc - branch_model_rmse_log2fc,
  branch_vs_metadata_q95_delta =
    metadata_model_q95_abs_error_log2fc -
    branch_model_q95_abs_error_log2fc,
  branch_vs_metadata_direction_accuracy_delta =
    branch_model_direction_accuracy -
    metadata_model_direction_accuracy
)]

gene_design <- model_dt[, {
  worst_idx <- which.max(branch_model_abs_error_log2fc)
  if (!length(worst_idx) || is.na(worst_idx)) worst_idx <- 1L
  data.table(
    branch_covariate_predictions = .N,
    original_mae_log2fc = mean_safe(prediction_abs_error_log2fc),
    metadata_model_mae_log2fc =
      mean_safe(metadata_model_abs_error_log2fc),
    branch_model_mae_log2fc =
      mean_safe(branch_model_abs_error_log2fc),
    original_rmse_log2fc = rmse_safe(prediction_error_log2fc),
    metadata_model_rmse_log2fc =
      rmse_safe(metadata_model_error_log2fc),
    branch_model_rmse_log2fc =
      rmse_safe(branch_model_error_log2fc),
    branch_vs_metadata_mae_delta =
      mean_safe(metadata_model_abs_error_log2fc) -
      mean_safe(branch_model_abs_error_log2fc),
    branch_vs_metadata_rmse_delta =
      rmse_safe(metadata_model_error_log2fc) -
      rmse_safe(branch_model_error_log2fc),
    branch_vs_metadata_q95_delta =
      quantile_safe(metadata_model_abs_error_log2fc, 0.95) -
      quantile_safe(branch_model_abs_error_log2fc, 0.95),
    branch_fraction_improved_vs_metadata =
      mean(branch_model_improved_vs_metadata, na.rm = TRUE),
    branch_covariate_any_training_rows =
      any(branch_covariate_any_training_rows > 0, na.rm = TRUE),
    branch_covariate_any_count_aware_support =
      any(branch_covariate_any_count_aware_support > 0, na.rm = TRUE),
    branch_covariate_max_abs_delta =
      max_safe(branch_covariate_max_abs_delta),
    branch_covariate_total_strict_rows =
      max_safe(branch_covariate_total_strict_rows),
    branch_covariate_total_evaluable_rows =
      max_safe(branch_covariate_total_evaluable_rows),
    worst_branch_model_context =
      paste(heldout_study[worst_idx], cell_lines[worst_idx],
            tissues[worst_idx], condition_tokens[worst_idx], sep = " | "),
    top_branch_improved_contexts = collapse_unique(
      paste(heldout_study, cell_lines, tissues, condition_tokens, sep = " | ")
        [order(-branch_vs_metadata_abs_error_delta_log2fc)],
      5L
    )
  )
}, by = .(gene_symbol, tx_id, design_family)]
gene_design[, branch_covariate_model_status := fifelse(
  branch_vs_metadata_rmse_delta > 0.03 &
    branch_vs_metadata_mae_delta >= -0.01,
  "branch_covariates_reduce_residual_error",
  fifelse(branch_vs_metadata_rmse_delta < -0.03 |
            branch_vs_metadata_mae_delta < -0.03,
          "branch_covariates_increase_residual_error",
          "branch_covariates_neutral")
)]
setorder(gene_design, -branch_vs_metadata_rmse_delta,
         -branch_vs_metadata_mae_delta, gene_symbol, design_family)

summary_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_summary.csv"
)
loo_file_out <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_loo.csv"
)
gene_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_gene_design.csv"
)
viral_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_viral_gene_design.csv"
)
coef_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_coefficients.csv"
)
lambda_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_lambda_selection.csv"
)
metadata_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_branch_covariate_model_metadata.csv"
)

fwrite(summary, summary_file)
fwrite(model_dt, loo_file_out)
fwrite(gene_design, gene_file)
fwrite(gene_design[design_family == "viral_infection"], viral_file)
fwrite(rbindlist(list(metadata_fit$coefficient, branch_fit$coefficient),
                 fill = TRUE),
       coef_file)
fwrite(rbindlist(list(metadata_fit$lambda_selection,
                     branch_fit$lambda_selection), fill = TRUE),
       lambda_file)
fwrite(rbindlist(list(metadata_fit$metadata, branch_fit$metadata),
                 fill = TRUE),
       metadata_file)

if (nrow(summary)) {
  plot_dt <- summary[grepl("^design_family:", scope)]
  plot_dt[, design_family := sub("^design_family:", "", scope)]
  plot_dt <- melt(
    plot_dt,
    id.vars = "design_family",
    measure.vars = c("original_rmse_log2fc",
                     "metadata_model_rmse_log2fc",
                     "branch_model_rmse_log2fc"),
    variable.name = "model",
    value.name = "rmse_log2fc"
  )
  plot_dt[, model := fifelse(
    model == "original_rmse_log2fc",
    "hierarchical",
    fifelse(model == "metadata_model_rmse_log2fc",
            "metadata residual model",
            "metadata + branch residual model")
  )]
  p <- ggplot(
    plot_dt,
    aes(x = reorder(design_family, rmse_log2fc, FUN = max),
        y = rmse_log2fc, fill = model,
        text = paste0(
          "Design: ", design_family,
          "<br>Model: ", model,
          "<br>RMSE: ", round(rmse_log2fc, 3)
        ))
  ) +
    geom_col(position = position_dodge(width = 0.76), width = 0.68) +
    coord_flip() +
    scale_fill_manual(values = c(
      "hierarchical" = "#6b7280",
      "metadata residual model" = "#047857",
      "metadata + branch residual model" = "#5b2a7a"
    )) +
    labs(
      title = "Do RDG Branch Covariates Improve CDS-Buffering Prediction?",
      subtitle = "Leave-one-study-out residual correction; branch covariates exclude held-out study",
      x = NULL, y = "RMSE (log2FC)", fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "branch_covariate_model_rmse.png"),
         p, width = 10.2, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "branch_covariate_model_rmse.pdf"),
         p, width = 10.2, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "branch_covariate_model_rmse.html"),
      title = "RDG branch covariate model RMSE"
    )
  }
}

viral_gene <- gene_design[design_family == "viral_infection"]
if (nrow(viral_gene)) {
  viral_plot <- viral_gene[
    order(-abs(branch_vs_metadata_rmse_delta))
  ][seq_len(min(.N, 30))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  viral_plot[, branch_covariate_direction_class := fifelse(
    branch_vs_metadata_rmse_delta >= 0.001,
    "small_improvement",
    fifelse(branch_vs_metadata_rmse_delta <= -0.001,
            "small_worse",
            "near_zero")
  )]
  p2 <- ggplot(
    viral_plot,
    aes(x = branch_vs_metadata_rmse_delta, y = gene_symbol,
        fill = branch_covariate_direction_class,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Status: ", branch_covariate_model_status,
          "<br>RMSE delta: ", round(branch_vs_metadata_rmse_delta, 3),
          "<br>MAE delta: ", round(branch_vs_metadata_mae_delta, 3),
          "<br>Branch training rows: ",
          branch_covariate_total_evaluable_rows,
          "<br>Branch max abs delta: ",
          round(branch_covariate_max_abs_delta, 3),
          "<br>Worst context: ", worst_branch_model_context
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(
      "small_improvement" = "#0f766e",
      "small_worse" = "#b91c1c",
      "near_zero" = "#6b7280"
    ), labels = c(
      "small_improvement" = "small improvement",
      "small_worse" = "small worse",
      "near_zero" = "near zero"
    )) +
    labs(
      title = "Viral-Infection Branch-Covariate Residual Model by Gene",
      subtitle = "Positive values mean branch covariates improve over metadata-only residual correction",
      x = "RMSE reduction (log2FC)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "branch_covariate_model_viral_genes.png"),
         p2, width = 9.6, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "branch_covariate_model_viral_genes.pdf"),
         p2, width = 9.6, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir, "branch_covariate_model_viral_genes.html"),
      title = "Viral branch-covariate residual model"
    )
  }
}

message("Saved hierarchical branch-covariate model outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
