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
output_dir <- file.path(analysis_dir, "dominant_rdg_mechanistic_motif_model")
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

finite_or_zero <- function(x) {
  x <- safe_num(x)
  x[!is.finite(x)] <- 0
  x
}

clip01 <- function(x) pmin(pmax(finite_or_zero(x), 0), 1)

scaled_log <- function(x, scale) {
  clip01(log1p(pmax(finite_or_zero(x), 0)) / log1p(scale))
}

bell_score <- function(x, center, width = 0.8) {
  x <- pmax(finite_or_zero(x), 0)
  exp(-((log1p(x) - log1p(center))^2) / (2 * width^2))
}

mean_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
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

clean_factor <- function(x, prefix = "unknown") {
  x <- tolower(as.character(x))
  x[is.na(x) | !nzchar(x)] <- prefix
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- prefix
  factor(x)
}

collapse_unique <- function(x, n = 6L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

metric_row <- function(dt, scope, model_prefixes) {
  out <- data.table(
    scope = scope,
    n_predictions = nrow(dt),
    n_genes = uniqueN(dt$gene_symbol),
    n_studies = uniqueN(dt$heldout_study),
    hierarchical_mae_log2fc = mean_safe(dt$hierarchical_abs_error_log2fc),
    hierarchical_rmse_log2fc = rmse_safe(dt$hierarchical_error_log2fc),
    hierarchical_q95_abs_error_log2fc =
      quantile_safe(dt$hierarchical_abs_error_log2fc, 0.95),
    hierarchical_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$predicted_log2fc),
           na.rm = TRUE)
  )
  for (prefix in model_prefixes) {
    pred_col <- paste0(prefix, "_predicted_log2fc")
    err_col <- paste0(prefix, "_error_log2fc")
    abs_col <- paste0(prefix, "_abs_error_log2fc")
    out[, (paste0(prefix, "_mae_log2fc")) := mean_safe(dt[[abs_col]])]
    out[, (paste0(prefix, "_rmse_log2fc")) := rmse_safe(dt[[err_col]])]
    out[, (paste0(prefix, "_q95_abs_error_log2fc")) :=
          quantile_safe(dt[[abs_col]], 0.95)]
    out[, (paste0(prefix, "_direction_accuracy")) :=
          mean(sign(dt$heldout_log2fc) == sign(dt[[pred_col]]),
               na.rm = TRUE)]
  }
  out
}

fit_blocked_direct_model <- function(model_dt, formula, label,
                                     block_col = "gene_symbol") {
  x_all <- Matrix::sparse.model.matrix(formula, data = model_dt)
  y_all <- finite_or_zero(model_dt$heldout_log2fc)
  blocks <- sort(unique(as.character(model_dt[[block_col]])))
  fold_id <- match(as.character(model_dt[[block_col]]), blocks)

  set.seed(3)
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
    pred <- rep(NA_real_, nrow(model_dt))
    train_rows <- rep(NA_integer_, nrow(model_dt))
    for (block_id in blocks) {
      test_idx <- which(as.character(model_dt[[block_col]]) == block_id)
      train_idx <- which(as.character(model_dt[[block_col]]) != block_id)
      if (!length(test_idx) || length(train_idx) < 100) next
      fit <- glmnet::glmnet(
        x_all[train_idx, , drop = FALSE],
        y_all[train_idx],
        family = "gaussian",
        alpha = 0,
        lambda = lambda_value,
        standardize = TRUE,
        intercept = TRUE
      )
      pred[test_idx] <- as.numeric(stats::predict(
        fit,
        newx = x_all[test_idx, , drop = FALSE],
        s = lambda_value
      ))
      train_rows[test_idx] <- length(train_idx)
    }
    pred[!is.finite(pred)] <- mean_safe(y_all)
    data.table(predicted_log2fc = pred, train_rows = train_rows)
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
    err <- y_all - pred$predicted_log2fc
    pred_list[[i]] <<- pred
    data.table(
      model = label,
      lambda_name = lambdas$lambda_name[[i]],
      lambda_value = lambdas$lambda_value[[i]],
      mae_log2fc = mean_safe(abs(err)),
      rmse_log2fc = rmse_safe(err),
      q95_abs_error_log2fc = quantile_safe(abs(err), 0.95),
      direction_accuracy = mean(sign(y_all) == sign(pred$predicted_log2fc),
                                na.rm = TRUE)
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
      blocked_cv = block_col,
      blocked_cv_groups = length(blocks)
    )
  )
}

message("RDG mechanistic motif model:")
message("  1. Collapse raw graph geometry into constrained motif scores")
message("  2. Fit leave-one-gene-out motif main-effect and motif-design models")
message("  3. Compare against context-only prediction without gene identity")

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
geometry_covariate_file <- file.path(
  analysis_dir, "dominant_rdg_graph_geometry_model",
  "dominant_rdg_graph_geometry_model_covariates.csv"
)
rescue_file <- file.path(
  analysis_dir, "dominant_downstream_tis_rescue",
  "downstream_tis_rescue_feature_summary.csv"
)

loo <- read_dt(loo_file, required = TRUE)
geometry_cov <- read_dt(geometry_covariate_file, required = TRUE)
rescue <- read_dt(rescue_file)

base_gene_cols <- c(
  "gene_symbol", "tx_id", "n_uorfs", "n_leader_uorfs",
  "n_overlapping_uorfs", "n_overlapping_or_nested_pairs",
  "n_abutting_pairs", "n_close_pairs_10bp", "n_close_pairs_30bp",
  "n_close_pairs_60bp", "leader_length", "clean_cds_bases",
  "uorf_union_bases", "leader_uorf_union_bases",
  "overlapping_uorf_union_bases", "uorf_union_fraction_of_leader",
  "leader_uorf_union_fraction_of_leader", "longest_uorf_bases",
  "mean_uorf_bases", "first_uorf_bases", "first_uorf_to_cds_gap",
  "min_leader_stop_to_clean_cds_gap",
  "n_leader_uorfs_stop_within_30bp_cds",
  "n_leader_uorfs_stop_within_60bp_cds",
  "n_leader_uorfs_stop_within_100bp_cds",
  "last_leader_stop_to_clean_cds_gap",
  "last_any_stop_to_clean_cds_gap", "max_cds_overlap_bases",
  "total_cds_overlap_bases_double_counted",
  "max_cds_overlap_fraction_of_clean_cds",
  "annotation_atg_uorfs", "annotation_near_cognate_uorfs",
  "annotation_strong_kozak_uorfs", "annotation_moderate_kozak_uorfs",
  "annotation_mean_kozak_score", "annotation_max_kozak_score",
  "annotation_mean_feature_bases", "annotation_max_feature_bases",
  "annotation_first_uorf_kozak_score", "annotation_first_uorf_is_atg",
  "annotation_first_uorf_is_near_cognate",
  "uorf_count_bin", "overlap_class", "spacing_class",
  "last_leader_gap_class", "first_uorf_length_class",
  "has_overlapping_uorf", "has_modeled_uorf"
)
geometry <- unique(geometry_cov[, intersect(base_gene_cols, names(geometry_cov)),
                                with = FALSE])
for (column in setdiff(names(geometry), c("gene_symbol", "tx_id",
                                          "uorf_count_bin", "overlap_class",
                                          "spacing_class",
                                          "last_leader_gap_class",
                                          "first_uorf_length_class",
                                          "has_overlapping_uorf",
                                          "has_modeled_uorf"))) {
  geometry[, (column) := finite_or_zero(get(column))]
}

if (nrow(rescue)) {
  rescue_summary <- rescue[, .(
    rescue_best_evidence_score =
      max(finite_or_zero(best_downstream_tis_evidence_score), na.rm = TRUE),
    rescue_best_review_priority =
      max(finite_or_zero(review_priority_score), na.rm = TRUE),
    rescue_frame_preserving_starts =
      max(finite_or_zero(n_frame_preserving_start_codons), na.rm = TRUE),
    rescue_sequence_shape_supported =
      max(finite_or_zero(n_sequence_and_shape_supported), na.rm = TRUE),
    rescue_best_distance_from_ouorf_stop =
      max(finite_or_zero(best_distance_from_ouorf_stop), na.rm = TRUE),
    rescue_best_post_vs_pre_log2 =
      max(finite_or_zero(best_total_post_vs_pre_density_log2), na.rm = TRUE)
  ), by = .(gene_symbol, tx_id)]
} else {
  rescue_summary <- data.table(gene_symbol = character(), tx_id = character())
}

motif <- merge(geometry, rescue_summary, by = c("gene_symbol", "tx_id"),
               all.x = TRUE, sort = FALSE)
rescue_cols <- setdiff(names(rescue_summary), c("gene_symbol", "tx_id"))
for (column in rescue_cols) {
  motif[is.na(get(column)) | !is.finite(finite_or_zero(get(column))),
        (column) := 0]
}

motif[, inhibitory_overlap_burden := clip01(
  0.30 * scaled_log(n_overlapping_uorfs, 4) +
    0.30 * clip01(max_cds_overlap_fraction_of_clean_cds / 0.35) +
    0.20 * scaled_log(max_cds_overlap_bases, 300) +
    0.20 * clip01(overlapping_uorf_union_bases /
                    pmax(clean_cds_bases, 1))
)]
motif[, reinitiation_distance_opportunity := fifelse(
  n_leader_uorfs > 0,
  clip01(
    0.45 * bell_score(last_leader_stop_to_clean_cds_gap, 120, 0.85) +
      0.25 * bell_score(min_leader_stop_to_clean_cds_gap, 80, 0.75) +
      0.30 * clip01(n_leader_uorfs_stop_within_100bp_cds /
                      pmax(n_leader_uorfs, 1))
  ),
  0
)]
motif[, dense_uorf_collision_burden := clip01(
  0.35 * scaled_log(n_uorfs, 10) +
    0.30 * scaled_log(n_close_pairs_30bp, 8) +
    0.20 * clip01(uorf_union_fraction_of_leader) +
    0.15 * scaled_log(n_overlapping_or_nested_pairs, 4)
)]
motif[, strong_start_leakage_burden := clip01(
  0.25 * scaled_log(annotation_atg_uorfs, 4) +
    0.25 * scaled_log(annotation_strong_kozak_uorfs, 4) +
    0.20 * clip01(annotation_mean_kozak_score / 2) +
    0.15 * clip01(annotation_first_uorf_kozak_score / 2) +
    0.15 * clip01(annotation_first_uorf_is_atg)
)]
motif[, downstream_rescue_opportunity := clip01(
  clip01(inhibitory_overlap_burden + 0.25) *
    (
      0.35 * clip01(rescue_best_evidence_score / 15) +
        0.25 * scaled_log(rescue_frame_preserving_starts, 50) +
        0.20 * scaled_log(rescue_sequence_shape_supported, 4) +
        0.20 * bell_score(rescue_best_distance_from_ouorf_stop, 180, 1.1)
    )
)]
motif[, overlap_without_rescue := clip01(
  inhibitory_overlap_burden * (1 - downstream_rescue_opportunity)
)]
motif[, leakage_with_reinitiation := clip01(
  strong_start_leakage_burden * reinitiation_distance_opportunity
)]
motif[, dense_short_reinitiation := clip01(
  dense_uorf_collision_burden * (1 - reinitiation_distance_opportunity)
)]
motif[, rescue_after_overlap := clip01(
  inhibitory_overlap_burden * downstream_rescue_opportunity
)]
motif[, motif_complexity_index := clip01(
  0.25 * inhibitory_overlap_burden +
    0.25 * dense_uorf_collision_burden +
    0.20 * strong_start_leakage_burden +
    0.15 * reinitiation_distance_opportunity +
    0.15 * downstream_rescue_opportunity
)]

motif_score_cols <- c(
  "inhibitory_overlap_burden",
  "reinitiation_distance_opportunity",
  "dense_uorf_collision_burden",
  "strong_start_leakage_burden",
  "downstream_rescue_opportunity",
  "overlap_without_rescue",
  "leakage_with_reinitiation",
  "dense_short_reinitiation",
  "rescue_after_overlap",
  "motif_complexity_index"
)
motif_file <- file.path(
  output_dir,
  "dominant_rdg_mechanistic_motif_scores.csv"
)
fwrite(motif[, c("gene_symbol", "tx_id", motif_score_cols,
                 intersect(c("n_uorfs", "n_overlapping_uorfs",
                             "uorf_count_bin", "overlap_class",
                             "spacing_class",
                             "last_leader_gap_class"),
                           names(motif))),
            with = FALSE],
       motif_file)

required <- c(
  "gene_symbol", "tx_id", "design_family", "heldout_study",
  "heldout_log2fc", "predicted_log2fc", "cell_lines", "tissues",
  "condition_tokens"
)
missing <- setdiff(required, names(loo))
if (length(missing)) {
  stop("LOO context table is missing required columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

model_dt <- merge(loo, motif[, c("gene_symbol", "tx_id", motif_score_cols),
                             with = FALSE],
                  by = c("gene_symbol", "tx_id"), all.x = TRUE,
                  sort = FALSE)
for (column in motif_score_cols) {
  model_dt[is.na(get(column)) | !is.finite(finite_or_zero(get(column))),
           (column) := 0]
}
for (column in c("gene_symbol", "design_family", "heldout_study",
                 "cell_lines", "tissues", "condition_tokens")) {
  model_dt[is.na(get(column)) | !nzchar(as.character(get(column))),
           (column) := "unknown"]
}
model_dt <- model_dt[
  is.finite(finite_or_zero(heldout_log2fc)) &
    is.finite(finite_or_zero(predicted_log2fc))
]
if (!nrow(model_dt)) stop("No rows with motif scores.", call. = FALSE)

model_dt[, `:=`(
  f_design = clean_factor(design_family, "design_unknown"),
  f_cell = clean_factor(cell_lines, "cell_unknown"),
  f_tissue = clean_factor(tissues, "tissue_unknown"),
  f_condition = clean_factor(condition_tokens, "condition_unknown"),
  log_matched_rows = log1p(pmax(finite_or_zero(matched_rows_in_context), 0)),
  log_strict_rows = log1p(pmax(finite_or_zero(strict_rows_in_context), 0)),
  log_case_counts = log1p(pmax(finite_or_zero(total_case_counts), 0)),
  log_control_counts = log1p(pmax(finite_or_zero(total_control_counts), 0)),
  hierarchical_error_log2fc = heldout_log2fc - predicted_log2fc,
  hierarchical_abs_error_log2fc = abs(heldout_log2fc - predicted_log2fc)
)]

base_numeric_cols <- c("log_matched_rows", "log_strict_rows",
                       "log_case_counts", "log_control_counts")
context_terms <- paste(
  "0 + f_design + f_cell + f_tissue + f_condition",
  paste(base_numeric_cols, collapse = " + "),
  sep = " + "
)
motif_terms <- paste(context_terms,
                     paste(motif_score_cols, collapse = " + "),
                     sep = " + ")
motif_design_terms <- paste(
  motif_terms,
  paste(paste0("f_design:", motif_score_cols), collapse = " + "),
  sep = " + "
)

context_fit <- fit_blocked_direct_model(
  model_dt, stats::as.formula(paste("~", context_terms)),
  "context_nogene"
)
motif_fit <- fit_blocked_direct_model(
  model_dt, stats::as.formula(paste("~", motif_terms)),
  "context_plus_motifs"
)
motif_design_fit <- fit_blocked_direct_model(
  model_dt, stats::as.formula(paste("~", motif_design_terms)),
  "context_plus_motif_design"
)

fit_list <- list(
  context_nogene = context_fit,
  motif_model = motif_fit,
  motif_design_model = motif_design_fit
)
for (prefix in names(fit_list)) {
  pred <- fit_list[[prefix]]$prediction
  model_dt[, (paste0(prefix, "_predicted_log2fc")) :=
             pred$predicted_log2fc]
  model_dt[, (paste0(prefix, "_train_rows")) := pred$train_rows]
  model_dt[, (paste0(prefix, "_error_log2fc")) :=
             heldout_log2fc - get(paste0(prefix, "_predicted_log2fc"))]
  model_dt[, (paste0(prefix, "_abs_error_log2fc")) :=
             abs(get(paste0(prefix, "_error_log2fc")))]
  model_dt[, (paste0(prefix, "_vs_context_abs_error_delta_log2fc")) :=
             context_nogene_abs_error_log2fc -
             get(paste0(prefix, "_abs_error_log2fc"))]
}

model_prefixes <- names(fit_list)
summary <- rbindlist(c(
  list(metric_row(model_dt, "all", model_prefixes)),
  lapply(split(model_dt, model_dt$design_family), function(x) {
    metric_row(x, paste0("design_family:", x$design_family[[1]]),
               model_prefixes)
  })
), fill = TRUE)
for (prefix in setdiff(model_prefixes, "context_nogene")) {
  summary[, (paste0(prefix, "_vs_context_mae_delta")) :=
            context_nogene_mae_log2fc - get(paste0(prefix, "_mae_log2fc"))]
  summary[, (paste0(prefix, "_vs_context_rmse_delta")) :=
            context_nogene_rmse_log2fc - get(paste0(prefix, "_rmse_log2fc"))]
  summary[, (paste0(prefix, "_vs_context_q95_delta")) :=
            context_nogene_q95_abs_error_log2fc -
            get(paste0(prefix, "_q95_abs_error_log2fc"))]
  summary[, (paste0(prefix, "_vs_context_direction_accuracy_delta")) :=
            get(paste0(prefix, "_direction_accuracy")) -
            context_nogene_direction_accuracy]
}

all_row <- summary[scope == "all"][1]
candidate_prefixes <- c("motif_model", "motif_design_model")
rmse_values <- unlist(all_row[, paste0(candidate_prefixes, "_rmse_log2fc"),
                              with = FALSE])
best_prefix <- candidate_prefixes[which.min(rmse_values)]
model_dt[, best_motif_model_prefix := best_prefix]
model_dt[, best_motif_predicted_log2fc := get(paste0(best_prefix,
                                                     "_predicted_log2fc"))]
model_dt[, best_motif_error_log2fc := get(paste0(best_prefix,
                                                 "_error_log2fc"))]
model_dt[, best_motif_abs_error_log2fc := get(paste0(best_prefix,
                                                     "_abs_error_log2fc"))]
model_dt[, best_motif_vs_context_abs_error_delta_log2fc :=
           context_nogene_abs_error_log2fc - best_motif_abs_error_log2fc]
model_dt[, best_motif_correction_class := fifelse(
  best_motif_vs_context_abs_error_delta_log2fc > 0.03,
  "motif_model_improves",
  fifelse(best_motif_vs_context_abs_error_delta_log2fc < -0.03,
          "motif_model_worse",
          "motif_model_neutral")
)]
summary[, best_motif_model_prefix := best_prefix]

gene_design <- model_dt[, {
  worst_idx <- which.max(best_motif_abs_error_log2fc)
  if (!length(worst_idx) || is.na(worst_idx)) worst_idx <- 1L
  data.table(
    motif_predictions = .N,
    n_studies = uniqueN(heldout_study),
    hierarchical_mae_log2fc = mean_safe(hierarchical_abs_error_log2fc),
    context_nogene_mae_log2fc =
      mean_safe(context_nogene_abs_error_log2fc),
    motif_model_mae_log2fc = mean_safe(motif_model_abs_error_log2fc),
    motif_design_model_mae_log2fc =
      mean_safe(motif_design_model_abs_error_log2fc),
    best_motif_mae_log2fc = mean_safe(best_motif_abs_error_log2fc),
    hierarchical_rmse_log2fc = rmse_safe(hierarchical_error_log2fc),
    context_nogene_rmse_log2fc =
      rmse_safe(context_nogene_error_log2fc),
    motif_model_rmse_log2fc = rmse_safe(motif_model_error_log2fc),
    motif_design_model_rmse_log2fc =
      rmse_safe(motif_design_model_error_log2fc),
    best_motif_rmse_log2fc = rmse_safe(best_motif_error_log2fc),
    best_motif_vs_context_mae_delta =
      mean_safe(context_nogene_abs_error_log2fc) -
      mean_safe(best_motif_abs_error_log2fc),
    best_motif_vs_context_rmse_delta =
      rmse_safe(context_nogene_error_log2fc) -
      rmse_safe(best_motif_error_log2fc),
    best_motif_fraction_improved_vs_context =
      mean(best_motif_vs_context_abs_error_delta_log2fc > 0,
           na.rm = TRUE),
    worst_best_motif_context =
      paste(heldout_study[worst_idx], cell_lines[worst_idx],
            tissues[worst_idx], condition_tokens[worst_idx], sep = " | "),
    top_motif_improved_contexts = collapse_unique(
      paste(heldout_study, cell_lines, tissues, condition_tokens, sep = " | ")
        [order(-best_motif_vs_context_abs_error_delta_log2fc)],
      5L
    )
  )
}, by = .(gene_symbol, tx_id, design_family)]
gene_design <- merge(gene_design,
                     motif[, c("gene_symbol", "tx_id", motif_score_cols,
                               "n_uorfs", "n_overlapping_uorfs",
                               "uorf_count_bin", "overlap_class"),
                           with = FALSE],
                     by = c("gene_symbol", "tx_id"),
                     all.x = TRUE, sort = FALSE)
gene_design[, best_motif_model_status := fifelse(
  best_motif_vs_context_rmse_delta > 0.03 &
    best_motif_vs_context_mae_delta >= -0.01,
  "motif_model_reduces_error",
  fifelse(best_motif_vs_context_rmse_delta < -0.03 |
            best_motif_vs_context_mae_delta < -0.03,
          "motif_model_increases_error",
          "motif_model_neutral")
)]
setorder(gene_design, -best_motif_vs_context_rmse_delta,
         -best_motif_vs_context_mae_delta, gene_symbol, design_family)

summary_file <- file.path(output_dir,
                          "dominant_rdg_mechanistic_motif_model_summary.csv")
loo_file_out <- file.path(output_dir,
                          "dominant_rdg_mechanistic_motif_model_loo.csv")
gene_file <- file.path(output_dir,
                       "dominant_rdg_mechanistic_motif_model_gene_design.csv")
viral_file <- file.path(output_dir,
                        "dominant_rdg_mechanistic_motif_model_viral_gene_design.csv")
coef_file <- file.path(output_dir,
                       "dominant_rdg_mechanistic_motif_model_coefficients.csv")
lambda_file <- file.path(output_dir,
                         "dominant_rdg_mechanistic_motif_model_lambda_selection.csv")
metadata_file <- file.path(output_dir,
                           "dominant_rdg_mechanistic_motif_model_metadata.csv")

fwrite(summary, summary_file)
fwrite(model_dt, loo_file_out)
fwrite(gene_design, gene_file)
fwrite(gene_design[design_family == "viral_infection"], viral_file)
fwrite(rbindlist(lapply(fit_list, `[[`, "coefficient"), fill = TRUE),
       coef_file)
fwrite(rbindlist(lapply(fit_list, `[[`, "lambda_selection"), fill = TRUE),
       lambda_file)
fwrite(rbindlist(lapply(fit_list, `[[`, "metadata"), fill = TRUE),
       metadata_file)

plot_dt <- summary[grepl("^design_family:", scope)]
if (nrow(plot_dt)) {
  plot_dt[, design_family := sub("^design_family:", "", scope)]
  plot_dt <- melt(
    plot_dt,
    id.vars = "design_family",
    measure.vars = c("context_nogene_rmse_log2fc",
                     "motif_model_rmse_log2fc",
                     "motif_design_model_rmse_log2fc"),
    variable.name = "model",
    value.name = "rmse_log2fc"
  )
  plot_dt[, model := fifelse(
    model == "context_nogene_rmse_log2fc",
    "context, no gene ID",
    fifelse(model == "motif_model_rmse_log2fc",
            "context + motif scores",
            "context + motif-design interactions")
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
    geom_col(position = position_dodge(width = 0.78), width = 0.68) +
    coord_flip() +
    scale_fill_manual(values = c(
      "context, no gene ID" = "#6b7280",
      "context + motif scores" = "#0f766e",
      "context + motif-design interactions" = "#5b2a7a"
    )) +
    labs(
      title = "Do Mechanistic RDG Motif Scores Generalize Across Genes?",
      subtitle = "Leave-one-gene-out prediction; no gene identity covariate",
      x = NULL, y = "RMSE (log2FC)", fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "mechanistic_motif_model_rmse.png"),
         p, width = 10.6, height = 7.4, dpi = 180)
  ggsave(file.path(figure_dir, "mechanistic_motif_model_rmse.pdf"),
         p, width = 10.6, height = 7.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "mechanistic_motif_model_rmse.html"),
      title = "RDG mechanistic motif model RMSE"
    )
  }
}

viral_gene <- gene_design[design_family == "viral_infection"]
if (nrow(viral_gene)) {
  viral_plot <- viral_gene[
    order(-abs(best_motif_vs_context_rmse_delta))
  ][seq_len(min(.N, 30))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  viral_plot[, motif_direction_class := fifelse(
    best_motif_vs_context_rmse_delta >= 0.01,
    "improvement",
    fifelse(best_motif_vs_context_rmse_delta <= -0.01,
            "worse",
            "near_zero")
  )]
  p2 <- ggplot(
    viral_plot,
    aes(x = best_motif_vs_context_rmse_delta, y = gene_symbol,
        fill = motif_direction_class,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Status: ", best_motif_model_status,
          "<br>RMSE delta: ", round(best_motif_vs_context_rmse_delta, 3),
          "<br>MAE delta: ", round(best_motif_vs_context_mae_delta, 3),
          "<br>overlap burden: ", round(inhibitory_overlap_burden, 2),
          "<br>reinitiation opportunity: ",
          round(reinitiation_distance_opportunity, 2),
          "<br>dense-uORF burden: ", round(dense_uorf_collision_burden, 2),
          "<br>rescue opportunity: ",
          round(downstream_rescue_opportunity, 2)
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(
      "improvement" = "#0f766e",
      "worse" = "#b91c1c",
      "near_zero" = "#6b7280"
    ), labels = c(
      "improvement" = "motifs improve",
      "worse" = "motifs worse",
      "near_zero" = "near zero"
    )) +
    labs(
      title = "Viral-Infection Mechanistic Motif Model by Gene",
      subtitle = "Positive values mean motif scores improve over metadata/context without gene ID",
      x = "RMSE reduction (log2FC)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "mechanistic_motif_model_viral_genes.png"),
         p2, width = 9.6, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "mechanistic_motif_model_viral_genes.pdf"),
         p2, width = 9.6, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir, "mechanistic_motif_model_viral_genes.html"),
      title = "Viral mechanistic RDG motif model"
    )
  }
}

coef_plot_dt <- rbindlist(lapply(fit_list, `[[`, "coefficient"),
                           fill = TRUE)
motif_feature_pattern <- paste0(
  "(^|:)(",
  paste(motif_score_cols, collapse = "|"),
  ")$"
)
coef_plot_dt <- coef_plot_dt[
  model %chin% c("context_plus_motifs", "context_plus_motif_design") &
    grepl(motif_feature_pattern, feature, perl = TRUE)
][seq_len(min(.N, 36))]
if (nrow(coef_plot_dt)) {
  coef_plot_dt[, feature := factor(feature, levels = rev(unique(feature)))]
  p3 <- ggplot(
    coef_plot_dt,
    aes(x = coefficient, y = feature,
        fill = coefficient > 0,
        text = paste0(
          "Model: ", model,
          "<br>Feature: ", feature,
          "<br>Coefficient: ", round(coefficient, 4)
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.68) +
    facet_wrap(~model, scales = "free_y") +
    scale_fill_manual(values = c("TRUE" = "#0f766e",
                                 "FALSE" = "#b91c1c"),
                      labels = c("TRUE" = "positive",
                                 "FALSE" = "negative")) +
    labs(
      title = "Mechanistic RDG Motif Coefficients",
      subtitle = "Ridge coefficients from leave-one-gene-out models",
      x = "Coefficient", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "mechanistic_motif_model_coefficients.png"),
         p3, width = 10.2, height = 7.6, dpi = 180)
  ggsave(file.path(figure_dir, "mechanistic_motif_model_coefficients.pdf"),
         p3, width = 10.2, height = 7.6)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p3,
      file.path(figure_dir, "mechanistic_motif_model_coefficients.html"),
      title = "Mechanistic RDG motif coefficients"
    )
  }
}

message("Saved mechanistic motif model outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
