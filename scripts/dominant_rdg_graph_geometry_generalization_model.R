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
output_dir <- file.path(analysis_dir, "dominant_rdg_graph_geometry_model")
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

kozak_score <- function(x) {
  x <- tolower(as.character(x))
  fifelse(x == "strong", 2,
          fifelse(x == "moderate", 1,
                  fifelse(x == "weak", 0, 0)))
}

max_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) 0 else max(x)
}

metric_row <- function(dt, scope) {
  data.table(
    scope = scope,
    n_predictions = nrow(dt),
    n_genes = uniqueN(dt$gene_symbol),
    n_studies = uniqueN(dt$heldout_study),
    hierarchical_mae_log2fc =
      mean_safe(abs(dt$heldout_log2fc - dt$predicted_log2fc)),
    context_nogene_mae_log2fc =
      mean_safe(dt$context_nogene_abs_error_log2fc),
    graph_geometry_mae_log2fc =
      mean_safe(dt$graph_geometry_abs_error_log2fc),
    hierarchical_rmse_log2fc =
      rmse_safe(dt$heldout_log2fc - dt$predicted_log2fc),
    context_nogene_rmse_log2fc =
      rmse_safe(dt$context_nogene_error_log2fc),
    graph_geometry_rmse_log2fc =
      rmse_safe(dt$graph_geometry_error_log2fc),
    hierarchical_q95_abs_error_log2fc =
      quantile_safe(abs(dt$heldout_log2fc - dt$predicted_log2fc), 0.95),
    context_nogene_q95_abs_error_log2fc =
      quantile_safe(dt$context_nogene_abs_error_log2fc, 0.95),
    graph_geometry_q95_abs_error_log2fc =
      quantile_safe(dt$graph_geometry_abs_error_log2fc, 0.95),
    hierarchical_direction_accuracy =
      mean(sign(dt$heldout_log2fc) == sign(dt$predicted_log2fc),
           na.rm = TRUE),
    context_nogene_direction_accuracy =
      mean(sign(dt$heldout_log2fc) ==
             sign(dt$context_nogene_predicted_log2fc), na.rm = TRUE),
    graph_geometry_direction_accuracy =
      mean(sign(dt$heldout_log2fc) ==
             sign(dt$graph_geometry_predicted_log2fc), na.rm = TRUE)
  )
}

fit_blocked_direct_model <- function(model_dt, formula, label,
                                     block_col = "gene_symbol") {
  x_all <- Matrix::sparse.model.matrix(formula, data = model_dt)
  y_all <- safe_num(model_dt$heldout_log2fc)
  blocks <- sort(unique(as.character(model_dt[[block_col]])))
  fold_id <- match(as.character(model_dt[[block_col]]), blocks)

  set.seed(2)
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

message("RDG graph-geometry generalization model:")
message("  1. Build transcript architecture covariates")
message("  2. Fit leave-one-gene-out context-only and context+geometry models")
message("  3. Test whether graph structure generalizes beyond gene identity")

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
geometry_file <- file.path(
  analysis_dir, "dominant_uorf_structure_rules",
  "uorf_structure_gene_geometry.csv"
)
feature_file <- file.path(
  analysis_dir, "dominant_rdg_outputs", "rdg_feature_annotation.csv"
)

loo <- read_dt(loo_file, required = TRUE)
geometry <- read_dt(geometry_file, required = TRUE)
features <- read_dt(feature_file, required = TRUE)

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

features[, `:=`(
  kozak_score = kozak_score(kozak_strength),
  feature_bases_num = safe_num(feature_bases),
  peptide_length_aa_num = safe_num(peptide_length_aa),
  cds_overlap_bases_num = safe_num(cds_overlap_bases),
  frame_relative_to_cds_num = safe_num(frame_relative_to_cds),
  feature_order_num = safe_num(feature_order)
)]
branch_features <- features[
  feature_type %chin% c("leader_uorf", "overlapping_uorf")
]
feature_summary <- branch_features[, {
  first_idx <- order(feature_order_num, tx_start, na.last = TRUE)[1]
  data.table(
    annotation_uorf_features = .N,
    annotation_leader_uorfs = sum(feature_type == "leader_uorf",
                                  na.rm = TRUE),
    annotation_overlapping_uorfs =
      sum(feature_type == "overlapping_uorf", na.rm = TRUE),
    annotation_atg_uorfs = sum(start_codon_class == "ATG", na.rm = TRUE),
    annotation_near_cognate_uorfs =
      sum(start_codon_class == "near_cognate", na.rm = TRUE),
    annotation_other_start_uorfs =
      sum(start_codon_class == "other", na.rm = TRUE),
    annotation_strong_kozak_uorfs =
      sum(kozak_strength == "strong", na.rm = TRUE),
    annotation_moderate_kozak_uorfs =
      sum(kozak_strength == "moderate", na.rm = TRUE),
    annotation_mean_kozak_score = mean_safe(kozak_score),
    annotation_max_kozak_score = max_safe(kozak_score),
    annotation_mean_feature_bases = mean_safe(feature_bases_num),
    annotation_max_feature_bases = max_safe(feature_bases_num),
    annotation_mean_peptide_length_aa = mean_safe(peptide_length_aa_num),
    annotation_max_cds_overlap_bases = max_safe(cds_overlap_bases_num),
    annotation_frame0_uorfs =
      sum(frame_relative_to_cds_num == 0, na.rm = TRUE),
    annotation_frame1_uorfs =
      sum(frame_relative_to_cds_num == 1, na.rm = TRUE),
    annotation_frame2_uorfs =
      sum(frame_relative_to_cds_num == 2, na.rm = TRUE),
    annotation_manual_supported_uorfs =
      sum(grepl("M", translon_sources_merged, fixed = TRUE), na.rm = TRUE),
    annotation_first_uorf_kozak_score = kozak_score[first_idx],
    annotation_first_uorf_is_atg =
      as.integer(start_codon_class[first_idx] == "ATG"),
    annotation_first_uorf_is_near_cognate =
      as.integer(start_codon_class[first_idx] == "near_cognate")
  )
}, by = .(gene_symbol, tx_id)]

architecture <- merge(
  geometry,
  feature_summary,
  by = c("gene_symbol", "tx_id"),
  all.x = TRUE,
  sort = FALSE
)

feature_cols <- setdiff(names(feature_summary), c("gene_symbol", "tx_id"))
for (column in feature_cols) {
  architecture[is.na(get(column)) | !is.finite(safe_num(get(column))),
               (column) := 0]
}

model_dt <- copy(loo)
model_dt <- merge(model_dt, architecture,
                  by = c("gene_symbol", "tx_id"),
                  all.x = TRUE, sort = FALSE)
for (column in c("gene_symbol", "design_family", "heldout_study",
                 "cell_lines", "tissues", "condition_tokens")) {
  model_dt[is.na(get(column)) | !nzchar(as.character(get(column))),
           (column) := "unknown"]
}

keep <- is.finite(safe_num(model_dt$heldout_log2fc)) &
  is.finite(safe_num(model_dt$predicted_log2fc)) &
  !is.na(model_dt$n_uorfs)
model_dt <- model_dt[keep]
if (!nrow(model_dt)) stop("No rows with graph geometry.", call. = FALSE)

model_dt[, `:=`(
  f_design = clean_factor(design_family, "design_unknown"),
  f_cell = clean_factor(cell_lines, "cell_unknown"),
  f_tissue = clean_factor(tissues, "tissue_unknown"),
  f_condition = clean_factor(condition_tokens, "condition_unknown"),
  f_uorf_count_bin = clean_factor(uorf_count_bin, "uorf_count_unknown"),
  f_overlap_class = clean_factor(overlap_class, "overlap_unknown"),
  f_spacing_class = clean_factor(spacing_class, "spacing_unknown"),
  f_last_leader_gap_class =
    clean_factor(last_leader_gap_class, "gap_unknown"),
  f_first_uorf_length_class =
    clean_factor(first_uorf_length_class, "first_length_unknown"),
  f_has_overlapping_uorf =
    factor(ifelse(has_overlapping_uorf == TRUE,
                  "has_overlapping_uorf", "no_overlapping_uorf")),
  f_has_modeled_uorf =
    factor(ifelse(has_modeled_uorf == TRUE,
                  "has_modeled_uorf", "no_modeled_uorf")),
  log_matched_rows = log1p(safe_num(matched_rows_in_context)),
  log_strict_rows = log1p(safe_num(strict_rows_in_context)),
  log_case_counts = log1p(safe_num(total_case_counts)),
  log_control_counts = log1p(safe_num(total_control_counts))
)]

base_numeric_cols <- c("log_matched_rows", "log_strict_rows",
                       "log_case_counts", "log_control_counts")
for (column in base_numeric_cols) {
  model_dt[!is.finite(get(column)), (column) := 0]
}

raw_geometry_numeric <- c(
  "n_uorfs", "n_leader_uorfs", "n_overlapping_uorfs",
  "n_uorfs_t_only", "n_uorfs_tc_only", "n_uorfs_m_only",
  "n_uorfs_t_tc_shared", "n_uorfs_with_tc_support",
  "n_uorfs_with_m_support", "n_overlapping_or_nested_pairs",
  "n_abutting_pairs", "n_close_pairs_10bp", "n_close_pairs_30bp",
  "n_close_pairs_60bp", "min_inter_uorf_gap", "median_inter_uorf_gap",
  "mean_inter_uorf_gap", "min_start_spacing", "leader_length",
  "clean_cds_bases", "total_uorf_bases_double_counted",
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
  feature_cols
)
raw_geometry_numeric <- intersect(raw_geometry_numeric, names(model_dt))
for (column in raw_geometry_numeric) {
  model_dt[, (column) := safe_num(get(column))]
  model_dt[!is.finite(get(column)), (column) := 0]
}

log_geometry_source <- intersect(c(
  "n_uorfs", "n_leader_uorfs", "n_overlapping_uorfs",
  "n_close_pairs_10bp", "n_close_pairs_30bp", "n_close_pairs_60bp",
  "min_inter_uorf_gap", "median_inter_uorf_gap", "mean_inter_uorf_gap",
  "min_start_spacing", "leader_length", "clean_cds_bases",
  "uorf_union_bases", "leader_uorf_union_bases",
  "overlapping_uorf_union_bases", "longest_uorf_bases",
  "mean_uorf_bases", "first_uorf_bases", "first_uorf_to_cds_gap",
  "min_leader_stop_to_clean_cds_gap",
  "last_leader_stop_to_clean_cds_gap",
  "last_any_stop_to_clean_cds_gap", "max_cds_overlap_bases",
  "annotation_uorf_features", "annotation_atg_uorfs",
  "annotation_near_cognate_uorfs", "annotation_strong_kozak_uorfs",
  "annotation_mean_feature_bases", "annotation_max_feature_bases",
  "annotation_mean_peptide_length_aa", "annotation_max_cds_overlap_bases"
), names(model_dt))
log_geometry_cols <- paste0("log1p_", log_geometry_source)
for (i in seq_along(log_geometry_source)) {
  source_col <- log_geometry_source[[i]]
  dest_col <- log_geometry_cols[[i]]
  model_dt[, (dest_col) := log1p(pmax(safe_num(get(source_col)), 0))]
}

geometry_numeric_cols <- unique(c(
  raw_geometry_numeric,
  log_geometry_cols
))
geometry_numeric_cols <- geometry_numeric_cols[
  vapply(model_dt[, ..geometry_numeric_cols], function(x) {
    uniqueN(safe_num(x)) > 1
  }, logical(1))
]

covariate_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_covariates.csv"
)
covariate_cols <- unique(c(
  "gene_symbol", "tx_id", "design_family", "heldout_study",
  "heldout_log2fc", "predicted_log2fc",
  base_numeric_cols, geometry_numeric_cols,
  "uorf_count_bin", "overlap_class", "spacing_class",
  "last_leader_gap_class", "first_uorf_length_class",
  "has_overlapping_uorf", "has_modeled_uorf"
))
fwrite(model_dt[, ..covariate_cols], covariate_file)

context_terms <- paste(
  "0 + f_design + f_cell + f_tissue + f_condition",
  paste(base_numeric_cols, collapse = " + "),
  sep = " + "
)
geometry_factor_terms <- c(
  "f_uorf_count_bin", "f_overlap_class", "f_spacing_class",
  "f_last_leader_gap_class", "f_first_uorf_length_class",
  "f_has_overlapping_uorf", "f_has_modeled_uorf"
)
geometry_factor_terms <- geometry_factor_terms[
  vapply(geometry_factor_terms, function(x) {
    x %in% names(model_dt) && length(levels(model_dt[[x]])) > 1
  }, logical(1))
]
geometry_interaction_terms <- c()
if ("f_overlap_class" %in% geometry_factor_terms &&
    length(levels(model_dt$f_design)) > 1) {
  geometry_interaction_terms <- c(geometry_interaction_terms,
                                  "f_design:f_overlap_class")
}
if ("f_uorf_count_bin" %in% geometry_factor_terms &&
    length(levels(model_dt$f_design)) > 1) {
  geometry_interaction_terms <- c(geometry_interaction_terms,
                                  "f_design:f_uorf_count_bin")
}
geometry_extra_terms <- c(
  geometry_numeric_cols,
  geometry_factor_terms,
  geometry_interaction_terms
)
geometry_terms <- paste(
  context_terms,
  paste(geometry_extra_terms, collapse = " + "),
  sep = " + "
)

context_formula <- stats::as.formula(paste("~", context_terms))
geometry_formula <- stats::as.formula(paste("~", geometry_terms))

context_fit <- fit_blocked_direct_model(model_dt, context_formula,
                                        "context_nogene")
geometry_fit <- fit_blocked_direct_model(model_dt, geometry_formula,
                                         "context_plus_graph_geometry")

model_dt[, `:=`(
  context_nogene_predicted_log2fc =
    context_fit$prediction$predicted_log2fc,
  context_nogene_train_rows = context_fit$prediction$train_rows,
  graph_geometry_predicted_log2fc =
    geometry_fit$prediction$predicted_log2fc,
  graph_geometry_train_rows = geometry_fit$prediction$train_rows
)]
model_dt[, `:=`(
  context_nogene_error_log2fc =
    heldout_log2fc - context_nogene_predicted_log2fc,
  graph_geometry_error_log2fc =
    heldout_log2fc - graph_geometry_predicted_log2fc
)]
model_dt[, `:=`(
  context_nogene_abs_error_log2fc = abs(context_nogene_error_log2fc),
  graph_geometry_abs_error_log2fc = abs(graph_geometry_error_log2fc),
  hierarchical_error_log2fc = heldout_log2fc - predicted_log2fc,
  hierarchical_abs_error_log2fc = abs(heldout_log2fc - predicted_log2fc)
)]
model_dt[, `:=`(
  graph_vs_context_abs_error_delta_log2fc =
    context_nogene_abs_error_log2fc - graph_geometry_abs_error_log2fc,
  graph_vs_context_direction_improved =
    (sign(heldout_log2fc) == sign(graph_geometry_predicted_log2fc)) &
    (sign(heldout_log2fc) != sign(context_nogene_predicted_log2fc))
)]
model_dt[, graph_geometry_correction_class := fifelse(
  graph_vs_context_abs_error_delta_log2fc > 0.03,
  "graph_geometry_improves",
  fifelse(graph_vs_context_abs_error_delta_log2fc < -0.03,
          "graph_geometry_worse",
          "graph_geometry_neutral")
)]

summary <- rbindlist(c(
  list(metric_row(model_dt, "all")),
  lapply(split(model_dt, model_dt$design_family), function(x) {
    metric_row(x, paste0("design_family:", x$design_family[[1]]))
  })
), fill = TRUE)
summary[, `:=`(
  context_vs_hierarchical_mae_delta =
    hierarchical_mae_log2fc - context_nogene_mae_log2fc,
  graph_vs_hierarchical_mae_delta =
    hierarchical_mae_log2fc - graph_geometry_mae_log2fc,
  graph_vs_context_mae_delta =
    context_nogene_mae_log2fc - graph_geometry_mae_log2fc,
  context_vs_hierarchical_rmse_delta =
    hierarchical_rmse_log2fc - context_nogene_rmse_log2fc,
  graph_vs_hierarchical_rmse_delta =
    hierarchical_rmse_log2fc - graph_geometry_rmse_log2fc,
  graph_vs_context_rmse_delta =
    context_nogene_rmse_log2fc - graph_geometry_rmse_log2fc,
  graph_vs_context_q95_delta =
    context_nogene_q95_abs_error_log2fc -
    graph_geometry_q95_abs_error_log2fc,
  graph_vs_context_direction_accuracy_delta =
    graph_geometry_direction_accuracy -
    context_nogene_direction_accuracy
)]

gene_design <- model_dt[, {
  worst_idx <- which.max(graph_geometry_abs_error_log2fc)
  if (!length(worst_idx) || is.na(worst_idx)) worst_idx <- 1L
  data.table(
    graph_geometry_predictions = .N,
    n_studies = uniqueN(heldout_study),
    hierarchical_mae_log2fc = mean_safe(hierarchical_abs_error_log2fc),
    context_nogene_mae_log2fc =
      mean_safe(context_nogene_abs_error_log2fc),
    graph_geometry_mae_log2fc =
      mean_safe(graph_geometry_abs_error_log2fc),
    hierarchical_rmse_log2fc = rmse_safe(hierarchical_error_log2fc),
    context_nogene_rmse_log2fc =
      rmse_safe(context_nogene_error_log2fc),
    graph_geometry_rmse_log2fc =
      rmse_safe(graph_geometry_error_log2fc),
    graph_vs_context_mae_delta =
      mean_safe(context_nogene_abs_error_log2fc) -
      mean_safe(graph_geometry_abs_error_log2fc),
    graph_vs_context_rmse_delta =
      rmse_safe(context_nogene_error_log2fc) -
      rmse_safe(graph_geometry_error_log2fc),
    graph_fraction_improved_vs_context =
      mean(graph_vs_context_abs_error_delta_log2fc > 0, na.rm = TRUE),
    worst_graph_geometry_context =
      paste(heldout_study[worst_idx], cell_lines[worst_idx],
            tissues[worst_idx], condition_tokens[worst_idx], sep = " | "),
    top_graph_improved_contexts = collapse_unique(
      paste(heldout_study, cell_lines, tissues, condition_tokens, sep = " | ")
        [order(-graph_vs_context_abs_error_delta_log2fc)],
      5L
    )
  )
}, by = .(gene_symbol, tx_id, design_family)]

gene_static_cols <- unique(c(
  "gene_symbol", "tx_id", "n_uorfs", "n_leader_uorfs",
  "n_overlapping_uorfs", "has_overlapping_uorf", "uorf_count_bin",
  "overlap_class", "spacing_class", "last_leader_gap_class",
  "first_uorf_to_cds_gap", "max_cds_overlap_fraction_of_clean_cds",
  "annotation_mean_kozak_score", "annotation_strong_kozak_uorfs",
  "annotation_atg_uorfs", "annotation_near_cognate_uorfs"
))
gene_static <- unique(model_dt[, ..gene_static_cols])
gene_design <- merge(gene_design, gene_static,
                     by = c("gene_symbol", "tx_id"),
                     all.x = TRUE, sort = FALSE)
gene_design[, graph_geometry_model_status := fifelse(
  graph_vs_context_rmse_delta > 0.03 &
    graph_vs_context_mae_delta >= -0.01,
  "graph_geometry_reduces_error",
  fifelse(graph_vs_context_rmse_delta < -0.03 |
            graph_vs_context_mae_delta < -0.03,
          "graph_geometry_increases_error",
          "graph_geometry_neutral")
)]
setorder(gene_design, -graph_vs_context_rmse_delta,
         -graph_vs_context_mae_delta, gene_symbol, design_family)

summary_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_summary.csv"
)
loo_file_out <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_loo.csv"
)
gene_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_gene_design.csv"
)
viral_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_viral_gene_design.csv"
)
coef_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_coefficients.csv"
)
lambda_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_lambda_selection.csv"
)
metadata_file <- file.path(
  output_dir,
  "dominant_rdg_graph_geometry_model_metadata.csv"
)

fwrite(summary, summary_file)
fwrite(model_dt, loo_file_out)
fwrite(gene_design, gene_file)
fwrite(gene_design[design_family == "viral_infection"], viral_file)
fwrite(rbindlist(list(context_fit$coefficient, geometry_fit$coefficient),
                 fill = TRUE),
       coef_file)
fwrite(rbindlist(list(context_fit$lambda_selection,
                     geometry_fit$lambda_selection), fill = TRUE),
       lambda_file)
fwrite(rbindlist(list(context_fit$metadata, geometry_fit$metadata),
                 fill = TRUE),
       metadata_file)

plot_dt <- summary[grepl("^design_family:", scope)]
if (nrow(plot_dt)) {
  plot_dt[, design_family := sub("^design_family:", "", scope)]
  plot_dt <- melt(
    plot_dt,
    id.vars = "design_family",
    measure.vars = c("context_nogene_rmse_log2fc",
                     "graph_geometry_rmse_log2fc"),
    variable.name = "model",
    value.name = "rmse_log2fc"
  )
  plot_dt[, model := fifelse(
    model == "context_nogene_rmse_log2fc",
    "context, no gene ID",
    "context + graph geometry"
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
    geom_col(position = position_dodge(width = 0.74), width = 0.66) +
    coord_flip() +
    scale_fill_manual(values = c(
      "context, no gene ID" = "#6b7280",
      "context + graph geometry" = "#0f766e"
    )) +
    labs(
      title = "Does RDG Geometry Generalize Across Genes?",
      subtitle = "Leave-one-gene-out prediction; no gene identity covariate",
      x = NULL, y = "RMSE (log2FC)", fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "graph_geometry_model_rmse.png"),
         p, width = 10.2, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "graph_geometry_model_rmse.pdf"),
         p, width = 10.2, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "graph_geometry_model_rmse.html"),
      title = "RDG graph geometry generalization RMSE"
    )
  }
}

viral_gene <- gene_design[design_family == "viral_infection"]
if (nrow(viral_gene)) {
  viral_plot <- viral_gene[
    order(-abs(graph_vs_context_rmse_delta))
  ][seq_len(min(.N, 30))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  viral_plot[, graph_geometry_direction_class := fifelse(
    graph_vs_context_rmse_delta >= 0.01,
    "improvement",
    fifelse(graph_vs_context_rmse_delta <= -0.01,
            "worse",
            "near_zero")
  )]
  p2 <- ggplot(
    viral_plot,
    aes(x = graph_vs_context_rmse_delta, y = gene_symbol,
        fill = graph_geometry_direction_class,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Status: ", graph_geometry_model_status,
          "<br>RMSE delta: ", round(graph_vs_context_rmse_delta, 3),
          "<br>MAE delta: ", round(graph_vs_context_mae_delta, 3),
          "<br>uORFs: ", n_uorfs,
          "<br>overlapping uORFs: ", n_overlapping_uorfs,
          "<br>overlap class: ", overlap_class,
          "<br>worst context: ", worst_graph_geometry_context
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(
      "improvement" = "#0f766e",
      "worse" = "#b91c1c",
      "near_zero" = "#6b7280"
    ), labels = c(
      "improvement" = "geometry improves",
      "worse" = "geometry worse",
      "near_zero" = "near zero"
    )) +
    labs(
      title = "Viral-Infection Leave-One-Gene-Out Geometry Model",
      subtitle = "Positive values mean graph geometry improves over metadata/context without gene ID",
      x = "RMSE reduction (log2FC)", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "graph_geometry_model_viral_genes.png"),
         p2, width = 9.6, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "graph_geometry_model_viral_genes.pdf"),
         p2, width = 9.6, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir, "graph_geometry_model_viral_genes.html"),
      title = "Viral graph geometry generalization"
    )
  }
}

coef_plot_dt <- rbindlist(list(context_fit$coefficient,
                               geometry_fit$coefficient), fill = TRUE)
coef_plot_dt <- coef_plot_dt[
  model == "context_plus_graph_geometry" &
    !grepl("^f_(cell|tissue|condition)", feature) &
    !grepl("^f_design", feature)
][seq_len(min(.N, 30))]
if (nrow(coef_plot_dt)) {
  coef_plot_dt[, feature := factor(feature, levels = rev(feature))]
  p3 <- ggplot(
    coef_plot_dt,
    aes(x = coefficient, y = feature,
        fill = coefficient > 0,
        text = paste0(
          "Feature: ", feature,
          "<br>Coefficient: ", round(coefficient, 4)
        ))
  ) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#0f766e",
                                 "FALSE" = "#b91c1c"),
                      labels = c("TRUE" = "positive",
                                 "FALSE" = "negative")) +
    labs(
      title = "Largest Graph-Geometry Ridge Coefficients",
      subtitle = "Leave-one-gene-out model coefficients; magnitudes are standardized by glmnet",
      x = "Coefficient", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "graph_geometry_model_coefficients.png"),
         p3, width = 9.6, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "graph_geometry_model_coefficients.pdf"),
         p3, width = 9.6, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p3,
      file.path(figure_dir, "graph_geometry_model_coefficients.html"),
      title = "RDG graph geometry coefficients"
    )
  }
}

message("Saved graph-geometry model outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
