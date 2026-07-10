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
output_dir <- file.path(analysis_dir, "dominant_rdg_dirichlet_multinomial")
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

max_safe <- function(x, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

median_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

collapse_unique <- function(x, n = 8L) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

q_component <- function(q) {
  q <- safe_num(q)
  out <- rep(0, length(q))
  finite <- is.finite(q) & q > 0
  out[finite] <- pmin(-log10(q[finite]) / 3, 1)
  out[is.finite(q) & q == 0] <- 1
  out
}

weighted_mean_safe <- function(x, w = NULL) {
  x <- safe_num(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- safe_num(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

random_effects_summary <- function(effect, se, weight = NULL) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- safe_num(weight)
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L, effect = NA_real_, se = NA_real_, lower = NA_real_,
      upper = NA_real_, p = NA_real_, tau2 = NA_real_, i2 = NA_real_,
      q_stat = NA_real_, direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  sei <- pmax(se[keep], 1e-6)
  wi_external <- weight[keep] / mean(weight[keep])
  vi <- sei^2
  wi <- wi_external / vi
  fixed <- sum(wi * yi) / sum(wi)
  q_stat <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q_stat - df) / c_val) else 0
  wi_re <- wi_external / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  z <- mu / se_mu
  p <- 2 * pnorm(-abs(z))
  sign_ref <- sign(mu)
  direction_fraction <- if (sign_ref == 0) mean(yi == 0) else mean(sign(yi) == sign_ref)
  i2 <- if (q_stat > 0 && df > 0) max(0, (q_stat - df) / q_stat) else 0
  data.table(
    n = length(yi), effect = mu, se = se_mu, lower = lower,
    upper = upper, p = p, tau2 = tau2, i2 = i2, q_stat = q_stat,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

beta_binom_two_sided_p <- function(k, n, alpha, beta, exact_limit = 650L) {
  k <- round(safe_num(k))
  n <- round(safe_num(n))
  alpha <- safe_num(alpha)
  beta <- safe_num(beta)
  if (!is.finite(k) || !is.finite(n) || !is.finite(alpha) ||
      !is.finite(beta) || n <= 0 || alpha <= 0 || beta <= 0 ||
      k < 0 || k > n) {
    return(NA_real_)
  }
  a0 <- alpha + beta
  mean_x <- n * alpha / a0
  var_x <- n * alpha * beta * (a0 + n) / (a0^2 * (a0 + 1))
  if (!is.finite(var_x) || var_x <= 0) return(NA_real_)
  if (n <= exact_limit) {
    xs <- 0:n
    logp <- lchoose(n, xs) + lbeta(xs + alpha, n - xs + beta) -
      lbeta(alpha, beta)
    logp <- logp - max(logp)
    p <- exp(logp)
    p <- p / sum(p)
    lower <- sum(p[xs <= k])
    upper <- sum(p[xs >= k])
    return(min(1, 2 * min(lower, upper)))
  }
  z <- (k - mean_x) / sqrt(var_x)
  min(1, 2 * pnorm(-abs(z)))
}

prefixed_sum <- function(dt, by_cols, count_cols, prefix) {
  out <- dt[, lapply(.SD, function(x) sum(finite_or_zero(x), na.rm = TRUE)),
            by = by_cols, .SDcols = count_cols]
  setnames(out, count_cols, paste0(prefix, count_cols))
  out[, (paste0(prefix, "total")) :=
        rowSums(as.matrix(.SD), na.rm = TRUE),
      .SDcols = paste0(prefix, count_cols)]
  out
}

message("Empirical-Bayes Dirichlet-multinomial RDG branch model:")
message("  1. Learn weak branch priors from other studies when possible")
message("  2. Use row control counts as posterior baseline")
message("  3. Test case branch counts with beta-binomial posterior prediction")

joint_file <- file.path(
  analysis_dir,
  "dominant_rdg_joint_branch_allocation",
  "dominant_rdg_joint_branch_allocation_rows.csv"
)
joint <- read_dt(joint_file, required = TRUE)

branch_levels <- c("leader_uORF", "overlapping_uORF", "clean_CDS", "other")
suffix <- c(
  leader_uORF = "leader_uORF",
  overlapping_uORF = "overlapping_uORF",
  clean_CDS = "clean_CDS",
  other = "other"
)
case_cols <- paste0("joint_case_count_", suffix[branch_levels])
control_cols <- paste0("joint_control_count_", suffix[branch_levels])
required <- c(
  "gene_symbol", "tx_id", "design_family", "study", "case_condition",
  "control_conditions", "CELL_LINE", "TISSUE", "match_scope",
  "sample_set_signature", "joint_evaluable", "joint_strict",
  "joint_review_class", "joint_review_priority",
  "case_branch_total_counts", "control_branch_total_counts",
  "n_case_runs_merged", "n_control_runs_merged", case_cols, control_cols
)
missing <- setdiff(required, names(joint))
if (length(missing)) {
  stop("Joint allocation row table is missing columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

joint[, dm_row_id := .I]
for (column in intersect(
  c(case_cols, control_cols, "case_branch_total_counts",
    "control_branch_total_counts", "n_case_runs_merged",
    "n_control_runs_merged", "joint_review_priority"),
  names(joint)
)) {
  joint[, (column) := finite_or_zero(get(column))]
}
for (column in intersect(c("joint_evaluable", "joint_strict"), names(joint))) {
  if (!is.logical(joint[[column]])) {
    joint[, (column) := tolower(as.character(get(column))) %chin%
            c("true", "t", "1", "yes")]
  }
  joint[is.na(get(column)), (column) := FALSE]
}

gene_totals <- prefixed_sum(joint, c("gene_symbol", "tx_id"),
                            control_cols, "gene_")
gene_study_totals <- prefixed_sum(joint, c("gene_symbol", "tx_id", "study"),
                                  control_cols, "gene_study_")
design_totals <- prefixed_sum(joint, "design_family", control_cols, "design_")
design_study_totals <- prefixed_sum(joint, c("design_family", "study"),
                                    control_cols, "design_study_")
global_counts <- colSums(as.matrix(joint[, ..control_cols]), na.rm = TRUE)
global_total <- sum(global_counts)
global_prop <- if (global_total > 0) global_counts / global_total else rep(1 / 4, 4)

joint <- merge(joint, gene_totals, by = c("gene_symbol", "tx_id"),
               all.x = TRUE, sort = FALSE)
joint <- merge(joint, gene_study_totals,
               by = c("gene_symbol", "tx_id", "study"),
               all.x = TRUE, sort = FALSE)
joint <- merge(joint, design_totals, by = "design_family",
               all.x = TRUE, sort = FALSE)
joint <- merge(joint, design_study_totals,
               by = c("design_family", "study"),
               all.x = TRUE, sort = FALSE)
setorder(joint, dm_row_id)

gene_cols <- paste0("gene_", control_cols)
gene_study_cols <- paste0("gene_study_", control_cols)
design_cols <- paste0("design_", control_cols)
design_study_cols <- paste0("design_study_", control_cols)
gene_mat <- as.matrix(joint[, ..gene_cols])
gene_study_mat <- as.matrix(joint[, ..gene_study_cols])
design_mat <- as.matrix(joint[, ..design_cols])
design_study_mat <- as.matrix(joint[, ..design_study_cols])
for (mat_name in c("gene_mat", "gene_study_mat", "design_mat",
                   "design_study_mat")) {
  m <- get(mat_name)
  m[!is.finite(m)] <- 0
  assign(mat_name, m)
}
gene_other <- pmax(gene_mat - gene_study_mat, 0)
design_other <- pmax(design_mat - design_study_mat, 0)
gene_total <- rowSums(gene_mat, na.rm = TRUE)
gene_other_total <- rowSums(gene_other, na.rm = TRUE)
design_other_total <- rowSums(design_other, na.rm = TRUE)

prior_counts <- matrix(0, nrow(joint), length(branch_levels))
prior_source <- rep("global_control_background", nrow(joint))
prior_strength <- rep(12, nrow(joint))

use_gene_other <- gene_other_total >= 100
prior_counts[use_gene_other, ] <- gene_other[use_gene_other, , drop = FALSE]
prior_source[use_gene_other] <- "gene_other_studies"
prior_strength[use_gene_other] <- pmin(120, 2 * sqrt(gene_other_total[use_gene_other]))

use_gene_all <- !use_gene_other & gene_total >= 100
prior_counts[use_gene_all, ] <- gene_mat[use_gene_all, , drop = FALSE]
prior_source[use_gene_all] <- "gene_all_studies"
prior_strength[use_gene_all] <- pmin(60, sqrt(gene_total[use_gene_all]))

use_design_other <- !use_gene_other & !use_gene_all & design_other_total >= 1000
prior_counts[use_design_other, ] <- design_other[use_design_other, , drop = FALSE]
prior_source[use_design_other] <- "design_other_studies"
prior_strength[use_design_other] <- pmin(35, sqrt(design_other_total[use_design_other]) / 4)

use_global <- rowSums(prior_counts) <= 0
if (any(use_global)) {
  prior_counts[use_global, ] <- matrix(
    global_prop, nrow = sum(use_global), ncol = length(branch_levels),
    byrow = TRUE
  )
  prior_strength[use_global] <- 12
}
prior_prop <- prior_counts / rowSums(prior_counts)
prior_prop[!is.finite(prior_prop)] <- 1 / length(branch_levels)
prior_alpha <- sweep(prior_prop, 1, prior_strength, "*") + 0.5
colnames(prior_alpha) <- branch_levels

case_mat <- as.matrix(joint[, ..case_cols])
control_mat <- as.matrix(joint[, ..control_cols])
case_mat[!is.finite(case_mat)] <- 0
control_mat[!is.finite(control_mat)] <- 0
case_mat <- pmax(case_mat, 0)
control_mat <- pmax(control_mat, 0)
case_total <- rowSums(case_mat)
control_total <- rowSums(control_mat)
posterior_alpha <- control_mat + prior_alpha
posterior_total <- rowSums(posterior_alpha)
posterior_prop <- posterior_alpha / posterior_total

long <- rbindlist(lapply(seq_along(branch_levels), function(i) {
  branch <- branch_levels[[i]]
  alpha_success <- posterior_alpha[, i]
  alpha_failure <- pmax(posterior_total - alpha_success, 1e-8)
  pred_mean_count <- case_total * alpha_success / posterior_total
  pred_var_count <- case_total * alpha_success * alpha_failure *
    (posterior_total + case_total) /
    (posterior_total^2 * (posterior_total + 1))
  pred_sd_prop <- sqrt(pmax(pred_var_count, 0)) / pmax(case_total, 1)
  observed_prop <- case_mat[, i] / pmax(case_total, 1)
  predicted_prop <- alpha_success / posterior_total
  p_value <- mapply(
    beta_binom_two_sided_p,
    k = case_mat[, i],
    n = case_total,
    alpha = alpha_success,
    beta = alpha_failure
  )
  data.table(
    dm_row_id = joint$dm_row_id,
    gene_symbol = joint$gene_symbol,
    tx_id = joint$tx_id,
    design_family = joint$design_family,
    study = joint$study,
    case_condition = joint$case_condition,
    control_conditions = joint$control_conditions,
    CELL_LINE = joint$CELL_LINE,
    TISSUE = joint$TISSUE,
    match_scope = joint$match_scope,
    sample_set_signature = joint$sample_set_signature,
    branch_class = branch,
    dm_case_count = case_mat[, i],
    dm_control_count = control_mat[, i],
    dm_case_total = case_total,
    dm_control_total = control_total,
    dm_prior_source = prior_source,
    dm_prior_strength = prior_strength,
    dm_prior_prop = prior_prop[, i],
    dm_posterior_baseline_prop = predicted_prop,
    dm_case_prop = observed_prop,
    dm_expected_case_count = pred_mean_count,
    dm_delta = observed_prop - predicted_prop,
    dm_delta_se = pmax(pred_sd_prop, 1e-8),
    dm_delta_z = (observed_prop - predicted_prop) / pmax(pred_sd_prop, 1e-8),
    dm_delta_p = p_value,
    dm_joint_review_class = joint$joint_review_class,
    dm_joint_review_priority = joint$joint_review_priority,
    dm_joint_evaluable = joint$joint_evaluable,
    dm_joint_strict = joint$joint_strict,
    dm_n_case_runs_merged = joint$n_case_runs_merged,
    dm_n_control_runs_merged = joint$n_control_runs_merged
  )
}), use.names = TRUE)

long[, dm_delta_q := p.adjust(dm_delta_p, method = "BH")]
long[, dm_min_runs_merged := pmin(dm_n_case_runs_merged,
                                  dm_n_control_runs_merged)]
long[, dm_min_total_counts := pmin(dm_case_total, dm_control_total)]
long[, dm_max_branch_count := pmax(dm_case_count, dm_control_count)]
long[, dm_min_branch_count := pmin(dm_case_count, dm_control_count)]
long[, dm_evaluable :=
       dm_joint_evaluable == TRUE &
       dm_min_runs_merged >= 2 &
       dm_min_total_counts >= 100 &
       dm_max_branch_count >= 30 &
       is.finite(dm_delta) &
       is.finite(dm_delta_se) &
       is.finite(dm_delta_p)]
long[, dm_strict :=
       dm_evaluable == TRUE &
       dm_joint_strict == TRUE &
       dm_min_branch_count >= 30]
long[, dm_row_priority := clip01(
  0.24 * pmin(abs(dm_delta) / 0.15, 1) +
    0.20 * q_component(dm_delta_q) +
    0.16 * pmin(log1p(dm_max_branch_count) / log1p(500), 1) +
    0.12 * pmin(dm_min_runs_merged / 3, 1) +
    0.10 * as.numeric(dm_strict == TRUE) +
    0.08 * as.numeric(dm_prior_source == "gene_other_studies") +
    0.06 * as.numeric(design_family == "viral_infection") +
    0.04 * pmin(dm_joint_review_priority, 1)
)]
long[dm_evaluable != TRUE, dm_row_priority := pmin(dm_row_priority, 0.25)]
long[, dm_row_support_class := fcase(
  !dm_evaluable,
  "dm_low_count_or_not_evaluable",
  dm_strict == TRUE &
    is.finite(dm_delta_q) & dm_delta_q <= 0.10 &
    abs(dm_delta) >= 0.04,
  "dm_strict_posterior_predictive_shift",
  dm_evaluable == TRUE &
    is.finite(dm_delta_p) & dm_delta_p <= 0.01 &
    abs(dm_delta) >= 0.05,
  "dm_posterior_predictive_review",
  default = "dm_evaluable_no_strong_shift"
)]
setorder(long, -dm_row_priority, gene_symbol, design_family, branch_class)
fwrite(long, file.path(output_dir,
                       "dominant_rdg_dirichlet_multinomial_branch_rows.csv"))

study_effects <- long[dm_evaluable == TRUE, {
  re <- random_effects_summary(dm_delta, dm_delta_se,
                               pmax(log1p(dm_max_branch_count), 1))
  top_i <- which.max(dm_row_priority)
  if (!length(top_i) || is.na(top_i)) top_i <- 1L
  data.table(
    dm_study_context_rows = .N,
    dm_study_strict_context_rows = sum(dm_strict == TRUE, na.rm = TRUE),
    dm_study_delta = re$effect,
    dm_study_delta_se = re$se,
    dm_study_delta_lower = re$lower,
    dm_study_delta_upper = re$upper,
    dm_study_delta_p = re$p,
    dm_study_i2 = re$i2,
    dm_study_direction_fraction = re$direction_fraction,
    dm_study_ci_excludes_zero = re$ci_excludes_zero,
    dm_study_min_branch_p = min(dm_delta_p, na.rm = TRUE),
    dm_study_max_priority = max_safe(dm_row_priority, 0),
    dm_study_top_context = paste(
      study[top_i], CELL_LINE[top_i], TISSUE[top_i],
      case_condition[top_i], paste("vs", control_conditions[top_i]),
      sep = " | "
    ),
    dm_study_top_prior_source = dm_prior_source[top_i],
    dm_study_top_match_scope = match_scope[top_i]
  )
}, by = .(gene_symbol, tx_id, design_family, branch_class, study)]
study_effects <- study_effects[is.finite(dm_study_delta) &
                                 is.finite(dm_study_delta_se) &
                                 dm_study_delta_se > 0]
fwrite(study_effects, file.path(
  output_dir,
  "dominant_rdg_dirichlet_multinomial_study_effects.csv"
))

branch_effects <- study_effects[, {
  re <- random_effects_summary(dm_study_delta, dm_study_delta_se)
  top_i <- which.max(dm_study_max_priority)
  if (!length(top_i) || is.na(top_i)) top_i <- 1L
  data.table(
    dm_studies = .N,
    dm_context_rows = sum(dm_study_context_rows, na.rm = TRUE),
    dm_strict_context_rows = sum(dm_study_strict_context_rows, na.rm = TRUE),
    dm_delta = re$effect,
    dm_delta_se = re$se,
    dm_delta_lower = re$lower,
    dm_delta_upper = re$upper,
    dm_delta_p = re$p,
    dm_i2 = re$i2,
    dm_q_stat = re$q_stat,
    dm_direction_fraction = re$direction_fraction,
    dm_ci_excludes_zero = re$ci_excludes_zero,
    dm_min_study_branch_p = min(dm_study_min_branch_p, na.rm = TRUE),
    dm_max_context_priority = max_safe(dm_study_max_priority, 0),
    dm_top_study = study[top_i],
    dm_top_context = dm_study_top_context[top_i],
    dm_top_prior_source = dm_study_top_prior_source[top_i],
    dm_top_match_scope = dm_study_top_match_scope[top_i],
    dm_study_list = collapse_unique(study, 8L)
  )
}, by = .(gene_symbol, tx_id, design_family, branch_class)]
branch_effects[, dm_delta_q := p.adjust(dm_delta_p, method = "BH")]
branch_effects[, dm_support_class := fcase(
  dm_studies >= 2 &
    dm_ci_excludes_zero == TRUE &
    is.finite(dm_delta_q) & dm_delta_q <= 0.20 &
    abs(dm_delta) >= 0.035 &
    dm_direction_fraction >= 0.70,
  "dm_replicated_shift",
  dm_studies >= 2 &
    (dm_ci_excludes_zero == TRUE |
       (is.finite(dm_delta_q) & dm_delta_q <= 0.30)) &
    abs(dm_delta) >= 0.025 &
    dm_direction_fraction >= 0.65,
  "dm_replicated_review",
  dm_studies == 1 &
    is.finite(dm_delta_p) & dm_delta_p <= 0.05 &
    abs(dm_delta) >= 0.08,
  "dm_single_study_strong",
  default = "dm_neutral"
)]
branch_effects[, dm_priority := clip01(
  0.24 * pmin(abs(dm_delta) / 0.15, 1) +
    0.18 * as.numeric(dm_ci_excludes_zero == TRUE) +
    0.16 * pmin(dm_studies / 3, 1) +
    0.13 * pmin(dm_direction_fraction, 1) +
    0.10 * (1 - pmin(dm_i2, 1)) +
    0.09 * q_component(dm_delta_q) +
    0.06 * pmin(dm_max_context_priority, 1) +
    0.04 * as.numeric(design_family == "viral_infection")
)]
branch_effects[dm_support_class == "dm_neutral" & dm_priority > 0.60,
               dm_priority := 0.60]
branch_effects[, dm_abs_delta := abs(dm_delta)]
setorder(branch_effects, -dm_priority, -dm_abs_delta,
         gene_symbol, design_family, branch_class)
branch_effects[, dm_abs_delta := NULL]
fwrite(branch_effects, file.path(
  output_dir,
  "dominant_rdg_dirichlet_multinomial_branch_effects.csv"
))

top_branch <- branch_effects[, .SD[1], by = .(gene_symbol, tx_id, design_family)]
gene_design <- top_branch[, .(
  gene_symbol, tx_id, design_family,
  dm_top_branch_class = as.character(branch_class),
  dm_top_delta = dm_delta,
  dm_top_delta_lower = dm_delta_lower,
  dm_top_delta_upper = dm_delta_upper,
  dm_top_delta_q = dm_delta_q,
  dm_top_support_class = dm_support_class,
  dm_top_priority = dm_priority,
  dm_top_studies = dm_studies,
  dm_top_context_rows = dm_context_rows,
  dm_top_direction_fraction = dm_direction_fraction,
  dm_top_i2 = dm_i2,
  dm_top_context = dm_top_context,
  dm_top_prior_source = dm_top_prior_source,
  dm_study_list = dm_study_list
)]

wide <- dcast(
  branch_effects,
  gene_symbol + tx_id + design_family ~ branch_class,
  value.var = c("dm_delta", "dm_delta_q", "dm_priority",
                "dm_support_class"),
  fill = NA
)
gene_design <- merge(gene_design, wide,
                     by = c("gene_symbol", "tx_id", "design_family"),
                     all.x = TRUE, sort = FALSE)
delta_cols <- intersect(paste0("dm_delta_", branch_levels),
                        names(gene_design))
gene_design[, dm_l1_effect := {
  m <- as.matrix(.SD)
  m[!is.finite(m)] <- 0
  rowSums(abs(m)) / 2
}, .SDcols = delta_cols]
support_cols <- intersect(paste0("dm_support_class_", branch_levels),
                          names(gene_design))
supported_classes <- c("dm_replicated_shift", "dm_replicated_review",
                       "dm_single_study_strong")
if (length(support_cols)) {
  support_matrix <- sapply(
    support_cols,
    function(col) gene_design[[col]] %chin% supported_classes
  )
  if (is.null(dim(support_matrix))) support_matrix <- matrix(support_matrix, ncol = 1)
  gene_design[, dm_supported_branch_count := rowSums(support_matrix,
                                                     na.rm = TRUE)]
} else {
  gene_design[, dm_supported_branch_count := 0L]
}
setorder(gene_design, -dm_top_priority, -dm_l1_effect,
         gene_symbol, design_family)
fwrite(gene_design, file.path(
  output_dir,
  "dominant_rdg_dirichlet_multinomial_gene_design_summary.csv"
))

review <- gene_design[
  dm_top_support_class != "dm_neutral" | dm_top_priority >= 0.60
]
setorder(review, -dm_top_priority, -dm_l1_effect,
         gene_symbol, design_family)
fwrite(review, file.path(
  output_dir,
  "dominant_rdg_dirichlet_multinomial_review.csv"
))

viral_review <- gene_design[design_family == "viral_infection"]
setorder(viral_review, -dm_top_priority, -dm_l1_effect, gene_symbol)
fwrite(viral_review, file.path(
  output_dir,
  "dominant_rdg_dirichlet_multinomial_viral_review.csv"
))

summary_metrics <- rbindlist(list(
  branch_effects[, .(
    branch_effect_rows = .N,
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    genes = uniqueN(gene_symbol),
    replicated_shift_branches = sum(dm_support_class == "dm_replicated_shift",
                                    na.rm = TRUE),
    replicated_review_branches = sum(dm_support_class == "dm_replicated_review",
                                     na.rm = TRUE),
    single_study_strong_branches = sum(dm_support_class == "dm_single_study_strong",
                                       na.rm = TRUE),
    median_abs_delta = median_safe(abs(dm_delta)),
    max_abs_delta = max_safe(abs(dm_delta), 0),
    median_i2 = median_safe(dm_i2)
  )][, scope := "all"],
  branch_effects[, .(
    branch_effect_rows = .N,
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    genes = uniqueN(gene_symbol),
    replicated_shift_branches = sum(dm_support_class == "dm_replicated_shift",
                                    na.rm = TRUE),
    replicated_review_branches = sum(dm_support_class == "dm_replicated_review",
                                     na.rm = TRUE),
    single_study_strong_branches = sum(dm_support_class == "dm_single_study_strong",
                                       na.rm = TRUE),
    median_abs_delta = median_safe(abs(dm_delta)),
    max_abs_delta = max_safe(abs(dm_delta), 0),
    median_i2 = median_safe(dm_i2)
  ), by = .(scope = paste0("design_family:", design_family))]
), fill = TRUE)
setcolorder(summary_metrics, c("scope", setdiff(names(summary_metrics), "scope")))
fwrite(summary_metrics, file.path(
  output_dir,
  "dominant_rdg_dirichlet_multinomial_summary_metrics.csv"
))

plot_dt <- head(review[order(-dm_top_priority)], 45)
if (nrow(plot_dt)) {
  plot_dt[, plot_label_text := paste0(gene_symbol, " / ", design_family)]
  plot_dt[, plot_label := factor(plot_label_text,
                                 levels = rev(unique(plot_label_text)))]
  p_top <- ggplot(
    plot_dt,
    aes(
      x = dm_top_delta,
      y = plot_label,
      fill = dm_top_branch_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Design: ", design_family,
        "<br>Top branch: ", dm_top_branch_class,
        "<br>DM delta: ", round(dm_top_delta, 3),
        "<br>95% interval: [", round(dm_top_delta_lower, 3), ", ",
        round(dm_top_delta_upper, 3), "]",
        "<br>q: ", signif(dm_top_delta_q, 3),
        "<br>Support: ", dm_top_support_class,
        "<br>Studies: ", dm_top_studies,
        "<br>Top context: ", dm_top_context
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "#666666", linewidth = 0.25) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31",
      other = "#667085"
    )) +
    labs(
      title = "Empirical-Bayes Dirichlet-Multinomial RDG Effects",
      subtitle = "Posterior-predictive branch deltas against control-informed baselines",
      x = "Branch allocation delta versus posterior baseline",
      y = NULL,
      fill = "Top branch"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "dirichlet_multinomial_top_effects.png"),
         p_top, width = 10.8, height = 8.4, dpi = 180)
  ggsave(file.path(figure_dir, "dirichlet_multinomial_top_effects.pdf"),
         p_top, width = 10.8, height = 8.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir, "dirichlet_multinomial_top_effects.html"),
      title = "Dirichlet-multinomial branch effects"
    )
  }
}

viral_plot <- head(viral_review[order(-dm_top_priority)], 35)
if (nrow(viral_plot)) {
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  viral_plot[, dm_support_group := fcase(
    dm_top_support_class %chin% c("dm_replicated_shift",
                                  "dm_replicated_review"),
    "replicated",
    dm_top_support_class == "dm_single_study_strong",
    "single study",
    default = "neutral / review"
  )]
  p_viral <- ggplot(
    viral_plot,
    aes(
      x = dm_top_priority,
      y = gene_symbol,
      fill = dm_top_branch_class,
      alpha = dm_support_group,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Top branch: ", dm_top_branch_class,
        "<br>DM delta: ", round(dm_top_delta, 3),
        "<br>Support: ", dm_top_support_class,
        "<br>Studies: ", dm_top_studies,
        "<br>Top context: ", dm_top_context
      )
    )
  ) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31",
      other = "#667085"
    )) +
    scale_alpha_manual(values = c(
      replicated = 1,
      `single study` = 0.78,
      `neutral / review` = 0.42
    )) +
    labs(
      title = "Viral Dirichlet-Multinomial RDG Review",
      subtitle = "Posterior-predictive branch outlier priority",
      x = "Dirichlet-multinomial priority",
      y = NULL,
      fill = "Top branch",
      alpha = "Support"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "dirichlet_multinomial_viral_review.png"),
         p_viral, width = 9.8, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "dirichlet_multinomial_viral_review.pdf"),
         p_viral, width = 9.8, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_viral,
      file.path(figure_dir, "dirichlet_multinomial_viral_review.html"),
      title = "Viral Dirichlet-multinomial review"
    )
  }
}

heat_dt <- head(gene_design[order(-dm_top_priority)], 40)
if (nrow(heat_dt)) {
  heat_dt[, plot_label := paste0(gene_symbol, " / ", design_family)]
  heat_long <- melt(
    heat_dt,
    id.vars = c("plot_label", "gene_symbol", "design_family",
                "dm_top_priority", "dm_top_support_class"),
    measure.vars = intersect(paste0("dm_delta_", branch_levels),
                             names(heat_dt)),
    variable.name = "branch",
    value.name = "dm_delta"
  )
  heat_long[, branch := sub("^dm_delta_", "", branch)]
  heat_long[, branch := factor(branch, levels = branch_levels)]
  heat_long[, plot_label := factor(plot_label,
                                   levels = rev(unique(heat_dt$plot_label)))]
  p_heat <- ggplot(
    heat_long,
    aes(
      x = branch,
      y = plot_label,
      fill = dm_delta,
      text = paste0(
        "Gene/design: ", plot_label,
        "<br>Branch: ", branch,
        "<br>DM delta: ", round(dm_delta, 3),
        "<br>Support: ", dm_top_support_class,
        "<br>Priority: ", round(dm_top_priority, 3)
      )
    )
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#1f5aa6", mid = "white", high = "#8b1a1a",
      midpoint = 0, limits = c(-0.25, 0.25), oob = scales::squish,
      na.value = "#f2f4f7"
    ) +
    labs(
      title = "Dirichlet-Multinomial Branch Delta Matrix",
      subtitle = "Study-aware branch deltas after control-informed posterior prediction",
      x = NULL,
      y = NULL,
      fill = "Delta"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(figure_dir, "dirichlet_multinomial_delta_heatmap.png"),
         p_heat, width = 9.6, height = 8.8, dpi = 180)
  ggsave(file.path(figure_dir, "dirichlet_multinomial_delta_heatmap.pdf"),
         p_heat, width = 9.6, height = 8.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_heat,
      file.path(figure_dir, "dirichlet_multinomial_delta_heatmap.html"),
      title = "Dirichlet-multinomial branch delta matrix"
    )
  }
}

message("Saved empirical-Bayes Dirichlet-multinomial outputs in: ", output_dir)
print(summary_metrics[scope %chin% c("all", "design_family:viral_infection")])
