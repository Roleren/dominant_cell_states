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
input_dir <- file.path(analysis_dir, "dominant_rdg_grouped_branch_usage")
output_dir <- file.path(analysis_dir, "dominant_rdg_branch_allocation")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

contrast_file <- file.path(
  input_dir,
  "dominant_rdg_grouped_branch_usage_contrast_table.csv"
)
if (!file.exists(contrast_file)) {
  stop("Missing grouped branch contrast table: ", contrast_file,
       call. = FALSE)
}

stage_message <- function(label) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", label)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

first_nonempty <- function(x, default = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else default
}

collapse_examples <- function(x, n = 6L) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (!length(x)) return("")
  paste(head(x, n), collapse = ";")
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

weighted_mean_safe <- function(x, w = NULL) {
  x <- safe_num(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- safe_num(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

beta_var <- function(a, b) {
  total <- a + b
  a * b / ((total^2) * (total + 1))
}

logit <- function(p) log(p / (1 - p))

random_effects_summary <- function(effect, se, weight = NULL) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- safe_num(weight)
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L,
      effect = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      z = NA_real_,
      p = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }

  yi <- effect[keep]
  vi <- se[keep]^2
  external_weight <- weight[keep] / mean(weight[keep])
  wi <- external_weight / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1L
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q - df) / c_val) else 0
  wi_re <- external_weight / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  z <- mu / se_mu
  p <- 2 * pnorm(-abs(z))
  sign_ref <- sign(mu)
  direction_fraction <- if (sign_ref == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign_ref)
  }
  i2 <- if (q > 0 && df > 0) max(0, (q - df) / q) else 0
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    z = z,
    p = p,
    tau2 = tau2,
    i2 = i2,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

message("Count-aware RDG branch-allocation model")
message("  Input: grouped same-condition replicate branch counts")
message("  Model: independent beta posterior usage deltas plus random-effects summary")
message("  Note: this is a screening layer, not a full hierarchical beta-binomial model")

stage_message("Loading grouped branch contrast table")
dt <- fread(contrast_file, showProgress = FALSE)
if (nrow(dt) == 0) {
  warning("Grouped branch contrast table is empty.")
  fwrite(data.table(), file.path(
    output_dir,
    "dominant_rdg_branch_allocation_posterior_contrasts.csv"
  ))
  fwrite(data.table(), file.path(
    output_dir,
    "dominant_rdg_branch_allocation_feature_summary.csv"
  ))
  fwrite(data.table(), file.path(
    output_dir,
    "dominant_rdg_branch_allocation_gene_design_summary.csv"
  ))
  fwrite(data.table(), file.path(
    output_dir,
    "dominant_rdg_branch_allocation_review.csv"
  ))
  fwrite(data.table(metric = "posterior_contrast_rows", value = 0),
         file.path(output_dir,
                   "dominant_rdg_branch_allocation_summary_metrics.csv"))
  quit(save = "no", status = 0)
}

numeric_cols <- c(
  "case_sum_raw_counts", "control_sum_raw_counts",
  "case_branch_total_counts", "control_branch_total_counts",
  "case_other_branch_counts", "control_other_branch_counts",
  "n_case_runs_merged", "n_control_runs_merged",
  "grouped_min_branch_total_counts", "grouped_max_feature_counts",
  "grouped_min_feature_counts", "grouped_row_weight"
)
for (column in intersect(numeric_cols, names(dt))) {
  dt[, (column) := safe_num(get(column))]
}

required_cols <- c(
  "gene_symbol", "tx_id", "design_family", "feature_id",
  "branch_feature_class", "case_sum_raw_counts", "control_sum_raw_counts",
  "case_branch_total_counts", "control_branch_total_counts"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols)) {
  stop("Grouped branch contrast table is missing columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

alpha <- 0.5
dt[, case_other_branch_counts :=
     pmax(case_branch_total_counts - case_sum_raw_counts, 0)]
dt[, control_other_branch_counts :=
     pmax(control_branch_total_counts - control_sum_raw_counts, 0)]

dt[, `:=`(
  case_beta_a = case_sum_raw_counts + alpha,
  case_beta_b = case_other_branch_counts + alpha,
  control_beta_a = control_sum_raw_counts + alpha,
  control_beta_b = control_other_branch_counts + alpha
)]
dt[, `:=`(
  case_beta_usage_mean = case_beta_a / (case_beta_a + case_beta_b),
  control_beta_usage_mean = control_beta_a /
    (control_beta_a + control_beta_b),
  case_beta_usage_var = beta_var(case_beta_a, case_beta_b),
  control_beta_usage_var = beta_var(control_beta_a, control_beta_b)
)]
dt[, `:=`(
  count_aware_usage_delta =
    case_beta_usage_mean - control_beta_usage_mean,
  count_aware_usage_delta_se =
    sqrt(case_beta_usage_var + control_beta_usage_var),
  count_aware_log_or =
    logit(case_beta_usage_mean) - logit(control_beta_usage_mean)
)]
dt[, `:=`(
  count_aware_usage_delta_lower =
    count_aware_usage_delta - 1.96 * count_aware_usage_delta_se,
  count_aware_usage_delta_upper =
    count_aware_usage_delta + 1.96 * count_aware_usage_delta_se,
  count_aware_usage_delta_z =
    count_aware_usage_delta / count_aware_usage_delta_se
)]
dt[, count_aware_usage_delta_p :=
     2 * pnorm(-abs(count_aware_usage_delta_z))]
dt[!is.finite(count_aware_usage_delta_p),
   count_aware_usage_delta_p := NA_real_]
dt[, count_aware_usage_delta_q :=
     p.adjust(count_aware_usage_delta_p, method = "BH")]

dt[, `:=`(
  count_aware_min_runs_merged =
    pmin(n_case_runs_merged, n_control_runs_merged),
  count_aware_min_branch_total_counts =
    pmin(case_branch_total_counts, control_branch_total_counts),
  count_aware_max_feature_counts =
    pmax(case_sum_raw_counts, control_sum_raw_counts),
  count_aware_min_feature_counts =
    pmin(case_sum_raw_counts, control_sum_raw_counts)
)]
dt[, `:=`(
  count_aware_passes_replicate_gate =
    count_aware_min_runs_merged >= 2,
  count_aware_passes_total_gate =
    count_aware_min_branch_total_counts >= 100,
  count_aware_passes_feature_gate_any =
    count_aware_max_feature_counts >= 30,
  count_aware_passes_feature_gate_both =
    count_aware_min_feature_counts >= 30
)]
dt[, count_aware_evaluable :=
     count_aware_passes_replicate_gate &
     count_aware_passes_total_gate &
     count_aware_passes_feature_gate_any &
     is.finite(count_aware_usage_delta) &
     is.finite(count_aware_usage_delta_se) &
     count_aware_usage_delta_se > 0]
dt[, count_aware_strict :=
     count_aware_evaluable &
     count_aware_passes_feature_gate_both &
     match_scope == "exact_context"]
dt[, count_aware_row_weight := pmax(
  1,
  log1p(count_aware_max_feature_counts) *
    sqrt(pmax(count_aware_min_branch_total_counts, 1)) *
    fifelse(match_scope == "exact_context", 1, 0.55)
)]
dt[!is.finite(count_aware_row_weight), count_aware_row_weight := 1]

posterior_file <- file.path(
  output_dir,
  "dominant_rdg_branch_allocation_posterior_contrasts.csv"
)
fwrite(dt, posterior_file)
stage_message("Wrote posterior contrast table")

feature_summary <- dt[, {
  re_strict <- random_effects_summary(
    count_aware_usage_delta[count_aware_strict == TRUE],
    count_aware_usage_delta_se[count_aware_strict == TRUE],
    count_aware_row_weight[count_aware_strict == TRUE]
  )
  re_all <- random_effects_summary(
    count_aware_usage_delta[count_aware_evaluable == TRUE],
    count_aware_usage_delta_se[count_aware_evaluable == TRUE],
    count_aware_row_weight[count_aware_evaluable == TRUE]
  )
  use_re <- if (re_strict$n[[1]] >= 2L) re_strict else re_all
  eval_rows <- count_aware_evaluable == TRUE
  strict_rows <- count_aware_strict == TRUE
  data.table(
    branch_feature_class = first_nonempty(branch_feature_class),
    feature_type = first_nonempty(feature_type),
    category = first_nonempty(category),
    n_count_aware_contrast_rows = .N,
    n_count_aware_evaluable_rows = sum(eval_rows, na.rm = TRUE),
    n_count_aware_strict_rows = sum(strict_rows, na.rm = TRUE),
    n_count_aware_strict_studies = uniqueN(study[strict_rows]),
    n_count_aware_strict_cell_lines = uniqueN(CELL_LINE[strict_rows]),
    n_count_aware_evaluable_studies = uniqueN(study[eval_rows]),
    n_count_aware_case_conditions = uniqueN(case_condition),
    total_case_branch_counts =
      sum(case_branch_total_counts[eval_rows], na.rm = TRUE),
    total_control_branch_counts =
      sum(control_branch_total_counts[eval_rows], na.rm = TRUE),
    max_case_feature_count = max_safe(case_sum_raw_counts),
    max_control_feature_count = max_safe(control_sum_raw_counts),
    count_aware_weighted_usage_delta = weighted_mean_safe(
      count_aware_usage_delta[eval_rows],
      count_aware_row_weight[eval_rows]
    ),
    count_aware_median_usage_delta =
      median_safe(count_aware_usage_delta[eval_rows]),
    count_aware_max_abs_usage_delta =
      max_safe(abs(count_aware_usage_delta[eval_rows])),
    count_aware_fraction_positive_usage_delta =
      mean(count_aware_usage_delta[eval_rows] > 0, na.rm = TRUE),
    count_aware_random_effects_delta = use_re$effect,
    count_aware_random_effects_delta_se = use_re$se,
    count_aware_random_effects_delta_lower = use_re$lower,
    count_aware_random_effects_delta_upper = use_re$upper,
    count_aware_random_effects_delta_p = use_re$p,
    count_aware_random_effects_delta_i2 = use_re$i2,
    count_aware_direction_fraction = use_re$direction_fraction,
    count_aware_interval_supported = use_re$ci_excludes_zero,
    count_aware_interval_source =
      if (re_strict$n[[1]] >= 2L) {
        "strict_exact_random_effects"
      } else {
        "all_evaluable_random_effects"
      },
    strongest_count_aware_case_condition =
      case_condition[which.max(abs(count_aware_usage_delta))][1],
    strongest_count_aware_control_conditions =
      control_conditions[which.max(abs(count_aware_usage_delta))][1],
    strongest_count_aware_study =
      study[which.max(abs(count_aware_usage_delta))][1],
    nearest_count_aware_contexts =
      collapse_examples(sample_set_signature, 5L)
  )
}, by = .(gene_symbol, tx_id, design_family, feature_id)]

feature_summary[, count_aware_random_effects_delta_q :=
                  p.adjust(count_aware_random_effects_delta_p, method = "BH")]
feature_summary[, count_aware_branch_allocation_score :=
                  0.22 * pmin(n_count_aware_strict_rows / 3, 1) +
                  0.18 * pmin(n_count_aware_strict_studies / 2, 1) +
                  0.18 * pmin(n_count_aware_evaluable_rows / 4, 1) +
                  0.20 * pmin(abs(count_aware_weighted_usage_delta) / 0.12, 1) +
                  0.16 * fifelse(
                    is.finite(count_aware_random_effects_delta_q),
                    pmin(-log10(pmax(
                      count_aware_random_effects_delta_q,
                      1e-300
                    )) / 4, 1),
                    0
                  ) +
                  0.04 * as.numeric(
                    count_aware_interval_supported == TRUE
                  ) +
                  0.02 * pmin(log1p(pmax(max_case_feature_count,
                                         max_control_feature_count)) /
                                log1p(1000), 1)]
feature_summary[!is.finite(count_aware_branch_allocation_score),
                count_aware_branch_allocation_score := 0]

feature_summary[, count_aware_support_class := fifelse(
  n_count_aware_strict_rows >= 2 &
    n_count_aware_strict_studies >= 2 &
    count_aware_interval_supported == TRUE &
    is.finite(count_aware_random_effects_delta_q) &
    count_aware_random_effects_delta_q <= 0.10 &
    abs(count_aware_weighted_usage_delta) >= 0.04 &
    count_aware_direction_fraction >= 0.70,
  "replicated_count_aware_branch_shift",
  fifelse(
    n_count_aware_strict_rows >= 1 &
      is.finite(count_aware_random_effects_delta_p) &
      count_aware_random_effects_delta_p <= 0.01 &
      abs(count_aware_weighted_usage_delta) >= 0.08,
    "single_count_aware_branch_shift",
    fifelse(
      n_count_aware_evaluable_rows >= 2 &
        count_aware_interval_supported == TRUE &
        is.finite(count_aware_random_effects_delta_q) &
        count_aware_random_effects_delta_q <= 0.10 &
        abs(count_aware_weighted_usage_delta) >= 0.04,
      "count_aware_branch_shift_review",
      fifelse(
        n_count_aware_evaluable_rows > 0,
        "count_aware_evaluable_no_strong_shift",
        "low_count_or_no_grouped_control"
      )
    )
  )
)]

setorder(feature_summary, -count_aware_branch_allocation_score,
         count_aware_random_effects_delta_q, gene_symbol, design_family,
         feature_id)
feature_summary_file <- file.path(
  output_dir,
  "dominant_rdg_branch_allocation_feature_summary.csv"
)
fwrite(feature_summary, feature_summary_file)
stage_message("Wrote feature summary")

top_feature <- feature_summary[, .SD[1],
                               by = .(gene_symbol, tx_id, design_family)]
class_wide <- dcast(
  feature_summary[
    branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                  "overlapping_uORF")
  ],
  gene_symbol + tx_id + design_family ~ branch_feature_class,
  value.var = c(
    "feature_id",
    "count_aware_weighted_usage_delta",
    "count_aware_random_effects_delta",
    "count_aware_random_effects_delta_q",
    "count_aware_branch_allocation_score"
  ),
  fun.aggregate = function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) {
      if (is.numeric(x)) return(NA_real_)
      return(NA_character_)
    }
    if (is.numeric(x)) x[which.max(abs(x))] else as.character(x[[1]])
  },
  fill = NA
)

gene_design_summary <- merge(
  top_feature[, .(
    gene_symbol, tx_id, design_family,
    count_aware_top_feature_id = feature_id,
    count_aware_top_feature_class = branch_feature_class,
    count_aware_top_feature_type = feature_type,
    count_aware_top_support_class = count_aware_support_class,
    count_aware_top_score = count_aware_branch_allocation_score,
    count_aware_top_weighted_usage_delta =
      count_aware_weighted_usage_delta,
    count_aware_top_random_effects_delta =
      count_aware_random_effects_delta,
    count_aware_top_delta_lower =
      count_aware_random_effects_delta_lower,
    count_aware_top_delta_upper =
      count_aware_random_effects_delta_upper,
    count_aware_top_q = count_aware_random_effects_delta_q,
    count_aware_top_i2 = count_aware_random_effects_delta_i2,
    count_aware_top_direction_fraction =
      count_aware_direction_fraction,
    count_aware_top_strict_rows = n_count_aware_strict_rows,
    count_aware_top_strict_studies = n_count_aware_strict_studies,
    count_aware_top_evaluable_rows = n_count_aware_evaluable_rows,
    count_aware_top_context = strongest_count_aware_case_condition,
    count_aware_top_control = strongest_count_aware_control_conditions,
    count_aware_top_study = strongest_count_aware_study
  )],
  class_wide,
  by = c("gene_symbol", "tx_id", "design_family"),
  all.x = TRUE,
  sort = FALSE
)
setorder(gene_design_summary, -count_aware_top_score,
         count_aware_top_q, gene_symbol, design_family)
gene_design_file <- file.path(
  output_dir,
  "dominant_rdg_branch_allocation_gene_design_summary.csv"
)
fwrite(gene_design_summary, gene_design_file)
stage_message("Wrote gene-design summary")

review <- gene_design_summary[
  count_aware_top_support_class != "low_count_or_no_grouped_control" &
    (count_aware_top_support_class !=
       "count_aware_evaluable_no_strong_shift" |
       count_aware_top_score >= 0.55)
]
setorder(review, -count_aware_top_score, count_aware_top_q,
         gene_symbol, design_family)
review_file <- file.path(
  output_dir,
  "dominant_rdg_branch_allocation_review.csv"
)
fwrite(review, review_file)
stage_message("Wrote review table")

summary_metrics <- data.table(
  metric = c(
    "posterior_contrast_rows",
    "posterior_evaluable_rows",
    "posterior_strict_rows",
    "feature_summary_rows",
    "gene_design_rows",
    "review_rows",
    "replicated_count_aware_branch_shift_rows",
    "single_count_aware_branch_shift_rows",
    "count_aware_branch_shift_review_rows",
    "viral_infection_review_rows",
    "viral_replicated_count_aware_rows",
    "genes",
    "design_families"
  ),
  value = c(
    nrow(dt),
    sum(dt$count_aware_evaluable, na.rm = TRUE),
    sum(dt$count_aware_strict, na.rm = TRUE),
    nrow(feature_summary),
    nrow(gene_design_summary),
    nrow(review),
    sum(feature_summary$count_aware_support_class ==
          "replicated_count_aware_branch_shift", na.rm = TRUE),
    sum(feature_summary$count_aware_support_class ==
          "single_count_aware_branch_shift", na.rm = TRUE),
    sum(feature_summary$count_aware_support_class ==
          "count_aware_branch_shift_review", na.rm = TRUE),
    nrow(review[design_family == "viral_infection"]),
    nrow(review[
      design_family == "viral_infection" &
        count_aware_top_support_class ==
        "replicated_count_aware_branch_shift"
    ]),
    uniqueN(dt$gene_symbol),
    uniqueN(dt$design_family)
  )
)
metrics_file <- file.path(
  output_dir,
  "dominant_rdg_branch_allocation_summary_metrics.csv"
)
fwrite(summary_metrics, metrics_file)

plot_dt <- feature_summary[
  count_aware_support_class != "low_count_or_no_grouped_control" &
    count_aware_branch_allocation_score >= 0.45
][order(-count_aware_branch_allocation_score)][seq_len(min(.N, 45))]
if (nrow(plot_dt) > 0) {
  plot_dt[, plot_label := paste0(gene_symbol, " / ", design_family,
                                 " / ", feature_id)]
  plot_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_top <- ggplot(
    plot_dt,
    aes(
      x = count_aware_weighted_usage_delta,
      y = plot_label,
      color = branch_feature_class,
      shape = count_aware_support_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Design family: ", design_family,
        "<br>Feature: ", feature_id,
        "<br>Class: ", branch_feature_class,
        "<br>Posterior usage delta: ",
        round(count_aware_weighted_usage_delta, 3),
        "<br>RE delta: ", round(count_aware_random_effects_delta, 3),
        "<br>q: ", signif(count_aware_random_effects_delta_q, 3),
        "<br>Score: ", round(count_aware_branch_allocation_score, 3),
        "<br>Strict rows: ", n_count_aware_strict_rows,
        "<br>Strict studies: ", n_count_aware_strict_studies
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "#7b8494", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = count_aware_random_effects_delta_lower,
          xmax = count_aware_random_effects_delta_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.35,
      alpha = 0.45,
      show.legend = FALSE
    ) +
    geom_point(size = 2.6, alpha = 0.9) +
    scale_color_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31"
    )) +
    scale_shape_manual(values = c(
      replicated_count_aware_branch_shift = 17,
      single_count_aware_branch_shift = 15,
      count_aware_branch_shift_review = 18,
      count_aware_evaluable_no_strong_shift = 16
    ), labels = c(
      replicated_count_aware_branch_shift = "replicated",
      single_count_aware_branch_shift = "single",
      count_aware_branch_shift_review = "review",
      count_aware_evaluable_no_strong_shift = "evaluable"
    )) +
    labs(
      title = "Count-aware RDG branch-allocation shifts",
      subtitle = "Beta-posterior case/control usage deltas over merged replicate groups",
      x = "Case-control branch usage delta",
      y = "",
      color = "Feature class",
      shape = "Support"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir,
                   "dominant_rdg_branch_allocation_top_shifts.png"),
         p_top, width = 11.2, height = 8.8, dpi = 180)
  ggsave(file.path(figure_dir,
                   "dominant_rdg_branch_allocation_top_shifts.pdf"),
         p_top, width = 11.2, height = 8.8)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir,
                "dominant_rdg_branch_allocation_top_shifts.html"),
      title = "Count-aware RDG branch-allocation shifts"
    )
  }
}

viral_dt <- feature_summary[
  design_family == "viral_infection" &
    n_count_aware_evaluable_rows > 0
][order(-count_aware_branch_allocation_score)][seq_len(min(.N, 45))]
if (nrow(viral_dt) > 0) {
  viral_dt[, plot_label := paste0(gene_symbol, " / ", feature_id)]
  viral_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_viral <- ggplot(
    viral_dt,
    aes(
      x = count_aware_weighted_usage_delta,
      y = plot_label,
      color = branch_feature_class,
      size = pmin(count_aware_branch_allocation_score, 1),
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Feature: ", feature_id,
        "<br>Class: ", branch_feature_class,
        "<br>Posterior usage delta: ",
        round(count_aware_weighted_usage_delta, 3),
        "<br>q: ", signif(count_aware_random_effects_delta_q, 3),
        "<br>Support: ", count_aware_support_class,
        "<br>Top context: ", strongest_count_aware_case_condition,
        "<br>Study: ", strongest_count_aware_study
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "#7b8494", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = count_aware_random_effects_delta_lower,
          xmax = count_aware_random_effects_delta_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.35,
      alpha = 0.45,
      show.legend = FALSE
    ) +
    geom_point(alpha = 0.9) +
    scale_size(range = c(1.8, 5), guide = "none") +
    scale_color_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31"
    )) +
    labs(
      title = "Viral-infection count-aware branch-allocation review",
      subtitle = "Positive x means higher feature share in infected/case groups",
      x = "Case-control branch usage delta",
      y = "",
      color = "Feature class"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir,
                   "dominant_rdg_branch_allocation_viral_review.png"),
         p_viral, width = 10.6, height = 8.2, dpi = 180)
  ggsave(file.path(figure_dir,
                   "dominant_rdg_branch_allocation_viral_review.pdf"),
         p_viral, width = 10.6, height = 8.2)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_viral,
      file.path(figure_dir,
                "dominant_rdg_branch_allocation_viral_review.html"),
      title = "Viral-infection count-aware branch-allocation review"
    )
  }
}

message("Saved posterior contrasts: ", posterior_file)
message("Saved feature summary: ", feature_summary_file)
message("Saved gene-design summary: ", gene_design_file)
message("Saved review table: ", review_file)
message("Saved summary metrics: ", metrics_file)
