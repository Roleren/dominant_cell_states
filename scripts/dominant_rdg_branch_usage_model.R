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
output_dir <- file.path(analysis_dir, "dominant_rdg_branch_usage")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

message("DOTSeq-inspired RDG branch-usage model")
message("  1. Reuse matched case/control RDG feature counts")
message("  2. Convert each feature into a feature-vs-other branch usage contrast")
message("  3. Aggregate contrasts by gene x design family x feature")
message("  4. Save review tables and compact figures")

input_file <- file.path(
  analysis_dir, "dominant_next_model",
  "dominant_next_model_matched_control_feature_contrasts.csv"
)
if (!file.exists(input_file)) {
  stop("Missing matched feature contrast input: ", input_file, call. = FALSE)
}

first_nonempty <- function(x, default = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else default
}

collapse_unique <- function(x, n = 5L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

weighted_mean_safe <- function(x, w = NULL) {
  x <- suppressWarnings(as.numeric(x))
  if (is.null(w)) w <- rep(1, length(x))
  w <- suppressWarnings(as.numeric(w))
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

median_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

max_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

bool_col <- function(dt, col, default = FALSE) {
  if (!col %in% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[col]]
  if (is.logical(x)) return(fifelse(is.na(x), default, x))
  if (is.numeric(x)) return(fifelse(is.na(x), default, x != 0))
  y <- tolower(as.character(x))
  fifelse(is.na(y), default, y %chin% c("true", "t", "1", "yes"))
}

random_effects_summary <- function(effect, se, weight = NULL) {
  effect <- suppressWarnings(as.numeric(effect))
  se <- suppressWarnings(as.numeric(se))
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- suppressWarnings(as.numeric(weight))
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
  df <- length(yi) - 1
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

branch_feature_class <- function(feature_type, category) {
  fifelse(
    feature_type == "clean_CDS" | category == "clean_CDS",
    "clean_CDS",
    fifelse(
      feature_type == "overlapping_uorf" |
        category %chin% c("uoORF", "ouORF"),
      "overlapping_uORF",
      fifelse(
        feature_type == "leader_uorf" | category == "uORF",
        "leader_uORF",
        fifelse(
          grepl("internal|iorf", feature_type, ignore.case = TRUE) |
            grepl("iORF", category, fixed = TRUE),
          "iORF",
          fifelse(
            grepl("nte|ntt", feature_type, ignore.case = TRUE) |
              category %chin% c("NTE", "NTT"),
            "NTE_NTT",
            "other_feature"
          )
        )
      )
    )
  )
}

logit <- function(p) log(p / (1 - p))

feature_dt <- fread(input_file, showProgress = FALSE)
if (!"design_family" %in% names(feature_dt)) {
  feature_dt[, design_family := perturbation_class]
}
feature_dt[is.na(design_family) | !nzchar(design_family),
           design_family := "unknown"]

numeric_cols <- c(
  "n_case_samples", "n_control_samples", "case_sum_raw_counts",
  "control_sum_raw_counts", "case_mean_relative_use",
  "control_mean_relative_use", "matched_relative_use_delta",
  "matched_relative_use_delta_se", "min_feature_raw_counts",
  "max_feature_raw_counts"
)
for (column in intersect(numeric_cols, names(feature_dt))) {
  feature_dt[, (column) := suppressWarnings(as.numeric(get(column)))]
}

feature_dt[, branch_feature_class :=
             branch_feature_class(feature_type, category)]
feature_dt <- feature_dt[
  branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                "overlapping_uORF", "iORF", "NTE_NTT")
]

design_cols <- intersect(c(
  "gene_symbol", "tx_id", "stratum_id", "study", "AUTHOR", "CELL_LINE",
  "TISSUE", "GENE", "case_condition", "control_conditions",
  "perturbation_class", "design_family", "sample_set_signature",
  "case_run_signature", "control_run_signature"
), names(feature_dt))

feature_dt[, contrast_id := do.call(
  paste,
  c(.SD, sep = " | ")
), .SDcols = design_cols]

totals <- feature_dt[, .(
  case_branch_total_counts = sum(case_sum_raw_counts, na.rm = TRUE),
  control_branch_total_counts = sum(control_sum_raw_counts, na.rm = TRUE),
  n_branch_features_in_contrast = .N,
  n_allocation_branch_features = sum(
    branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                  "overlapping_uORF"),
    na.rm = TRUE
  )
), by = contrast_id]
feature_dt <- merge(feature_dt, totals, by = "contrast_id",
                    all.x = TRUE, sort = FALSE)

feature_dt[, case_other_branch_counts :=
             pmax(case_branch_total_counts - case_sum_raw_counts, 0)]
feature_dt[, control_other_branch_counts :=
             pmax(control_branch_total_counts - control_sum_raw_counts, 0)]

alpha <- 0.5
feature_dt[, case_branch_usage_proportion :=
             (case_sum_raw_counts + alpha) /
             (case_branch_total_counts + 2 * alpha)]
feature_dt[, control_branch_usage_proportion :=
             (control_sum_raw_counts + alpha) /
             (control_branch_total_counts + 2 * alpha)]
feature_dt[, branch_usage_delta_count_model :=
             case_branch_usage_proportion -
             control_branch_usage_proportion]
feature_dt[, branch_usage_log_odds_ratio :=
             logit(case_branch_usage_proportion) -
             logit(control_branch_usage_proportion)]
feature_dt[, branch_usage_log_or_se := sqrt(
  1 / (case_sum_raw_counts + alpha) +
    1 / (case_other_branch_counts + alpha) +
    1 / (control_sum_raw_counts + alpha) +
    1 / (control_other_branch_counts + alpha)
)]
feature_dt[, branch_usage_z :=
             branch_usage_log_odds_ratio / branch_usage_log_or_se]
feature_dt[, branch_usage_p := 2 * pnorm(-abs(branch_usage_z))]
feature_dt[, branch_usage_q := p.adjust(branch_usage_p, method = "BH")]

feature_dt[, branch_usage_passes_sample_gate :=
             n_case_samples >= 2 & n_control_samples >= 2]
feature_dt[, branch_usage_passes_total_count_gate :=
             pmin(case_branch_total_counts,
                  control_branch_total_counts,
                  na.rm = TRUE) >= 100]
feature_dt[, branch_usage_passes_feature_count_gate_any :=
             pmax(case_sum_raw_counts,
                  control_sum_raw_counts,
                  na.rm = TRUE) >= 30]
feature_dt[, branch_usage_passes_feature_count_gate_both :=
             pmin(case_sum_raw_counts,
                  control_sum_raw_counts,
                  na.rm = TRUE) >= 30]
feature_dt[, branch_usage_evaluable :=
             branch_usage_passes_sample_gate &
             branch_usage_passes_total_count_gate &
             branch_usage_passes_feature_count_gate_any &
             is.finite(branch_usage_log_odds_ratio) &
             is.finite(branch_usage_log_or_se)]
feature_dt[, branch_usage_strict :=
             branch_usage_evaluable &
             branch_usage_passes_feature_count_gate_both]
feature_dt[, branch_usage_row_weight := pmax(
  1,
  log1p(pmax(case_sum_raw_counts, control_sum_raw_counts, na.rm = TRUE)) *
    pmin(n_case_samples, n_control_samples, na.rm = TRUE)
)]
feature_dt[!is.finite(branch_usage_row_weight),
           branch_usage_row_weight := 1]

fwrite(
  feature_dt,
  file.path(output_dir, "dominant_rdg_branch_usage_contrast_table.csv")
)

feature_summary <- feature_dt[, {
  re_all <- random_effects_summary(
    branch_usage_log_odds_ratio,
    branch_usage_log_or_se,
    branch_usage_row_weight
  )
  re_strict <- random_effects_summary(
    branch_usage_log_odds_ratio[branch_usage_strict == TRUE],
    branch_usage_log_or_se[branch_usage_strict == TRUE],
    branch_usage_row_weight[branch_usage_strict == TRUE]
  )
  use_re <- if (re_strict$n[[1]] >= 2L) re_strict else re_all
  data.table(
    branch_feature_class = first_nonempty(branch_feature_class),
    feature_type = first_nonempty(feature_type),
    category = first_nonempty(category),
    n_design_rows = .N,
    n_studies = uniqueN(study),
    n_cell_lines = uniqueN(CELL_LINE),
    n_tissues = uniqueN(TISSUE),
    n_case_conditions = uniqueN(case_condition),
    branch_usage_evaluable_rows = sum(branch_usage_evaluable, na.rm = TRUE),
    branch_usage_strict_rows = sum(branch_usage_strict, na.rm = TRUE),
    branch_usage_strict_studies = uniqueN(study[branch_usage_strict == TRUE]),
    branch_usage_strict_cell_lines =
      uniqueN(CELL_LINE[branch_usage_strict == TRUE]),
    total_case_branch_counts = sum(case_branch_total_counts, na.rm = TRUE),
    total_control_branch_counts =
      sum(control_branch_total_counts, na.rm = TRUE),
    max_case_feature_count = max_safe(case_sum_raw_counts),
    max_control_feature_count = max_safe(control_sum_raw_counts),
    weighted_usage_delta = weighted_mean_safe(
      branch_usage_delta_count_model,
      branch_usage_row_weight
    ),
    median_usage_delta = median_safe(branch_usage_delta_count_model),
    max_abs_usage_delta = max_safe(abs(branch_usage_delta_count_model)),
    fraction_positive_usage_delta =
      mean(branch_usage_delta_count_model > 0, na.rm = TRUE),
    branch_usage_log_or = use_re$effect,
    branch_usage_log_or_se = use_re$se,
    branch_usage_log_or_lower = use_re$lower,
    branch_usage_log_or_upper = use_re$upper,
    branch_usage_log_or_p = use_re$p,
    branch_usage_log_or_i2 = use_re$i2,
    branch_usage_direction_fraction = use_re$direction_fraction,
    branch_usage_interval_supported = use_re$ci_excludes_zero,
    branch_usage_interval_source =
      if (re_strict$n[[1]] >= 2L) "strict_random_effects" else "all_rows_random_effects",
    strongest_case_condition =
      case_condition[which.max(abs(branch_usage_delta_count_model))][1],
    strongest_control_conditions =
      control_conditions[which.max(abs(branch_usage_delta_count_model))][1],
    strongest_study = study[which.max(abs(branch_usage_delta_count_model))][1],
    nearest_contexts = collapse_unique(sample_set_signature, n = 5L)
  )
}, by = .(gene_symbol, tx_id, design_family, feature_id)]

feature_summary[, branch_usage_log_or_q :=
                  p.adjust(branch_usage_log_or_p, method = "BH")]
feature_summary[, branch_usage_score :=
                  0.24 * pmin(branch_usage_strict_rows / 3, 1) +
                  0.18 * pmin(branch_usage_strict_studies / 2, 1) +
                  0.18 * pmin(abs(weighted_usage_delta) / 0.20, 1) +
                  0.18 * fifelse(is.finite(branch_usage_log_or_q),
                                  pmin(-log10(pmax(branch_usage_log_or_q,
                                                   1e-300)) / 4, 1),
                                  0) +
                  0.12 * as.numeric(branch_usage_interval_supported == TRUE) +
                  0.10 * pmin(log1p(pmax(max_case_feature_count,
                                         max_control_feature_count,
                                         na.rm = TRUE)) / log1p(1000), 1)]
feature_summary[!is.finite(branch_usage_score), branch_usage_score := 0]

feature_summary[, branch_usage_support_class := fifelse(
  branch_usage_strict_rows >= 2 &
    branch_usage_strict_studies >= 2 &
    branch_usage_interval_supported == TRUE &
    is.finite(branch_usage_log_or_q) &
    branch_usage_log_or_q <= 0.10 &
    abs(weighted_usage_delta) >= 0.05 &
    branch_usage_direction_fraction >= 0.70,
  "replicated_branch_usage_shift",
  fifelse(
    branch_usage_strict_rows >= 1 &
      is.finite(branch_usage_log_or_p) &
      branch_usage_log_or_p <= 0.01 &
      abs(weighted_usage_delta) >= 0.10,
    "single_design_branch_usage_shift",
    fifelse(
      branch_usage_evaluable_rows > 0,
      "evaluable_no_strong_shift",
      "low_count_or_unbalanced"
    )
  )
)]
setorder(feature_summary, -branch_usage_score, branch_usage_log_or_q,
         gene_symbol, design_family, feature_id)

fwrite(
  feature_summary,
  file.path(output_dir, "dominant_rdg_branch_usage_feature_summary.csv")
)

top_feature <- feature_summary[, .SD[1], by = .(gene_symbol, tx_id, design_family)]
class_wide <- dcast(
  feature_summary[
    branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                  "overlapping_uORF")
  ],
  gene_symbol + tx_id + design_family ~ branch_feature_class,
  value.var = c("feature_id", "weighted_usage_delta",
                "branch_usage_log_or", "branch_usage_log_or_q",
                "branch_usage_score"),
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
    branch_usage_top_feature_id = feature_id,
    branch_usage_top_feature_class = branch_feature_class,
    branch_usage_top_feature_type = feature_type,
    branch_usage_top_support_class = branch_usage_support_class,
    branch_usage_top_score = branch_usage_score,
    branch_usage_top_weighted_usage_delta = weighted_usage_delta,
    branch_usage_top_log_or = branch_usage_log_or,
    branch_usage_top_q = branch_usage_log_or_q,
    branch_usage_top_i2 = branch_usage_log_or_i2,
    branch_usage_top_direction_fraction = branch_usage_direction_fraction,
    branch_usage_top_strict_rows = branch_usage_strict_rows,
    branch_usage_top_strict_studies = branch_usage_strict_studies,
    branch_usage_top_context = strongest_case_condition,
    branch_usage_top_control = strongest_control_conditions,
    branch_usage_top_study = strongest_study
  )],
  class_wide,
  by = c("gene_symbol", "tx_id", "design_family"),
  all.x = TRUE,
  sort = FALSE
)
setorder(gene_design_summary, -branch_usage_top_score,
         branch_usage_top_q, gene_symbol, design_family)

fwrite(
  gene_design_summary,
  file.path(output_dir, "dominant_rdg_branch_usage_gene_design_summary.csv")
)

review <- gene_design_summary[
  branch_usage_top_support_class != "low_count_or_unbalanced" &
    (branch_usage_top_support_class != "evaluable_no_strong_shift" |
       branch_usage_top_score >= 0.55)
]
setorder(review, -branch_usage_top_score, branch_usage_top_q,
         gene_symbol, design_family)
fwrite(
  review,
  file.path(output_dir, "dominant_rdg_branch_usage_review.csv")
)

summary_metrics <- data.table(
  metric = c(
    "contrast_rows",
    "feature_summary_rows",
    "gene_design_rows",
    "review_rows",
    "replicated_branch_usage_shift_rows",
    "single_design_branch_usage_shift_rows",
    "viral_infection_review_rows",
    "genes",
    "design_families"
  ),
  value = c(
    nrow(feature_dt),
    nrow(feature_summary),
    nrow(gene_design_summary),
    nrow(review),
    sum(feature_summary$branch_usage_support_class ==
          "replicated_branch_usage_shift", na.rm = TRUE),
    sum(feature_summary$branch_usage_support_class ==
          "single_design_branch_usage_shift", na.rm = TRUE),
    nrow(review[design_family == "viral_infection"]),
    uniqueN(feature_dt$gene_symbol),
    uniqueN(feature_dt$design_family)
  )
)
fwrite(
  summary_metrics,
  file.path(output_dir, "dominant_rdg_branch_usage_summary_metrics.csv")
)

if (nrow(review) > 0) {
  plot_dt <- feature_summary[
    branch_usage_support_class != "low_count_or_unbalanced" &
      branch_usage_score >= 0.45
  ][order(-branch_usage_score)][seq_len(min(.N, 45))]
  if (nrow(plot_dt) > 0) {
    plot_dt[, plot_label := paste0(gene_symbol, " / ", design_family,
                                   " / ", feature_id)]
    plot_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
    plot_dt[, support_label := fifelse(
      branch_usage_support_class == "replicated_branch_usage_shift",
      "replicated",
      fifelse(
        branch_usage_support_class == "single_design_branch_usage_shift",
        "single design",
        "evaluable"
      )
    )]
    p_top <- ggplot(
      plot_dt,
      aes(
        x = weighted_usage_delta,
        y = plot_label,
        color = branch_feature_class,
        shape = support_label,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Design family: ", design_family,
          "<br>Feature: ", feature_id,
          "<br>Class: ", branch_feature_class,
          "<br>Usage delta: ", round(weighted_usage_delta, 3),
          "<br>logOR: ", round(branch_usage_log_or, 3),
          "<br>q: ", signif(branch_usage_log_or_q, 3),
          "<br>Score: ", round(branch_usage_score, 3),
          "<br>Strict rows: ", branch_usage_strict_rows,
          "<br>Studies: ", branch_usage_strict_studies
        )
      )
    ) +
      geom_vline(xintercept = 0, color = "#7b8494", linewidth = 0.35) +
      geom_point(size = 2.6, alpha = 0.9) +
      scale_color_manual(values = c(
        clean_CDS = "#1f5aa6",
        leader_uORF = "#c57b00",
        overlapping_uORF = "#9b1c31",
        iORF = "#5f3dc4",
        NTE_NTT = "#008060",
        other_feature = "#596273"
      )) +
      labs(
        title = "DOTSeq-like RDG branch-usage shifts",
        subtitle = "Feature-vs-other branch usage from matched Ribo-seq contrasts",
        x = "Weighted case-control branch usage delta",
        y = "",
        color = "Feature class",
        shape = "Support"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
    ggsave(
      file.path(figure_dir, "dominant_rdg_branch_usage_top_shifts.png"),
      p_top, width = 10.8, height = 8.8, dpi = 180
    )
    ggsave(
      file.path(figure_dir, "dominant_rdg_branch_usage_top_shifts.pdf"),
      p_top, width = 10.8, height = 8.8
    )
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_top,
        file.path(figure_dir, "dominant_rdg_branch_usage_top_shifts.html"),
        title = "DOTSeq-like RDG branch-usage shifts"
      )
    }
  }

  viral_dt <- feature_summary[
    design_family == "viral_infection" &
      branch_usage_evaluable_rows > 0
  ][order(-branch_usage_score)][seq_len(min(.N, 40))]
  if (nrow(viral_dt) > 0) {
    viral_dt[, plot_label := paste0(gene_symbol, " / ", feature_id)]
    viral_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
    p_viral <- ggplot(
      viral_dt,
      aes(
        x = weighted_usage_delta,
        y = plot_label,
        color = branch_feature_class,
        size = pmin(branch_usage_score, 1),
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Feature: ", feature_id,
          "<br>Class: ", branch_feature_class,
          "<br>Usage delta: ", round(weighted_usage_delta, 3),
          "<br>logOR q: ", signif(branch_usage_log_or_q, 3),
          "<br>Support: ", branch_usage_support_class,
          "<br>Top context: ", strongest_case_condition,
          "<br>Study: ", strongest_study
        )
      )
    ) +
      geom_vline(xintercept = 0, color = "#7b8494", linewidth = 0.35) +
      geom_point(alpha = 0.9) +
      scale_size(range = c(1.8, 5), guide = "none") +
      scale_color_manual(values = c(
        clean_CDS = "#1f5aa6",
        leader_uORF = "#c57b00",
        overlapping_uORF = "#9b1c31",
        iORF = "#5f3dc4",
        NTE_NTT = "#008060",
        other_feature = "#596273"
      )) +
      labs(
        title = "Viral-infection RDG branch-usage review",
        subtitle = "DOTSeq-like feature-vs-other branch shifts; positive means higher feature share in infected/case samples",
        x = "Weighted case-control branch usage delta",
        y = "",
        color = "Feature class"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
    ggsave(
      file.path(figure_dir, "dominant_rdg_branch_usage_viral_review.png"),
      p_viral, width = 10.4, height = 8.2, dpi = 180
    )
    ggsave(
      file.path(figure_dir, "dominant_rdg_branch_usage_viral_review.pdf"),
      p_viral, width = 10.4, height = 8.2
    )
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_viral,
        file.path(figure_dir, "dominant_rdg_branch_usage_viral_review.html"),
        title = "Viral-infection RDG branch-usage review"
      )
    }
  }
}

message("Saved branch-usage contrast table: ",
        file.path(output_dir, "dominant_rdg_branch_usage_contrast_table.csv"))
message("Saved branch-usage feature summary: ",
        file.path(output_dir, "dominant_rdg_branch_usage_feature_summary.csv"))
message("Saved branch-usage gene-design summary: ",
        file.path(output_dir, "dominant_rdg_branch_usage_gene_design_summary.csv"))
message("Saved branch-usage review table: ",
        file.path(output_dir, "dominant_rdg_branch_usage_review.csv"))
message("Saved branch-usage metrics: ",
        file.path(output_dir, "dominant_rdg_branch_usage_summary_metrics.csv"))
