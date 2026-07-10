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
output_dir <- file.path(analysis_dir, "dominant_rdg_context_partial_pooling")
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

clip01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- 0
  pmin(pmax(x, 0), 1)
}

num_col <- function(dt, column, default = NA_real_) {
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

bool_col <- function(dt, column, default = FALSE) {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[column]]
  if (is.logical(x)) {
    x[is.na(x)] <- default
    return(x)
  }
  if (is.numeric(x)) {
    out <- x != 0
    out[is.na(out)] <- default
    return(out)
  }
  y <- tolower(as.character(x))
  out <- y %chin% c("true", "t", "1", "yes")
  out[is.na(out)] <- default
  out
}

safe_max <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

safe_mean <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else mean(x)
}

collapse_unique <- function(x, n = 4L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

wrap_text <- function(x, width = 94L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

partial_pool_group <- function(dt) {
  y <- suppressWarnings(as.numeric(dt$interaction_delta))
  se <- suppressWarnings(as.numeric(dt$interaction_delta_se))
  keep <- is.finite(y) & is.finite(se) & se > 0
  out <- copy(dt)
  if (!any(keep)) {
    out[, `:=`(
      pooling_group_contexts = 0L,
      pooling_group_mu = NA_real_,
      pooling_group_mu_se = NA_real_,
      pooling_group_tau2 = NA_real_,
      partial_pooling_weight = 0,
      pooled_interaction_delta = NA_real_,
      pooled_interaction_se = NA_real_
    )]
    return(out)
  }

  yk <- y[keep]
  sek <- se[keep]
  w <- 1 / (sek^2)
  sum_w <- sum(w)
  mu <- sum(w * yk) / sum_w
  mu_se2 <- 1 / sum_w
  k <- length(yk)
  if (k > 1L) {
    q <- sum(w * (yk - mu)^2)
    denom <- sum_w - sum(w^2) / sum_w
    tau2 <- if (denom > 0) max(0, (q - (k - 1)) / denom) else 0
  } else {
    tau2 <- 0
  }
  shrink_w <- ifelse(is.finite(se) & se > 0,
                     tau2 / (tau2 + se^2),
                     0)
  shrink_w[!is.finite(shrink_w)] <- 0
  pooled <- mu + shrink_w * (y - mu)
  pooled_se <- sqrt((shrink_w^2) * se^2 + ((1 - shrink_w)^2) * mu_se2)
  pooled_se[!is.finite(pooled_se) | pooled_se <= 0] <- sqrt(mu_se2)

  out[, `:=`(
    pooling_group_contexts = k,
    pooling_group_mu = mu,
    pooling_group_mu_se = sqrt(mu_se2),
    pooling_group_tau2 = tau2,
    partial_pooling_weight = shrink_w,
    pooled_interaction_delta = pooled,
    pooled_interaction_se = pooled_se
  )]
  out
}

context_file <- file.path(
  analysis_dir,
  "dominant_rdg_context_interactions",
  "dominant_rdg_context_interaction_branch_context_effects.csv"
)
context_rows <- read_dt(context_file, required = TRUE)
if (!nrow(context_rows)) {
  stop("Context-interaction branch-context table is empty: ",
       context_file, call. = FALSE)
}

context_rows[, gene_symbol := as.character(gene_symbol)]
context_rows[, tx_id := as.character(tx_id)]
context_rows[, design_family := as.character(design_family)]
context_rows[, branch_class := as.character(branch_class)]
context_rows[, context_axis := as.character(context_axis)]
context_rows[, context_value := as.character(context_value)]
context_interaction_evaluable_vec <- bool_col(
  context_rows, "context_interaction_evaluable"
)
context_opposite_complement_vec <- bool_col(
  context_rows, "context_opposite_complement_direction"
)
interaction_delta_vec <- num_col(context_rows, "interaction_delta")
interaction_delta_se_vec <- num_col(context_rows, "interaction_delta_se")
context_rows[, `:=`(
  context_interaction_evaluable = context_interaction_evaluable_vec,
  context_opposite_complement_direction = context_opposite_complement_vec,
  interaction_delta = interaction_delta_vec,
  interaction_delta_se = interaction_delta_se_vec
)]
context_rows[, raw_context_supported := context_interaction_class %chin% c(
  "context_direction_reversal",
  "context_specific_supported_shift",
  "context_magnitude_modulated"
)]

pool_input <- context_rows[
  context_interaction_evaluable == TRUE &
    is.finite(interaction_delta) &
    is.finite(interaction_delta_se) &
    interaction_delta_se > 0
]
if (!nrow(pool_input)) {
  stop("No evaluable context-interaction rows have finite interaction SE.",
       call. = FALSE)
}

pooled_rows <- pool_input[
  ,
  partial_pool_group(.SD),
  by = .(gene_symbol, tx_id, design_family, branch_class, context_axis)
]
pooled_rows[, pooled_interaction_lower :=
              pooled_interaction_delta - 1.96 * pooled_interaction_se]
pooled_rows[, pooled_interaction_upper :=
              pooled_interaction_delta + 1.96 * pooled_interaction_se]
pooled_rows[, pooled_interaction_p := 2 * pnorm(
  -abs(pooled_interaction_delta / pooled_interaction_se)
)]
pooled_rows[, pooled_interaction_q := p.adjust(pooled_interaction_p, "BH")]
pooled_rows[, pooled_ci_excludes_zero :=
              pooled_interaction_lower > 0 | pooled_interaction_upper < 0]
pooled_rows[, pooled_abs_delta := abs(pooled_interaction_delta)]
pooled_rows[, shrinkage_abs_loss := abs(interaction_delta) -
              abs(pooled_interaction_delta)]
pooled_rows[, shrinkage_fraction := fifelse(
  abs(interaction_delta) > 0,
  pmax(0, shrinkage_abs_loss) / abs(interaction_delta),
  0
)]
pooled_rows[, pooled_same_direction_as_raw := sign(pooled_interaction_delta) ==
              sign(interaction_delta)]
pooled_rows[, partial_pooled_context_class := fcase(
  pooling_group_contexts < 2 &
    raw_context_supported == TRUE &
    pooled_ci_excludes_zero &
    pooled_abs_delta >= 0.06,
  "single_context_supported_not_pooled",
  pooling_group_contexts < 2 &
    (raw_context_supported == TRUE |
       (pooled_interaction_q <= 0.20 & pooled_abs_delta >= 0.04)),
  "single_context_review_not_pooled",
  pooling_group_contexts < 2,
  "single_context_neutral_not_pooled",
  pooled_ci_excludes_zero &
    pooled_interaction_q <= 0.05 &
    pooled_abs_delta >= 0.08 &
    context_opposite_complement_direction == TRUE,
  "partial_pooled_direction_reversal",
  pooled_ci_excludes_zero &
    pooled_interaction_q <= 0.05 &
    pooled_abs_delta >= 0.06,
  "partial_pooled_supported_context_shift",
  pooled_interaction_q <= 0.20 &
    pooled_abs_delta >= 0.04,
  "partial_pooled_review",
  raw_context_supported == TRUE &
    (pooled_interaction_q > 0.20 | shrinkage_fraction >= 0.50),
  "raw_signal_shrunk_fragile",
  default = "partial_pooled_neutral"
)]
pooled_rows[, partial_pooling_survival_class := fcase(
  raw_context_supported == TRUE &
    partial_pooled_context_class %chin% c(
      "single_context_supported_not_pooled",
      "single_context_review_not_pooled"
    ),
  "single_context_raw_signal_not_pooled",
  raw_context_supported == TRUE &
    partial_pooled_context_class %chin% c(
      "partial_pooled_direction_reversal",
      "partial_pooled_supported_context_shift",
      "partial_pooled_review"
    ),
  "raw_signal_survives_partial_pooling",
  raw_context_supported == TRUE,
  "raw_signal_shrunk_by_partial_pooling",
  raw_context_supported == FALSE &
    partial_pooled_context_class %chin% c(
      "partial_pooled_direction_reversal",
      "partial_pooled_supported_context_shift",
      "partial_pooled_review"
    ),
  "pooled_context_emerges",
  default = "no_context_signal"
)]
pooled_rows[, partial_pooling_review_priority := clip01(
  0.28 * clip01(pooled_abs_delta / 0.12) +
    0.18 * clip01(-log10(pooled_interaction_q + 1e-12) / 4) +
    0.14 * clip01(context_interaction_score) +
    0.14 * clip01(model_agreement_score) +
    0.10 * clip01(transport_stability_score) +
    0.08 * clip01((context_studies + complement_studies) / 8) +
    0.06 * as.numeric(pooled_ci_excludes_zero) +
    0.04 * as.numeric(raw_context_supported) -
    0.10 * clip01(shrinkage_fraction)
)]
setorder(pooled_rows, -partial_pooling_review_priority,
         gene_symbol, design_family, branch_class, context_axis)
pooled_rows[, partial_pooling_review_rank := seq_len(.N)]

gene_design_summary <- pooled_rows[
  ,
  {
    idx <- which.max(partial_pooling_review_priority)
    if (!length(idx)) idx <- 1L
    .(
      n_partial_pooling_rows = .N,
      n_pooled_supported_rows = sum(partial_pooled_context_class %chin% c(
        "partial_pooled_direction_reversal",
        "partial_pooled_supported_context_shift"
      )),
      n_pooled_review_rows = sum(partial_pooled_context_class ==
                                   "partial_pooled_review"),
      n_raw_supported_rows = sum(raw_context_supported),
      n_raw_survived_rows = sum(
        partial_pooling_survival_class ==
          "raw_signal_survives_partial_pooling"
      ),
      n_raw_shrunk_rows = sum(
        partial_pooling_survival_class ==
          "raw_signal_shrunk_by_partial_pooling"
      ),
      n_single_context_raw_signal_rows = sum(
        partial_pooling_survival_class ==
          "single_context_raw_signal_not_pooled"
      ),
      n_pooled_emergent_rows = sum(
        partial_pooling_survival_class == "pooled_context_emerges"
      ),
      top_partial_pooling_class = partial_pooled_context_class[idx],
      top_partial_pooling_survival_class =
        partial_pooling_survival_class[idx],
      top_context_axis = context_axis[idx],
      top_context_value = context_value[idx],
      top_branch_class = branch_class[idx],
      top_pooled_interaction_delta = pooled_interaction_delta[idx],
      top_pooled_interaction_lower = pooled_interaction_lower[idx],
      top_pooled_interaction_upper = pooled_interaction_upper[idx],
      top_pooled_interaction_q = pooled_interaction_q[idx],
      top_raw_interaction_delta = interaction_delta[idx],
      top_partial_pooling_weight = partial_pooling_weight[idx],
      top_shrinkage_fraction = shrinkage_fraction[idx],
      top_review_priority = partial_pooling_review_priority[idx],
      top_context_summary = collapse_unique(context_top_labels[idx], 1L),
      pooled_context_values = collapse_unique(context_value, 6L)
    )
  },
  by = .(gene_symbol, tx_id, design_family)
]
gene_design_summary[, partial_pooling_gene_class := fcase(
  n_pooled_supported_rows > 0,
  "partial_pooled_context_supported",
  n_single_context_raw_signal_rows > 0,
  "single_context_raw_signal_not_pooled",
  n_pooled_review_rows > 0,
  "partial_pooled_context_review",
  n_raw_supported_rows > 0 &
    n_raw_shrunk_rows == n_raw_supported_rows,
  "raw_context_signal_fragile_after_pooling",
  default = "partial_pooled_context_neutral"
)]
setorder(gene_design_summary, -top_review_priority,
         gene_symbol, design_family)
gene_design_summary[, partial_pooling_gene_rank := seq_len(.N)]

viral_review <- gene_design_summary[
  design_family == "viral_infection" &
    partial_pooling_gene_class != "partial_pooled_context_neutral"
]
setorder(viral_review, -top_review_priority)

comparison <- pooled_rows[
  raw_context_supported == TRUE |
    partial_pooled_context_class != "partial_pooled_neutral"
][
  ,
  .(
    n_rows = .N,
    median_raw_abs_delta = median(abs(interaction_delta), na.rm = TRUE),
    median_pooled_abs_delta = median(pooled_abs_delta, na.rm = TRUE),
    mean_shrinkage_fraction = safe_mean(shrinkage_fraction, 0),
    n_raw_supported = sum(raw_context_supported),
    n_raw_survived = sum(
      partial_pooling_survival_class ==
        "raw_signal_survives_partial_pooling"
    ),
    n_raw_shrunk = sum(
      partial_pooling_survival_class ==
        "raw_signal_shrunk_by_partial_pooling"
    ),
    n_single_context_raw_signal = sum(
      partial_pooling_survival_class ==
        "single_context_raw_signal_not_pooled"
    ),
    n_pooled_emergent = sum(
      partial_pooling_survival_class == "pooled_context_emerges"
    )
  ),
  by = .(design_family, context_axis, branch_class)
]
setorder(comparison, -n_raw_supported, -n_raw_survived)

summary_metrics <- data.table(
  metric = c(
    "n_evaluable_context_rows",
    "n_gene_designs",
    "n_pooled_supported_rows",
    "n_pooled_review_rows",
    "n_raw_supported_rows",
    "n_raw_survived_rows",
    "n_raw_shrunk_rows",
    "n_single_context_raw_signal_rows",
    "n_pooled_emergent_rows",
    "top_gene",
    "top_design_family",
    "top_context",
    "top_branch",
    "top_pooled_interaction_delta",
    "top_pooled_interaction_q"
  ),
  value = c(
    as.character(nrow(pooled_rows)),
    as.character(nrow(gene_design_summary)),
    as.character(sum(pooled_rows$partial_pooled_context_class %chin% c(
      "partial_pooled_direction_reversal",
      "partial_pooled_supported_context_shift"
    ))),
    as.character(sum(
      pooled_rows$partial_pooled_context_class == "partial_pooled_review"
    )),
    as.character(sum(pooled_rows$raw_context_supported)),
    as.character(sum(
      pooled_rows$partial_pooling_survival_class ==
        "raw_signal_survives_partial_pooling"
    )),
    as.character(sum(
      pooled_rows$partial_pooling_survival_class ==
        "raw_signal_shrunk_by_partial_pooling"
    )),
    as.character(sum(
      pooled_rows$partial_pooling_survival_class ==
        "single_context_raw_signal_not_pooled"
    )),
    as.character(sum(
      pooled_rows$partial_pooling_survival_class == "pooled_context_emerges"
    )),
    pooled_rows$gene_symbol[[1]],
    pooled_rows$design_family[[1]],
    paste0(pooled_rows$context_axis[[1]], "=",
           pooled_rows$context_value[[1]]),
    pooled_rows$branch_class[[1]],
    sprintf("%.4f", pooled_rows$pooled_interaction_delta[[1]]),
    sprintf("%.3g", pooled_rows$pooled_interaction_q[[1]])
  )
)

fwrite(
  pooled_rows,
  file.path(output_dir, "dominant_rdg_context_partial_pooling_rows.csv")
)
fwrite(
  gene_design_summary,
  file.path(output_dir,
            "dominant_rdg_context_partial_pooling_gene_design_summary.csv")
)
fwrite(
  viral_review,
  file.path(output_dir, "dominant_rdg_context_partial_pooling_viral_review.csv")
)
fwrite(
  comparison,
  file.path(output_dir,
            "dominant_rdg_context_partial_pooling_raw_vs_pooled.csv")
)
fwrite(
  summary_metrics,
  file.path(output_dir,
            "dominant_rdg_context_partial_pooling_summary_metrics.csv")
)

plot_top <- head(copy(viral_review), 24L)
if (nrow(plot_top)) {
  plot_top[, label := paste(
    gene_symbol, top_branch_class,
    paste0(top_context_axis, "=", top_context_value),
    sep = " | "
  )]
  plot_top[, label := factor(label, levels = rev(label))]
  plot_top[, display_class := fcase(
    partial_pooling_gene_class == "partial_pooled_context_supported",
    "pooled supported",
    partial_pooling_gene_class == "partial_pooled_context_review",
    "pooled review",
    partial_pooling_gene_class == "single_context_raw_signal_not_pooled",
    "single context only",
    partial_pooling_gene_class == "raw_context_signal_fragile_after_pooling",
    "fragile after pooling",
    default = partial_pooling_gene_class
  )]
  p_top <- ggplot(
    plot_top,
    aes(
      x = top_pooled_interaction_delta,
      y = label,
      color = display_class
    )
  ) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_errorbar(
      aes(
        xmin = top_pooled_interaction_lower,
        xmax = top_pooled_interaction_upper
      ),
      orientation = "y",
      width = 0.18,
      alpha = 0.7
    ) +
    geom_point(aes(size = top_review_priority), alpha = 0.9) +
    scale_size_continuous(limits = c(0, 1), range = c(2, 5.5)) +
    scale_color_manual(values = c(
      "pooled supported" = "#76A80D",
      "pooled review" = "#F8766D",
      "fragile after pooling" = "#00BFC4",
      "single context only" = "#C77CFF"
    )) +
    guides(
      size = "none",
      color = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    labs(
      title = "Partial-Pooled Viral Context Effects",
      subtitle = wrap_text(
        "Empirical-Bayes shrinkage across context values within each gene/design/branch/axis. Intervals show pooled interaction effects."
      ),
      x = "partial-pooled interaction delta",
      y = NULL,
      color = "evidence"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  ggsave(
    file.path(figure_dir,
              "context_partial_pooling_viral_top_effects.png"),
    p_top,
    width = 11,
    height = 7.5,
    dpi = 180
  )
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir,
                "context_partial_pooling_viral_top_effects.html"),
      title = "RDG partial-pooled viral context effects"
    )
  }
}

ifih1 <- pooled_rows[
  gene_symbol == "IFIH1" &
    design_family == "viral_infection" &
    context_axis == "tissue_context" &
    context_value == "lung" &
    branch_class %in% c("clean_CDS", "leader_uORF")
]
if (nrow(ifih1)) {
  ifih1_long <- melt(
    ifih1[
      ,
      .(
        branch_class,
        raw = interaction_delta,
        partial_pooled = pooled_interaction_delta
      )
    ],
    id.vars = "branch_class",
    variable.name = "estimate_type",
    value.name = "interaction_delta"
  )
  p_ifih1 <- ggplot(
    ifih1_long,
    aes(x = branch_class, y = interaction_delta, fill = estimate_type)
  ) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    labs(
      title = "IFIH1 Lung Viral Context: Raw vs Shrinkage-Aware Estimate",
      subtitle = wrap_text(
        "No sibling tissue contexts were available for shrinkage; the estimate equals raw and is labeled single-context evidence.",
        72L
      ),
      x = NULL,
      y = "interaction delta",
      fill = "estimate"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
  ggsave(
    file.path(figure_dir, "context_partial_pooling_ifih1_lung.png"),
    p_ifih1,
    width = 7.5,
    height = 4.8,
    dpi = 180
  )
}

message("Saved partial-pooled context outputs to: ", output_dir)
