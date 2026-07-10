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
output_dir <- file.path(analysis_dir, "dominant_rdg_hierarchical_joint_allocation")
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

weighted_mean_safe <- function(x, w = NULL) {
  x <- safe_num(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- safe_num(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
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
      p = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      q_stat = NA_real_,
      direction_fraction = NA_real_,
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
  direction_fraction <- if (sign_ref == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign_ref)
  }
  i2 <- if (q_stat > 0 && df > 0) max(0, (q_stat - df) / q_stat) else 0
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    p = p,
    tau2 = tau2,
    i2 = i2,
    q_stat = q_stat,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

message("Hierarchical joint RDG allocation model:")
message("  1. Use row-level joint Dirichlet branch deltas")
message("  2. Collapse contexts within study")
message("  3. Estimate study-aware random-effects branch deltas per gene/design")

class_file <- file.path(
  analysis_dir,
  "dominant_rdg_joint_branch_allocation",
  "dominant_rdg_joint_branch_allocation_class_rows.csv"
)
dt <- read_dt(class_file, required = TRUE)

required <- c(
  "gene_symbol", "tx_id", "design_family", "study", "case_condition",
  "control_conditions", "CELL_LINE", "TISSUE", "match_scope",
  "sample_set_signature", "branch_class", "joint_delta", "joint_delta_se",
  "joint_case_count", "joint_control_count", "joint_evaluable",
  "joint_strict", "joint_l1_allocation_shift", "joint_review_class",
  "joint_review_priority"
)
missing <- setdiff(required, names(dt))
if (length(missing)) {
  stop("Joint class rows are missing columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

for (column in intersect(
  c("joint_delta", "joint_delta_se", "joint_delta_p", "joint_delta_q",
    "joint_case_count", "joint_control_count", "joint_l1_allocation_shift",
    "joint_js_divergence", "joint_chisq_p", "joint_chisq_q",
    "joint_review_priority"),
  names(dt)
)) {
  dt[, (column) := safe_num(get(column))]
}
for (column in intersect(c("joint_evaluable", "joint_strict"), names(dt))) {
  if (!is.logical(dt[[column]])) {
    dt[, (column) := tolower(as.character(get(column))) %chin%
         c("true", "t", "1", "yes")]
  }
  dt[is.na(get(column)), (column) := FALSE]
}

branch_levels <- c("leader_uORF", "overlapping_uORF", "clean_CDS", "other")
dt <- dt[
  branch_class %chin% branch_levels &
    joint_evaluable == TRUE &
    is.finite(joint_delta) &
    is.finite(joint_delta_se) &
    joint_delta_se > 0
]

if (!nrow(dt)) {
  empty <- data.table(metric = "hierarchical_joint_branch_rows", value = 0)
  fwrite(empty, file.path(
    output_dir,
    "dominant_rdg_hierarchical_joint_allocation_summary_metrics.csv"
  ))
  quit(save = "no", status = 0)
}

dt[, branch_class := factor(branch_class, levels = branch_levels)]
dt[, context_label := paste(
  study, CELL_LINE, TISSUE, case_condition,
  paste("vs", control_conditions),
  sep = " | "
)]
dt[, context_weight := pmax(
  0.25,
  fifelse(match_scope == "exact_context", 1, 0.55) *
    pmin(log1p(pmax(joint_case_count, joint_control_count, na.rm = TRUE)) /
           log1p(300), 1.5)
)]
dt[!is.finite(context_weight), context_weight := 0.25]

study_effects <- dt[, {
  re <- random_effects_summary(joint_delta, joint_delta_se, context_weight)
  top_i <- which.max(joint_review_priority)
  if (!length(top_i) || is.na(top_i)) top_i <- 1L
  data.table(
    study_context_rows = .N,
    study_strict_context_rows = sum(joint_strict == TRUE, na.rm = TRUE),
    study_delta = re$effect,
    study_delta_se = re$se,
    study_delta_lower = re$lower,
    study_delta_upper = re$upper,
    study_delta_p = re$p,
    study_tau2 = re$tau2,
    study_i2 = re$i2,
    study_direction_fraction = re$direction_fraction,
    study_ci_excludes_zero = re$ci_excludes_zero,
    study_max_abs_context_delta = max_safe(abs(joint_delta), 0),
    study_max_context_priority = max_safe(joint_review_priority, 0),
    study_max_l1_shift = max_safe(joint_l1_allocation_shift, 0),
    study_top_context = context_label[top_i],
    study_top_case_condition = case_condition[top_i],
    study_top_control_conditions = control_conditions[top_i],
    study_top_cell_line = CELL_LINE[top_i],
    study_top_tissue = TISSUE[top_i],
    study_top_match_scope = match_scope[top_i]
  )
}, by = .(gene_symbol, tx_id, design_family, branch_class, study)]
study_effects <- study_effects[is.finite(study_delta) &
                                 is.finite(study_delta_se) &
                                 study_delta_se > 0]
study_effects[, hierarchical_study_delta_q := p.adjust(study_delta_p,
                                                       method = "BH")]
study_effects[, study_abs_delta := abs(study_delta)]
setorder(study_effects, gene_symbol, tx_id, design_family, branch_class,
         -study_max_context_priority, -study_abs_delta)
study_effects[, study_abs_delta := NULL]
fwrite(study_effects, file.path(
  output_dir,
  "dominant_rdg_hierarchical_joint_allocation_study_effects.csv"
))

branch_effects <- study_effects[, {
  re <- random_effects_summary(study_delta, study_delta_se)
  top_i <- which.max(study_max_context_priority)
  if (!length(top_i) || is.na(top_i)) top_i <- 1L
  data.table(
    hierarchical_joint_studies = .N,
    hierarchical_joint_context_rows = sum(study_context_rows, na.rm = TRUE),
    hierarchical_joint_strict_context_rows =
      sum(study_strict_context_rows, na.rm = TRUE),
    hierarchical_joint_delta = re$effect,
    hierarchical_joint_delta_se = re$se,
    hierarchical_joint_delta_lower = re$lower,
    hierarchical_joint_delta_upper = re$upper,
    hierarchical_joint_delta_p = re$p,
    hierarchical_joint_tau2 = re$tau2,
    hierarchical_joint_i2 = re$i2,
    hierarchical_joint_q_stat = re$q_stat,
    hierarchical_joint_direction_fraction = re$direction_fraction,
    hierarchical_joint_ci_excludes_zero = re$ci_excludes_zero,
    hierarchical_joint_weighted_study_delta =
      weighted_mean_safe(study_delta, 1 / pmax(study_delta_se^2, 1e-8)),
    hierarchical_joint_median_study_delta = median_safe(study_delta),
    hierarchical_joint_max_abs_study_delta = max_safe(abs(study_delta), 0),
    hierarchical_joint_max_context_priority =
      max_safe(study_max_context_priority, 0),
    hierarchical_joint_max_l1_shift = max_safe(study_max_l1_shift, 0),
    hierarchical_joint_top_study = study[top_i],
    hierarchical_joint_top_context = study_top_context[top_i],
    hierarchical_joint_top_case_condition = study_top_case_condition[top_i],
    hierarchical_joint_top_control_conditions =
      study_top_control_conditions[top_i],
    hierarchical_joint_top_cell_line = study_top_cell_line[top_i],
    hierarchical_joint_top_tissue = study_top_tissue[top_i],
    hierarchical_joint_top_match_scope = study_top_match_scope[top_i],
    hierarchical_joint_study_list = collapse_unique(study, 8L)
  )
}, by = .(gene_symbol, tx_id, design_family, branch_class)]

branch_effects[, hierarchical_joint_delta_q :=
                 p.adjust(hierarchical_joint_delta_p, method = "BH")]
branch_effects[, hierarchical_joint_support_class := fcase(
  hierarchical_joint_studies >= 2 &
    hierarchical_joint_ci_excludes_zero == TRUE &
    abs(hierarchical_joint_delta) >= 0.04 &
    hierarchical_joint_direction_fraction >= 0.70,
  "hierarchical_joint_replicated_shift",
  hierarchical_joint_studies >= 2 &
    is.finite(hierarchical_joint_delta_q) &
    hierarchical_joint_delta_q <= 0.20 &
    abs(hierarchical_joint_delta) >= 0.03 &
    hierarchical_joint_direction_fraction >= 0.65,
  "hierarchical_joint_replicated_review",
  hierarchical_joint_studies == 1 &
    abs(hierarchical_joint_delta) >= 0.08 &
    is.finite(hierarchical_joint_delta_p) &
    hierarchical_joint_delta_p <= 0.05,
  "hierarchical_joint_single_study_strong",
  hierarchical_joint_context_rows < 2,
  "hierarchical_joint_low_context_support",
  default = "hierarchical_joint_neutral"
)]
branch_effects[, hierarchical_joint_priority := clip01(
  0.24 * pmin(abs(hierarchical_joint_delta) / 0.15, 1) +
    0.18 * as.numeric(hierarchical_joint_ci_excludes_zero == TRUE) +
    0.16 * pmin(hierarchical_joint_studies / 3, 1) +
    0.13 * pmin(hierarchical_joint_direction_fraction, 1) +
    0.10 * (1 - pmin(hierarchical_joint_i2, 1)) +
    0.09 * q_component(hierarchical_joint_delta_q) +
    0.06 * pmin(hierarchical_joint_max_context_priority, 1) +
    0.04 * as.numeric(design_family == "viral_infection")
)]
branch_effects[hierarchical_joint_support_class ==
                 "hierarchical_joint_low_context_support",
               hierarchical_joint_priority :=
                 pmin(hierarchical_joint_priority, 0.25)]
branch_effects[, hierarchical_joint_abs_delta := abs(hierarchical_joint_delta)]
setorder(branch_effects, -hierarchical_joint_priority,
         -hierarchical_joint_abs_delta, gene_symbol, design_family,
         branch_class)
branch_effects[, hierarchical_joint_abs_delta := NULL]

branch_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_joint_allocation_branch_effects.csv"
)
fwrite(branch_effects, branch_file)

top_branch <- branch_effects[, .SD[1],
                             by = .(gene_symbol, tx_id, design_family)]
gene_design <- top_branch[, .(
  gene_symbol, tx_id, design_family,
  hierarchical_joint_top_branch_class = as.character(branch_class),
  hierarchical_joint_top_delta = hierarchical_joint_delta,
  hierarchical_joint_top_delta_lower = hierarchical_joint_delta_lower,
  hierarchical_joint_top_delta_upper = hierarchical_joint_delta_upper,
  hierarchical_joint_top_delta_q = hierarchical_joint_delta_q,
  hierarchical_joint_top_support_class = hierarchical_joint_support_class,
  hierarchical_joint_top_priority = hierarchical_joint_priority,
  hierarchical_joint_top_studies = hierarchical_joint_studies,
  hierarchical_joint_top_context_rows = hierarchical_joint_context_rows,
  hierarchical_joint_top_direction_fraction =
    hierarchical_joint_direction_fraction,
  hierarchical_joint_top_i2 = hierarchical_joint_i2,
  hierarchical_joint_top_context = hierarchical_joint_top_context,
  hierarchical_joint_top_study = hierarchical_joint_top_study,
  hierarchical_joint_top_match_scope = hierarchical_joint_top_match_scope,
  hierarchical_joint_study_list = hierarchical_joint_study_list
)]

wide <- dcast(
  branch_effects,
  gene_symbol + tx_id + design_family ~ branch_class,
  value.var = c(
    "hierarchical_joint_delta",
    "hierarchical_joint_delta_lower",
    "hierarchical_joint_delta_upper",
    "hierarchical_joint_delta_q",
    "hierarchical_joint_priority",
    "hierarchical_joint_support_class"
  ),
  fill = NA
)
gene_design <- merge(gene_design, wide,
                     by = c("gene_symbol", "tx_id", "design_family"),
                     all.x = TRUE, sort = FALSE)
delta_cols <- intersect(
  paste0("hierarchical_joint_delta_", branch_levels),
  names(gene_design)
)
gene_design[, hierarchical_joint_l1_effect := {
  m <- as.matrix(.SD)
  m[!is.finite(m)] <- 0
  rowSums(abs(m)) / 2
}, .SDcols = delta_cols]
support_cols <- intersect(
  paste0("hierarchical_joint_support_class_", branch_levels),
  names(gene_design)
)
supported_classes <- c(
  "hierarchical_joint_replicated_shift",
  "hierarchical_joint_replicated_review",
  "hierarchical_joint_single_study_strong"
)
if (length(support_cols)) {
  support_matrix <- sapply(
    support_cols,
    function(col) gene_design[[col]] %chin% supported_classes
  )
  if (is.null(dim(support_matrix))) {
    support_matrix <- matrix(support_matrix, ncol = 1)
  }
  gene_design[, hierarchical_joint_supported_branch_count :=
                rowSums(support_matrix, na.rm = TRUE)]
} else {
  gene_design[, hierarchical_joint_supported_branch_count := 0L]
}
setorder(gene_design, -hierarchical_joint_top_priority,
         -hierarchical_joint_l1_effect, gene_symbol, design_family)

gene_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv"
)
fwrite(gene_design, gene_file)

review <- gene_design[
  hierarchical_joint_top_support_class != "hierarchical_joint_neutral" |
    hierarchical_joint_top_priority >= 0.55
]
setorder(review, -hierarchical_joint_top_priority,
         -hierarchical_joint_l1_effect, gene_symbol, design_family)
review_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_joint_allocation_review.csv"
)
fwrite(review, review_file)

viral_review <- gene_design[design_family == "viral_infection"]
setorder(viral_review, -hierarchical_joint_top_priority,
         -hierarchical_joint_l1_effect, gene_symbol)
viral_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_joint_allocation_viral_review.csv"
)
fwrite(viral_review, viral_file)

summary_metrics <- rbindlist(list(
  branch_effects[, .(
    branch_effect_rows = .N,
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    genes = uniqueN(gene_symbol),
    replicated_shift_branches =
      sum(hierarchical_joint_support_class ==
            "hierarchical_joint_replicated_shift", na.rm = TRUE),
    replicated_review_branches =
      sum(hierarchical_joint_support_class ==
            "hierarchical_joint_replicated_review", na.rm = TRUE),
    single_study_strong_branches =
      sum(hierarchical_joint_support_class ==
            "hierarchical_joint_single_study_strong", na.rm = TRUE),
    median_abs_delta = median_safe(abs(hierarchical_joint_delta)),
    max_abs_delta = max_safe(abs(hierarchical_joint_delta), 0),
    median_i2 = median_safe(hierarchical_joint_i2)
  )][, scope := "all"],
  branch_effects[, .(
    branch_effect_rows = .N,
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    genes = uniqueN(gene_symbol),
    replicated_shift_branches =
      sum(hierarchical_joint_support_class ==
            "hierarchical_joint_replicated_shift", na.rm = TRUE),
    replicated_review_branches =
      sum(hierarchical_joint_support_class ==
            "hierarchical_joint_replicated_review", na.rm = TRUE),
    single_study_strong_branches =
      sum(hierarchical_joint_support_class ==
            "hierarchical_joint_single_study_strong", na.rm = TRUE),
    median_abs_delta = median_safe(abs(hierarchical_joint_delta)),
    max_abs_delta = max_safe(abs(hierarchical_joint_delta), 0),
    median_i2 = median_safe(hierarchical_joint_i2)
  ), by = .(scope = paste0("design_family:", design_family))]
), fill = TRUE)
setcolorder(summary_metrics, c("scope", setdiff(names(summary_metrics), "scope")))
summary_file <- file.path(
  output_dir,
  "dominant_rdg_hierarchical_joint_allocation_summary_metrics.csv"
)
fwrite(summary_metrics, summary_file)

plot_dt <- head(review[order(-hierarchical_joint_top_priority)], 45)
if (nrow(plot_dt)) {
  plot_dt[, plot_label := paste0(gene_symbol, " / ", design_family)]
  plot_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_top <- ggplot(
    plot_dt,
    aes(
      x = hierarchical_joint_top_delta,
      y = plot_label,
      fill = hierarchical_joint_top_branch_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Design: ", design_family,
        "<br>Top branch: ", hierarchical_joint_top_branch_class,
        "<br>RE delta: ", round(hierarchical_joint_top_delta, 3),
        "<br>95% interval: [",
        round(hierarchical_joint_top_delta_lower, 3), ", ",
        round(hierarchical_joint_top_delta_upper, 3), "]",
        "<br>q: ", signif(hierarchical_joint_top_delta_q, 3),
        "<br>Support: ", hierarchical_joint_top_support_class,
        "<br>Studies: ", hierarchical_joint_top_studies,
        "<br>Top context: ", hierarchical_joint_top_context
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
      title = "Hierarchical Joint RDG Allocation Effects",
      subtitle = "Study-aware random-effects estimates of branch allocation delta",
      x = "Random-effects branch allocation delta",
      y = NULL,
      fill = "Top branch"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_joint_allocation_top_effects.png"),
         p_top, width = 10.8, height = 8.4, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_joint_allocation_top_effects.pdf"),
         p_top, width = 10.8, height = 8.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir, "hierarchical_joint_allocation_top_effects.html"),
      title = "Hierarchical joint allocation effects"
    )
  }
}

viral_plot <- head(viral_review[order(-hierarchical_joint_top_priority)], 35)
if (nrow(viral_plot)) {
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p_viral <- ggplot(
    viral_plot,
    aes(
      x = hierarchical_joint_top_priority,
      y = gene_symbol,
      fill = hierarchical_joint_top_branch_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Top branch: ", hierarchical_joint_top_branch_class,
        "<br>RE delta: ", round(hierarchical_joint_top_delta, 3),
        "<br>Support: ", hierarchical_joint_top_support_class,
        "<br>Studies: ", hierarchical_joint_top_studies,
        "<br>Top context: ", hierarchical_joint_top_context
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
    labs(
      title = "Viral Hierarchical Joint RDG Allocation Review",
      subtitle = "Study-aware branch redistribution priority",
      x = "Hierarchical joint allocation priority",
      y = NULL,
      fill = "Top branch"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_joint_allocation_viral_review.png"),
         p_viral, width = 9.8, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_joint_allocation_viral_review.pdf"),
         p_viral, width = 9.8, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_viral,
      file.path(figure_dir, "hierarchical_joint_allocation_viral_review.html"),
      title = "Viral hierarchical joint allocation review"
    )
  }
}

heat_dt <- head(gene_design[order(-hierarchical_joint_top_priority)], 40)
if (nrow(heat_dt)) {
  heat_dt[, plot_label := paste0(gene_symbol, " / ", design_family)]
  heat_long <- melt(
    heat_dt,
    id.vars = c("plot_label", "gene_symbol", "design_family",
                "hierarchical_joint_top_priority",
                "hierarchical_joint_top_support_class"),
    measure.vars = intersect(
      paste0("hierarchical_joint_delta_", branch_levels),
      names(heat_dt)
    ),
    variable.name = "branch",
    value.name = "hierarchical_delta"
  )
  heat_long[, branch := sub("^hierarchical_joint_delta_", "", branch)]
  heat_long[, branch := factor(branch, levels = branch_levels)]
  heat_long[, plot_label := factor(plot_label, levels = rev(unique(heat_dt$plot_label)))]
  p_heat <- ggplot(
    heat_long,
    aes(
      x = branch,
      y = plot_label,
      fill = hierarchical_delta,
      text = paste0(
        "Gene/design: ", plot_label,
        "<br>Branch: ", branch,
        "<br>RE delta: ", round(hierarchical_delta, 3),
        "<br>Support: ", hierarchical_joint_top_support_class,
        "<br>Priority: ", round(hierarchical_joint_top_priority, 3)
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
      title = "Hierarchical Joint Allocation Delta Matrix",
      subtitle = "Study-aware random-effects branch deltas",
      x = NULL,
      y = NULL,
      fill = "Delta"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(figure_dir, "hierarchical_joint_allocation_delta_heatmap.png"),
         p_heat, width = 9.6, height = 8.8, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_joint_allocation_delta_heatmap.pdf"),
         p_heat, width = 9.6, height = 8.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_heat,
      file.path(figure_dir, "hierarchical_joint_allocation_delta_heatmap.html"),
      title = "Hierarchical joint allocation delta matrix"
    )
  }
}

message("Saved hierarchical joint allocation outputs in: ", output_dir)
print(summary_metrics[scope %chin% c("all", "design_family:viral_infection")])
