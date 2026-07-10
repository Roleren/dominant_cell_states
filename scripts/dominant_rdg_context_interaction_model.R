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
output_dir <- file.path(analysis_dir, "dominant_rdg_context_interactions")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

read_dt <- function(path, required = FALSE, ...) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing required input: ", path, call. = FALSE)
    return(data.table())
  }
  fread(path, showProgress = FALSE, ...)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

finite_or_zero <- function(x) {
  x <- safe_num(x)
  x[!is.finite(x)] <- 0
  x
}

clip01 <- function(x) pmin(pmax(finite_or_zero(x), 0), 1)

min_safe <- function(x, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else min(x)
}

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
  x <- as.character(x[!is.na(x) & nzchar(as.character(x))])
  x <- trimws(unlist(strsplit(x, ";\\s*", perl = TRUE), use.names = FALSE))
  x <- sort(unique(x[nzchar(x)]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

normalize_level <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | tolower(x) %chin%
      c("na", "n/a", "unknown", "missing", "none")] <- NA_character_
  x
}

bool_col <- function(x) {
  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  }
  out <- tolower(as.character(x)) %chin% c("true", "t", "1", "yes")
  out[is.na(out)] <- FALSE
  out
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

support_call <- function(n, effect, p, ci, direction_fraction) {
  if (!is.finite(n) || n < 2 || !is.finite(effect) ||
      !is.finite(p) || !is.finite(direction_fraction)) {
    return("not_evaluable")
  }
  fcase(
    ci & abs(effect) >= 0.04 & direction_fraction >= 0.70,
    "replicated_shift",
    p <= 0.20 & abs(effect) >= 0.03 & direction_fraction >= 0.65,
    "replicated_review",
    default = "neutral"
  )
}

context_class <- function(evaluable, context_support, complement_support,
                          context_effect, complement_effect,
                          interaction_delta, interaction_p) {
  if (!isTRUE(evaluable)) return("not_evaluable")
  same_direction <- is.finite(context_effect) & is.finite(complement_effect) &
    sign(context_effect) == sign(complement_effect) & sign(context_effect) != 0
  opposite_direction <- is.finite(context_effect) & is.finite(complement_effect) &
    sign(context_effect) != sign(complement_effect) & sign(context_effect) != 0 &
    sign(complement_effect) != 0
  interaction_supported <- is.finite(interaction_delta) &
    abs(interaction_delta) >= 0.045 &
    (is.finite(interaction_p) & interaction_p <= 0.25)
  context_has_support <- context_support %chin%
    c("replicated_shift", "replicated_review")
  complement_has_support <- complement_support %chin%
    c("replicated_shift", "replicated_review")
  if (opposite_direction && abs(interaction_delta) >= 0.04) {
    return("context_direction_reversal")
  }
  if (context_has_support && !complement_has_support &&
      abs(interaction_delta) >= 0.035) {
    return("context_specific_supported_shift")
  }
  if (interaction_supported) return("context_magnitude_modulated")
  if (context_has_support && complement_has_support && same_direction) {
    return("broad_context_consistent")
  }
  "weak_context_difference"
}

message("RDG context-interaction model:")
message("  1. Use joint branch-allocation class rows")
message("  2. Collapse rows within study and tissue/cell-line context")
message("  3. Compare context-specific viral branch effects to complements")

class_file <- file.path(
  analysis_dir,
  "dominant_rdg_joint_branch_allocation",
  "dominant_rdg_joint_branch_allocation_class_rows.csv"
)
dt <- read_dt(class_file, required = TRUE)

agreement <- read_dt(file.path(
  analysis_dir,
  "dominant_rdg_model_agreement_transport",
  "dominant_rdg_model_agreement_transport_gene_design.csv"
))

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
  dt[, (column) := bool_col(get(column))]
}

branch_levels <- c("leader_uORF", "overlapping_uORF", "clean_CDS", "other")
dt <- dt[
  design_family == "viral_infection" &
    branch_class %chin% branch_levels &
    joint_evaluable == TRUE &
    is.finite(joint_delta) &
    is.finite(joint_delta_se) &
    joint_delta_se > 0
]
if (!nrow(dt)) {
  empty <- data.table(metric = "context_interaction_rows", value = 0)
  fwrite(empty, file.path(
    output_dir,
    "dominant_rdg_context_interaction_summary_metrics.csv"
  ))
  quit(save = "no", status = 0)
}

dt[, branch_class := factor(branch_class, levels = branch_levels)]
dt[, TISSUE := normalize_level(TISSUE)]
dt[, CELL_LINE := normalize_level(CELL_LINE)]
dt[, context_row_weight := pmax(
  0.25,
  fifelse(match_scope == "exact_context", 1, 0.55) *
    pmin(log1p(pmax(joint_case_count, joint_control_count, na.rm = TRUE)) /
           log1p(300), 1.5)
)]
dt[!is.finite(context_row_weight), context_row_weight := 0.25]
dt[, context_label := paste(
  study,
  fifelse(is.na(CELL_LINE), "unknown_cell", CELL_LINE),
  fifelse(is.na(TISSUE), "unknown_tissue", TISSUE),
  case_condition,
  paste("vs", control_conditions),
  sep = " | "
)]

context_long <- rbindlist(list(
  dt[!is.na(TISSUE)][, .SD][, `:=`(
    context_axis = "tissue_context",
    context_value = TISSUE
  )],
  dt[!is.na(CELL_LINE)][, .SD][, `:=`(
    context_axis = "cell_line_context",
    context_value = CELL_LINE
  )]
), use.names = TRUE, fill = TRUE)

context_long <- context_long[!is.na(context_value) & nzchar(context_value)]

study_context_effects <- context_long[, {
  re <- random_effects_summary(joint_delta, joint_delta_se, context_row_weight)
  top_i <- which.max(joint_review_priority)
  if (!length(top_i) || is.na(top_i)) top_i <- 1L
  data.table(
    study_context_rows = .N,
    study_context_strict_rows = sum(joint_strict == TRUE, na.rm = TRUE),
    study_context_delta = re$effect,
    study_context_delta_se = re$se,
    study_context_delta_lower = re$lower,
    study_context_delta_upper = re$upper,
    study_context_delta_p = re$p,
    study_context_i2 = re$i2,
    study_context_direction_fraction = re$direction_fraction,
    study_context_ci_excludes_zero = re$ci_excludes_zero,
    study_context_max_l1_shift = max_safe(joint_l1_allocation_shift, 0),
    study_context_max_priority = max_safe(joint_review_priority, 0),
    study_context_case_count = sum(finite_or_zero(joint_case_count)),
    study_context_control_count = sum(finite_or_zero(joint_control_count)),
    study_context_top_label = context_label[top_i],
    study_context_top_case_condition = case_condition[top_i],
    study_context_top_control_conditions = control_conditions[top_i],
    study_context_top_match_scope = match_scope[top_i]
  )
}, by = .(gene_symbol, tx_id, design_family, branch_class,
          context_axis, context_value, study)]
study_context_effects <- study_context_effects[
  is.finite(study_context_delta) &
    is.finite(study_context_delta_se) &
    study_context_delta_se > 0
]

branch_key <- c("gene_symbol", "tx_id", "design_family", "branch_class")
axis_key <- c(branch_key, "context_axis")

context_effects <- study_context_effects[, {
  axis_dt <- .SD
  all_re <- random_effects_summary(
    axis_dt$study_context_delta,
    axis_dt$study_context_delta_se,
    pmax(axis_dt$study_context_case_count +
           axis_dt$study_context_control_count, 1)
  )
  rbindlist(lapply(sort(unique(axis_dt$context_value)), function(current_value) {
    context_dt <- axis_dt[context_value == current_value]
    complement_dt <- axis_dt[context_value != current_value]
    context_re <- random_effects_summary(
      context_dt$study_context_delta,
      context_dt$study_context_delta_se,
      pmax(context_dt$study_context_case_count +
             context_dt$study_context_control_count, 1)
    )
    complement_re <- random_effects_summary(
      complement_dt$study_context_delta,
      complement_dt$study_context_delta_se,
      pmax(complement_dt$study_context_case_count +
             complement_dt$study_context_control_count, 1)
    )
    interaction_delta <- context_re$effect - complement_re$effect
    interaction_se <- sqrt(context_re$se^2 + complement_re$se^2)
    interaction_z <- interaction_delta / interaction_se
    interaction_p <- 2 * pnorm(-abs(interaction_z))
    context_support <- support_call(
      context_re$n, context_re$effect, context_re$p,
      context_re$ci_excludes_zero, context_re$direction_fraction
    )
    complement_support <- support_call(
      complement_re$n, complement_re$effect, complement_re$p,
      complement_re$ci_excludes_zero, complement_re$direction_fraction
    )
    evaluable <- context_re$n >= 2 && complement_re$n >= 2 &&
      is.finite(interaction_delta) && is.finite(interaction_se) &&
      interaction_se > 0
    class <- context_class(
      evaluable, context_support, complement_support,
      context_re$effect, complement_re$effect,
      interaction_delta, interaction_p
    )
    opposite <- is.finite(context_re$effect) & is.finite(complement_re$effect) &
      sign(context_re$effect) != sign(complement_re$effect) &
      sign(context_re$effect) != 0 & sign(complement_re$effect) != 0
    data.table(
      context_value = current_value,
      context_studies = context_re$n,
      complement_studies = complement_re$n,
      total_axis_studies = uniqueN(axis_dt$study),
      context_rows = nrow(context_dt),
      complement_rows = nrow(complement_dt),
      context_delta = context_re$effect,
      context_delta_se = context_re$se,
      context_delta_lower = context_re$lower,
      context_delta_upper = context_re$upper,
      context_delta_p = context_re$p,
      context_direction_fraction = context_re$direction_fraction,
      context_i2 = context_re$i2,
      context_ci_excludes_zero = context_re$ci_excludes_zero,
      context_support_class = context_support,
      complement_delta = complement_re$effect,
      complement_delta_se = complement_re$se,
      complement_delta_lower = complement_re$lower,
      complement_delta_upper = complement_re$upper,
      complement_delta_p = complement_re$p,
      complement_direction_fraction = complement_re$direction_fraction,
      complement_i2 = complement_re$i2,
      complement_ci_excludes_zero = complement_re$ci_excludes_zero,
      complement_support_class = complement_support,
      global_axis_delta = all_re$effect,
      global_axis_delta_se = all_re$se,
      global_axis_delta_p = all_re$p,
      interaction_delta = interaction_delta,
      interaction_delta_se = interaction_se,
      interaction_delta_lower = interaction_delta - 1.96 * interaction_se,
      interaction_delta_upper = interaction_delta + 1.96 * interaction_se,
      interaction_delta_p = interaction_p,
      context_opposite_complement_direction = opposite,
      context_interaction_evaluable = evaluable,
      context_interaction_class = class,
      context_study_names = collapse_unique(context_dt$study, 10L),
      complement_study_names = collapse_unique(complement_dt$study, 10L),
      context_top_labels =
        collapse_unique(context_dt$study_context_top_label, 5L),
      complement_top_labels =
        collapse_unique(complement_dt$study_context_top_label, 5L),
      context_case_counts = sum(context_dt$study_context_case_count),
      context_control_counts = sum(context_dt$study_context_control_count),
      complement_case_counts = sum(complement_dt$study_context_case_count),
      complement_control_counts =
        sum(complement_dt$study_context_control_count)
    )
  }), use.names = TRUE)
}, by = axis_key]

context_effects[, interaction_delta_q := p.adjust(interaction_delta_p,
                                                  method = "BH")]
context_effects[, context_interaction_score := clip01(
  0.25 * as.numeric(context_interaction_evaluable) +
    0.18 * pmin(context_studies / 4, 1) +
    0.12 * pmin(complement_studies / 4, 1) +
    0.16 * q_component(interaction_delta_p) +
    0.12 * pmin(abs(interaction_delta) / 0.12, 1) +
    0.10 * as.numeric(context_support_class %chin%
                        c("replicated_shift", "replicated_review")) +
    0.07 * as.numeric(context_opposite_complement_direction)
)]
context_effects[
  context_interaction_class == "not_evaluable",
  context_interaction_score := 0
]
setorder(context_effects, -context_interaction_score,
         gene_symbol, branch_class, context_axis, context_value)

if (nrow(agreement)) {
  keep <- intersect(names(agreement), c(
    "gene_symbol", "tx_id", "design_family",
    "model_agreement_tier", "model_agreement_score",
    "transport_consensus_class", "transport_stability_score",
    "transport_agreement_class", "transport_worst_omitted_groups"
  ))
  context_effects <- merge(
    context_effects,
    agreement[, ..keep],
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE,
    sort = FALSE
  )
}
for (column in c("model_agreement_tier", "transport_consensus_class",
                 "transport_agreement_class",
                 "transport_worst_omitted_groups")) {
  if (!column %in% names(context_effects)) context_effects[, (column) := ""]
  context_effects[is.na(get(column)), (column) := ""]
}
for (column in c("model_agreement_score", "transport_stability_score")) {
  if (!column %in% names(context_effects)) context_effects[, (column) := 0]
  context_effects[!is.finite(get(column)), (column) := 0]
}
context_effects[, context_review_priority := pmin(
  1.5,
  context_interaction_score +
    0.20 * as.numeric(model_agreement_tier %chin% c(
      "A_browser_validated_consensus",
      "B_replicated_branch_consensus",
      "C_multi_model_consensus"
    )) +
    0.15 * as.numeric(transport_consensus_class %chin% c(
      "consensus_transport_fragile",
      "consensus_transport_direction_robust"
    )) +
    0.10 * pmin(model_agreement_score, 1) +
    0.08 * pmin(transport_stability_score, 1)
)]
setorder(context_effects, -context_review_priority,
         -context_interaction_score)

review_classes <- c(
  "context_direction_reversal",
  "context_specific_supported_shift",
  "context_magnitude_modulated"
)
context_review <- context_effects[
  context_interaction_class %chin% review_classes |
    context_review_priority >= 0.55
]

gene_design_summary <- context_effects[, {
  review_i <- which(context_interaction_class %chin% review_classes)
  top_i <- if (length(review_i)) {
    review_i[which.max(context_review_priority[review_i])]
  } else {
    which.max(context_review_priority)
  }
  if (!length(top_i) || is.na(top_i)) top_i <- 1L
  data.table(
    context_interaction_rows = .N,
    context_interaction_evaluable_rows =
      sum(context_interaction_evaluable == TRUE, na.rm = TRUE),
    context_direction_reversal_rows =
      sum(context_interaction_class == "context_direction_reversal"),
    context_specific_supported_rows =
      sum(context_interaction_class == "context_specific_supported_shift"),
    context_magnitude_modulated_rows =
      sum(context_interaction_class == "context_magnitude_modulated"),
    context_interaction_top_class = context_interaction_class[top_i],
    context_interaction_top_score = context_interaction_score[top_i],
    context_interaction_review_priority = context_review_priority[top_i],
    context_interaction_top_axis = context_axis[top_i],
    context_interaction_top_value = context_value[top_i],
    context_interaction_top_branch_class = as.character(branch_class[top_i]),
    context_interaction_context_delta = context_delta[top_i],
    context_interaction_context_delta_lower = context_delta_lower[top_i],
    context_interaction_context_delta_upper = context_delta_upper[top_i],
    context_interaction_complement_delta = complement_delta[top_i],
    context_interaction_complement_delta_lower = complement_delta_lower[top_i],
    context_interaction_complement_delta_upper = complement_delta_upper[top_i],
    context_interaction_delta = interaction_delta[top_i],
    context_interaction_delta_q = interaction_delta_q[top_i],
    context_interaction_context_studies = context_studies[top_i],
    context_interaction_complement_studies = complement_studies[top_i],
    context_interaction_context_study_names = context_study_names[top_i],
    context_interaction_complement_study_names = complement_study_names[top_i],
    context_interaction_summary = paste0(
      context_axis[top_i], "=", context_value[top_i],
      " / ", as.character(branch_class[top_i]),
      " context delta ", round(context_delta[top_i], 3),
      " vs complement ", round(complement_delta[top_i], 3),
      " (interaction ", round(interaction_delta[top_i], 3), ")"
    )
  )
}, by = .(gene_symbol, tx_id, design_family)]
gene_design_summary[, context_interaction_gene_class := fcase(
  context_direction_reversal_rows > 0, "context_direction_reversal",
  context_specific_supported_rows > 0, "context_specific_supported_shift",
  context_magnitude_modulated_rows > 0, "context_magnitude_modulated",
  context_interaction_evaluable_rows > 0, "context_evaluable_weak",
  default = "context_not_evaluable"
)]
setorder(gene_design_summary, -context_interaction_review_priority,
         gene_symbol)

viral_consensus <- gene_design_summary[
  gene_symbol %chin% c("DDIT3", "PPP1R15A", "ATF4", "VEGFA", "IFIH1")
]

summary_metrics <- data.table(
  metric = c(
    "context_interaction_study_context_rows",
    "context_interaction_branch_context_rows",
    "context_interaction_review_rows",
    "context_interaction_gene_designs",
    "context_interaction_evaluable_gene_designs",
    "context_direction_reversal_gene_designs",
    "context_specific_supported_gene_designs",
    "context_magnitude_modulated_gene_designs",
    "viral_consensus_rows",
    "viral_consensus_context_direction_reversal_rows",
    "median_context_interaction_score_review",
    "median_context_interaction_context_studies_review",
    "median_context_interaction_complement_studies_review"
  ),
  value = c(
    nrow(study_context_effects),
    nrow(context_effects),
    nrow(context_review),
    nrow(gene_design_summary),
    sum(gene_design_summary$context_interaction_evaluable_rows > 0),
    sum(gene_design_summary$context_direction_reversal_rows > 0),
    sum(gene_design_summary$context_specific_supported_rows > 0),
    sum(gene_design_summary$context_magnitude_modulated_rows > 0),
    nrow(viral_consensus),
    sum(viral_consensus$context_direction_reversal_rows > 0),
    median_safe(context_review$context_interaction_score),
    median_safe(context_review$context_studies),
    median_safe(context_review$complement_studies)
  )
)

fwrite(study_context_effects, file.path(
  output_dir,
  "dominant_rdg_context_interaction_study_context_effects.csv"
))
fwrite(context_effects, file.path(
  output_dir,
  "dominant_rdg_context_interaction_branch_context_effects.csv"
))
fwrite(context_review, file.path(
  output_dir,
  "dominant_rdg_context_interaction_review.csv"
))
fwrite(gene_design_summary, file.path(
  output_dir,
  "dominant_rdg_context_interaction_gene_design_summary.csv"
))
fwrite(viral_consensus, file.path(
  output_dir,
  "dominant_rdg_context_interaction_viral_consensus_summary.csv"
))
fwrite(summary_metrics, file.path(
  output_dir,
  "dominant_rdg_context_interaction_summary_metrics.csv"
))

top_plot <- context_review[
  context_interaction_evaluable == TRUE
][seq_len(min(.N, 35))]
if (nrow(top_plot)) {
  top_plot[, plot_label := paste(
    gene_symbol,
    as.character(branch_class),
    paste0(context_axis, "=", context_value),
    sep = " / "
  )]
  top_plot[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_top <- ggplot(
    top_plot,
    aes(
      x = interaction_delta,
      y = plot_label,
      color = context_interaction_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Branch: ", branch_class,
        "<br>Context: ", context_axis, "=", context_value,
        "<br>Context delta: ", round(context_delta, 3),
        "<br>Complement delta: ", round(complement_delta, 3),
        "<br>Interaction delta: ", round(interaction_delta, 3),
        "<br>Interaction q: ", round(interaction_delta_q, 4),
        "<br>Context studies: ", context_studies,
        "<br>Complement studies: ", complement_studies
      )
    )
  ) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_errorbar(
      aes(xmin = interaction_delta_lower, xmax = interaction_delta_upper),
      orientation = "y",
      height = 0.16,
      alpha = 0.6
    ) +
    geom_point(aes(size = context_review_priority), alpha = 0.85) +
    labs(
      title = "Viral RDG Context-Interaction Review",
      subtitle = "Interaction delta = context branch effect minus complement branch effect",
      x = "Context minus complement branch delta",
      y = NULL,
      color = "Context class",
      size = "Review priority"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(figure_dir, "context_interaction_viral_top_contexts.png"),
         p_top, width = 11, height = 8.2, dpi = 180)
  ggsave(file.path(figure_dir, "context_interaction_viral_top_contexts.pdf"),
         p_top, width = 11, height = 8.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir, "context_interaction_viral_top_contexts.html"),
      title = "Viral RDG context-interaction review"
    )
  }
}

ifih1_plot <- context_effects[
  gene_symbol == "IFIH1" &
    context_axis == "tissue_context" &
    branch_class %chin% c("clean_CDS", "leader_uORF")
]
if (nrow(ifih1_plot)) {
  ifih1_long <- melt(
    ifih1_plot,
    id.vars = c("gene_symbol", "branch_class", "context_axis",
                "context_value", "context_interaction_class"),
    measure.vars = list(
      delta = c("context_delta", "complement_delta"),
      lower = c("context_delta_lower", "complement_delta_lower"),
      upper = c("context_delta_upper", "complement_delta_upper")
    ),
    variable.name = "estimate_type"
  )
  ifih1_long[, estimate_type := factor(
    estimate_type,
    levels = c(1, 2),
    labels = c("context", "complement")
  )]
  ifih1_long[, branch_class := factor(
    as.character(branch_class),
    levels = c("clean_CDS", "leader_uORF")
  )]
  p_ifih1 <- ggplot(
    ifih1_long,
    aes(
      x = context_value,
      y = delta,
      ymin = lower,
      ymax = upper,
      color = estimate_type,
      group = estimate_type
    )
  ) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_pointrange(
      position = position_dodge(width = 0.45),
      linewidth = 0.45
    ) +
    facet_wrap(~ branch_class, ncol = 1, scales = "free_y") +
    labs(
      title = "IFIH1 Viral Branch Effects By Tissue Context",
      subtitle = "Context estimates are compared with the remaining non-matching tissues",
      x = NULL,
      y = "Branch allocation delta",
      color = "Estimate"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "context_interaction_ifih1_tissue.png"),
         p_ifih1, width = 10.5, height = 6.2, dpi = 180)
  ggsave(file.path(figure_dir, "context_interaction_ifih1_tissue.pdf"),
         p_ifih1, width = 10.5, height = 6.2)
}

consensus_heat <- context_effects[
  gene_symbol %chin% c("DDIT3", "PPP1R15A", "ATF4", "VEGFA", "IFIH1") &
    context_interaction_evaluable == TRUE &
    context_axis == "tissue_context"
]
if (nrow(consensus_heat)) {
  consensus_heat[, plot_row := paste(gene_symbol, as.character(branch_class),
                                     sep = " / ")]
  row_order <- consensus_heat[
    , .(score = max(context_review_priority, na.rm = TRUE)),
    by = plot_row
  ][order(score), plot_row]
  consensus_heat[, plot_row := factor(plot_row, levels = row_order)]
  p_heat <- ggplot(
    consensus_heat,
    aes(
      x = context_value,
      y = plot_row,
      fill = pmax(pmin(interaction_delta, 0.2), -0.2),
      text = paste0(
        "Gene/branch: ", gene_symbol, " / ", branch_class,
        "<br>Context: ", context_value,
        "<br>Class: ", context_interaction_class,
        "<br>Context delta: ", round(context_delta, 3),
        "<br>Complement delta: ", round(complement_delta, 3),
        "<br>Interaction delta: ", round(interaction_delta, 3),
        "<br>Context studies: ", context_studies,
        "<br>Complement studies: ", complement_studies
      )
    )
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", interaction_delta)), size = 3) +
    scale_fill_gradient2(
      low = "#245b9e",
      mid = "white",
      high = "#8b0000",
      midpoint = 0,
      limits = c(-0.2, 0.2),
      name = "Context -\ncomplement"
    ) +
    labs(
      title = "Viral Consensus Tissue Interaction Matrix",
      subtitle = "Positive values mean the branch delta is stronger in that tissue than in the complement",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  ggsave(file.path(figure_dir,
                   "context_interaction_viral_consensus_tissue_heatmap.png"),
         p_heat, width = 10.8, height = 6.4, dpi = 180)
  ggsave(file.path(figure_dir,
                   "context_interaction_viral_consensus_tissue_heatmap.pdf"),
         p_heat, width = 10.8, height = 6.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_heat,
      file.path(figure_dir,
                "context_interaction_viral_consensus_tissue_heatmap.html"),
      title = "Viral consensus tissue interaction matrix"
    )
  }
}

message("Saved RDG context-interaction outputs in: ", output_dir)
print(summary_metrics)
