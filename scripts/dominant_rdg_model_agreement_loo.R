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
output_dir <- file.path(analysis_dir, "dominant_rdg_model_agreement_loo")
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
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

weighted_mean_safe <- function(x, w = NULL) {
  x <- safe_num(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- safe_num(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

random_effects_summary <- function(effect, se) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  keep <- is.finite(effect) & is.finite(se) & se > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L, effect = NA_real_, se = NA_real_, lower = NA_real_,
      upper = NA_real_, p = NA_real_, tau2 = NA_real_, i2 = NA_real_,
      direction_fraction = NA_real_, ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  vi <- pmax(se[keep], 1e-6)^2
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  q_stat <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q_stat - df) / c_val) else 0
  wi_re <- 1 / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  direction_fraction <- if (sign(mu) == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign(mu))
  }
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    p = 2 * pnorm(-abs(mu / se_mu)),
    tau2 = tau2,
    i2 = if (q_stat > 0 && df > 0) max(0, (q_stat - df) / q_stat) else 0,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

classify_loo_support <- function(model, n, effect, p, ci, direction_fraction) {
  if (!is.finite(effect) || !is.finite(p) || !is.finite(direction_fraction)) {
    return("not_evaluable")
  }
  if (model == "hierarchical_joint") {
    return(fcase(
      n >= 2 & ci & abs(effect) >= 0.04 & direction_fraction >= 0.70,
      "replicated_shift",
      n >= 2 & p <= 0.20 & abs(effect) >= 0.03 &
        direction_fraction >= 0.65,
      "replicated_review",
      n == 1 & p <= 0.05 & abs(effect) >= 0.08,
      "single_study_strong",
      default = "neutral"
    ))
  }
  fcase(
    n >= 2 & ci & p <= 0.20 & abs(effect) >= 0.035 &
      direction_fraction >= 0.70,
    "replicated_shift",
    n >= 2 & (ci | p <= 0.30) & abs(effect) >= 0.025 &
      direction_fraction >= 0.65,
    "replicated_review",
    n == 1 & p <= 0.05 & abs(effect) >= 0.08,
    "single_study_strong",
    default = "neutral"
  )
}

standardize_hierarchical <- function(dt) {
  out <- copy(dt)
  setnames(out, c(
    "study_delta", "study_delta_se", "study_delta_p",
    "study_context_rows", "study_strict_context_rows",
    "study_max_context_priority", "study_top_context",
    "study_top_case_condition", "study_top_control_conditions",
    "study_top_cell_line", "study_top_tissue", "study_top_match_scope"
  ), c(
    "model_delta", "model_delta_se", "model_delta_p",
    "model_context_rows", "model_strict_context_rows",
    "model_max_context_priority", "model_top_context",
    "model_top_case_condition", "model_top_control_conditions",
    "model_top_cell_line", "model_top_tissue", "model_top_match_scope"
  ), skip_absent = TRUE)
  out[, model := "hierarchical_joint"]
  out
}

standardize_dm <- function(dt) {
  out <- copy(dt)
  setnames(out, c(
    "dm_study_delta", "dm_study_delta_se", "dm_study_delta_p",
    "dm_study_context_rows", "dm_study_strict_context_rows",
    "dm_study_max_priority", "dm_study_top_context",
    "dm_study_top_match_scope", "dm_study_top_prior_source"
  ), c(
    "model_delta", "model_delta_se", "model_delta_p",
    "model_context_rows", "model_strict_context_rows",
    "model_max_context_priority", "model_top_context",
    "model_top_match_scope", "model_top_prior_source"
  ), skip_absent = TRUE)
  out[, model := "dirichlet_multinomial"]
  out
}

message("Leave-one-study-out RDG model-agreement analysis:")
message("  1. Preserve study-level branch evidence")
message("  2. Refit each supported branch after omitting one study")
message("  3. Distinguish direction-robust from support-robust consensus")

agreement <- read_dt(
  file.path(analysis_dir, "dominant_rdg_model_agreement",
            "dominant_rdg_model_agreement_gene_design.csv"),
  required = TRUE
)
branch_comparison <- read_dt(
  file.path(analysis_dir, "dominant_rdg_model_agreement",
            "dominant_rdg_model_agreement_branch_comparison.csv"),
  required = TRUE
)
hierarchical_studies <- standardize_hierarchical(read_dt(
  file.path(analysis_dir, "dominant_rdg_hierarchical_joint_allocation",
            "dominant_rdg_hierarchical_joint_allocation_study_effects.csv"),
  required = TRUE
))
dm_studies <- standardize_dm(read_dt(
  file.path(analysis_dir, "dominant_rdg_dirichlet_multinomial",
            "dominant_rdg_dirichlet_multinomial_study_effects.csv"),
  required = TRUE
))

key <- c("gene_symbol", "tx_id", "design_family")
branch_key <- c(key, "branch_class")
study_key <- c(branch_key, "study")
required_study <- c(
  study_key, "model", "model_delta", "model_delta_se",
  "model_context_rows", "model_strict_context_rows",
  "model_max_context_priority", "model_top_context"
)
studies <- rbindlist(list(hierarchical_studies, dm_studies),
                     use.names = TRUE, fill = TRUE)
missing <- setdiff(required_study, names(studies))
if (length(missing)) {
  stop("Study-effect inputs are missing columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}
for (column in c(
  "model_delta", "model_delta_se", "model_delta_p",
  "model_context_rows", "model_strict_context_rows",
  "model_max_context_priority"
)) {
  studies[, (column) := safe_num(get(column))]
}
studies <- studies[
  is.finite(model_delta) & is.finite(model_delta_se) & model_delta_se > 0
]
studies[, study_model_direction := fcase(
  model_delta > 0, "positive",
  model_delta < 0, "negative",
  default = "zero"
)]
setorder(studies, gene_symbol, tx_id, design_family, study, branch_class,
         model)

study_wide <- dcast(
  studies,
  gene_symbol + tx_id + design_family + branch_class + study ~ model,
  value.var = c(
    "model_delta", "model_delta_se", "model_delta_p",
    "model_context_rows", "model_strict_context_rows",
    "model_max_context_priority", "model_top_context",
    "model_top_match_scope", "model_top_prior_source"
  ),
  fill = NA
)
for (column in c(
  "model_delta_hierarchical_joint",
  "model_delta_dirichlet_multinomial",
  "model_max_context_priority_hierarchical_joint",
  "model_max_context_priority_dirichlet_multinomial"
)) {
  if (!column %in% names(study_wide)) study_wide[, (column) := NA_real_]
  study_wide[, (column) := safe_num(get(column))]
}
study_wide[, both_models_evaluable :=
  is.finite(model_delta_hierarchical_joint) &
    is.finite(model_delta_dirichlet_multinomial)]
study_wide[, study_model_direction_agree :=
  both_models_evaluable &
    sign(model_delta_hierarchical_joint) ==
      sign(model_delta_dirichlet_multinomial) &
    sign(model_delta_hierarchical_joint) != 0]
study_wide[, study_model_direction_opposite :=
  both_models_evaluable &
    sign(model_delta_hierarchical_joint) !=
      sign(model_delta_dirichlet_multinomial) &
    sign(model_delta_hierarchical_joint) != 0 &
    sign(model_delta_dirichlet_multinomial) != 0]
study_wide[, study_evidence_class := fcase(
  study_model_direction_agree, "both_models_same_direction",
  study_model_direction_opposite, "both_models_opposite_direction",
  is.finite(model_delta_hierarchical_joint) &
    !is.finite(model_delta_dirichlet_multinomial),
  "hierarchical_joint_only",
  !is.finite(model_delta_hierarchical_joint) &
    is.finite(model_delta_dirichlet_multinomial),
  "dirichlet_multinomial_only",
  default = "not_evaluable"
)]
study_wide[, study_mean_delta := rowMeans(cbind(
  model_delta_hierarchical_joint,
  model_delta_dirichlet_multinomial
), na.rm = TRUE)]
study_wide[!is.finite(study_mean_delta), study_mean_delta := NA_real_]
study_wide[, study_max_priority := pmax(
  finite_or_zero(model_max_context_priority_hierarchical_joint),
  finite_or_zero(model_max_context_priority_dirichlet_multinomial)
)]
study_wide[, study_evidence_priority := clip01(
  0.45 * pmin(abs(study_mean_delta) / 0.15, 1) +
    0.25 * study_max_priority +
    0.20 * as.numeric(study_model_direction_agree) +
    0.10 * as.numeric(study_model_direction_opposite)
)]
study_wide[, study_top_context := fcase(
  nzchar(as.character(model_top_context_hierarchical_joint)),
  as.character(model_top_context_hierarchical_joint),
  nzchar(as.character(model_top_context_dirichlet_multinomial)),
  as.character(model_top_context_dirichlet_multinomial),
  default = ""
)]
setorder(study_wide, gene_symbol, tx_id, design_family, study,
         -study_evidence_priority, branch_class)

study_compact <- study_wide[, {
  top <- .SD[1]
  data.table(
    study_evaluable_branches = .N,
    study_both_model_branches = sum(both_models_evaluable),
    study_same_direction_branches = sum(study_model_direction_agree),
    study_opposite_direction_branches = sum(study_model_direction_opposite),
    study_top_branch_class = top$branch_class,
    study_top_evidence_class = top$study_evidence_class,
    study_top_mean_delta = top$study_mean_delta,
    study_top_hierarchical_delta = top$model_delta_hierarchical_joint,
    study_top_dm_delta = top$model_delta_dirichlet_multinomial,
    study_top_priority = top$study_evidence_priority,
    study_top_context = top$study_top_context,
    study_branch_summary = paste0(
      branch_class, "=",
      round(model_delta_hierarchical_joint, 3), "/",
      round(model_delta_dirichlet_multinomial, 3),
      collapse = "; "
    )
  )
}, by = c(key, "study")]
setorder(study_compact, gene_symbol, design_family, -study_top_priority, study)

run_loo <- function(dt) {
  dt[, {
    full <- random_effects_summary(model_delta, model_delta_se)
    if (.N < 2) {
      data.table(
        omitted_study = study[1],
        model_studies = .N,
        remaining_studies = 0L,
        full_delta = full$effect,
        loo_delta = NA_real_,
        loo_delta_se = NA_real_,
        loo_delta_lower = NA_real_,
        loo_delta_upper = NA_real_,
        loo_delta_p = NA_real_,
        loo_direction_fraction = NA_real_,
        loo_ci_excludes_zero = FALSE,
        loo_same_direction = FALSE,
        loo_abs_change = NA_real_,
        omitted_study_delta = model_delta[1],
        omitted_study_context = model_top_context[1],
        loo_support_class = "insufficient_studies",
        loo_support_retained = FALSE
      )
    } else {
      rbindlist(lapply(seq_len(.N), function(i) {
        remaining <- .SD[-i]
        re <- random_effects_summary(
          remaining$model_delta,
          remaining$model_delta_se
        )
        support_class <- classify_loo_support(
          model[1], re$n, re$effect, re$p, re$ci_excludes_zero,
          re$direction_fraction
        )
        data.table(
          omitted_study = study[i],
          model_studies = .N,
          remaining_studies = re$n,
          full_delta = full$effect,
          loo_delta = re$effect,
          loo_delta_se = re$se,
          loo_delta_lower = re$lower,
          loo_delta_upper = re$upper,
          loo_delta_p = re$p,
          loo_direction_fraction = re$direction_fraction,
          loo_ci_excludes_zero = re$ci_excludes_zero,
          loo_same_direction =
            is.finite(re$effect) & sign(re$effect) == sign(full$effect) &
              sign(full$effect) != 0,
          loo_abs_change = abs(re$effect - full$effect),
          omitted_study_delta = model_delta[i],
          omitted_study_context = model_top_context[i],
          loo_support_class = support_class,
          loo_support_retained =
            support_class %chin% c(
              "replicated_shift", "replicated_review", "single_study_strong"
            )
        )
      }), use.names = TRUE, fill = TRUE)
    }
  }, by = c(branch_key, "model")]
}

loo_rows <- run_loo(studies)
setorder(loo_rows, gene_symbol, design_family, branch_class, model,
         -loo_abs_change)
loo_summary <- loo_rows[, .(
  model_studies = max_safe(model_studies, 0),
  loo_iterations = .N,
  full_delta = first(full_delta),
  loo_fraction_same_direction = mean(loo_same_direction, na.rm = TRUE),
  loo_fraction_support_retained = mean(loo_support_retained, na.rm = TRUE),
  loo_min_delta = min_safe(loo_delta),
  loo_max_delta = max_safe(loo_delta),
  loo_max_abs_change = max_safe(loo_abs_change, 0),
  loo_worst_omitted_study = omitted_study[which.max(loo_abs_change)][1],
  loo_worst_omitted_context =
    omitted_study_context[which.max(loo_abs_change)][1],
  loo_direction_robust = mean(loo_same_direction, na.rm = TRUE) >= 0.80,
  loo_support_robust = mean(loo_support_retained, na.rm = TRUE) >= 0.67
), by = c(branch_key, "model")]
loo_summary[, loo_model_stability_class := fcase(
  model_studies < 2, "insufficient_studies",
  loo_direction_robust & loo_support_robust, "loo_support_robust",
  loo_direction_robust, "loo_direction_robust_support_fragile",
  default = "loo_direction_fragile"
)]

loo_branch <- dcast(
  loo_summary,
  gene_symbol + tx_id + design_family + branch_class ~ model,
  value.var = c(
    "model_studies",
    "loo_fraction_same_direction",
    "loo_fraction_support_retained",
    "loo_max_abs_change",
    "loo_worst_omitted_study",
    "loo_worst_omitted_context",
    "loo_direction_robust",
    "loo_support_robust",
    "loo_model_stability_class"
  ),
  fill = NA
)
loo_branch <- merge(
  branch_comparison,
  loo_branch,
  by = branch_key,
  all.x = TRUE,
  sort = FALSE
)
for (model in c("hierarchical_joint", "dirichlet_multinomial")) {
  dir_col <- paste0("loo_fraction_same_direction_", model)
  support_col <- paste0("loo_fraction_support_retained_", model)
  study_col <- paste0("model_studies_", model)
  if (!dir_col %in% names(loo_branch)) loo_branch[, (dir_col) := NA_real_]
  if (!support_col %in% names(loo_branch)) {
    loo_branch[, (support_col) := NA_real_]
  }
  if (!study_col %in% names(loo_branch)) loo_branch[, (study_col) := NA_real_]
}
loo_branch[, loo_min_direction_fraction := pmin(
  fifelse(hierarchical_joint_supported,
          loo_fraction_same_direction_hierarchical_joint, Inf),
  fifelse(dm_supported,
          loo_fraction_same_direction_dirichlet_multinomial, Inf),
  na.rm = TRUE
)]
loo_branch[!is.finite(loo_min_direction_fraction),
           loo_min_direction_fraction := NA_real_]
loo_branch[, loo_min_support_retention := pmin(
  fifelse(hierarchical_joint_supported,
          loo_fraction_support_retained_hierarchical_joint, Inf),
  fifelse(dm_supported,
          loo_fraction_support_retained_dirichlet_multinomial, Inf),
  na.rm = TRUE
)]
loo_branch[!is.finite(loo_min_support_retention),
           loo_min_support_retention := NA_real_]
loo_branch[, loo_min_model_studies := pmin(
  fifelse(hierarchical_joint_supported,
          model_studies_hierarchical_joint, Inf),
  fifelse(dm_supported,
          model_studies_dirichlet_multinomial, Inf),
  na.rm = TRUE
)]
loo_branch[!is.finite(loo_min_model_studies), loo_min_model_studies := NA_real_]
loo_branch[, relevant_supported_branch :=
  hierarchical_joint_supported | dm_supported]
loo_branch[, loo_branch_stability_class := fcase(
  !relevant_supported_branch, "not_supported",
  !is.finite(loo_min_model_studies) | loo_min_model_studies < 2,
  "loo_not_evaluable",
  is.finite(loo_min_direction_fraction) & loo_min_direction_fraction >= 0.80 &
    is.finite(loo_min_support_retention) & loo_min_support_retention >= 0.67,
  "loo_support_robust",
  is.finite(loo_min_direction_fraction) & loo_min_direction_fraction >= 0.80,
  "loo_direction_robust_support_fragile",
  default = "loo_direction_fragile"
)]
setorder(loo_branch, gene_symbol, design_family, branch_class)

gene_loo <- loo_branch[relevant_supported_branch == TRUE, .(
  loo_relevant_supported_branches = .N,
  loo_support_robust_branches =
    sum(loo_branch_stability_class == "loo_support_robust"),
  loo_direction_robust_support_fragile_branches =
    sum(loo_branch_stability_class ==
          "loo_direction_robust_support_fragile"),
  loo_direction_fragile_branches =
    sum(loo_branch_stability_class == "loo_direction_fragile"),
  loo_not_evaluable_branches =
    sum(loo_branch_stability_class == "loo_not_evaluable"),
  loo_min_direction_fraction = min_safe(loo_min_direction_fraction, 0),
  loo_min_support_retention = min_safe(loo_min_support_retention, 0),
  loo_min_model_studies = min_safe(loo_min_model_studies, 0),
  loo_support_robust_branch_names =
    collapse_unique(branch_class[
      loo_branch_stability_class == "loo_support_robust"
    ], 3L),
  loo_fragile_branch_names =
    collapse_unique(branch_class[
      loo_branch_stability_class != "loo_support_robust"
    ], 3L),
  loo_worst_omitted_studies = collapse_unique(c(
    loo_worst_omitted_study_hierarchical_joint,
    loo_worst_omitted_study_dirichlet_multinomial
  ), 6L)
), by = key]
agreement_loo <- merge(agreement, gene_loo, by = key, all.x = TRUE,
                       sort = FALSE)
for (column in c(
  "loo_relevant_supported_branches", "loo_support_robust_branches",
  "loo_direction_robust_support_fragile_branches",
  "loo_direction_fragile_branches", "loo_not_evaluable_branches",
  "loo_min_direction_fraction",
  "loo_min_support_retention", "loo_min_model_studies"
)) {
  agreement_loo[is.na(get(column)), (column) := 0]
}
for (column in c(
  "loo_support_robust_branch_names", "loo_fragile_branch_names",
  "loo_worst_omitted_studies"
)) {
  agreement_loo[is.na(get(column)), (column) := ""]
}
agreement_loo[, loo_agreement_class := fcase(
  loo_relevant_supported_branches == 0, "no_supported_branch",
  loo_direction_fragile_branches > 0, "loo_direction_fragile",
  loo_not_evaluable_branches > 0, "loo_not_evaluable",
  loo_support_robust_branches == loo_relevant_supported_branches,
  "loo_support_robust",
  default = "loo_direction_robust_support_fragile"
)]
agreement_loo[, loo_consensus_class := fcase(
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  ) & loo_agreement_class == "loo_support_robust",
  "consensus_loo_support_robust",
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  ) & loo_agreement_class == "loo_direction_robust_support_fragile",
  "consensus_loo_direction_robust",
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  ) & loo_agreement_class == "loo_not_evaluable",
  "consensus_loo_not_evaluable",
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  ),
  "consensus_loo_fragile",
  model_agreement_tier == "D_model_disagreement_review",
  "disagreement_review",
  default = "non_consensus_review"
)]
agreement_loo[, loo_stability_score := clip01(
  0.55 * loo_min_direction_fraction +
    0.30 * loo_min_support_retention +
    0.15 * pmin(loo_min_model_studies / 3, 1)
)]
agreement_loo[loo_agreement_class == "loo_not_evaluable",
              loo_stability_score := 0]
agreement_loo[, loo_review_priority := pmin(
  1.5,
  model_agreement_review_priority +
    0.12 * loo_stability_score +
    0.10 * as.numeric(loo_consensus_class == "consensus_loo_fragile") +
    0.06 * as.numeric(loo_consensus_class ==
                        "consensus_loo_not_evaluable") +
    0.06 * as.numeric(loo_agreement_class ==
                        "loo_direction_robust_support_fragile")
)]
setorder(agreement_loo, -loo_review_priority, -loo_stability_score,
         gene_symbol, design_family)

viral_review <- agreement_loo[design_family == "viral_infection"]
setorder(viral_review, -loo_review_priority, -loo_stability_score, gene_symbol)
fragile_review <- agreement_loo[
  loo_agreement_class %chin% c(
    "loo_direction_fragile", "loo_direction_robust_support_fragile",
    "loo_not_evaluable"
  )
]
setorder(fragile_review, -loo_review_priority, loo_stability_score)

summarize_scope <- function(dt, scope) {
  consensus <- dt$model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  )
  evaluable <- dt$loo_min_model_studies >= 2
  data.table(
    scope = scope,
    loo_support_rule = "within_branch_p_proxy_no_global_BH",
    gene_designs = nrow(dt),
    consensus_gene_designs = sum(consensus),
    loo_evaluable_gene_designs = sum(evaluable),
    consensus_loo_evaluable_gene_designs = sum(consensus & evaluable),
    consensus_loo_support_robust =
      sum(dt$loo_consensus_class == "consensus_loo_support_robust"),
    consensus_loo_direction_robust =
      sum(dt$loo_consensus_class == "consensus_loo_direction_robust"),
    consensus_loo_not_evaluable =
      sum(dt$loo_consensus_class == "consensus_loo_not_evaluable"),
    consensus_loo_fragile =
      sum(dt$loo_consensus_class == "consensus_loo_fragile"),
    all_loo_support_robust =
      sum(dt$loo_agreement_class == "loo_support_robust"),
    all_loo_direction_robust_support_fragile =
      sum(dt$loo_agreement_class == "loo_direction_robust_support_fragile"),
    all_loo_not_evaluable =
      sum(dt$loo_agreement_class == "loo_not_evaluable"),
    all_loo_direction_fragile =
      sum(dt$loo_agreement_class == "loo_direction_fragile"),
    median_loo_stability_score = median_safe(dt$loo_stability_score),
    median_loo_min_direction_fraction =
      median_safe(dt$loo_min_direction_fraction),
    median_loo_min_support_retention =
      median_safe(dt$loo_min_support_retention),
    median_evaluable_loo_stability_score =
      median_safe(dt$loo_stability_score[evaluable]),
    median_evaluable_loo_min_direction_fraction =
      median_safe(dt$loo_min_direction_fraction[evaluable]),
    median_evaluable_loo_min_support_retention =
      median_safe(dt$loo_min_support_retention[evaluable])
  )
}
summary_metrics <- rbindlist(list(
  summarize_scope(agreement_loo, "all"),
  summarize_scope(viral_review, "design_family:viral_infection")
))

fwrite(study_wide, file.path(
  output_dir, "dominant_rdg_model_agreement_study_branch_evidence.csv"
))
fwrite(study_compact, file.path(
  output_dir, "dominant_rdg_model_agreement_study_evidence.csv"
))
fwrite(loo_rows, file.path(
  output_dir, "dominant_rdg_model_agreement_loo_rows.csv"
))
fwrite(loo_branch, file.path(
  output_dir, "dominant_rdg_model_agreement_loo_branch_summary.csv"
))
fwrite(agreement_loo, file.path(
  output_dir, "dominant_rdg_model_agreement_loo_gene_design.csv"
))
fwrite(viral_review, file.path(
  output_dir, "dominant_rdg_model_agreement_loo_viral_review.csv"
))
fwrite(fragile_review, file.path(
  output_dir, "dominant_rdg_model_agreement_loo_fragile_review.csv"
))
fwrite(summary_metrics, file.path(
  output_dir, "dominant_rdg_model_agreement_loo_summary_metrics.csv"
))

viral_plot <- viral_review[
  loo_min_model_studies >= 2 & (
    model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
    ) | loo_review_priority >= 0.45
  )
][1:min(.N, 30)]
if (nrow(viral_plot)) {
  viral_plot[, plot_label := paste0(gene_symbol, " / ", loo_consensus_class)]
  viral_plot[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_viral <- ggplot(
    viral_plot,
    aes(
      x = loo_stability_score,
      y = plot_label,
      fill = loo_consensus_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Agreement tier: ", model_agreement_tier,
        "<br>LOO class: ", loo_consensus_class,
        "<br>Direction retention: ", round(loo_min_direction_fraction, 3),
        "<br>Support retention: ", round(loo_min_support_retention, 3),
        "<br>Minimum model studies: ", loo_min_model_studies,
        "<br>Worst omitted studies: ", loo_worst_omitted_studies
      )
    )
  ) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(
      consensus_loo_support_robust = "#245b9e",
      consensus_loo_direction_robust = "#7ba7cf",
      consensus_loo_not_evaluable = "#7d8795",
      consensus_loo_fragile = "#a84232",
      disagreement_review = "#c47a28",
      non_consensus_review = "#9aa4b2"
    )) +
    labs(
      title = "Viral RDG Leave-One-Study-Out Stability",
      subtitle = "Only study-evaluable rows shown; consensus is stronger when direction and support survive omission",
      x = "Leave-one-study-out stability score",
      y = NULL,
      fill = "LOO class"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "model_agreement_loo_viral_stability.png"),
         p_viral, width = 11, height = 7.8, dpi = 180)
  ggsave(file.path(figure_dir, "model_agreement_loo_viral_stability.pdf"),
         p_viral, width = 11, height = 7.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_viral,
      file.path(figure_dir, "model_agreement_loo_viral_stability.html"),
      title = "Viral RDG leave-one-study-out stability"
    )
  }
}

consensus_plot <- agreement_loo[
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus",
    "D_model_disagreement_review"
  )
]
if (nrow(consensus_plot)) {
  label_dt <- consensus_plot[
    design_family == "viral_infection" |
      loo_consensus_class == "consensus_loo_fragile"
  ][order(-loo_review_priority)][1:min(.N, 24)]
  p_consensus <- ggplot(
    consensus_plot,
    aes(
      x = loo_min_direction_fraction,
      y = loo_min_support_retention,
      color = loo_consensus_class,
      size = pmax(model_agreement_score, 0.05),
      text = paste0(
        "Gene/design: ", gene_symbol, " / ", design_family,
        "<br>Agreement tier: ", model_agreement_tier,
        "<br>LOO class: ", loo_consensus_class,
        "<br>Direction retention: ", round(loo_min_direction_fraction, 3),
        "<br>Support retention: ", round(loo_min_support_retention, 3),
        "<br>Minimum model studies: ", loo_min_model_studies
      )
    )
  ) +
    geom_vline(xintercept = 0.80, color = "#98a2b3", linetype = "dashed") +
    geom_hline(yintercept = 0.67, color = "#98a2b3", linetype = "dashed") +
    geom_point(alpha = 0.78) +
    geom_text(
      data = label_dt,
      aes(label = gene_symbol),
      size = 2.5, vjust = -0.8, check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      consensus_loo_support_robust = "#245b9e",
      consensus_loo_direction_robust = "#7ba7cf",
      consensus_loo_fragile = "#a84232",
      disagreement_review = "#c47a28",
      non_consensus_review = "#9aa4b2"
    )) +
    labs(
      title = "RDG Consensus Stability Under Study Omission",
      subtitle = "Dashed lines mark direction-robust and support-robust thresholds",
      x = "Minimum fraction retaining branch-effect direction",
      y = "Minimum fraction retaining branch support",
      color = "LOO class",
      size = "Agreement score"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "model_agreement_loo_consensus_stability.png"),
         p_consensus, width = 10.8, height = 7.8, dpi = 180)
  ggsave(file.path(figure_dir, "model_agreement_loo_consensus_stability.pdf"),
         p_consensus, width = 10.8, height = 7.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_consensus,
      file.path(figure_dir, "model_agreement_loo_consensus_stability.html"),
      title = "RDG consensus stability under study omission"
    )
  }
}

message("Saved leave-one-study-out model-agreement outputs in: ", output_dir)
print(summary_metrics)
