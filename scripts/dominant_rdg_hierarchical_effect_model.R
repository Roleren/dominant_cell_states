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
output_dir <- file.path(analysis_dir, "dominant_rdg_hierarchical")
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

safe_num <- function(x) suppressWarnings(as.numeric(x))

pct_from_log2 <- function(x) {
  x <- safe_num(x)
  out <- (2^x - 1) * 100
  out[!is.finite(out)] <- NA_real_
  out
}

bool_col <- function(dt, col, default = FALSE) {
  if (!col %in% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[col]]
  if (is.logical(x)) return(fifelse(is.na(x), default, x))
  if (is.numeric(x)) return(fifelse(is.na(x), default, x != 0))
  y <- tolower(as.character(x))
  fifelse(is.na(y), default, y %chin% c("true", "t", "1", "yes"))
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

sum_safe <- function(x) {
  x <- safe_num(x)
  sum(x[is.finite(x)], na.rm = TRUE)
}

weighted_mean_safe <- function(x, w = NULL) {
  x <- safe_num(x)
  if (is.null(w)) w <- rep(1, length(x))
  w <- safe_num(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

quantile_safe <- function(x, prob = 0.95, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 8,
                             na.rm = TRUE))
}

effect_se_from_rows <- function(dt) {
  se <- rep(NA_real_, nrow(dt))
  if ("bootstrap_fpkm_log2_fc_sd" %in% names(dt)) {
    se <- safe_num(dt$bootstrap_fpkm_log2_fc_sd)
  }
  if (all(!is.finite(se)) &&
      all(c("bootstrap_fpkm_log2_fc_lower",
            "bootstrap_fpkm_log2_fc_upper") %in% names(dt))) {
    lower <- safe_num(dt$bootstrap_fpkm_log2_fc_lower)
    upper <- safe_num(dt$bootstrap_fpkm_log2_fc_upper)
    se <- abs(upper - lower) / (2 * 1.96)
  }
  if (all(c("case_sum_raw_counts", "control_sum_raw_counts") %in% names(dt))) {
    fallback <- 1 / sqrt(pmax(
      1,
      safe_num(dt$case_sum_raw_counts) + safe_num(dt$control_sum_raw_counts)
    ))
    fallback <- pmin(pmax(fallback, 0.05), 1)
    replace <- !is.finite(se) | se <= 0
    se[replace] <- fallback[replace]
  }
  se[!is.finite(se) | se <= 0] <- 0.5
  pmin(pmax(se, 0.04), 2)
}

random_effects_summary <- function(effect, se, row_weight = NULL) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  if (is.null(row_weight)) row_weight <- rep(1, length(effect))
  row_weight <- safe_num(row_weight)
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(row_weight) & row_weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L,
      effect = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  vi <- se[keep]^2
  external_weight <- row_weight[keep] / mean(row_weight[keep])
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
    tau2 = tau2,
    i2 = i2,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

posterior_from_priors <- function(direct_effect, direct_se,
                                  gene_effect, gene_se,
                                  design_effect, design_se,
                                  global_effect, global_se,
                                  prior_strength = 0.35) {
  direct_info <- if (is.finite(direct_effect) && is.finite(direct_se) &&
                     direct_se > 0) {
    1 / direct_se^2
  } else {
    0
  }
  prior_effects <- c(gene_effect, design_effect, global_effect)
  prior_ses <- c(gene_se, design_se, global_se)
  prior_info <- ifelse(is.finite(prior_effects) & is.finite(prior_ses) &
                         prior_ses > 0,
                       1 / (prior_ses^2 + 0.05^2), 0)
  prior_info <- pmin(prior_info, c(35, 35, 15))
  if (!sum(prior_info) > 0) {
    prior_effect <- NA_real_
    prior_se <- NA_real_
    prior_total_info <- 0
  } else {
    prior_effect <- sum(prior_effects * prior_info, na.rm = TRUE) /
      sum(prior_info)
    prior_total_info <- sum(prior_info)
    prior_se <- sqrt(1 / prior_total_info)
  }
  if (direct_info > 0) {
    total_info <- direct_info + prior_strength * prior_total_info
    posterior <- (
      direct_effect * direct_info +
        prior_effect * prior_strength * prior_total_info
    ) / total_info
    posterior_se <- sqrt(1 / total_info)
    evidence_weight <- direct_info / total_info
  } else if (prior_total_info > 0) {
    posterior <- prior_effect
    posterior_se <- prior_se
    evidence_weight <- 0
  } else {
    posterior <- NA_real_
    posterior_se <- NA_real_
    evidence_weight <- NA_real_
  }
  posterior_se <- pmax(posterior_se, 0.04)
  data.table(
    hierarchical_prior_log2fc = prior_effect,
    hierarchical_prior_se = prior_se,
    hierarchical_posterior_log2fc = posterior,
    hierarchical_posterior_se = posterior_se,
    hierarchical_posterior_log2fc_lower = posterior - 1.96 * posterior_se,
    hierarchical_posterior_log2fc_upper = posterior + 1.96 * posterior_se,
    hierarchical_evidence_weight = evidence_weight,
    hierarchical_prior_info = prior_total_info,
    hierarchical_direct_info = direct_info
  )
}

summarise_groups <- function(dt, by_cols) {
  dt[, {
    all_re <- random_effects_summary(effect_log2fc, effect_se, row_weight)
    strict_re <- random_effects_summary(
      effect_log2fc[strict_support == TRUE],
      effect_se[strict_support == TRUE],
      row_weight[strict_support == TRUE]
    )
    supported_re <- random_effects_summary(
      effect_log2fc[supported_support == TRUE],
      effect_se[supported_support == TRUE],
      row_weight[supported_support == TRUE]
    )
    use_re <- if (strict_re$n[[1]] >= 2L) strict_re else
      if (supported_re$n[[1]] >= 2L) supported_re else all_re
    data.table(
      n_rows = .N,
      n_supported_rows = sum(supported_support, na.rm = TRUE),
      n_strict_rows = sum(strict_support, na.rm = TRUE),
      n_studies = uniqueN(study),
      n_strict_studies = uniqueN(study[strict_support == TRUE]),
      n_designs = uniqueN(matched_design_signature),
      n_cell_lines = uniqueN(CELL_LINE),
      n_tissues = uniqueN(TISSUE),
      total_case_samples = sum_safe(n_case_samples),
      total_control_samples = sum_safe(n_control_samples),
      total_case_counts = sum_safe(case_sum_raw_counts),
      total_control_counts = sum_safe(control_sum_raw_counts),
      median_effect_log2fc = median_safe(effect_log2fc),
      weighted_effect_log2fc = weighted_mean_safe(effect_log2fc, row_weight),
      fraction_positive_rows = mean(effect_log2fc > 0, na.rm = TRUE),
      strongest_condition = first_nonempty(
        case_condition[order(-abs(effect_log2fc))]
      ),
      strongest_control = first_nonempty(
        control_conditions[order(-abs(effect_log2fc))]
      ),
      strongest_study = first_nonempty(study[order(-abs(effect_log2fc))]),
      nearest_designs = collapse_unique(matched_design_signature, 5L),
      random_effects_log2fc = all_re$effect,
      random_effects_se = all_re$se,
      random_effects_log2fc_lower = all_re$lower,
      random_effects_log2fc_upper = all_re$upper,
      random_effects_i2 = all_re$i2,
      random_effects_direction_fraction = all_re$direction_fraction,
      random_effects_ci_excludes_zero = all_re$ci_excludes_zero,
      strict_random_effects_log2fc = strict_re$effect,
      strict_random_effects_se = strict_re$se,
      strict_random_effects_log2fc_lower = strict_re$lower,
      strict_random_effects_log2fc_upper = strict_re$upper,
      strict_random_effects_i2 = strict_re$i2,
      strict_random_effects_direction_fraction = strict_re$direction_fraction,
      strict_random_effects_ci_excludes_zero = strict_re$ci_excludes_zero,
      supported_random_effects_log2fc = supported_re$effect,
      supported_random_effects_se = supported_re$se,
      supported_random_effects_log2fc_lower = supported_re$lower,
      supported_random_effects_log2fc_upper = supported_re$upper,
      preferred_direct_source = if (strict_re$n[[1]] >= 2L) {
        "strict_random_effects"
      } else if (supported_re$n[[1]] >= 2L) {
        "supported_random_effects"
      } else {
        "all_rows_random_effects"
      },
      direct_log2fc = use_re$effect,
      direct_se = use_re$se,
      direct_log2fc_lower = use_re$lower,
      direct_log2fc_upper = use_re$upper,
      direct_i2 = use_re$i2,
      direct_direction_fraction = use_re$direction_fraction,
      direct_ci_excludes_zero = use_re$ci_excludes_zero
    )
  }, by = by_cols]
}

message("RDG hierarchical clean-CDS effect model:")
message("  1. Summarize matched clean-CDS effects by gene x design family")
message("  2. Build empirical gene/design/global priors")
message("  3. Shrink weak matched rows toward those priors")
message("  4. Run leave-one-study-out calibration")
message("  5. Widen posterior intervals with residual-calibrated conformal bands")
message("  6. Save atlas inputs")

cds_file <- file.path(
  analysis_dir, "dominant_next_model",
  "dominant_next_model_matched_control_cds_buffering.csv"
)
branch_file <- file.path(
  analysis_dir, "dominant_rdg_branch_usage",
  "dominant_rdg_branch_usage_gene_design_summary.csv"
)
count_aware_file <- file.path(
  analysis_dir, "dominant_rdg_branch_allocation",
  "dominant_rdg_branch_allocation_gene_design_summary.csv"
)
postviral_file <- file.path(
  analysis_dir, "postviral_fatigue_outputs",
  "postviral_fatigue_matched_gene_design_evidence.csv"
)
candidate_file <- file.path(
  analysis_dir, "dominant_next_model",
  "dominant_next_model_candidate_rankings_v4.csv"
)

cds <- read_dt(cds_file, required = TRUE)
branch_usage <- read_dt(branch_file)
count_aware_branch_allocation <- read_dt(count_aware_file)
postviral <- read_dt(postviral_file)
candidate <- read_dt(candidate_file)

if (!"design_family" %in% names(cds)) cds[, design_family := perturbation_class]
cds[is.na(design_family) | !nzchar(design_family), design_family := "unknown"]
cds <- cds[feature_id == "clean_CDS" | feature_type == "clean_CDS"]

cds[, effect_log2fc := fifelse(
  is.finite(safe_num(shrunk_matched_fpkm_log2_fc)),
  safe_num(shrunk_matched_fpkm_log2_fc),
  safe_num(matched_fpkm_log2_fc_pseudocount1)
)]
cds[, effect_se := effect_se_from_rows(.SD)]
cds[, strict_support := bool_col(.SD, "matched_feature_passes_sample_gate") &
      bool_col(.SD, "matched_feature_passes_count_gate_both") &
      safe_num(n_case_samples) >= 2 &
      safe_num(n_control_samples) >= 2 &
      safe_num(min_feature_raw_counts) >= 30]
cds[, supported_support := strict_support |
      (bool_col(.SD, "matched_feature_passes_count_gate_any") &
         safe_num(min_matched_samples) >= 2 &
         safe_num(max_feature_raw_counts) >= 30)]
cds[, row_weight := pmax(
  1,
  log1p(pmin(safe_num(case_sum_raw_counts),
             safe_num(control_sum_raw_counts),
             na.rm = TRUE))
)]
if ("matched_bootstrap_certainty_score" %in% names(cds)) {
  cds[, row_weight := row_weight *
        pmax(0.25, safe_num(matched_bootstrap_certainty_score))]
}
cds[!is.finite(row_weight), row_weight := 1]

model_rows <- cds[is.finite(effect_log2fc) & is.finite(effect_se)]
if (!nrow(model_rows)) stop("No matched clean-CDS rows were modelable.")

gene_design <- summarise_groups(
  model_rows,
  c("gene_symbol", "tx_id", "design_family")
)
gene_prior <- summarise_groups(model_rows, c("gene_symbol", "tx_id"))
setnames(
  gene_prior,
  setdiff(names(gene_prior), c("gene_symbol", "tx_id")),
  paste0("gene_prior_", setdiff(names(gene_prior), c("gene_symbol", "tx_id")))
)
design_prior <- summarise_groups(model_rows, c("design_family"))
setnames(
  design_prior,
  setdiff(names(design_prior), "design_family"),
  paste0("design_prior_", setdiff(names(design_prior), "design_family"))
)
global_re <- random_effects_summary(
  model_rows$effect_log2fc,
  model_rows$effect_se,
  model_rows$row_weight
)

gene_design <- merge(gene_design, gene_prior,
                     by = c("gene_symbol", "tx_id"), all.x = TRUE)
gene_design <- merge(gene_design, design_prior,
                     by = "design_family", all.x = TRUE)

posterior <- rbindlist(lapply(seq_len(nrow(gene_design)), function(i) {
  posterior_from_priors(
    direct_effect = gene_design$direct_log2fc[[i]],
    direct_se = gene_design$direct_se[[i]],
    gene_effect = gene_design$gene_prior_direct_log2fc[[i]],
    gene_se = gene_design$gene_prior_direct_se[[i]],
    design_effect = gene_design$design_prior_direct_log2fc[[i]],
    design_se = gene_design$design_prior_direct_se[[i]],
    global_effect = global_re$effect[[1]],
    global_se = global_re$se[[1]]
  )
}))
gene_design <- cbind(gene_design, posterior)

gene_design[, hierarchical_clean_cds_pct_change :=
              pct_from_log2(hierarchical_posterior_log2fc)]
gene_design[, hierarchical_clean_cds_pct_lower :=
              pct_from_log2(hierarchical_posterior_log2fc_lower)]
gene_design[, hierarchical_clean_cds_pct_upper :=
              pct_from_log2(hierarchical_posterior_log2fc_upper)]
gene_design[, hierarchical_ci_excludes_zero :=
              hierarchical_posterior_log2fc_lower > 0 |
              hierarchical_posterior_log2fc_upper < 0]
gene_design[, hierarchical_direction :=
              fifelse(hierarchical_posterior_log2fc > 0, "clean_CDS_up",
                      fifelse(hierarchical_posterior_log2fc < 0,
                              "clean_CDS_down", "flat_or_unknown"))]
gene_design[, hierarchical_direction_fraction :=
              fifelse(n_strict_rows >= 2 & is.finite(strict_random_effects_direction_fraction),
                      strict_random_effects_direction_fraction,
                      direct_direction_fraction)]
gene_design[, hierarchical_direction_stable :=
              is.finite(hierarchical_direction_fraction) &
              hierarchical_direction_fraction >= 0.75]
gene_design[, hierarchical_prior_source := fifelse(
  hierarchical_evidence_weight >= 0.75, "direct_dominant",
  fifelse(hierarchical_evidence_weight >= 0.35,
          "direct_plus_gene_design_prior", "gene_design_prior_dominant")
)]

if (nrow(branch_usage)) {
  branch_keep <- intersect(names(branch_usage), c(
    "gene_symbol", "tx_id", "design_family",
    "branch_usage_top_feature_id",
    "branch_usage_top_feature_class",
    "branch_usage_top_support_class",
    "branch_usage_top_score",
    "branch_usage_top_weighted_usage_delta",
    "branch_usage_top_strict_studies"
  ))
  gene_design <- merge(
    gene_design,
    branch_usage[, ..branch_keep],
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE
  )
}
if (!"branch_usage_top_score" %in% names(gene_design)) {
  gene_design[, branch_usage_top_score := NA_real_]
}

if (nrow(count_aware_branch_allocation)) {
  count_aware_keep <- intersect(names(count_aware_branch_allocation), c(
    "gene_symbol", "tx_id", "design_family",
    "count_aware_top_feature_id",
    "count_aware_top_feature_class",
    "count_aware_top_feature_type",
    "count_aware_top_support_class",
    "count_aware_top_score",
    "count_aware_top_weighted_usage_delta",
    "count_aware_top_random_effects_delta",
    "count_aware_top_delta_lower",
    "count_aware_top_delta_upper",
    "count_aware_top_q",
    "count_aware_top_i2",
    "count_aware_top_direction_fraction",
    "count_aware_top_strict_rows",
    "count_aware_top_strict_studies",
    "count_aware_top_evaluable_rows",
    "count_aware_top_context",
    "count_aware_top_control",
    "count_aware_top_study",
    "feature_id_clean_CDS",
    "feature_id_leader_uORF",
    "feature_id_overlapping_uORF",
    "count_aware_weighted_usage_delta_clean_CDS",
    "count_aware_weighted_usage_delta_leader_uORF",
    "count_aware_weighted_usage_delta_overlapping_uORF",
    "count_aware_random_effects_delta_clean_CDS",
    "count_aware_random_effects_delta_leader_uORF",
    "count_aware_random_effects_delta_overlapping_uORF",
    "count_aware_random_effects_delta_q_clean_CDS",
    "count_aware_random_effects_delta_q_leader_uORF",
    "count_aware_random_effects_delta_q_overlapping_uORF",
    "count_aware_branch_allocation_score_clean_CDS",
    "count_aware_branch_allocation_score_leader_uORF",
    "count_aware_branch_allocation_score_overlapping_uORF"
  ))
  count_aware_small <- count_aware_branch_allocation[, ..count_aware_keep]
  feature_id_cols <- intersect(
    c("feature_id_clean_CDS", "feature_id_leader_uORF",
      "feature_id_overlapping_uORF"),
    names(count_aware_small)
  )
  if (length(feature_id_cols)) {
    setnames(count_aware_small, feature_id_cols,
             paste0("count_aware_", feature_id_cols))
  }
  gene_design <- merge(
    gene_design,
    count_aware_small,
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE
  )
}
if (!"count_aware_top_score" %in% names(gene_design)) {
  gene_design[, count_aware_top_score := NA_real_]
}
if (!"count_aware_top_support_class" %in% names(gene_design)) {
  gene_design[, count_aware_top_support_class := NA_character_]
}

if (nrow(postviral)) {
  post_keep <- intersect(names(postviral), c(
    "gene_symbol", "tx_id", "design_family",
    "matched_gene_design_score",
    "gene_design_evidence_class",
    "strict_random_effects_interval_supported",
    "loo_effect_direction_stable",
    "study_loo_effect_direction_stable"
  ))
  post_small <- postviral[, ..post_keep]
  setnames(
    post_small,
    setdiff(names(post_small), c("gene_symbol", "tx_id", "design_family")),
    paste0("postviral_", setdiff(names(post_small),
                                 c("gene_symbol", "tx_id", "design_family")))
  )
  gene_design <- merge(
    gene_design, post_small,
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = TRUE
  )
}
if (!"postviral_matched_gene_design_score" %in% names(gene_design)) {
  gene_design[, postviral_matched_gene_design_score := NA_real_]
}

if (nrow(candidate)) {
  candidate_keep <- intersect(names(candidate), c(
    "gene_symbol", "tx_id", "candidate_rank_v4", "novelty_class",
    "known_uorf_regulation", "manual_m_source_gene", "hidden_uorf_candidate",
    "has_overlapping_uorf", "n_uorfs", "n_overlapping_uorfs",
    "n_uorfs_with_m_support", "next_iteration_candidate_score_v4"
  ))
  gene_design <- merge(
    gene_design,
    unique(candidate[, ..candidate_keep]),
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE
  )
}

gene_design[, hierarchical_support_component :=
              0.35 * pmin(n_strict_studies / 2, 1) +
              0.25 * pmin(n_studies / 3, 1) +
              0.20 * pmin(n_strict_rows / 4, 1) +
              0.20 * pmin(n_supported_rows / 8, 1)]
gene_design[, hierarchical_effect_component :=
              pmin(abs(hierarchical_clean_cds_pct_change) / 100, 1)]
gene_design[, hierarchical_branch_component := fifelse(
  is.finite(branch_usage_top_score), pmin(branch_usage_top_score, 1), 0
)]
gene_design[, hierarchical_count_aware_branch_component := fifelse(
  grepl("count_aware_branch_shift", count_aware_top_support_class,
        fixed = TRUE) &
    is.finite(count_aware_top_score),
  pmin(count_aware_top_score, 1),
  0
)]
gene_design[!is.finite(hierarchical_count_aware_branch_component),
            hierarchical_count_aware_branch_component := 0]
gene_design[, hierarchical_mechanistic_support_class := fifelse(
  count_aware_top_support_class == "replicated_count_aware_branch_shift",
  "count_aware_replicated_branch_support",
  fifelse(
    grepl("count_aware_branch_shift", count_aware_top_support_class,
          fixed = TRUE),
    "count_aware_branch_support",
    fifelse(
      grepl("branch_usage_shift", branch_usage_top_support_class,
            fixed = TRUE),
      "branch_usage_support",
      "no_mechanistic_branch_support"
    )
  )
)]
gene_design[
  is.na(hierarchical_mechanistic_support_class) |
    !nzchar(hierarchical_mechanistic_support_class),
  hierarchical_mechanistic_support_class := "no_mechanistic_branch_support"
]
gene_design[, hierarchical_postviral_component := fifelse(
  is.finite(postviral_matched_gene_design_score),
  pmin(postviral_matched_gene_design_score / 1.5, 1),
  0
)]
gene_design[, hierarchical_interval_component :=
              as.numeric(hierarchical_ci_excludes_zero == TRUE)]
gene_design[, hierarchical_stability_component :=
              as.numeric(hierarchical_direction_stable == TRUE)]
gene_design[, hierarchical_priority_score := pmin(
  1.5,
  0.29 * hierarchical_support_component +
    0.22 * hierarchical_interval_component +
    0.17 * pmin(hierarchical_evidence_weight, 1) +
    0.12 * hierarchical_stability_component +
    0.10 * hierarchical_effect_component +
    0.04 * hierarchical_branch_component +
    0.04 * hierarchical_count_aware_branch_component +
    0.03 * hierarchical_postviral_component
)]
gene_design[, hierarchical_evidence_class := fifelse(
  n_strict_studies >= 2 & n_strict_rows >= 2 &
    hierarchical_ci_excludes_zero == TRUE &
    hierarchical_direction_stable == TRUE,
  "hierarchical_atlas_ready",
  fifelse(n_studies >= 2 & hierarchical_ci_excludes_zero == TRUE,
          "hierarchical_replicated_interval_review",
          fifelse(hierarchical_ci_excludes_zero == TRUE,
                  "hierarchical_single_design_interval_review",
                  fifelse(n_studies >= 2,
                          "hierarchical_replicated_no_interval",
                          "hierarchical_low_support"))))]
gene_design[, abs_hierarchical_clean_cds_pct_change :=
              abs(hierarchical_clean_cds_pct_change)]

setorder(gene_design, -hierarchical_priority_score,
         -abs_hierarchical_clean_cds_pct_change, gene_symbol, design_family)

loo_rows <- list()
loo_i <- 0L
eligible <- gene_design[n_studies >= 2, .(gene_symbol, tx_id, design_family)]
if (nrow(eligible)) {
  for (i in seq_len(nrow(eligible))) {
    key <- eligible[i]
    rows <- model_rows[
      gene_symbol == key$gene_symbol &
        tx_id == key$tx_id &
        design_family == key$design_family
    ]
    studies <- unique(rows$study[!is.na(rows$study) & nzchar(rows$study)])
    if (length(studies) < 2) next
    for (study_id in studies) {
      train <- rows[study != study_id]
      test <- rows[study == study_id]
      if (!nrow(train) || !nrow(test)) next
      train_re <- random_effects_summary(
        train$effect_log2fc, train$effect_se, train$row_weight
      )
      test_re <- random_effects_summary(
        test$effect_log2fc, test$effect_se, test$row_weight
      )
      gd <- gene_design[
        gene_symbol == key$gene_symbol &
          tx_id == key$tx_id &
          design_family == key$design_family
      ][1]
      pred <- posterior_from_priors(
        direct_effect = train_re$effect[[1]],
        direct_se = train_re$se[[1]],
        gene_effect = gd$gene_prior_direct_log2fc[[1]],
        gene_se = gd$gene_prior_direct_se[[1]],
        design_effect = gd$design_prior_direct_log2fc[[1]],
        design_se = gd$design_prior_direct_se[[1]],
        global_effect = global_re$effect[[1]],
        global_se = global_re$se[[1]]
      )
      loo_i <- loo_i + 1L
      loo_rows[[loo_i]] <- data.table(
        gene_symbol = key$gene_symbol,
        tx_id = key$tx_id,
        design_family = key$design_family,
        heldout_study = study_id,
        train_studies = uniqueN(train$study),
        heldout_rows = nrow(test),
        heldout_log2fc = test_re$effect,
        predicted_log2fc = pred$hierarchical_posterior_log2fc,
        predicted_log2fc_lower = pred$hierarchical_posterior_log2fc -
          1.96 * pred$hierarchical_posterior_se,
        predicted_log2fc_upper = pred$hierarchical_posterior_log2fc +
          1.96 * pred$hierarchical_posterior_se,
        prediction_error_log2fc =
          test_re$effect - pred$hierarchical_posterior_log2fc,
        prediction_abs_error_log2fc =
          abs(test_re$effect - pred$hierarchical_posterior_log2fc),
        prediction_interval_contains_heldout =
          test_re$effect >= pred$hierarchical_posterior_log2fc -
            1.96 * pred$hierarchical_posterior_se &
          test_re$effect <= pred$hierarchical_posterior_log2fc +
            1.96 * pred$hierarchical_posterior_se
      )
    }
  }
}
loo <- if (length(loo_rows)) rbindlist(loo_rows, fill = TRUE) else data.table()

if (nrow(loo)) {
  loo[, prediction_interval_half_width_log2fc :=
        pmax(abs(predicted_log2fc_upper - predicted_log2fc) ,
             abs(predicted_log2fc - predicted_log2fc_lower),
             na.rm = TRUE)]
  loo[!is.finite(prediction_interval_half_width_log2fc),
      prediction_interval_half_width_log2fc := 0]
  global_q90 <- quantile_safe(loo$prediction_abs_error_log2fc, 0.90,
                              default = 0.5)
  global_q95 <- quantile_safe(loo$prediction_abs_error_log2fc, 0.95,
                              default = 0.75)
  design_calibration <- loo[, .(
    hierarchical_calibration_predictions = .N,
    hierarchical_calibration_mae_log2fc =
      mean(prediction_abs_error_log2fc, na.rm = TRUE),
    hierarchical_conformal_q90_log2fc =
      quantile_safe(prediction_abs_error_log2fc, 0.90,
                    default = global_q90),
    hierarchical_conformal_q95_log2fc =
      quantile_safe(prediction_abs_error_log2fc, 0.95,
                    default = global_q95)
  ), by = design_family]
  design_calibration[
    hierarchical_calibration_predictions < 50 |
      !is.finite(hierarchical_conformal_q95_log2fc),
    `:=`(
      hierarchical_conformal_q90_log2fc = global_q90,
      hierarchical_conformal_q95_log2fc = global_q95
    )
  ]
  design_calibration[, hierarchical_calibration_scope := fifelse(
    hierarchical_calibration_predictions >= 50,
    "design_family_residual_quantile",
    "global_residual_quantile"
  )]
  loo <- merge(
    loo,
    design_calibration[, .(
      design_family,
      hierarchical_calibration_predictions,
      hierarchical_conformal_q90_log2fc,
      hierarchical_conformal_q95_log2fc,
      hierarchical_calibration_scope
    )],
    by = "design_family",
    all.x = TRUE
  )
  loo[!is.finite(hierarchical_conformal_q95_log2fc),
      `:=`(
        hierarchical_conformal_q90_log2fc = global_q90,
        hierarchical_conformal_q95_log2fc = global_q95,
        hierarchical_calibration_scope = "global_residual_quantile"
      )]
  loo[, calibrated_prediction_half_width_log2fc :=
        pmax(prediction_interval_half_width_log2fc,
             hierarchical_conformal_q95_log2fc,
             na.rm = TRUE)]
  loo[, calibrated_predicted_log2fc_lower :=
        predicted_log2fc - calibrated_prediction_half_width_log2fc]
  loo[, calibrated_predicted_log2fc_upper :=
        predicted_log2fc + calibrated_prediction_half_width_log2fc]
  loo[, calibrated_prediction_interval_contains_heldout :=
        heldout_log2fc >= calibrated_predicted_log2fc_lower &
        heldout_log2fc <= calibrated_predicted_log2fc_upper]

  loo_gene_design <- loo[, .(
    hierarchical_loo_rows = .N,
    hierarchical_loo_mae_log2fc = mean(prediction_abs_error_log2fc,
                                       na.rm = TRUE),
    hierarchical_loo_rmse_log2fc = sqrt(mean(prediction_error_log2fc^2,
                                             na.rm = TRUE)),
    hierarchical_loo_interval_coverage =
      mean(prediction_interval_contains_heldout, na.rm = TRUE),
    hierarchical_loo_calibrated_interval_coverage =
      mean(calibrated_prediction_interval_contains_heldout, na.rm = TRUE),
    hierarchical_loo_direction_accuracy = mean(
      sign(heldout_log2fc) == sign(predicted_log2fc),
      na.rm = TRUE
    )
  ), by = .(gene_symbol, tx_id, design_family)]
  gene_design <- merge(
    gene_design, loo_gene_design,
    by = c("gene_symbol", "tx_id", "design_family"), all.x = TRUE
  )
  gene_design <- merge(
    gene_design,
    design_calibration,
    by = "design_family",
    all.x = TRUE
  )
} else {
  gene_design[, `:=`(
    hierarchical_loo_rows = 0L,
    hierarchical_loo_mae_log2fc = NA_real_,
    hierarchical_loo_rmse_log2fc = NA_real_,
    hierarchical_loo_interval_coverage = NA_real_,
    hierarchical_loo_calibrated_interval_coverage = NA_real_,
    hierarchical_loo_direction_accuracy = NA_real_
  )]
  gene_design[, `:=`(
    hierarchical_calibration_predictions = 0L,
    hierarchical_calibration_mae_log2fc = NA_real_,
    hierarchical_conformal_q90_log2fc = NA_real_,
    hierarchical_conformal_q95_log2fc = NA_real_,
    hierarchical_calibration_scope = "no_loo_calibration"
  )]
}

if (!"hierarchical_conformal_q95_log2fc" %in% names(gene_design)) {
  gene_design[, `:=`(
    hierarchical_calibration_predictions = 0L,
    hierarchical_calibration_mae_log2fc = NA_real_,
    hierarchical_conformal_q90_log2fc = NA_real_,
    hierarchical_conformal_q95_log2fc = NA_real_,
    hierarchical_calibration_scope = "no_loo_calibration"
  )]
}
if (nrow(loo)) {
  gene_design[
    !is.finite(hierarchical_conformal_q95_log2fc),
    `:=`(
      hierarchical_calibration_predictions = nrow(loo),
      hierarchical_calibration_mae_log2fc =
        mean(loo$prediction_abs_error_log2fc, na.rm = TRUE),
      hierarchical_conformal_q90_log2fc = global_q90,
      hierarchical_conformal_q95_log2fc = global_q95,
      hierarchical_calibration_scope =
        "global_residual_quantile_no_family_loo"
    )
  ]
}
gene_design[, hierarchical_uncalibrated_evidence_class :=
              hierarchical_evidence_class]
gene_design[, hierarchical_uncalibrated_priority_score :=
              hierarchical_priority_score]
gene_design[, hierarchical_uncalibrated_ci_excludes_zero :=
              hierarchical_ci_excludes_zero]
gene_design[, hierarchical_normal_half_width_log2fc :=
              pmax(abs(hierarchical_posterior_log2fc_upper -
                         hierarchical_posterior_log2fc),
                   abs(hierarchical_posterior_log2fc -
                         hierarchical_posterior_log2fc_lower),
                   na.rm = TRUE)]
gene_design[!is.finite(hierarchical_normal_half_width_log2fc),
            hierarchical_normal_half_width_log2fc :=
              1.96 * hierarchical_posterior_se]
gene_design[, hierarchical_calibrated_half_width_log2fc := fifelse(
  is.finite(hierarchical_conformal_q95_log2fc),
  pmax(hierarchical_normal_half_width_log2fc,
       hierarchical_conformal_q95_log2fc,
       na.rm = TRUE),
  hierarchical_normal_half_width_log2fc
)]
gene_design[!is.finite(hierarchical_calibrated_half_width_log2fc),
            hierarchical_calibrated_half_width_log2fc :=
              hierarchical_normal_half_width_log2fc]
gene_design[, hierarchical_calibrated_posterior_log2fc_lower :=
              hierarchical_posterior_log2fc -
              hierarchical_calibrated_half_width_log2fc]
gene_design[, hierarchical_calibrated_posterior_log2fc_upper :=
              hierarchical_posterior_log2fc +
              hierarchical_calibrated_half_width_log2fc]
gene_design[, hierarchical_calibrated_clean_cds_pct_lower :=
              pct_from_log2(hierarchical_calibrated_posterior_log2fc_lower)]
gene_design[, hierarchical_calibrated_clean_cds_pct_upper :=
              pct_from_log2(hierarchical_calibrated_posterior_log2fc_upper)]
gene_design[, hierarchical_calibrated_ci_excludes_zero :=
              hierarchical_calibrated_posterior_log2fc_lower > 0 |
              hierarchical_calibrated_posterior_log2fc_upper < 0]
gene_design[, hierarchical_interval_component :=
              as.numeric(hierarchical_calibrated_ci_excludes_zero == TRUE)]
gene_design[, hierarchical_priority_score := pmin(
  1.5,
  0.29 * hierarchical_support_component +
    0.22 * hierarchical_interval_component +
    0.17 * pmin(hierarchical_evidence_weight, 1) +
    0.12 * hierarchical_stability_component +
    0.10 * hierarchical_effect_component +
    0.04 * hierarchical_branch_component +
    0.04 * hierarchical_count_aware_branch_component +
    0.03 * hierarchical_postviral_component
)]
gene_design[, hierarchical_evidence_class := fifelse(
  n_strict_studies >= 2 & n_strict_rows >= 2 &
    hierarchical_calibrated_ci_excludes_zero == TRUE &
    hierarchical_direction_stable == TRUE,
  "hierarchical_calibrated_atlas_ready",
  fifelse(n_studies >= 2 & hierarchical_calibrated_ci_excludes_zero == TRUE,
          "hierarchical_calibrated_replicated_interval_review",
          fifelse(hierarchical_calibrated_ci_excludes_zero == TRUE,
                  "hierarchical_calibrated_single_design_interval_review",
                  fifelse(n_studies >= 2,
                          "hierarchical_replicated_no_calibrated_interval",
                          "hierarchical_low_support"))))]
gene_design[, abs_hierarchical_clean_cds_pct_change :=
              abs(hierarchical_clean_cds_pct_change)]
setorder(gene_design, -hierarchical_priority_score,
         -abs_hierarchical_clean_cds_pct_change, gene_symbol, design_family)

genes <- unique(model_rows[, .(gene_symbol, tx_id)])
families <- sort(unique(model_rows$design_family))
counterfactual <- CJ(
  gene_symbol = genes$gene_symbol,
  design_family = families,
  unique = TRUE
)
counterfactual <- merge(counterfactual, genes, by = "gene_symbol", all.x = TRUE,
                        allow.cartesian = TRUE)
observed_keys <- gene_design[, .(gene_symbol, tx_id, design_family)]
counterfactual <- fsetdiff(
  unique(counterfactual[, .(gene_symbol, tx_id, design_family)]),
  observed_keys
)
if (nrow(counterfactual)) {
  counterfactual <- merge(counterfactual, gene_prior,
                          by = c("gene_symbol", "tx_id"), all.x = TRUE)
  counterfactual <- merge(counterfactual, design_prior,
                          by = "design_family", all.x = TRUE)
  counter_posterior <- rbindlist(lapply(seq_len(nrow(counterfactual)), function(i) {
    posterior_from_priors(
      direct_effect = NA_real_,
      direct_se = NA_real_,
      gene_effect = counterfactual$gene_prior_direct_log2fc[[i]],
      gene_se = counterfactual$gene_prior_direct_se[[i]],
      design_effect = counterfactual$design_prior_direct_log2fc[[i]],
      design_se = counterfactual$design_prior_direct_se[[i]],
      global_effect = global_re$effect[[1]],
      global_se = global_re$se[[1]]
    )
  }))
  counterfactual <- cbind(counterfactual, counter_posterior)
  counterfactual[, hierarchical_clean_cds_pct_change :=
                   pct_from_log2(hierarchical_posterior_log2fc)]
  counterfactual[, hierarchical_clean_cds_pct_lower :=
                   pct_from_log2(hierarchical_posterior_log2fc_lower)]
  counterfactual[, hierarchical_clean_cds_pct_upper :=
                   pct_from_log2(hierarchical_posterior_log2fc_upper)]
  counterfactual[, counterfactual_prediction_class := fifelse(
    gene_prior_n_studies >= 2 & design_prior_n_studies >= 2,
    "gene_and_design_prior",
    fifelse(gene_prior_n_studies >= 2, "gene_prior_only",
            fifelse(design_prior_n_studies >= 2, "design_prior_only",
                    "global_prior_only"))
  )]
  counterfactual[, abs_hierarchical_clean_cds_pct_change :=
                   abs(hierarchical_clean_cds_pct_change)]
  setorder(counterfactual, -abs_hierarchical_clean_cds_pct_change,
           gene_symbol, design_family)
}

preferred_gene_design <- c(
  "gene_symbol", "tx_id", "design_family", "hierarchical_evidence_class",
  "hierarchical_uncalibrated_evidence_class",
  "hierarchical_priority_score", "hierarchical_direction",
  "hierarchical_uncalibrated_priority_score",
  "hierarchical_clean_cds_pct_change", "hierarchical_clean_cds_pct_lower",
  "hierarchical_clean_cds_pct_upper", "hierarchical_posterior_log2fc",
  "hierarchical_posterior_log2fc_lower",
  "hierarchical_posterior_log2fc_upper",
  "hierarchical_posterior_se", "hierarchical_ci_excludes_zero",
  "hierarchical_uncalibrated_ci_excludes_zero",
  "hierarchical_calibrated_posterior_log2fc_lower",
  "hierarchical_calibrated_posterior_log2fc_upper",
  "hierarchical_calibrated_clean_cds_pct_lower",
  "hierarchical_calibrated_clean_cds_pct_upper",
  "hierarchical_calibrated_ci_excludes_zero",
  "hierarchical_calibrated_half_width_log2fc",
  "hierarchical_normal_half_width_log2fc",
  "hierarchical_conformal_q90_log2fc",
  "hierarchical_conformal_q95_log2fc",
  "hierarchical_calibration_scope",
  "hierarchical_calibration_predictions",
  "hierarchical_calibration_mae_log2fc",
  "hierarchical_direction_fraction", "hierarchical_direction_stable",
  "hierarchical_evidence_weight", "hierarchical_prior_source",
  "preferred_direct_source", "direct_log2fc", "direct_se",
  "direct_log2fc_lower", "direct_log2fc_upper", "direct_i2",
  "n_rows", "n_supported_rows", "n_strict_rows", "n_studies",
  "n_strict_studies", "n_designs", "n_cell_lines", "n_tissues",
  "total_case_samples", "total_control_samples", "total_case_counts",
  "total_control_counts", "fraction_positive_rows", "strongest_condition",
  "strongest_control", "strongest_study", "nearest_designs",
  "gene_prior_direct_log2fc", "gene_prior_direct_se",
  "design_prior_direct_log2fc", "design_prior_direct_se",
  "branch_usage_top_feature_id", "branch_usage_top_feature_class",
  "branch_usage_top_support_class", "branch_usage_top_score",
  "count_aware_top_feature_id", "count_aware_top_feature_class",
  "count_aware_top_support_class", "count_aware_top_score",
  "count_aware_top_weighted_usage_delta",
  "count_aware_top_random_effects_delta",
  "count_aware_top_delta_lower", "count_aware_top_delta_upper",
  "count_aware_top_q", "count_aware_top_i2",
  "count_aware_top_direction_fraction",
  "count_aware_top_strict_rows", "count_aware_top_strict_studies",
  "count_aware_top_evaluable_rows", "count_aware_top_context",
  "count_aware_top_study",
  "hierarchical_count_aware_branch_component",
  "hierarchical_mechanistic_support_class",
  "postviral_matched_gene_design_score",
  "postviral_gene_design_evidence_class",
  "candidate_rank_v4", "novelty_class", "known_uorf_regulation",
  "manual_m_source_gene", "hidden_uorf_candidate", "has_overlapping_uorf",
  "n_uorfs", "n_overlapping_uorfs", "n_uorfs_with_m_support",
  "hierarchical_loo_rows", "hierarchical_loo_mae_log2fc",
  "hierarchical_loo_rmse_log2fc", "hierarchical_loo_interval_coverage",
  "hierarchical_loo_calibrated_interval_coverage",
  "hierarchical_loo_direction_accuracy"
)
gene_design <- gene_design[
  ,
  c(intersect(preferred_gene_design, names(gene_design)),
    setdiff(names(gene_design), preferred_gene_design)),
  with = FALSE
]

gene_design_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_gene_design_effects.csv"
)
counterfactual_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_counterfactual_predictions.csv"
)
loo_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_leave_one_study_calibration.csv"
)
loo_summary_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_calibration_summary.csv"
)
interval_calibration_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_interval_calibration.csv"
)
summary_file <- file.path(
  output_dir, "dominant_rdg_hierarchical_summary_metrics.csv"
)

fwrite(gene_design, gene_design_file)
fwrite(counterfactual, counterfactual_file)
fwrite(loo, loo_file)

if (nrow(loo) && !"calibrated_prediction_interval_contains_heldout" %in%
    names(loo)) {
  loo[, calibrated_prediction_interval_contains_heldout := NA]
}
loo_summary <- if (nrow(loo)) {
  rbindlist(list(
    loo[, .(
      scope = "all",
      n_predictions = .N,
      mae_log2fc = mean(prediction_abs_error_log2fc, na.rm = TRUE),
      rmse_log2fc = sqrt(mean(prediction_error_log2fc^2, na.rm = TRUE)),
      interval_coverage = mean(prediction_interval_contains_heldout,
                               na.rm = TRUE),
      calibrated_interval_coverage =
        mean(calibrated_prediction_interval_contains_heldout, na.rm = TRUE),
      conformal_q90_log2fc =
        quantile_safe(prediction_abs_error_log2fc, 0.90),
      conformal_q95_log2fc =
        quantile_safe(prediction_abs_error_log2fc, 0.95),
      direction_accuracy = mean(sign(heldout_log2fc) ==
                                  sign(predicted_log2fc), na.rm = TRUE)
    )],
    loo[, .(
      scope = paste0("design_family:", design_family),
      n_predictions = .N,
      mae_log2fc = mean(prediction_abs_error_log2fc, na.rm = TRUE),
      rmse_log2fc = sqrt(mean(prediction_error_log2fc^2, na.rm = TRUE)),
      interval_coverage = mean(prediction_interval_contains_heldout,
                               na.rm = TRUE),
      calibrated_interval_coverage =
        mean(calibrated_prediction_interval_contains_heldout, na.rm = TRUE),
      conformal_q90_log2fc =
        quantile_safe(prediction_abs_error_log2fc, 0.90),
      conformal_q95_log2fc =
        quantile_safe(prediction_abs_error_log2fc, 0.95),
      direction_accuracy = mean(sign(heldout_log2fc) ==
                                  sign(predicted_log2fc), na.rm = TRUE)
    ), by = design_family][, design_family := NULL]
  ), fill = TRUE)
} else {
  data.table(
    scope = "all",
    n_predictions = 0L,
    mae_log2fc = NA_real_,
    rmse_log2fc = NA_real_,
    interval_coverage = NA_real_,
    calibrated_interval_coverage = NA_real_,
    conformal_q90_log2fc = NA_real_,
    conformal_q95_log2fc = NA_real_,
    direction_accuracy = NA_real_
  )
}
fwrite(loo_summary, loo_summary_file)

interval_calibration <- if (nrow(gene_design)) {
  unique(gene_design[, .(
    design_family,
    hierarchical_calibration_scope,
    hierarchical_calibration_predictions,
    hierarchical_calibration_mae_log2fc,
    hierarchical_conformal_q90_log2fc,
    hierarchical_conformal_q95_log2fc
  )])
} else {
  data.table()
}
if (nrow(loo_summary)) {
  interval_calibration <- rbindlist(list(
    data.table(
      design_family = "all",
      hierarchical_calibration_scope = "global_residual_quantile",
      hierarchical_calibration_predictions =
        loo_summary[scope == "all", n_predictions][1],
      hierarchical_calibration_mae_log2fc =
        loo_summary[scope == "all", mae_log2fc][1],
      hierarchical_conformal_q90_log2fc =
        loo_summary[scope == "all", conformal_q90_log2fc][1],
      hierarchical_conformal_q95_log2fc =
        loo_summary[scope == "all", conformal_q95_log2fc][1]
    ),
    interval_calibration
  ), fill = TRUE)
}
fwrite(interval_calibration, interval_calibration_file)

summary_metrics <- data.table(
  metric = c(
    "model_rows",
    "gene_design_rows",
    "genes",
    "design_families",
    "hierarchical_atlas_ready_rows",
    "hierarchical_replicated_interval_review_rows",
    "hierarchical_interval_rows",
    "hierarchical_uncalibrated_atlas_ready_rows",
    "hierarchical_uncalibrated_interval_rows",
    "hierarchical_count_aware_supported_rows",
    "hierarchical_count_aware_replicated_rows",
    "viral_hierarchical_count_aware_supported_rows",
    "viral_hierarchical_count_aware_replicated_rows",
    "median_hierarchical_priority_score",
    "loo_predictions",
    "loo_mae_log2fc",
    "loo_rmse_log2fc",
    "loo_interval_coverage",
    "loo_calibrated_interval_coverage",
    "loo_conformal_q95_log2fc",
    "counterfactual_predictions"
  ),
  value = c(
    nrow(model_rows),
    nrow(gene_design),
    uniqueN(gene_design$gene_symbol),
    uniqueN(gene_design$design_family),
    sum(gene_design$hierarchical_evidence_class ==
          "hierarchical_calibrated_atlas_ready", na.rm = TRUE),
    sum(gene_design$hierarchical_evidence_class ==
          "hierarchical_calibrated_replicated_interval_review", na.rm = TRUE),
    sum(gene_design$hierarchical_calibrated_ci_excludes_zero == TRUE,
        na.rm = TRUE),
    sum(gene_design$hierarchical_uncalibrated_evidence_class ==
          "hierarchical_atlas_ready", na.rm = TRUE),
    sum(gene_design$hierarchical_uncalibrated_ci_excludes_zero == TRUE,
        na.rm = TRUE),
    sum(grepl("count_aware_branch_shift",
              gene_design$count_aware_top_support_class,
              fixed = TRUE), na.rm = TRUE),
    sum(gene_design$count_aware_top_support_class ==
          "replicated_count_aware_branch_shift", na.rm = TRUE),
    sum(gene_design$design_family == "viral_infection" &
          grepl("count_aware_branch_shift",
                gene_design$count_aware_top_support_class,
                fixed = TRUE), na.rm = TRUE),
    sum(gene_design$design_family == "viral_infection" &
          gene_design$count_aware_top_support_class ==
          "replicated_count_aware_branch_shift", na.rm = TRUE),
    median_safe(gene_design$hierarchical_priority_score),
    if (nrow(loo)) nrow(loo) else 0,
    loo_summary[scope == "all", mae_log2fc][1],
    loo_summary[scope == "all", rmse_log2fc][1],
    loo_summary[scope == "all", interval_coverage][1],
    loo_summary[scope == "all", calibrated_interval_coverage][1],
    loo_summary[scope == "all", conformal_q95_log2fc][1],
    if (exists("counterfactual")) nrow(counterfactual) else 0
  )
)
fwrite(summary_metrics, summary_file)

plot_dt <- gene_design[
  hierarchical_priority_score >= 0.25 |
    hierarchical_evidence_class != "hierarchical_low_support"
]
if (nrow(plot_dt)) {
  keep_genes <- plot_dt[
    , .(score = max(hierarchical_priority_score *
                      pmin(abs(hierarchical_clean_cds_pct_change) / 100, 2),
                    na.rm = TRUE)),
    by = gene_symbol
  ][order(-score)][seq_len(min(.N, 40)), gene_symbol]
  plot_dt <- plot_dt[gene_symbol %chin% keep_genes]
  plot_dt[, gene_symbol := factor(gene_symbol, levels = rev(keep_genes))]
  family_order <- plot_dt[
    , .(score = max(hierarchical_priority_score, na.rm = TRUE)),
    by = design_family
  ][order(-score), design_family]
  plot_dt[, design_family := factor(design_family, levels = family_order)]
  plot_dt[, capped_pct := pmax(pmin(hierarchical_clean_cds_pct_change, 300), -100)]

  p_heat <- ggplot(
    plot_dt,
    aes(x = design_family, y = gene_symbol, fill = capped_pct,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Design: ", design_family,
          "<br>Hierarchical clean-CDS: ",
          round(hierarchical_clean_cds_pct_change, 1), "%",
          "<br>Calibrated interval: [",
          round(hierarchical_calibrated_clean_cds_pct_lower, 1),
          ", ", round(hierarchical_calibrated_clean_cds_pct_upper, 1), "]",
          "<br>Conformal q95 log2FC: ",
          round(hierarchical_conformal_q95_log2fc, 3),
          "<br>Class: ", hierarchical_evidence_class,
          "<br>Priority: ", round(hierarchical_priority_score, 2),
          "<br>Evidence weight: ", round(hierarchical_evidence_weight, 2),
          "<br>Mechanistic support: ",
          hierarchical_mechanistic_support_class,
          "<br>Count-aware support: ", count_aware_top_support_class
        ))
  ) +
    geom_tile(color = "grey88", linewidth = 0.15) +
    geom_point(aes(size = hierarchical_priority_score),
               shape = 21, fill = "white", color = "grey20",
               stroke = 0.2, alpha = 0.65) +
    geom_point(
      data = plot_dt[
        hierarchical_mechanistic_support_class %chin%
          c("count_aware_replicated_branch_support",
            "count_aware_branch_support")
      ],
      shape = 23,
      size = 2.1,
      fill = "#ead7f6",
      color = "#5b2a7a",
      stroke = 0.45,
      inherit.aes = TRUE
    ) +
    scale_fill_gradient2(
      low = "#2b6cb0", mid = "white", high = "#8b0000",
      midpoint = 0, limits = c(-100, 300), oob = scales::squish,
      name = "Clean-CDS\nchange (%)"
    ) +
    scale_size_continuous(range = c(0.4, 2.2), name = "Priority") +
    labs(
      title = "Hierarchical RDG Clean-CDS Effects",
      subtitle = "Purple diamond marks count-aware branch-allocation support",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  ggsave(file.path(figure_dir, "hierarchical_gene_design_heatmap.png"),
         p_heat, width = 12.5, height = 8.5, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_gene_design_heatmap.pdf"),
         p_heat, width = 12.5, height = 8.5)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_heat,
      file.path(figure_dir, "hierarchical_gene_design_heatmap.html"),
      title = "Hierarchical gene-design RDG effects"
    )
  }
}

if (nrow(loo)) {
  p_cal <- ggplot(
    loo,
    aes(x = predicted_log2fc, y = heldout_log2fc, color = design_family,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Design: ", design_family,
          "<br>Held-out study: ", heldout_study,
          "<br>Predicted log2FC: ", round(predicted_log2fc, 3),
          "<br>Held-out log2FC: ", round(heldout_log2fc, 3),
          "<br>Error: ", round(prediction_error_log2fc, 3)
        ))
  ) +
    geom_hline(yintercept = 0, color = "grey78", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "grey78", linewidth = 0.25) +
    geom_abline(slope = 1, intercept = 0, color = "grey35",
                linewidth = 0.35, linetype = "dashed") +
    geom_point(alpha = 0.62, size = 1.8) +
    labs(
      title = "Hierarchical Model Leave-One-Study-Out Calibration",
      subtitle = "Prediction from remaining studies versus held-out study effect",
      x = "Predicted clean-CDS log2 fold-change",
      y = "Held-out clean-CDS log2 fold-change"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "hierarchical_leave_one_study_calibration.png"),
         p_cal, width = 9.5, height = 6.4, dpi = 180)
  ggsave(file.path(figure_dir, "hierarchical_leave_one_study_calibration.pdf"),
         p_cal, width = 9.5, height = 6.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_cal,
      file.path(figure_dir, "hierarchical_leave_one_study_calibration.html"),
      title = "Hierarchical leave-one-study-out calibration"
    )
  }
}

message("Saved hierarchical outputs in: ", output_dir)
print(summary_metrics)
