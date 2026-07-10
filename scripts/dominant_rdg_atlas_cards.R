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
results_dir <- Sys.getenv(
  "DOMINANT_CELL_STATES_RESULTS_DIR",
  unset = file.path(analysis_dir, "results")
)
output_dir <- file.path(analysis_dir, "dominant_rdg_atlas")
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

collapse_unique <- function(x, n = 4L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

collapse_top <- function(x, n = 5L) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & !x %chin% c("MISSING", "NA", "none")]
  if (!length(x)) return("")
  tab <- sort(table(x), decreasing = TRUE)
  paste(names(tab)[seq_len(min(length(tab), n))], collapse = "; ")
}

norm_label <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[is.na(y)] <- ""
  y <- gsub("[^a-z0-9]+", "_", y)
  gsub("^_+|_+$", "", y)
}

control_label_like <- function(x) {
  y <- norm_label(x)
  grepl("(^|_)(wt|wildtype|wild_type|mock|control|ctrl|vehicle|dmso|untreated)($|_)",
        y)
}

protocol_bias_class <- function(inhibitor, fraction = "", condition = "") {
  text <- norm_label(paste(inhibitor, fraction, condition, sep = " "))
  fcase(
    grepl("(^|_)(ltm|lactimidomycin)($|_)", text),
    "start_site_enriched_ltm",
    grepl("harringtonine|harringtonine_treatment|(^|_)harr($|_)", text),
    "start_site_enriched_harringtonine",
    grepl("(^|_)(chx|cycloheximide)($|_)", text),
    "elongation_freeze_chx",
    grepl("(^|_)(puro|puromycin|emetine|anisomycin)($|_)", text),
    "translation_drug_stress",
    grepl("polysome|monosome|ribosome|fraction|sucrose", text),
    "fractionated_ribosome_pool",
    grepl("frozen|none|missing|untreated|dmso|vehicle", text),
    "no_specific_inhibitor_or_unknown",
    default = "other_protocol"
  )
}

count_informative_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & !x %chin% c("MISSING", "NA", "none")]
  uniqueN(x)
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

min_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

bool_col <- function(dt, col, default = FALSE) {
  if (!col %in% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[col]]
  if (is.logical(x)) return(fifelse(is.na(x), default, x))
  if (is.numeric(x)) return(fifelse(is.na(x), default, x != 0))
  y <- tolower(as.character(x))
  fifelse(is.na(y), default, y %chin% c("true", "t", "1", "yes"))
}

pct_from_log2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- (2^x - 1) * 100
  out[!is.finite(out)] <- NA_real_
  out
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
      lower = NA_real_,
      upper = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      q = NA_real_,
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
    lower = lower,
    upper = upper,
    tau2 = tau2,
    i2 = i2,
    q = q,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

effect_se_from_rows <- function(dt) {
  se <- rep(NA_real_, nrow(dt))
  if ("bootstrap_fpkm_log2_fc_sd" %in% names(dt)) {
    se <- suppressWarnings(as.numeric(dt$bootstrap_fpkm_log2_fc_sd))
  }
  if (all(!is.finite(se)) &&
      all(c("bootstrap_fpkm_log2_fc_lower",
            "bootstrap_fpkm_log2_fc_upper") %in% names(dt))) {
    lower <- suppressWarnings(as.numeric(dt$bootstrap_fpkm_log2_fc_lower))
    upper <- suppressWarnings(as.numeric(dt$bootstrap_fpkm_log2_fc_upper))
    se <- abs(upper - lower) / (2 * 1.96)
  }
  if ("matched_fpkm_log2_fc_pseudocount1" %in% names(dt)) {
    fallback <- 1 / sqrt(pmax(
      1,
      suppressWarnings(as.numeric(dt$case_sum_raw_counts)) +
        suppressWarnings(as.numeric(dt$control_sum_raw_counts))
    ))
    fallback <- pmax(fallback, 0.05)
    se[!is.finite(se) | se <= 0] <- fallback[!is.finite(se) | se <= 0]
  }
  se
}

source_file <- function(...) file.path(analysis_dir, ...)
first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) paths[1] else hit
}

cds <- read_dt(
  source_file("dominant_next_model",
              "dominant_next_model_matched_control_cds_buffering.csv"),
  required = TRUE
)
features <- read_dt(
  source_file("dominant_next_model",
              "dominant_next_model_matched_control_feature_contrasts.csv"),
  required = TRUE
)
postviral <- read_dt(
  source_file("postviral_fatigue_outputs",
              "postviral_fatigue_matched_gene_design_evidence.csv")
)
sparse <- read_dt(
  source_file("dominant_ouorf_cds_coupling",
              "dominant_ouorf_cds_coupling_sparse_onoff_candidates.csv")
)
feature_conf <- read_dt(
  source_file("dominant_next_model",
              "dominant_next_model_feature_confidence.csv")
)
candidate <- read_dt(
  source_file("dominant_next_model",
              "dominant_next_model_candidate_rankings_v4.csv")
)
browser <- read_dt(
  source_file("dominant_next_model",
              "dominant_next_model_browser_review_handoff.csv")
)
downstream <- read_dt(
  source_file("dominant_downstream_tis_rescue",
              "downstream_tis_rescue_best_candidates.csv")
)
branch_usage <- read_dt(
  source_file("dominant_rdg_branch_usage",
              "dominant_rdg_branch_usage_gene_design_summary.csv")
)
grouped_branch_usage <- read_dt(
  source_file("dominant_rdg_grouped_branch_usage",
              "dominant_rdg_grouped_branch_usage_gene_design_summary.csv")
)
count_aware_branch_allocation <- read_dt(
  source_file("dominant_rdg_branch_allocation",
              "dominant_rdg_branch_allocation_gene_design_summary.csv")
)
joint_branch_allocation <- read_dt(
  source_file("dominant_rdg_joint_branch_allocation",
              "dominant_rdg_joint_branch_allocation_gene_design_summary.csv")
)
joint_branch_rows <- read_dt(
  source_file("dominant_rdg_joint_branch_allocation",
              "dominant_rdg_joint_branch_allocation_rows.csv")
)
hierarchical_joint_allocation <- read_dt(
  source_file("dominant_rdg_hierarchical_joint_allocation",
              "dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv")
)
dirichlet_multinomial_allocation <- read_dt(
  source_file("dominant_rdg_dirichlet_multinomial",
              "dominant_rdg_dirichlet_multinomial_gene_design_summary.csv")
)
model_agreement <- read_dt(
  source_file("dominant_rdg_model_agreement",
              "dominant_rdg_model_agreement_gene_design.csv")
)
model_agreement_loo <- read_dt(
  source_file("dominant_rdg_model_agreement_loo",
              "dominant_rdg_model_agreement_loo_gene_design.csv")
)
model_agreement_transport <- read_dt(
  source_file("dominant_rdg_model_agreement_transport",
              "dominant_rdg_model_agreement_transport_gene_design.csv")
)
context_interactions <- read_dt(
  source_file("dominant_rdg_context_interactions",
              "dominant_rdg_context_interaction_gene_design_summary.csv")
)
context_partial_pooling <- read_dt(
  source_file("dominant_rdg_context_partial_pooling",
              "dominant_rdg_context_partial_pooling_gene_design_summary.csv")
)
hierarchical <- read_dt(
  source_file("dominant_rdg_hierarchical",
              "dominant_rdg_hierarchical_gene_design_effects.csv")
)
hierarchical_residual <- read_dt(
  source_file("dominant_rdg_hierarchical_residual_audit",
              "dominant_rdg_hierarchical_residual_audit_gene_design.csv")
)
hierarchical_context <- read_dt(
  source_file("dominant_rdg_hierarchical_context_calibration",
              "dominant_rdg_hierarchical_context_calibration_gene_design.csv")
)
hierarchical_context_model <- read_dt(
  source_file("dominant_rdg_hierarchical_context_model",
              "dominant_rdg_hierarchical_context_model_gene_design.csv")
)
browser_validation <- read_dt(
  first_existing_file(c(
    source_file("dominant_rdg_atlas", "browser_validation_notes.csv"),
    file.path(results_dir, "curated_inputs", "browser_validation_notes.csv")
  ))
)

if (!"design_family" %in% names(cds)) cds[, design_family := perturbation_class]
cds[is.na(design_family) | !nzchar(design_family), design_family := "unknown"]
cds <- cds[feature_id == "clean_CDS" | feature_type == "clean_CDS"]

cds[, clean_cds_log2fc := fifelse(
  is.finite(suppressWarnings(as.numeric(shrunk_matched_fpkm_log2_fc))),
  suppressWarnings(as.numeric(shrunk_matched_fpkm_log2_fc)),
  suppressWarnings(as.numeric(matched_fpkm_log2_fc_pseudocount1))
)]
cds[, clean_cds_pct_change := pct_from_log2(clean_cds_log2fc)]
cds[, clean_cds_se := effect_se_from_rows(.SD)]
cds[, strict_support := bool_col(.SD, "matched_feature_passes_sample_gate") &
      bool_col(.SD, "matched_feature_passes_count_gate_both") &
      suppressWarnings(as.numeric(n_case_samples)) >= 2 &
      suppressWarnings(as.numeric(n_control_samples)) >= 2 &
      suppressWarnings(as.numeric(min_feature_raw_counts)) >= 30]
cds[, supported_support := strict_support |
      (bool_col(.SD, "matched_feature_passes_count_gate_any") &
         suppressWarnings(as.numeric(min_matched_samples)) >= 2)]
cds[, row_weight := pmax(1, log1p(pmin(
  suppressWarnings(as.numeric(case_sum_raw_counts)),
  suppressWarnings(as.numeric(control_sum_raw_counts)),
  na.rm = TRUE
)))]
if ("matched_bootstrap_certainty_score" %in% names(cds)) {
  cds[, row_weight := row_weight *
        pmax(0.25, suppressWarnings(as.numeric(matched_bootstrap_certainty_score)))]
}
cds[!is.finite(row_weight), row_weight := 1]

cards <- cds[, {
  re_all <- random_effects_summary(clean_cds_log2fc, clean_cds_se, row_weight)
  re_strict <- random_effects_summary(
    clean_cds_log2fc[strict_support],
    clean_cds_se[strict_support],
    row_weight[strict_support]
  )
  use_re <- if (re_strict$n[[1]] >= 2L) re_strict else re_all
  data.table(
    n_design_rows = .N,
    n_designs = uniqueN(matched_design_signature),
    n_studies = uniqueN(study),
    n_cell_lines = uniqueN(CELL_LINE),
    n_tissues = uniqueN(TISSUE),
    n_case_conditions = uniqueN(case_condition),
    total_case_samples = sum(suppressWarnings(as.numeric(n_case_samples)),
                             na.rm = TRUE),
    total_control_samples = sum(suppressWarnings(as.numeric(n_control_samples)),
                                na.rm = TRUE),
    strict_supported_rows = sum(strict_support, na.rm = TRUE),
    strict_supported_designs = uniqueN(matched_design_signature[strict_support]),
    strict_supported_studies = uniqueN(study[strict_support]),
    supported_rows = sum(supported_support, na.rm = TRUE),
    supported_designs = uniqueN(matched_design_signature[supported_support]),
    supported_studies = uniqueN(study[supported_support]),
    weighted_clean_cds_log2fc = weighted_mean_safe(clean_cds_log2fc, row_weight),
    weighted_clean_cds_pct_change =
      pct_from_log2(weighted_mean_safe(clean_cds_log2fc, row_weight)),
    median_clean_cds_pct_change = median_safe(clean_cds_pct_change),
    median_abs_clean_cds_pct_change = median_safe(abs(clean_cds_pct_change)),
    fraction_positive_effects = mean(clean_cds_log2fc > 0, na.rm = TRUE),
    max_abs_clean_cds_pct_change = max_safe(abs(clean_cds_pct_change)),
    best_case_condition =
      case_condition[which.max(abs(clean_cds_pct_change))][1],
    best_control_conditions =
      control_conditions[which.max(abs(clean_cds_pct_change))][1],
    best_study = study[which.max(abs(clean_cds_pct_change))][1],
    nearest_designs = collapse_unique(matched_design_signature, n = 5L),
    random_effects_rows = re_all$n,
    random_effects_log2fc = re_all$effect,
    random_effects_log2fc_lower = re_all$lower,
    random_effects_log2fc_upper = re_all$upper,
    random_effects_pct_change = pct_from_log2(re_all$effect),
    random_effects_pct_change_lower = pct_from_log2(re_all$lower),
    random_effects_pct_change_upper = pct_from_log2(re_all$upper),
    random_effects_tau2 = re_all$tau2,
    random_effects_i2 = re_all$i2,
    random_effects_direction_fraction = re_all$direction_fraction,
    random_effects_ci_excludes_zero = re_all$ci_excludes_zero,
    strict_random_effects_rows = re_strict$n,
    strict_random_effects_log2fc = re_strict$effect,
    strict_random_effects_log2fc_lower = re_strict$lower,
    strict_random_effects_log2fc_upper = re_strict$upper,
    strict_random_effects_pct_change = pct_from_log2(re_strict$effect),
    strict_random_effects_pct_change_lower = pct_from_log2(re_strict$lower),
    strict_random_effects_pct_change_upper = pct_from_log2(re_strict$upper),
    strict_random_effects_i2 = re_strict$i2,
    strict_random_effects_direction_fraction = re_strict$direction_fraction,
    strict_random_effects_ci_excludes_zero = re_strict$ci_excludes_zero,
    atlas_clean_cds_log2fc = use_re$effect,
    atlas_clean_cds_pct_change = pct_from_log2(use_re$effect),
    atlas_clean_cds_pct_lower = pct_from_log2(use_re$lower),
    atlas_clean_cds_pct_upper = pct_from_log2(use_re$upper),
    atlas_interval_source =
      if (re_strict$n[[1]] >= 2L) "strict_random_effects" else "all_rows_random_effects",
    median_bootstrap_certainty_score =
      median_safe(matched_bootstrap_certainty_score),
    max_bootstrap_certainty_score =
      max_safe(matched_bootstrap_certainty_score),
    strongest_matched_support = first_nonempty(
      matched_cds_support[order(-abs(clean_cds_pct_change))]
    ),
    strongest_buffering_interpretation = first_nonempty(
      matched_buffering_interpretation[order(-abs(clean_cds_pct_change))]
    )
  )
}, by = .(gene_symbol, tx_id, design_family)]

features[is.na(design_family) | !nzchar(design_family),
         design_family := perturbation_class]
features[is.na(design_family) | !nzchar(design_family),
         design_family := "unknown"]
features <- features[feature_id != "clean_CDS" & feature_type != "clean_CDS"]
if (nrow(features)) {
  features[, strict_support := bool_col(.SD, "matched_feature_passes_sample_gate") &
             bool_col(.SD, "matched_feature_passes_count_gate_both") &
             suppressWarnings(as.numeric(n_case_samples)) >= 2 &
             suppressWarnings(as.numeric(n_control_samples)) >= 2 &
             suppressWarnings(as.numeric(min_feature_raw_counts)) >= 30]
  features[, feature_family := fifelse(
    grepl("internal|iorf", feature_type, ignore.case = TRUE) |
      grepl("iORF", category, ignore.case = FALSE),
    "iORF",
    fifelse(grepl("nte|ntt", feature_type, ignore.case = TRUE) |
              category %chin% c("NTE", "NTT"),
            "NTE_NTT",
            fifelse(grepl("uorf|uoorf|ouorf|overlap", feature_type,
                          ignore.case = TRUE) |
                      category %chin% c("uORF", "uoORF", "ouORF"),
                    "uORF_ouORF", "other_feature"))
  )]
  features[, feature_row_weight := pmax(1, log1p(pmax(
    suppressWarnings(as.numeric(case_sum_raw_counts)),
    suppressWarnings(as.numeric(control_sum_raw_counts)),
    na.rm = TRUE
  )))]
  features[!is.finite(feature_row_weight), feature_row_weight := 1]

  feature_summary <- features[, .(
    feature_type = first_nonempty(feature_type),
    category = first_nonempty(category),
    feature_family = first_nonempty(feature_family),
    feature_rows = .N,
    feature_designs = uniqueN(sample_set_signature),
    feature_studies = uniqueN(study),
    feature_strict_rows = sum(strict_support, na.rm = TRUE),
    weighted_relative_use_delta =
      weighted_mean_safe(matched_relative_use_delta, feature_row_weight),
    median_relative_use_delta = median_safe(matched_relative_use_delta),
    max_abs_relative_use_delta = max_safe(abs(matched_relative_use_delta)),
    weighted_feature_log2fc =
      weighted_mean_safe(matched_fpkm_log2_fc_pseudocount1, feature_row_weight),
    max_feature_raw_counts =
      max_safe(pmax(suppressWarnings(as.numeric(case_sum_raw_counts)),
                    suppressWarnings(as.numeric(control_sum_raw_counts)),
                    na.rm = TRUE)),
    feature_example_condition = first_nonempty(case_condition),
    feature_example_study = first_nonempty(study)
  ), by = .(gene_symbol, tx_id, design_family, feature_id)]
  feature_summary[, feature_score :=
                    pmin(feature_strict_rows / 2, 1) +
                    pmin(feature_studies / 2, 1) +
                    pmin(abs(weighted_relative_use_delta) / 0.25, 1) +
                    pmin(log1p(max_feature_raw_counts) / log1p(1000), 1)]
  setorder(feature_summary, gene_symbol, tx_id, design_family,
           -feature_score, -max_abs_relative_use_delta)
  top_feature <- feature_summary[, .SD[1], by = .(gene_symbol, tx_id, design_family)]
  setnames(
    top_feature,
    old = setdiff(names(top_feature), c("gene_symbol", "tx_id", "design_family")),
    new = paste0("top_shifted_", setdiff(names(top_feature),
                                         c("gene_symbol", "tx_id", "design_family")))
  )

  top_by_family <- feature_summary[
    feature_family %chin% c("uORF_ouORF", "iORF", "NTE_NTT"),
    .SD[1],
    by = .(gene_symbol, tx_id, design_family, feature_family)
  ]
  family_wide <- dcast(
    top_by_family,
    gene_symbol + tx_id + design_family ~ feature_family,
    value.var = c("feature_id", "weighted_relative_use_delta",
                  "weighted_feature_log2fc", "feature_score"),
    fill = NA
  )
} else {
  top_feature <- data.table(gene_symbol = character(), tx_id = character(),
                            design_family = character())
  family_wide <- copy(top_feature)
}

cards <- merge(cards, top_feature,
               by = c("gene_symbol", "tx_id", "design_family"), all.x = TRUE)
cards <- merge(cards, family_wide,
               by = c("gene_symbol", "tx_id", "design_family"), all.x = TRUE)

if (nrow(sparse)) {
  sparse[, design_family := fifelse(
    !is.na(onoff_design_family) & nzchar(onoff_design_family),
    onoff_design_family,
    "unknown"
  )]
  setorder(sparse, gene_symbol, tx_id, design_family, -sparse_onoff_score)
  sparse_top <- sparse[, .SD[1], by = .(gene_symbol, tx_id, design_family)]
  sparse_top <- sparse_top[, .(
    gene_symbol, tx_id, design_family,
    sparse_feature_id = feature_id,
    sparse_onoff_pattern = onoff_pattern,
    sparse_onoff_class = sparse_onoff_class,
    sparse_onoff_score,
    sparse_n_onoff_rows = n_onoff_rows,
    sparse_n_studies = n_studies,
    sparse_example_studies = example_studies,
    sparse_example_conditions = example_case_conditions,
    sparse_strongest_context = strongest_onoff_context
  )]
  cards <- merge(cards, sparse_top,
                 by = c("gene_symbol", "tx_id", "design_family"), all.x = TRUE)
}

if (nrow(postviral)) {
  postviral_keep <- intersect(names(postviral), c(
    "gene_symbol", "tx_id", "design_family", "context_priority",
    "gene_design_evidence_class", "matched_gene_design_score",
    "atlas_condition_evidence_gate", "atlas_condition_evidence_score",
    "loo_effect_direction_stable", "study_loo_effect_direction_stable",
    "loo_fraction_same_effect_direction", "study_loo_fraction_same_effect_direction",
    "strict_random_effects_interval_supported",
    "random_effects_strict_pct_change",
    "random_effects_strict_pct_change_lower",
    "random_effects_strict_pct_change_upper",
    "best_condition", "best_control", "best_study"
  ))
  postviral_small <- postviral[, ..postviral_keep]
  setnames(
    postviral_small,
    setdiff(names(postviral_small), c("gene_symbol", "tx_id", "design_family")),
    paste0("postviral_", setdiff(names(postviral_small),
                                 c("gene_symbol", "tx_id", "design_family")))
  )
  cards <- merge(cards, postviral_small,
                 by = c("gene_symbol", "tx_id", "design_family"), all.x = TRUE)
}

if (nrow(feature_conf)) {
  gene_conf <- feature_conf[, .(
    n_rdg_features = .N,
    n_measured_features = sum(runs_with_feature_counts > 0, na.rm = TRUE),
    transcript_qc_class = first_nonempty(transcript_qc_class),
    cds_exon_qc_status = first_nonempty(cds_exon_qc_status),
    cds_exon_qc_min_ratio = median_safe(cds_exon_qc_min_ratio),
    n_features_with_signal_bump = sum(cds_signal_bump_supported == TRUE,
                                      na.rm = TRUE),
    n_internal_review_features =
      sum(grepl("review", internal_compendium_bump_class,
                ignore.case = TRUE), na.rm = TRUE),
    has_measured_overlapping_uorf =
      any(feature_type == "overlapping_uorf" | category == "uoORF",
          na.rm = TRUE),
    n_measured_uorf_ouorf =
      sum(feature_type %chin% c("leader_uorf", "overlapping_uorf") |
            category %chin% c("uORF", "uoORF", "ouORF"), na.rm = TRUE),
    n_measured_internal_orf =
      sum(grepl("internal", feature_type, ignore.case = TRUE), na.rm = TRUE),
    n_measured_nte_ntt =
      sum(grepl("nte|ntt", feature_type, ignore.case = TRUE) |
            category %chin% c("NTE", "NTT"), na.rm = TRUE)
  ), by = .(gene_symbol, tx_id)]
  cards <- merge(cards, gene_conf, by = c("gene_symbol", "tx_id"), all.x = TRUE)
}

if (nrow(candidate)) {
  candidate_keep <- intersect(names(candidate), c(
    "gene_symbol", "tx_id", "candidate_rank_v4", "novelty_class",
    "known_uorf_regulation", "manual_m_source_gene", "hidden_uorf_candidate",
    "has_overlapping_uorf", "n_uorfs", "n_overlapping_uorfs",
    "n_uorfs_with_m_support", "primary_clean_cds_usage",
    "next_iteration_candidate_score_v4", "postviral_validation_score",
    "top_matched_postviral_family", "top_matched_postviral_condition",
    "top_matched_postviral_cds_pct", "postviral_design_evidence_class"
  ))
  cards <- merge(cards, unique(candidate[, ..candidate_keep]),
                 by = c("gene_symbol", "tx_id"), all.x = TRUE)
}

if (nrow(browser)) {
  browser[, browser_order := frank(browser_review_rank, ties.method = "first")]
  browser_top <- browser[order(browser_order), .SD[1],
                         by = .(gene_symbol, tx_id)]
  browser_top <- browser_top[, .(
    gene_symbol, tx_id,
    browser_review_rank,
    browser_review_feature = review_feature_id,
    browser_review_reason = review_reason,
    browser_review_question,
    rdg_png_file,
    state_rdg_png_file,
    relative_usage_matrix_file,
    matched_control_table_file,
    browser_review_priority_score
  )]
  cards <- merge(cards, browser_top, by = c("gene_symbol", "tx_id"),
                 all.x = TRUE)
}

if (nrow(downstream)) {
  setorder(downstream, gene_symbol, tx_id, -downstream_tis_evidence_score)
  downstream_top <- downstream[, .SD[1], by = .(gene_symbol, tx_id)]
  downstream_top <- downstream_top[, .(
    gene_symbol, tx_id,
    downstream_tis_feature = feature_id,
    downstream_tis_start = candidate_tis_tx_start,
    downstream_tis_start_codon = candidate_start_codon,
    downstream_tis_kozak_strength = candidate_kozak_strength,
    downstream_tis_review_class = rescue_review_class,
    downstream_tis_evidence_score
  )]
  cards <- merge(cards, downstream_top, by = c("gene_symbol", "tx_id"),
                 all.x = TRUE)
}

if (nrow(branch_usage)) {
  branch_keep <- intersect(names(branch_usage), c(
    "gene_symbol", "tx_id", "design_family",
    "branch_usage_top_feature_id",
    "branch_usage_top_feature_class",
    "branch_usage_top_feature_type",
    "branch_usage_top_support_class",
    "branch_usage_top_score",
    "branch_usage_top_weighted_usage_delta",
    "branch_usage_top_log_or",
    "branch_usage_top_q",
    "branch_usage_top_i2",
    "branch_usage_top_direction_fraction",
    "branch_usage_top_strict_rows",
    "branch_usage_top_strict_studies",
    "branch_usage_top_context",
    "branch_usage_top_control",
    "branch_usage_top_study",
    "feature_id_clean_CDS",
    "weighted_usage_delta_clean_CDS",
    "branch_usage_log_or_q_clean_CDS",
    "branch_usage_score_clean_CDS",
    "feature_id_leader_uORF",
    "weighted_usage_delta_leader_uORF",
    "branch_usage_log_or_q_leader_uORF",
    "branch_usage_score_leader_uORF",
    "feature_id_overlapping_uORF",
    "weighted_usage_delta_overlapping_uORF",
    "branch_usage_log_or_q_overlapping_uORF",
    "branch_usage_score_overlapping_uORF"
  ))
  branch_small <- branch_usage[, ..branch_keep]
  cards <- merge(cards, branch_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(grouped_branch_usage)) {
  grouped_keep <- intersect(names(grouped_branch_usage), c(
    "gene_symbol", "tx_id", "design_family",
    "grouped_branch_usage_top_feature_id",
    "grouped_branch_usage_top_feature_class",
    "grouped_branch_usage_top_feature_type",
    "grouped_branch_usage_top_support_class",
    "grouped_branch_usage_top_score",
    "grouped_branch_usage_top_weighted_usage_delta",
    "grouped_branch_usage_top_log_or",
    "grouped_branch_usage_top_q",
    "grouped_branch_usage_top_i2",
    "grouped_branch_usage_top_direction_fraction",
    "grouped_branch_usage_top_strict_rows",
    "grouped_branch_usage_top_strict_studies",
    "grouped_branch_usage_top_evaluable_rows",
    "grouped_branch_usage_top_context",
    "grouped_branch_usage_top_control",
    "grouped_branch_usage_top_study",
    "grouped_weighted_usage_delta_clean_CDS",
    "grouped_weighted_usage_delta_leader_uORF",
    "grouped_weighted_usage_delta_overlapping_uORF",
    "grouped_branch_usage_log_or_q_clean_CDS",
    "grouped_branch_usage_log_or_q_leader_uORF",
    "grouped_branch_usage_log_or_q_overlapping_uORF"
  ))
  grouped_small <- grouped_branch_usage[, ..grouped_keep]
  cards <- merge(cards, grouped_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
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
  cards <- merge(cards, count_aware_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(joint_branch_allocation)) {
  joint_keep <- intersect(names(joint_branch_allocation), c(
    "gene_symbol", "tx_id", "design_family",
    "joint_rows",
    "joint_evaluable_rows",
    "joint_strict_rows",
    "joint_strict_studies",
    "joint_evaluable_studies",
    "joint_weighted_l1_shift",
    "joint_median_l1_shift",
    "joint_max_l1_shift",
    "joint_weighted_clean_CDS_delta",
    "joint_weighted_leader_uORF_delta",
    "joint_weighted_overlapping_uORF_delta",
    "joint_weighted_other_delta",
    "joint_top_branch_class",
    "joint_top_branch_delta",
    "joint_top_branch_q",
    "joint_top_l1_shift",
    "joint_top_js_divergence",
    "joint_top_review_class",
    "joint_top_review_priority",
    "joint_top_motif_alignment",
    "joint_top_motif_prior_signed",
    "joint_top_motif_prior_strength",
    "joint_top_context",
    "joint_top_studies",
    "joint_top_conditions"
  ))
  joint_small <- joint_branch_allocation[, ..joint_keep]
  cards <- merge(cards, joint_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(hierarchical_joint_allocation)) {
  hierarchical_joint_keep <- intersect(names(hierarchical_joint_allocation), c(
    "gene_symbol", "tx_id", "design_family",
    "hierarchical_joint_top_branch_class",
    "hierarchical_joint_top_delta",
    "hierarchical_joint_top_delta_lower",
    "hierarchical_joint_top_delta_upper",
    "hierarchical_joint_top_delta_q",
    "hierarchical_joint_top_support_class",
    "hierarchical_joint_top_priority",
    "hierarchical_joint_top_studies",
    "hierarchical_joint_top_context_rows",
    "hierarchical_joint_top_direction_fraction",
    "hierarchical_joint_top_i2",
    "hierarchical_joint_top_context",
    "hierarchical_joint_top_study",
    "hierarchical_joint_top_match_scope",
    "hierarchical_joint_l1_effect",
    "hierarchical_joint_supported_branch_count",
    "hierarchical_joint_delta_clean_CDS",
    "hierarchical_joint_delta_leader_uORF",
    "hierarchical_joint_delta_overlapping_uORF",
    "hierarchical_joint_delta_other",
    "hierarchical_joint_priority_clean_CDS",
    "hierarchical_joint_priority_leader_uORF",
    "hierarchical_joint_priority_overlapping_uORF",
    "hierarchical_joint_priority_other",
    "hierarchical_joint_support_class_clean_CDS",
    "hierarchical_joint_support_class_leader_uORF",
    "hierarchical_joint_support_class_overlapping_uORF",
    "hierarchical_joint_support_class_other"
  ))
  hierarchical_joint_small <- hierarchical_joint_allocation[
    , ..hierarchical_joint_keep
  ]
  cards <- merge(cards, hierarchical_joint_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(dirichlet_multinomial_allocation)) {
  dm_keep <- intersect(names(dirichlet_multinomial_allocation), c(
    "gene_symbol", "tx_id", "design_family",
    "dm_top_branch_class",
    "dm_top_delta",
    "dm_top_delta_lower",
    "dm_top_delta_upper",
    "dm_top_delta_q",
    "dm_top_support_class",
    "dm_top_priority",
    "dm_top_studies",
    "dm_top_context_rows",
    "dm_top_direction_fraction",
    "dm_top_i2",
    "dm_top_context",
    "dm_top_prior_source",
    "dm_l1_effect",
    "dm_supported_branch_count",
    "dm_delta_clean_CDS",
    "dm_delta_leader_uORF",
    "dm_delta_overlapping_uORF",
    "dm_delta_other",
    "dm_priority_clean_CDS",
    "dm_priority_leader_uORF",
    "dm_priority_overlapping_uORF",
    "dm_priority_other",
    "dm_support_class_clean_CDS",
    "dm_support_class_leader_uORF",
    "dm_support_class_overlapping_uORF",
    "dm_support_class_other"
  ))
  dm_small <- dirichlet_multinomial_allocation[, ..dm_keep]
  cards <- merge(cards, dm_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(model_agreement)) {
  model_agreement_keep <- intersect(names(model_agreement), c(
    "gene_symbol", "tx_id", "design_family",
    "model_agreement_class",
    "model_agreement_tier",
    "model_agreement_score",
    "model_agreement_review_priority",
    "model_support_count",
    "replicated_model_count",
    "model_support_signature",
    "model_disagreement_flag",
    "model_disagreement_flags",
    "branch_model_consensus",
    "branch_model_direction_disagreement",
    "different_branch_support",
    "allocation_output_disagreement",
    "high_branch_heterogeneity",
    "common_supported_branches",
    "agreeing_common_branches",
    "disagreeing_common_branches",
    "common_supported_branch_names",
    "agreeing_common_branch_names",
    "disagreeing_common_branch_names",
    "clean_cds_output_allocation_alignment",
    "top_model_disagreement_branch",
    "top_model_disagreement_abs_difference",
    "model_agreement_top_context"
  ))
  model_agreement_small <- model_agreement[, ..model_agreement_keep]
  cards <- merge(cards, model_agreement_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(model_agreement_loo)) {
  model_agreement_loo_keep <- intersect(names(model_agreement_loo), c(
    "gene_symbol", "tx_id", "design_family",
    "loo_agreement_class",
    "loo_consensus_class",
    "loo_stability_score",
    "loo_review_priority",
    "loo_relevant_supported_branches",
    "loo_support_robust_branches",
    "loo_direction_robust_support_fragile_branches",
    "loo_direction_fragile_branches",
    "loo_not_evaluable_branches",
    "loo_min_direction_fraction",
    "loo_min_support_retention",
    "loo_min_model_studies",
    "loo_support_robust_branch_names",
    "loo_fragile_branch_names",
    "loo_worst_omitted_studies"
  ))
  model_agreement_loo_small <- model_agreement_loo[
    , ..model_agreement_loo_keep
  ]
  cards <- merge(cards, model_agreement_loo_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(model_agreement_transport)) {
  model_agreement_transport_keep <- intersect(
    names(model_agreement_transport),
    c(
      "gene_symbol", "tx_id", "design_family",
      "transport_agreement_class",
      "transport_consensus_class",
      "transport_stability_score",
      "transport_review_priority",
      "transport_relevant_supported_branches",
      "transport_evaluable_supported_branches",
      "transport_evaluable_axis_tests",
      "transport_possible_axis_tests",
      "transport_min_evaluable_axes",
      "transport_support_robust_branches",
      "transport_direction_robust_support_fragile_branches",
      "transport_direction_fragile_branches",
      "transport_not_evaluable_branches",
      "transport_min_direction_fraction",
      "transport_min_support_retention",
      "transport_min_remaining_studies",
      "transport_coverage_fraction",
      "transport_axis_coverage_fraction",
      "transport_evaluable_axis_names",
      "transport_support_robust_axis_names",
      "transport_fragile_axis_names",
      "transport_worst_omitted_groups",
      "transport_worst_omitted_studies"
    )
  )
  model_agreement_transport_small <- model_agreement_transport[
    , ..model_agreement_transport_keep
  ]
  cards <- merge(cards, model_agreement_transport_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(context_interactions)) {
  context_interactions_keep <- intersect(names(context_interactions), c(
    "gene_symbol", "tx_id", "design_family",
    "context_interaction_rows",
    "context_interaction_evaluable_rows",
    "context_direction_reversal_rows",
    "context_specific_supported_rows",
    "context_magnitude_modulated_rows",
    "context_interaction_top_class",
    "context_interaction_gene_class",
    "context_interaction_top_score",
    "context_interaction_review_priority",
    "context_interaction_top_axis",
    "context_interaction_top_value",
    "context_interaction_top_branch_class",
    "context_interaction_context_delta",
    "context_interaction_context_delta_lower",
    "context_interaction_context_delta_upper",
    "context_interaction_complement_delta",
    "context_interaction_complement_delta_lower",
    "context_interaction_complement_delta_upper",
    "context_interaction_delta",
    "context_interaction_delta_q",
    "context_interaction_context_studies",
    "context_interaction_complement_studies",
    "context_interaction_context_study_names",
    "context_interaction_complement_study_names",
    "context_interaction_summary"
  ))
  context_interactions_small <- context_interactions[
    , ..context_interactions_keep
  ]
  cards <- merge(cards, context_interactions_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(context_partial_pooling)) {
  context_partial_pooling_keep <- intersect(
    names(context_partial_pooling),
    c(
      "gene_symbol", "tx_id", "design_family",
      "n_partial_pooling_rows",
      "n_pooled_supported_rows",
      "n_pooled_review_rows",
      "n_raw_supported_rows",
      "n_raw_survived_rows",
      "n_raw_shrunk_rows",
      "n_single_context_raw_signal_rows",
      "n_pooled_emergent_rows",
      "top_partial_pooling_class",
      "top_partial_pooling_survival_class",
      "top_context_axis",
      "top_context_value",
      "top_branch_class",
      "top_pooled_interaction_delta",
      "top_pooled_interaction_lower",
      "top_pooled_interaction_upper",
      "top_pooled_interaction_q",
      "top_raw_interaction_delta",
      "top_partial_pooling_weight",
      "top_shrinkage_fraction",
      "top_review_priority",
      "top_context_summary",
      "pooled_context_values",
      "partial_pooling_gene_class",
      "partial_pooling_gene_rank"
    )
  )
  context_partial_pooling_small <- context_partial_pooling[
    , ..context_partial_pooling_keep
  ]
  setnames(
    context_partial_pooling_small,
    old = intersect(
      names(context_partial_pooling_small),
      c(
        "n_partial_pooling_rows",
        "n_pooled_supported_rows",
        "n_pooled_review_rows",
        "n_raw_supported_rows",
        "n_raw_survived_rows",
        "n_raw_shrunk_rows",
        "n_single_context_raw_signal_rows",
        "n_pooled_emergent_rows",
        "top_partial_pooling_class",
        "top_partial_pooling_survival_class",
        "top_context_axis",
        "top_context_value",
        "top_branch_class",
        "top_pooled_interaction_delta",
        "top_pooled_interaction_lower",
        "top_pooled_interaction_upper",
        "top_pooled_interaction_q",
        "top_raw_interaction_delta",
        "top_partial_pooling_weight",
        "top_shrinkage_fraction",
        "top_review_priority",
        "top_context_summary",
        "pooled_context_values",
        "partial_pooling_gene_rank"
      )
    ),
    new = c(
      "partial_pooling_rows",
      "partial_pooling_pooled_supported_rows",
      "partial_pooling_pooled_review_rows",
      "partial_pooling_raw_supported_rows",
      "partial_pooling_raw_survived_rows",
      "partial_pooling_raw_shrunk_rows",
      "partial_pooling_single_context_raw_signal_rows",
      "partial_pooling_pooled_emergent_rows",
      "partial_pooling_top_class",
      "partial_pooling_top_survival_class",
      "partial_pooling_top_axis",
      "partial_pooling_top_value",
      "partial_pooling_top_branch_class",
      "partial_pooling_top_delta",
      "partial_pooling_top_lower",
      "partial_pooling_top_upper",
      "partial_pooling_top_q",
      "partial_pooling_top_raw_delta",
      "partial_pooling_top_weight",
      "partial_pooling_top_shrinkage_fraction",
      "partial_pooling_review_priority",
      "partial_pooling_top_context_summary",
      "partial_pooling_context_values",
      "partial_pooling_gene_rank"
    )
  )
  cards <- merge(cards, context_partial_pooling_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(hierarchical)) {
  hierarchical_keep <- intersect(names(hierarchical), c(
    "gene_symbol", "tx_id", "design_family",
    "hierarchical_evidence_class",
    "hierarchical_uncalibrated_evidence_class",
    "hierarchical_priority_score",
    "hierarchical_uncalibrated_priority_score",
    "hierarchical_direction",
    "hierarchical_clean_cds_pct_change",
    "hierarchical_clean_cds_pct_lower",
    "hierarchical_clean_cds_pct_upper",
    "hierarchical_calibrated_clean_cds_pct_lower",
    "hierarchical_calibrated_clean_cds_pct_upper",
    "hierarchical_posterior_log2fc",
    "hierarchical_posterior_log2fc_lower",
    "hierarchical_posterior_log2fc_upper",
    "hierarchical_calibrated_posterior_log2fc_lower",
    "hierarchical_calibrated_posterior_log2fc_upper",
    "hierarchical_posterior_se",
    "hierarchical_ci_excludes_zero",
    "hierarchical_uncalibrated_ci_excludes_zero",
    "hierarchical_calibrated_ci_excludes_zero",
    "hierarchical_calibrated_half_width_log2fc",
    "hierarchical_normal_half_width_log2fc",
    "hierarchical_conformal_q90_log2fc",
    "hierarchical_conformal_q95_log2fc",
    "hierarchical_calibration_scope",
    "hierarchical_calibration_predictions",
    "hierarchical_calibration_mae_log2fc",
    "hierarchical_direction_fraction",
    "hierarchical_direction_stable",
    "hierarchical_evidence_weight",
    "hierarchical_prior_source",
    "hierarchical_loo_rows",
    "hierarchical_loo_mae_log2fc",
    "hierarchical_loo_rmse_log2fc",
    "hierarchical_loo_interval_coverage",
    "hierarchical_loo_calibrated_interval_coverage",
    "hierarchical_loo_direction_accuracy"
  ))
  hierarchical_small <- hierarchical[, ..hierarchical_keep]
  cards <- merge(cards, hierarchical_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(hierarchical_residual)) {
  residual_keep <- intersect(names(hierarchical_residual), c(
    "gene_symbol", "tx_id", "design_family",
    "residual_predictions",
    "residual_mae_log2fc",
    "residual_rmse_log2fc",
    "residual_median_abs_error_log2fc",
    "residual_q90_abs_error_log2fc",
    "residual_q95_abs_error_log2fc",
    "residual_mean_signed_error_log2fc",
    "residual_fraction_over_design_q90",
    "residual_fraction_over_design_q95",
    "residual_direction_accuracy",
    "residual_interval_coverage",
    "residual_calibrated_interval_coverage",
    "residual_worst_abs_error_log2fc",
    "residual_worst_signed_error_log2fc",
    "residual_worst_study",
    "residual_worst_cell_line",
    "residual_worst_tissue",
    "residual_worst_case_condition",
    "residual_worst_control_condition",
    "residual_worst_context_key",
    "residual_top_studies",
    "residual_top_cell_lines",
    "residual_top_tissues",
    "residual_top_conditions",
    "residual_relative_to_family_q95",
    "residual_heterogeneity_score",
    "residual_risk_class"
  ))
  residual_small <- hierarchical_residual[, ..residual_keep]
  cards <- merge(cards, residual_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(hierarchical_context)) {
  context_keep <- intersect(names(hierarchical_context), c(
    "gene_symbol", "tx_id", "design_family",
    "context_calibration_predictions",
    "context_calibration_peer_corrected_fraction",
    "context_calibration_original_mae_log2fc",
    "context_calibration_context_mae_log2fc",
    "context_calibration_original_rmse_log2fc",
    "context_calibration_context_rmse_log2fc",
    "context_calibration_original_q95_log2fc",
    "context_calibration_context_q95_log2fc",
    "context_calibration_mae_delta_log2fc",
    "context_calibration_rmse_delta_log2fc",
    "context_calibration_q95_delta_log2fc",
    "context_calibration_fraction_improved",
    "context_calibration_top_offset_level",
    "context_calibration_top_worse_context",
    "context_calibration_top_worse_abs_error_log2fc",
    "context_calibration_top_improved_contexts",
    "context_calibration_status"
  ))
  context_small <- hierarchical_context[, ..context_keep]
  cards <- merge(cards, context_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(hierarchical_context_model)) {
  context_model_keep <- intersect(names(hierarchical_context_model), c(
    "gene_symbol", "tx_id", "design_family",
    "context_model_predictions",
    "context_model_original_mae_log2fc",
    "context_model_mae_log2fc",
    "context_model_original_rmse_log2fc",
    "context_model_rmse_log2fc",
    "context_model_original_q95_log2fc",
    "context_model_q95_log2fc",
    "context_model_mae_delta_log2fc",
    "context_model_rmse_delta_log2fc",
    "context_model_q95_delta_log2fc",
    "context_model_fraction_improved",
    "context_model_mean_abs_offset_log2fc",
    "context_model_q95_abs_offset_log2fc",
    "context_model_worst_context",
    "context_model_worst_abs_error_log2fc",
    "context_model_top_improved_contexts",
    "context_model_status"
  ))
  context_model_small <- hierarchical_context_model[, ..context_model_keep]
  cards <- merge(cards, context_model_small,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

ensure_column <- function(dt, column, value = NA) {
  if (!column %in% names(dt)) dt[, (column) := value]
  invisible(dt)
}

if (nrow(browser_validation)) {
  for (column in c("gene_symbol", "tx_id", "design_family",
                   "validation_context", "browser_validation_status",
                   "browser_validation_label", "validated_by",
                   "validation_date", "validation_notes",
                   "validation_score")) {
    ensure_column(browser_validation, column)
  }
  browser_validation[, validation_score :=
                       suppressWarnings(as.numeric(validation_score))]
  browser_validation[!is.finite(validation_score),
                     validation_score := fifelse(
                       grepl("support|confirm|valid|interesting",
                             paste(browser_validation_status,
                                   browser_validation_label,
                                   validation_notes),
                             ignore.case = TRUE),
                       1, 0.5
                     )]
  browser_validation[
    is.na(design_family) | !nzchar(design_family),
    design_family := "unknown"
  ]
  setorder(browser_validation, gene_symbol, tx_id, design_family,
           -validation_score, validation_date)
  browser_validation_top <- browser_validation[
    , .SD[1],
    by = .(gene_symbol, tx_id, design_family)
  ][, .(
    gene_symbol, tx_id, design_family,
    browser_validation_status,
    browser_validation_label,
    browser_validation_context = validation_context,
    browser_validation_by = validated_by,
    browser_validation_date = validation_date,
    browser_validation_notes = validation_notes,
    browser_validation_score = validation_score
  )]
  cards <- merge(cards, browser_validation_top,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

if (nrow(joint_branch_rows)) {
  pooling_rows <- copy(joint_branch_rows)
  for (column in c(
    "gene_symbol", "tx_id", "design_family", "study", "case_condition",
    "control_conditions", "TISSUE", "CELL_LINE", "INHIBITOR", "FRACTION",
    "match_scope", "n_case_runs_merged", "n_control_runs_merged",
    "case_branch_total_counts", "control_branch_total_counts",
    "joint_top_branch_delta", "joint_l1_allocation_shift",
    "joint_top_branch_q", "joint_evaluable", "joint_strict"
  )) {
    if (!column %in% names(pooling_rows)) pooling_rows[, (column) := NA]
  }
  for (column in c(
    "n_case_runs_merged", "n_control_runs_merged",
    "case_branch_total_counts", "control_branch_total_counts",
    "joint_top_branch_delta", "joint_l1_allocation_shift",
    "joint_top_branch_q"
  )) {
    pooling_rows[, (column) := suppressWarnings(as.numeric(get(column)))]
  }
  pooling_rows[, match_scope := as.character(match_scope)]
  pooling_rows[, exact_context := match_scope == "exact_context"]
  pooling_rows[is.na(exact_context), exact_context := FALSE]
  pooling_rows[, joint_evaluable_flag := bool_col(.SD, "joint_evaluable")]
  pooling_rows[, joint_strict_flag := bool_col(.SD, "joint_strict")]
  pooling_rows[, control_is_wt_like := control_label_like(control_conditions)]
  pooling_rows[, pooling_protocol_class :=
                 protocol_bias_class(INHIBITOR, FRACTION, case_condition)]
  pooling_rows[, start_site_protocol :=
                 pooling_protocol_class %chin%
                 c("start_site_enriched_ltm",
                   "start_site_enriched_harringtonine")]
  pooling_rows[, pooling_strong_signal :=
                 joint_evaluable_flag == TRUE &
                 is.finite(joint_top_branch_q) &
                 joint_top_branch_q <= 0.10 &
                 is.finite(joint_top_branch_delta) &
                 abs(joint_top_branch_delta) >= 0.05]
  pooling_context <- pooling_rows[, .(
    pooling_context_rows = .N,
    pooling_exact_context_rows = sum(exact_context == TRUE, na.rm = TRUE),
    pooling_relaxed_context_rows = sum(exact_context != TRUE, na.rm = TRUE),
    pooling_joint_strict_rows = sum(joint_strict_flag == TRUE, na.rm = TRUE),
    pooling_joint_evaluable_rows =
      sum(joint_evaluable_flag == TRUE, na.rm = TRUE),
    pooling_context_wt_like_control_rows =
      sum(control_is_wt_like == TRUE, na.rm = TRUE),
    pooling_context_strong_signal_rows =
      sum(pooling_strong_signal == TRUE, na.rm = TRUE),
    pooling_context_studies = count_informative_unique(study),
    pooling_protocol_classes =
      count_informative_unique(pooling_protocol_class),
    pooling_inhibitors = count_informative_unique(INHIBITOR),
    pooling_tissues = count_informative_unique(TISSUE),
    pooling_cell_lines = count_informative_unique(CELL_LINE),
    pooling_case_conditions = count_informative_unique(case_condition),
    pooling_median_case_runs = median_safe(n_case_runs_merged),
    pooling_median_control_runs = median_safe(n_control_runs_merged),
    pooling_median_min_branch_total_counts =
      median_safe(pmin(case_branch_total_counts,
                       control_branch_total_counts, na.rm = TRUE)),
    pooling_median_abs_branch_delta =
      median_safe(abs(joint_top_branch_delta)),
    pooling_median_l1_allocation_shift =
      median_safe(joint_l1_allocation_shift),
    pooling_top_protocol_classes = collapse_top(pooling_protocol_class),
    pooling_top_inhibitors = collapse_top(INHIBITOR),
    pooling_top_tissues = collapse_top(TISSUE),
    pooling_top_cell_lines = collapse_top(CELL_LINE),
    pooling_top_case_conditions = collapse_top(case_condition),
    pooling_top_controls = collapse_top(control_conditions),
    pooling_any_start_site_enriched_protocol =
      any(start_site_protocol == TRUE, na.rm = TRUE),
    pooling_any_relaxed_context = any(exact_context != TRUE, na.rm = TRUE)
  ), by = .(gene_symbol, tx_id, design_family)]
  pooling_context[, pooling_exact_context_fraction := fifelse(
    pooling_context_rows > 0,
    pooling_exact_context_rows / pooling_context_rows,
    0
  )]
  pooling_context[, pooling_joint_strict_fraction := fifelse(
    pooling_context_rows > 0,
    pooling_joint_strict_rows / pooling_context_rows,
    0
  )]
  pooling_context[, pooling_joint_evaluable_fraction := fifelse(
    pooling_context_rows > 0,
    pooling_joint_evaluable_rows / pooling_context_rows,
    0
  )]
  cards <- merge(cards, pooling_context,
                 by = c("gene_symbol", "tx_id", "design_family"),
                 all.x = TRUE)
}

for (column in c(
  "postviral_atlas_condition_evidence_score",
  "postviral_atlas_condition_evidence_gate",
  "postviral_gene_design_evidence_class",
  "postviral_loo_effect_direction_stable",
  "postviral_study_loo_effect_direction_stable",
  "sparse_onoff_score",
  "sparse_feature_id",
  "sparse_onoff_pattern",
  "sparse_onoff_class",
  "browser_validation_status",
  "browser_validation_label",
  "browser_validation_context",
  "browser_validation_by",
  "browser_validation_date",
  "browser_validation_notes",
  "browser_validation_score",
  "branch_usage_top_support_class",
  "branch_usage_top_feature_id",
  "branch_usage_top_control",
  "branch_usage_top_score",
  "joint_top_review_class",
  "joint_top_review_priority",
  "joint_top_branch_class",
  "joint_top_l1_shift",
  "grouped_branch_usage_top_evaluable_rows",
  "grouped_branch_usage_top_control",
  "count_aware_top_evaluable_rows",
  "count_aware_top_control",
  "hierarchical_joint_top_support_class",
  "hierarchical_joint_top_priority",
  "hierarchical_joint_top_branch_class",
  "hierarchical_joint_l1_effect",
  "dm_top_support_class",
  "dm_top_priority",
  "dm_top_branch_class",
  "dm_l1_effect",
  "model_agreement_class",
  "model_agreement_tier",
  "model_agreement_score",
  "model_agreement_review_priority",
  "model_support_count",
  "replicated_model_count",
  "model_support_signature",
  "model_disagreement_flag",
  "model_disagreement_flags",
  "branch_model_consensus",
  "branch_model_direction_disagreement",
  "allocation_output_disagreement",
  "model_agreement_top_context",
  "loo_agreement_class",
  "loo_consensus_class",
  "loo_stability_score",
  "loo_review_priority",
  "loo_relevant_supported_branches",
  "loo_support_robust_branches",
  "loo_direction_robust_support_fragile_branches",
  "loo_direction_fragile_branches",
  "loo_not_evaluable_branches",
  "loo_min_direction_fraction",
  "loo_min_support_retention",
  "loo_min_model_studies",
  "loo_support_robust_branch_names",
  "loo_fragile_branch_names",
  "loo_worst_omitted_studies",
  "transport_agreement_class",
  "transport_consensus_class",
  "transport_stability_score",
  "transport_review_priority",
  "transport_relevant_supported_branches",
  "transport_evaluable_supported_branches",
  "transport_evaluable_axis_tests",
  "transport_possible_axis_tests",
  "transport_min_evaluable_axes",
  "transport_support_robust_branches",
  "transport_direction_robust_support_fragile_branches",
  "transport_direction_fragile_branches",
  "transport_not_evaluable_branches",
  "transport_min_direction_fraction",
  "transport_min_support_retention",
  "transport_min_remaining_studies",
  "transport_coverage_fraction",
  "transport_axis_coverage_fraction",
  "transport_evaluable_axis_names",
  "transport_support_robust_axis_names",
  "transport_fragile_axis_names",
  "transport_worst_omitted_groups",
  "transport_worst_omitted_studies",
  "context_interaction_rows",
  "context_interaction_evaluable_rows",
  "context_direction_reversal_rows",
  "context_specific_supported_rows",
  "context_magnitude_modulated_rows",
  "context_interaction_top_class",
  "context_interaction_gene_class",
  "context_interaction_top_score",
  "context_interaction_review_priority",
  "context_interaction_top_axis",
  "context_interaction_top_value",
  "context_interaction_top_branch_class",
  "context_interaction_context_delta",
  "context_interaction_context_delta_lower",
  "context_interaction_context_delta_upper",
  "context_interaction_complement_delta",
  "context_interaction_complement_delta_lower",
  "context_interaction_complement_delta_upper",
  "context_interaction_delta",
  "context_interaction_delta_q",
	  "context_interaction_context_studies",
	  "context_interaction_complement_studies",
	  "context_interaction_context_study_names",
	  "context_interaction_complement_study_names",
	  "context_interaction_summary",
	  "partial_pooling_rows",
	  "partial_pooling_pooled_supported_rows",
	  "partial_pooling_pooled_review_rows",
	  "partial_pooling_raw_supported_rows",
	  "partial_pooling_raw_survived_rows",
	  "partial_pooling_raw_shrunk_rows",
	  "partial_pooling_single_context_raw_signal_rows",
	  "partial_pooling_pooled_emergent_rows",
	  "partial_pooling_top_class",
	  "partial_pooling_top_survival_class",
	  "partial_pooling_gene_class",
	  "partial_pooling_top_axis",
	  "partial_pooling_top_value",
	  "partial_pooling_top_branch_class",
	  "partial_pooling_top_delta",
	  "partial_pooling_top_lower",
	  "partial_pooling_top_upper",
	  "partial_pooling_top_q",
	  "partial_pooling_top_raw_delta",
	  "partial_pooling_top_weight",
	  "partial_pooling_top_shrinkage_fraction",
	  "partial_pooling_review_priority",
	  "partial_pooling_top_context_summary",
	  "partial_pooling_context_values",
	  "partial_pooling_gene_rank",
	  "pooling_context_rows",
	  "pooling_exact_context_rows",
	  "pooling_relaxed_context_rows",
	  "pooling_joint_strict_rows",
	  "pooling_joint_evaluable_rows",
	  "pooling_context_wt_like_control_rows",
	  "pooling_context_strong_signal_rows",
	  "pooling_context_studies",
	  "pooling_protocol_classes",
	  "pooling_inhibitors",
	  "pooling_tissues",
	  "pooling_cell_lines",
	  "pooling_case_conditions",
	  "pooling_median_case_runs",
	  "pooling_median_control_runs",
	  "pooling_median_min_branch_total_counts",
	  "pooling_median_abs_branch_delta",
	  "pooling_median_l1_allocation_shift",
	  "pooling_top_protocol_classes",
	  "pooling_top_inhibitors",
	  "pooling_top_tissues",
	  "pooling_top_cell_lines",
	  "pooling_top_case_conditions",
	  "pooling_top_controls",
	  "pooling_any_start_site_enriched_protocol",
	  "pooling_any_relaxed_context",
	  "pooling_exact_context_fraction",
	  "pooling_joint_strict_fraction",
	  "pooling_joint_evaluable_fraction",
	  "joint_weighted_clean_CDS_delta",
	  "joint_weighted_leader_uORF_delta",
  "joint_weighted_overlapping_uORF_delta",
  "hierarchical_evidence_class",
  "hierarchical_priority_score",
  "hierarchical_clean_cds_pct_change",
  "hierarchical_ci_excludes_zero",
  "hierarchical_calibrated_ci_excludes_zero",
  "hierarchical_direction_stable",
  "hierarchical_evidence_weight",
  "hierarchical_prior_source",
  "residual_risk_class",
  "residual_heterogeneity_score",
  "residual_rmse_log2fc",
  "residual_worst_context_key",
  "context_calibration_status",
  "context_calibration_rmse_delta_log2fc",
  "context_calibration_fraction_improved",
  "context_model_status",
  "context_model_rmse_delta_log2fc",
  "context_model_mae_delta_log2fc",
  "context_model_q95_delta_log2fc",
  "context_model_fraction_improved",
  "known_uorf_regulation",
  "manual_m_source_gene",
  "has_overlapping_uorf",
  "cds_exon_qc_status"
)) {
  ensure_column(cards, column)
}
for (column in intersect(
  c("joint_top_review_class", "joint_top_branch_class",
    "joint_top_motif_alignment", "hierarchical_joint_top_support_class",
    "hierarchical_joint_top_branch_class", "dm_top_support_class",
    "dm_top_branch_class", "dm_top_prior_source",
    "branch_usage_top_control", "grouped_branch_usage_top_control",
    "count_aware_top_control",
    "model_agreement_class", "model_agreement_tier",
    "model_support_signature", "model_disagreement_flags",
    "model_agreement_top_context", "loo_agreement_class",
    "loo_consensus_class", "loo_support_robust_branch_names",
    "loo_fragile_branch_names", "loo_worst_omitted_studies",
    "transport_agreement_class", "transport_consensus_class",
    "transport_evaluable_axis_names", "transport_support_robust_axis_names",
    "transport_fragile_axis_names", "transport_worst_omitted_groups",
    "transport_worst_omitted_studies",
    "context_interaction_top_class", "context_interaction_gene_class",
    "context_interaction_top_axis", "context_interaction_top_value",
	    "context_interaction_top_branch_class",
	    "context_interaction_context_study_names",
	    "context_interaction_complement_study_names",
	    "context_interaction_summary",
	    "partial_pooling_top_class",
	    "partial_pooling_top_survival_class",
	    "partial_pooling_gene_class",
	    "partial_pooling_top_axis",
	    "partial_pooling_top_value",
	    "partial_pooling_top_branch_class",
	    "partial_pooling_top_context_summary",
	    "partial_pooling_context_values",
	    "pooling_top_protocol_classes",
	    "pooling_top_inhibitors",
	    "pooling_top_tissues",
	    "pooling_top_cell_lines",
	    "pooling_top_case_conditions",
	    "pooling_top_controls"),
	  names(cards)
	)) {
  cards[, (column) := as.character(get(column))]
}
for (column in intersect(
  c("joint_top_review_priority", "joint_top_l1_shift",
    "grouped_branch_usage_top_evaluable_rows",
    "count_aware_top_evaluable_rows",
    "hierarchical_joint_top_priority", "hierarchical_joint_l1_effect",
    "dm_top_priority", "dm_l1_effect",
    "model_agreement_score", "model_agreement_review_priority",
    "model_support_count", "replicated_model_count",
    "loo_stability_score", "loo_review_priority",
    "loo_relevant_supported_branches", "loo_support_robust_branches",
    "loo_direction_robust_support_fragile_branches",
    "loo_direction_fragile_branches", "loo_not_evaluable_branches",
    "loo_min_direction_fraction",
    "loo_min_support_retention", "loo_min_model_studies",
    "transport_stability_score", "transport_review_priority",
    "transport_relevant_supported_branches",
    "transport_evaluable_supported_branches",
    "transport_evaluable_axis_tests",
    "transport_possible_axis_tests",
    "transport_min_evaluable_axes",
    "transport_support_robust_branches",
    "transport_direction_robust_support_fragile_branches",
    "transport_direction_fragile_branches",
    "transport_not_evaluable_branches",
    "transport_min_direction_fraction",
    "transport_min_support_retention",
    "transport_min_remaining_studies",
    "transport_coverage_fraction",
    "transport_axis_coverage_fraction",
    "context_interaction_rows",
    "context_interaction_evaluable_rows",
    "context_direction_reversal_rows",
    "context_specific_supported_rows",
    "context_magnitude_modulated_rows",
    "context_interaction_top_score",
    "context_interaction_review_priority",
    "context_interaction_context_delta",
    "context_interaction_context_delta_lower",
    "context_interaction_context_delta_upper",
    "context_interaction_complement_delta",
    "context_interaction_complement_delta_lower",
    "context_interaction_complement_delta_upper",
    "context_interaction_delta",
    "context_interaction_delta_q",
	    "context_interaction_context_studies",
	    "context_interaction_complement_studies",
	    "partial_pooling_rows",
	    "partial_pooling_pooled_supported_rows",
	    "partial_pooling_pooled_review_rows",
	    "partial_pooling_raw_supported_rows",
	    "partial_pooling_raw_survived_rows",
	    "partial_pooling_raw_shrunk_rows",
	    "partial_pooling_single_context_raw_signal_rows",
	    "partial_pooling_pooled_emergent_rows",
	    "partial_pooling_top_delta",
	    "partial_pooling_top_lower",
	    "partial_pooling_top_upper",
	    "partial_pooling_top_q",
	    "partial_pooling_top_raw_delta",
	    "partial_pooling_top_weight",
	    "partial_pooling_top_shrinkage_fraction",
	    "partial_pooling_review_priority",
	    "partial_pooling_gene_rank",
	    "pooling_context_rows",
	    "pooling_exact_context_rows",
	    "pooling_relaxed_context_rows",
	    "pooling_joint_strict_rows",
	    "pooling_joint_evaluable_rows",
	    "pooling_context_wt_like_control_rows",
	    "pooling_context_strong_signal_rows",
	    "pooling_context_studies",
	    "pooling_protocol_classes",
	    "pooling_inhibitors",
	    "pooling_tissues",
	    "pooling_cell_lines",
	    "pooling_case_conditions",
	    "pooling_median_case_runs",
	    "pooling_median_control_runs",
	    "pooling_median_min_branch_total_counts",
	    "pooling_median_abs_branch_delta",
	    "pooling_median_l1_allocation_shift",
	    "pooling_exact_context_fraction",
	    "pooling_joint_strict_fraction",
	    "pooling_joint_evaluable_fraction",
    "joint_weighted_clean_CDS_delta", "joint_weighted_leader_uORF_delta",
    "joint_weighted_overlapping_uORF_delta"),
  names(cards)
)) {
  cards[, (column) := suppressWarnings(as.numeric(get(column)))]
}

pooling_numeric_defaults <- intersect(c(
  "pooling_context_rows", "pooling_exact_context_rows",
  "pooling_relaxed_context_rows", "pooling_joint_strict_rows",
  "pooling_joint_evaluable_rows", "pooling_context_wt_like_control_rows",
  "pooling_context_strong_signal_rows", "pooling_context_studies",
  "pooling_protocol_classes", "pooling_inhibitors", "pooling_tissues",
  "pooling_cell_lines", "pooling_case_conditions",
  "pooling_median_case_runs", "pooling_median_control_runs",
  "pooling_median_min_branch_total_counts",
  "pooling_median_abs_branch_delta",
  "pooling_median_l1_allocation_shift",
  "pooling_exact_context_fraction", "pooling_joint_strict_fraction",
  "pooling_joint_evaluable_fraction"
), names(cards))
for (column in pooling_numeric_defaults) {
  cards[!is.finite(get(column)), (column) := 0]
}
pooling_text_defaults <- intersect(c(
  "pooling_top_protocol_classes", "pooling_top_inhibitors",
  "pooling_top_tissues", "pooling_top_cell_lines",
  "pooling_top_case_conditions", "pooling_top_controls"
), names(cards))
for (column in pooling_text_defaults) {
  cards[is.na(get(column)), (column) := ""]
}
for (column in intersect(c("pooling_any_start_site_enriched_protocol",
                           "pooling_any_relaxed_context"), names(cards))) {
  cards[, (column) := bool_col(.SD, column)]
}

cards[, replication_component :=
        0.5 * pmin(strict_supported_designs / 2, 1) +
        0.5 * pmin(strict_supported_studies / 2, 1)]
cards[, support_component :=
        0.5 * pmin(strict_supported_rows / 4, 1) +
        0.5 * pmin(supported_rows / 8, 1)]
cards[, interval_component := fifelse(
  strict_random_effects_ci_excludes_zero == TRUE |
    random_effects_ci_excludes_zero == TRUE,
  1, 0
)]
cards[, certainty_component := fifelse(
  is.finite(max_bootstrap_certainty_score),
  pmin(max_bootstrap_certainty_score, 1),
  0
)]
cards[, postviral_component := fifelse(
  is.finite(postviral_atlas_condition_evidence_score),
  pmin(postviral_atlas_condition_evidence_score / 1.5, 1),
  0
)]
cards[, sparse_component := fifelse(
  is.finite(sparse_onoff_score), pmin(sparse_onoff_score, 1), 0
)]
cards[, browser_validation_component := fifelse(
  is.finite(browser_validation_score),
  pmin(browser_validation_score, 1),
  0
)]
cards[, branch_usage_component := fifelse(
  is.finite(branch_usage_top_score),
  pmin(branch_usage_top_score, 1),
  0
)]
if (!"grouped_branch_usage_top_score" %in% names(cards)) {
  cards[, grouped_branch_usage_top_score := NA_real_]
}
if (!"grouped_branch_usage_top_support_class" %in% names(cards)) {
  cards[, grouped_branch_usage_top_support_class := NA_character_]
}
cards[, grouped_branch_usage_component := fifelse(
  is.finite(grouped_branch_usage_top_score),
  pmin(grouped_branch_usage_top_score, 1),
  0
)]
if (!"count_aware_top_score" %in% names(cards)) {
  cards[, count_aware_top_score := NA_real_]
}
if (!"count_aware_top_support_class" %in% names(cards)) {
  cards[, count_aware_top_support_class := NA_character_]
}
cards[, count_aware_branch_allocation_component := fifelse(
  is.finite(count_aware_top_score),
  pmin(count_aware_top_score, 1),
  0
)]
cards[, joint_branch_allocation_component := fifelse(
  joint_top_review_class %chin% c(
    "strict_joint_branch_shift",
    "joint_branch_shift_review",
    "motif_joint_discordant_review"
  ) &
    is.finite(joint_top_review_priority),
  pmin(joint_top_review_priority, 1),
  fifelse(
    joint_top_review_class == "high_motif_prior_weak_joint_evidence" &
      is.finite(joint_top_review_priority),
    0.25 * pmin(joint_top_review_priority, 1),
    0
  )
)]
cards[, hierarchical_joint_allocation_component := fifelse(
  hierarchical_joint_top_support_class %chin% c(
    "hierarchical_joint_replicated_shift",
    "hierarchical_joint_replicated_review",
    "hierarchical_joint_single_study_strong"
  ) &
    is.finite(hierarchical_joint_top_priority),
  pmin(hierarchical_joint_top_priority, 1),
  fifelse(
    is.finite(hierarchical_joint_top_priority),
    0.25 * pmin(hierarchical_joint_top_priority, 1),
    0
  )
)]
cards[, dirichlet_multinomial_component := fifelse(
  dm_top_support_class %chin% c(
    "dm_replicated_shift",
    "dm_replicated_review",
    "dm_single_study_strong"
  ) &
    is.finite(dm_top_priority),
  pmin(dm_top_priority, 1),
  fifelse(
    is.finite(dm_top_priority),
    0.20 * pmin(dm_top_priority, 1),
    0
  )
)]
cards[, model_agreement_component := fifelse(
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  ) &
    is.finite(model_agreement_score),
  pmin(model_agreement_score, 1),
  0
)]
cards[, model_disagreement_component := fifelse(
  model_agreement_tier == "D_model_disagreement_review", 1, 0
)]
cards[, loo_stability_component := fifelse(
  loo_consensus_class == "consensus_loo_support_robust" &
    is.finite(loo_stability_score),
  pmin(loo_stability_score, 1),
  fifelse(
    loo_consensus_class == "consensus_loo_direction_robust" &
      is.finite(loo_stability_score),
    0.5 * pmin(loo_stability_score, 1),
    0
  )
)]
cards[, loo_fragility_component := fifelse(
  loo_consensus_class == "consensus_loo_fragile", 1,
  fifelse(
    loo_consensus_class == "consensus_loo_not_evaluable", 0.5,
    fifelse(
      loo_agreement_class == "loo_direction_robust_support_fragile", 0.5, 0
    )
  )
)]
cards[, transport_stability_component := fifelse(
  transport_consensus_class == "consensus_transport_support_robust" &
    is.finite(transport_stability_score),
  pmin(transport_stability_score, 1),
  fifelse(
    transport_consensus_class == "consensus_transport_direction_robust" &
      is.finite(transport_stability_score),
    0.5 * pmin(transport_stability_score, 1),
    0
  )
)]
cards[, transport_fragility_component := fifelse(
  transport_consensus_class == "consensus_transport_fragile", 1,
  fifelse(
    transport_consensus_class == "consensus_transport_partially_evaluable",
    0.75,
    fifelse(
      transport_consensus_class == "consensus_transport_not_evaluable", 0.5,
      fifelse(
        transport_agreement_class ==
          "transport_direction_robust_support_fragile",
        0.5, 0
      )
    )
  )
)]
cards[, context_interaction_component := fifelse(
  context_interaction_gene_class %chin% c(
    "context_direction_reversal",
    "context_specific_supported_shift",
    "context_magnitude_modulated"
  ) &
    is.finite(context_interaction_top_score),
  pmin(context_interaction_top_score, 1),
  0
)]
cards[, partial_pooling_component := fifelse(
  partial_pooling_gene_class == "partial_pooled_context_supported" &
    is.finite(partial_pooling_review_priority),
  pmin(partial_pooling_review_priority, 1),
  fifelse(
    partial_pooling_gene_class == "partial_pooled_context_review" &
      is.finite(partial_pooling_review_priority),
    0.65 * pmin(partial_pooling_review_priority, 1),
    fifelse(
      partial_pooling_gene_class == "single_context_raw_signal_not_pooled" &
        is.finite(partial_pooling_review_priority),
      0.45 * pmin(partial_pooling_review_priority, 1),
      fifelse(
        partial_pooling_gene_class ==
          "raw_context_signal_fragile_after_pooling" &
          is.finite(partial_pooling_review_priority),
        0.20 * pmin(partial_pooling_review_priority, 1),
        0
      )
    )
  )
)]
cards[, hierarchical_component := fifelse(
  is.finite(hierarchical_priority_score),
  pmin(hierarchical_priority_score, 1),
  0
)]
cards[, residual_stability_component := fifelse(
  is.finite(residual_heterogeneity_score),
  1 - pmin(residual_heterogeneity_score, 1),
  0.5
)]
cards[, context_calibration_component := fifelse(
  !is.na(context_calibration_status) &
    context_calibration_status == "context_reduces_residual_error",
  pmin(pmax(context_calibration_rmse_delta_log2fc, 0) / 0.25, 1),
  0
)]
cards[, context_model_component := fifelse(
  !is.na(context_model_status) &
    context_model_status == "context_model_reduces_residual_error",
  pmin(pmax(context_model_rmse_delta_log2fc, 0) / 0.25, 1),
  0
)]
cards[, structure_component := as.numeric(
  known_uorf_regulation == TRUE |
    manual_m_source_gene == TRUE |
    has_overlapping_uorf == TRUE
)]
cards[!is.finite(structure_component), structure_component := 0]
cards[, atlas_i2 := fifelse(
  strict_random_effects_rows >= 2 & is.finite(strict_random_effects_i2),
  strict_random_effects_i2,
  random_effects_i2
)]
cards[, atlas_direction_fraction := fifelse(
  strict_random_effects_rows >= 2 &
    is.finite(strict_random_effects_direction_fraction),
  strict_random_effects_direction_fraction,
  random_effects_direction_fraction
)]
cards[, atlas_direction_stable :=
        is.finite(atlas_direction_fraction) & atlas_direction_fraction >= 0.75]
cards[, atlas_interval_supported :=
        strict_random_effects_ci_excludes_zero == TRUE |
        random_effects_ci_excludes_zero == TRUE]
cards[, pooling_has_wt_like_control :=
        pooling_context_wt_like_control_rows > 0 |
        control_label_like(paste(best_control_conditions,
                                 branch_usage_top_control,
                                 grouped_branch_usage_top_control,
                                 count_aware_top_control,
                                 pooling_top_controls,
                                 sep = " "))]
cards[, pooling_proof_tier := fcase(
  pooling_exact_context_rows > 0 &
    pooling_joint_strict_rows > 0 &
    pooling_has_wt_like_control == TRUE,
  "A_exact_strict_project_matched_control",
  pooling_exact_context_rows > 0 &
    (pooling_joint_evaluable_rows > 0 |
       grouped_branch_usage_top_evaluable_rows > 0 |
       count_aware_top_evaluable_rows > 0),
  "B_exact_or_replicated_branch_model",
  model_agreement_tier %chin% c(
    "A_browser_validated_consensus",
    "B_replicated_branch_consensus",
    "C_multi_model_consensus"
  ) |
    loo_consensus_class %chin% c(
      "consensus_loo_support_robust",
      "consensus_loo_direction_robust"
    ) |
    transport_consensus_class %chin% c(
      "consensus_transport_support_robust",
      "consensus_transport_direction_robust"
    ) |
    pooling_context_strong_signal_rows > 0,
  "C_pooled_recurrence_or_transport",
  context_interaction_gene_class %chin% c(
    "context_direction_reversal",
    "context_specific_supported_shift",
    "context_magnitude_modulated"
  ) |
    partial_pooling_gene_class %chin% c(
      "partial_pooled_context_supported",
      "partial_pooled_context_review",
      "single_context_raw_signal_not_pooled",
      "raw_context_signal_fragile_after_pooling"
    ),
  "D_context_specific_or_partial_pooled",
  default = "E_low_count_or_background"
)]
cards[, pooling_proof_component := fcase(
  pooling_proof_tier == "A_exact_strict_project_matched_control", 1,
  pooling_proof_tier == "B_exact_or_replicated_branch_model", 0.75,
  pooling_proof_tier == "C_pooled_recurrence_or_transport", 0.55,
  pooling_proof_tier == "D_context_specific_or_partial_pooled", 0.40,
  default = 0.10
)]
cards[, pooling_claim_language := fcase(
  pooling_proof_tier == "A_exact_strict_project_matched_control",
  "project-specific matched-control RDG claim candidate",
  pooling_proof_tier == "B_exact_or_replicated_branch_model",
  "exact or replicated branch-model support",
  pooling_proof_tier == "C_pooled_recurrence_or_transport",
  "pooled recurrent signal requiring local-control confirmation",
  pooling_proof_tier == "D_context_specific_or_partial_pooled",
  "context-specific review lead; do not generalize yet",
  default = "low-count or background screen"
)]
cards[, pooling_relaxed_fraction := fifelse(
  pooling_context_rows > 0,
  pooling_relaxed_context_rows / pooling_context_rows,
  0
)]
cards[, pooling_hazard_score := pmin(
  1,
  0.18 * pooling_relaxed_fraction +
    0.14 * pmin(pmax(pooling_protocol_classes - 1, 0) / 3, 1) +
    0.14 * pmin(pmax(pooling_tissues - 1, 0) / 6, 1) +
    0.14 * pmin(pmax(pooling_cell_lines - 1, 0) / 6, 1) +
    0.12 * pmin(pmax(pooling_case_conditions - 1, 0) / 8, 1) +
    0.12 * pmax(1 - pooling_joint_strict_fraction, 0) +
    0.08 * as.numeric(pooling_any_start_site_enriched_protocol == TRUE) +
    0.06 * as.numeric(design_family %chin%
                        c("tumor_cancer_context", "other_non_control")) +
    0.06 * as.numeric(partial_pooling_gene_class %chin%
                        c("single_context_raw_signal_not_pooled",
                          "raw_context_signal_fragile_after_pooling")) +
    0.06 * pmin(pmax(atlas_i2, 0), 1) +
    0.05 * as.numeric(high_branch_heterogeneity == TRUE)
)]
cards[!is.finite(pooling_hazard_score), pooling_hazard_score := 0]
cards[, pooling_hazard_class := fcase(
  pooling_hazard_score >= 0.65, "high_pooling_confounding_hazard",
  pooling_hazard_score >= 0.40, "moderate_pooling_confounding_hazard",
  default = "lower_pooling_confounding_hazard"
)]
cards[, pooling_hazard_sources := {
  tag_list <- vector("list", .N)
  for (i in seq_len(.N)) {
    tags <- character()
    if (pooling_relaxed_fraction[i] >= 0.25) {
      tags <- c(tags, "relaxed_context_pooling")
    }
    if (pooling_protocol_classes[i] > 1) tags <- c(tags, "multi_protocol_pool")
    if (pooling_tissues[i] > 1) tags <- c(tags, "multi_tissue_pool")
    if (pooling_cell_lines[i] > 1) tags <- c(tags, "multi_cell_line_pool")
    if (pooling_case_conditions[i] > 2) tags <- c(tags, "multi_condition_pool")
    if (isTRUE(pooling_any_start_site_enriched_protocol[i])) {
      tags <- c(tags, "start_site_protocol_present")
    }
    if (design_family[i] %chin% c("tumor_cancer_context",
                                  "other_non_control")) {
      tags <- c(tags, "broad_design_family")
    }
    if (partial_pooling_gene_class[i] %chin%
        c("single_context_raw_signal_not_pooled",
          "raw_context_signal_fragile_after_pooling")) {
      tags <- c(tags, "context_signal_pooling_sensitive")
    }
    if (isTRUE(atlas_i2[i] > 0.75) ||
        isTRUE(high_branch_heterogeneity[i] == TRUE)) {
      tags <- c(tags, "heterogeneous_effects")
    }
    if (!length(tags)) tags <- "none"
    tag_list[[i]] <- paste(unique(tags), collapse = "; ")
  }
  unlist(tag_list, use.names = FALSE)
}]
cards[, atlas_confidence_score :=
        0.30 * replication_component +
        0.18 * support_component +
        0.16 * interval_component +
        0.12 * certainty_component +
        0.08 * postviral_component +
        0.05 * sparse_component +
        0.03 * branch_usage_component +
        0.05 * grouped_branch_usage_component +
        0.04 * count_aware_branch_allocation_component +
        0.04 * joint_branch_allocation_component +
        0.03 * hierarchical_joint_allocation_component +
        0.02 * dirichlet_multinomial_component +
        0.03 * model_agreement_component +
        0.04 * loo_stability_component +
        0.04 * transport_stability_component +
        0.02 * context_interaction_component +
        0.02 * partial_pooling_component +
        0.03 * pooling_proof_component +
        0.09 * hierarchical_component +
        0.02 * residual_stability_component]
cards[is.finite(atlas_i2) & atlas_i2 > 0.75,
      atlas_confidence_score :=
        atlas_confidence_score *
        (1 - pmin((atlas_i2 - 0.75) / 0.25, 1) * 0.15)]
cards[design_family %chin% c("other_non_control", "tumor_cancer_context"),
      atlas_confidence_score := atlas_confidence_score * 0.9]
cards[strict_supported_rows == 0 & is.finite(sparse_onoff_score),
      atlas_confidence_score := pmax(atlas_confidence_score,
                                     0.25 + 0.35 * pmin(sparse_onoff_score, 1))]
cards[, atlas_confidence_score := pmin(pmax(atlas_confidence_score, 0), 1)]
cards[, atlas_browser_review_priority_score :=
        pmin(
          1.5,
          atlas_confidence_score +
            0.20 * browser_validation_component +
            0.10 * as.numeric(design_family == "viral_infection") +
            0.07 * postviral_component +
            0.04 * branch_usage_component +
            0.08 * grouped_branch_usage_component +
            0.08 * count_aware_branch_allocation_component +
            0.08 * joint_branch_allocation_component +
            0.07 * hierarchical_joint_allocation_component +
            0.05 * dirichlet_multinomial_component +
            0.07 * model_agreement_component +
            0.06 * model_disagreement_component +
            0.06 * loo_stability_component +
            0.08 * loo_fragility_component +
            0.06 * transport_stability_component +
            0.08 * transport_fragility_component +
            0.08 * context_interaction_component +
            0.07 * partial_pooling_component +
            0.05 * pooling_proof_component +
            0.05 * pooling_hazard_score +
            0.07 * hierarchical_component +
            0.03 * as.numeric(residual_risk_class ==
                                "high_residual_heterogeneity") +
            0.03 * context_calibration_component +
            0.04 * context_model_component +
            0.05 * structure_component +
            0.05 * pmin(abs(atlas_clean_cds_pct_change) / 100, 1) +
            0.03 * sparse_component
        )]
cards[!is.finite(atlas_browser_review_priority_score),
      atlas_browser_review_priority_score :=
        fifelse(is.finite(atlas_confidence_score), atlas_confidence_score, 0)]

cards[, atlas_evidence_class := fifelse(
  strict_supported_designs >= 2 & strict_supported_studies >= 2 &
    atlas_interval_supported == TRUE & atlas_direction_stable == TRUE,
  "atlas_ready_replicated_interval",
  fifelse(strict_supported_designs >= 2 & strict_supported_studies >= 2,
          "replicated_strict_review",
          fifelse(strict_supported_rows > 0,
                  "single_design_strict_review",
                  fifelse(is.finite(sparse_onoff_score),
                          "sparse_onoff_review",
                          "exploratory_low_support"))))]

cards[, atlas_direction := fifelse(
  atlas_clean_cds_log2fc > 0, "clean_CDS_up",
  fifelse(atlas_clean_cds_log2fc < 0, "clean_CDS_down", "flat_or_unknown")
)]

cards[, warnings := {
  warning_list <- vector("list", .N)
  for (i in seq_len(.N)) {
    w <- character()
    if (strict_supported_studies[i] < 2) w <- c(w, "low_study_replication")
    if (strict_supported_designs[i] < 2) w <- c(w, "low_design_replication")
    if (design_family[i] %chin% c("other_non_control", "tumor_cancer_context")) {
      w <- c(w, "broad_or_ambiguous_design_family")
    }
    if (isTRUE(atlas_i2[i] > 0.75)) w <- c(w, "high_between_design_heterogeneity")
    if (!isTRUE(atlas_direction_stable[i])) w <- c(w, "direction_not_stable")
    if (isTRUE(cds_exon_qc_status[i] != "pass")) w <- c(w, "cds_exon_qc_not_pass")
    if (isTRUE(browser_validation_score[i] > 0)) w <- c(w, "browser_supported")
    if (isTRUE(grepl("branch_usage_shift", branch_usage_top_support_class[i],
                     fixed = TRUE))) {
      w <- c(w, "branch_usage_shift_supported")
    }
    if (isTRUE(grepl("group_branch_shift",
                     grouped_branch_usage_top_support_class[i],
                     fixed = TRUE))) {
      w <- c(w, "grouped_branch_usage_shift_supported")
    }
    if (isTRUE(grepl("count_aware_branch_shift",
                     count_aware_top_support_class[i],
                     fixed = TRUE))) {
      w <- c(w, "count_aware_branch_allocation_supported")
    }
    if (isTRUE(joint_top_review_class[i] %chin%
               c("strict_joint_branch_shift",
                 "joint_branch_shift_review",
                 "motif_joint_discordant_review"))) {
      w <- c(w, "joint_branch_allocation_supported")
    }
    if (isTRUE(hierarchical_joint_top_support_class[i] %chin%
               c("hierarchical_joint_replicated_shift",
                 "hierarchical_joint_replicated_review",
                 "hierarchical_joint_single_study_strong"))) {
      w <- c(w, "hierarchical_joint_allocation_supported")
    }
    if (isTRUE(dm_top_support_class[i] %chin%
               c("dm_replicated_shift",
                 "dm_replicated_review",
                 "dm_single_study_strong"))) {
      w <- c(w, "dirichlet_multinomial_supported")
    }
    if (isTRUE(model_agreement_tier[i] %chin%
               c("A_browser_validated_consensus",
                 "B_replicated_branch_consensus",
                 "C_multi_model_consensus"))) {
      w <- c(w, "multi_model_consensus")
    }
    if (isTRUE(model_agreement_tier[i] ==
               "D_model_disagreement_review")) {
      w <- c(w, "model_disagreement_review")
    }
    if (isTRUE(loo_consensus_class[i] ==
               "consensus_loo_support_robust")) {
      w <- c(w, "loo_consensus_support_robust")
    } else if (isTRUE(loo_consensus_class[i] ==
                      "consensus_loo_direction_robust")) {
      w <- c(w, "loo_consensus_direction_only")
    } else if (isTRUE(loo_consensus_class[i] ==
                      "consensus_loo_fragile")) {
      w <- c(w, "loo_consensus_fragile")
    } else if (isTRUE(loo_consensus_class[i] ==
                      "consensus_loo_not_evaluable")) {
      w <- c(w, "loo_consensus_not_evaluable")
    }
    if (isTRUE(transport_consensus_class[i] ==
               "consensus_transport_support_robust")) {
      w <- c(w, "transport_consensus_support_robust")
    } else if (isTRUE(transport_consensus_class[i] ==
                      "consensus_transport_direction_robust")) {
      w <- c(w, "transport_consensus_direction_only")
    } else if (isTRUE(transport_consensus_class[i] ==
                      "consensus_transport_fragile")) {
      w <- c(w, "transport_consensus_fragile")
    } else if (isTRUE(transport_consensus_class[i] ==
                      "consensus_transport_partially_evaluable")) {
      w <- c(w, "transport_consensus_partially_evaluable")
    } else if (isTRUE(transport_consensus_class[i] ==
                      "consensus_transport_not_evaluable")) {
      w <- c(w, "transport_consensus_not_evaluable")
    }
    if (isTRUE(context_interaction_gene_class[i] ==
               "context_direction_reversal")) {
      w <- c(w, "context_direction_reversal")
    } else if (isTRUE(context_interaction_gene_class[i] ==
                      "context_specific_supported_shift")) {
      w <- c(w, "context_specific_branch_shift")
	    } else if (isTRUE(context_interaction_gene_class[i] ==
	                      "context_magnitude_modulated")) {
	      w <- c(w, "context_magnitude_modulated")
	    }
    if (isTRUE(partial_pooling_gene_class[i] ==
	               "partial_pooled_context_supported")) {
	      w <- c(w, "partial_pooled_context_supported")
	    } else if (isTRUE(partial_pooling_gene_class[i] ==
	                      "partial_pooled_context_review")) {
	      w <- c(w, "partial_pooled_context_review")
	    } else if (isTRUE(partial_pooling_gene_class[i] ==
	                      "single_context_raw_signal_not_pooled")) {
	      w <- c(w, "single_context_context_signal")
	    } else if (isTRUE(partial_pooling_gene_class[i] ==
	                      "raw_context_signal_fragile_after_pooling")) {
	      w <- c(w, "context_signal_fragile_after_pooling")
	    }
    if (isTRUE(pooling_proof_tier[i] ==
               "A_exact_strict_project_matched_control")) {
      w <- c(w, "pooling_exact_project_matched_control")
    } else if (isTRUE(pooling_proof_tier[i] ==
                      "C_pooled_recurrence_or_transport")) {
      w <- c(w, "pooled_recurrence_not_project_specific")
    }
    if (isTRUE(pooling_hazard_class[i] ==
               "high_pooling_confounding_hazard")) {
      w <- c(w, "pooling_high_confounding_hazard")
    } else if (isTRUE(pooling_hazard_class[i] ==
                      "moderate_pooling_confounding_hazard")) {
      w <- c(w, "pooling_moderate_confounding_hazard")
    }
    if (isTRUE(pooling_any_start_site_enriched_protocol[i])) {
      w <- c(w, "start_site_protocol_bias_present")
    }
	    if (isTRUE(hierarchical_calibrated_ci_excludes_zero[i] == TRUE)) {
	      w <- c(w, "hierarchical_calibrated_interval_supported")
	    } else if (isTRUE(hierarchical_ci_excludes_zero[i] == TRUE)) {
      w <- c(w, "hierarchical_interval_lost_after_calibration")
    }
    if (isTRUE(hierarchical_prior_source[i] == "gene_design_prior_dominant")) {
      w <- c(w, "hierarchical_prior_dominant")
    }
    if (isTRUE(residual_risk_class[i] == "high_residual_heterogeneity")) {
      w <- c(w, "high_hierarchical_residual_heterogeneity")
    } else if (isTRUE(residual_risk_class[i] ==
                      "moderate_residual_heterogeneity")) {
      w <- c(w, "moderate_hierarchical_residual_heterogeneity")
    }
    if (isTRUE(context_calibration_status[i] ==
               "context_reduces_residual_error")) {
      w <- c(w, "context_calibration_reduces_residual_error")
    } else if (isTRUE(context_calibration_status[i] ==
                      "context_increases_residual_error")) {
      w <- c(w, "context_calibration_increases_residual_error")
    }
    if (isTRUE(context_model_status[i] ==
               "context_model_reduces_residual_error")) {
      w <- c(w, "context_model_reduces_residual_error")
    } else if (isTRUE(context_model_status[i] ==
                      "context_model_increases_residual_error")) {
      w <- c(w, "context_model_increases_residual_error")
    }
    if (strict_supported_rows[i] == 0 && is.finite(sparse_onoff_score[i])) {
      w <- c(w, "sparse_onoff_not_fold_change")
    }
    if (!length(w)) w <- "none"
    warning_list[[i]] <- paste(unique(w), collapse = "; ")
  }
  unlist(warning_list, use.names = FALSE)
}]

cards[, abs_atlas_clean_cds_pct_change := abs(atlas_clean_cds_pct_change)]
setorder(cards, -atlas_confidence_score, -abs_atlas_clean_cds_pct_change,
         gene_symbol, design_family)

preferred_cols <- c(
  "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
  "atlas_confidence_score", "atlas_browser_review_priority_score",
  "atlas_direction",
  "atlas_clean_cds_pct_change", "atlas_clean_cds_pct_lower",
  "atlas_clean_cds_pct_upper", "abs_atlas_clean_cds_pct_change",
  "atlas_interval_source",
  "weighted_clean_cds_pct_change", "median_clean_cds_pct_change",
  "fraction_positive_effects", "atlas_direction_fraction", "atlas_i2",
  "atlas_interval_supported", "atlas_direction_stable",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_proof_component", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "pooling_context_rows", "pooling_exact_context_rows",
  "pooling_relaxed_context_rows", "pooling_joint_strict_rows",
  "pooling_joint_evaluable_rows", "pooling_exact_context_fraction",
  "pooling_joint_strict_fraction", "pooling_joint_evaluable_fraction",
  "pooling_context_wt_like_control_rows",
  "pooling_context_strong_signal_rows", "pooling_context_studies",
  "pooling_protocol_classes", "pooling_inhibitors", "pooling_tissues",
  "pooling_cell_lines", "pooling_case_conditions",
  "pooling_median_case_runs", "pooling_median_control_runs",
  "pooling_median_min_branch_total_counts",
  "pooling_median_abs_branch_delta",
  "pooling_median_l1_allocation_shift",
  "pooling_top_protocol_classes", "pooling_top_inhibitors",
  "pooling_top_tissues", "pooling_top_cell_lines",
  "pooling_top_case_conditions", "pooling_top_controls",
  "pooling_any_start_site_enriched_protocol",
  "pooling_any_relaxed_context",
  "n_designs", "n_studies", "n_cell_lines",
  "n_tissues", "strict_supported_rows", "strict_supported_designs",
  "strict_supported_studies", "supported_rows", "supported_designs",
  "supported_studies", "best_case_condition", "best_control_conditions",
  "best_study", "top_shifted_feature_id", "top_shifted_feature_type",
  "top_shifted_feature_family", "top_shifted_weighted_relative_use_delta",
  "feature_id_uORF_ouORF", "weighted_relative_use_delta_uORF_ouORF",
  "feature_id_iORF", "weighted_relative_use_delta_iORF",
  "feature_id_NTE_NTT", "weighted_relative_use_delta_NTE_NTT",
  "sparse_feature_id", "sparse_onoff_pattern", "sparse_onoff_class",
  "sparse_onoff_score", "postviral_atlas_condition_evidence_gate",
  "browser_validation_status", "browser_validation_label",
  "browser_validation_context", "browser_validation_by",
  "browser_validation_date", "browser_validation_notes",
  "browser_validation_score",
  "browser_validation_component",
  "branch_usage_top_feature_id", "branch_usage_top_feature_class",
  "branch_usage_top_support_class", "branch_usage_top_score",
  "branch_usage_top_weighted_usage_delta", "branch_usage_top_log_or",
  "branch_usage_top_q", "branch_usage_component", "branch_usage_top_i2",
  "branch_usage_top_direction_fraction", "branch_usage_top_strict_rows",
  "branch_usage_top_strict_studies", "branch_usage_top_context",
  "branch_usage_top_study", "weighted_usage_delta_clean_CDS",
  "weighted_usage_delta_leader_uORF",
  "weighted_usage_delta_overlapping_uORF",
  "grouped_branch_usage_top_feature_id",
  "grouped_branch_usage_top_feature_class",
  "grouped_branch_usage_top_support_class",
  "grouped_branch_usage_top_score",
  "grouped_branch_usage_top_weighted_usage_delta",
  "grouped_branch_usage_top_log_or",
  "grouped_branch_usage_top_q", "grouped_branch_usage_component",
  "grouped_branch_usage_top_i2",
  "grouped_branch_usage_top_direction_fraction",
  "grouped_branch_usage_top_strict_rows",
  "grouped_branch_usage_top_strict_studies",
  "grouped_branch_usage_top_evaluable_rows",
  "grouped_branch_usage_top_context",
  "grouped_branch_usage_top_study",
  "grouped_weighted_usage_delta_clean_CDS",
  "grouped_weighted_usage_delta_leader_uORF",
  "grouped_weighted_usage_delta_overlapping_uORF",
  "count_aware_top_feature_id",
  "count_aware_top_feature_class",
  "count_aware_top_support_class",
  "count_aware_top_score",
  "count_aware_top_weighted_usage_delta",
  "count_aware_top_random_effects_delta",
  "count_aware_top_delta_lower",
  "count_aware_top_delta_upper",
  "count_aware_top_q",
  "count_aware_branch_allocation_component",
  "count_aware_top_i2",
  "count_aware_top_direction_fraction",
  "count_aware_top_strict_rows",
  "count_aware_top_strict_studies",
  "count_aware_top_evaluable_rows",
  "count_aware_top_context",
  "count_aware_top_study",
  "count_aware_weighted_usage_delta_clean_CDS",
  "count_aware_weighted_usage_delta_leader_uORF",
  "count_aware_weighted_usage_delta_overlapping_uORF",
  "joint_top_review_class",
  "joint_top_review_priority",
  "joint_branch_allocation_component",
  "joint_top_branch_class",
  "joint_top_branch_delta",
  "joint_top_l1_shift",
  "joint_top_branch_q",
  "joint_top_motif_alignment",
  "joint_top_motif_prior_strength",
  "joint_evaluable_rows",
  "joint_strict_rows",
  "joint_strict_studies",
  "joint_weighted_l1_shift",
  "joint_weighted_clean_CDS_delta",
  "joint_weighted_leader_uORF_delta",
  "joint_weighted_overlapping_uORF_delta",
  "joint_weighted_other_delta",
  "joint_top_context",
  "hierarchical_joint_top_support_class",
  "hierarchical_joint_top_priority",
  "hierarchical_joint_allocation_component",
  "hierarchical_joint_top_branch_class",
  "hierarchical_joint_top_delta",
  "hierarchical_joint_top_delta_lower",
  "hierarchical_joint_top_delta_upper",
  "hierarchical_joint_top_delta_q",
  "hierarchical_joint_l1_effect",
  "hierarchical_joint_supported_branch_count",
  "hierarchical_joint_top_studies",
  "hierarchical_joint_top_context_rows",
  "hierarchical_joint_top_direction_fraction",
  "hierarchical_joint_top_i2",
  "hierarchical_joint_top_context",
  "hierarchical_joint_delta_clean_CDS",
  "hierarchical_joint_delta_leader_uORF",
  "hierarchical_joint_delta_overlapping_uORF",
  "hierarchical_joint_delta_other",
  "dm_top_support_class",
  "dm_top_priority",
  "dirichlet_multinomial_component",
  "dm_top_branch_class",
  "dm_top_delta",
  "dm_top_delta_lower",
  "dm_top_delta_upper",
  "dm_top_delta_q",
  "dm_l1_effect",
  "dm_supported_branch_count",
  "dm_top_studies",
  "dm_top_context_rows",
  "dm_top_direction_fraction",
  "dm_top_i2",
  "dm_top_prior_source",
  "dm_top_context",
  "dm_delta_clean_CDS",
  "dm_delta_leader_uORF",
  "dm_delta_overlapping_uORF",
  "dm_delta_other",
  "model_agreement_class",
  "model_agreement_tier",
  "model_agreement_score",
  "model_agreement_review_priority",
  "model_agreement_component",
  "model_disagreement_component",
  "model_support_count",
  "replicated_model_count",
  "model_support_signature",
  "model_disagreement_flag",
  "model_disagreement_flags",
  "branch_model_consensus",
  "branch_model_direction_disagreement",
  "different_branch_support",
  "allocation_output_disagreement",
  "high_branch_heterogeneity",
  "common_supported_branches",
  "agreeing_common_branches",
  "disagreeing_common_branches",
  "common_supported_branch_names",
  "agreeing_common_branch_names",
  "disagreeing_common_branch_names",
  "clean_cds_output_allocation_alignment",
  "top_model_disagreement_branch",
  "top_model_disagreement_abs_difference",
  "model_agreement_top_context",
  "loo_agreement_class",
  "loo_consensus_class",
  "loo_stability_score",
  "loo_review_priority",
  "loo_stability_component",
  "loo_fragility_component",
  "loo_relevant_supported_branches",
  "loo_support_robust_branches",
  "loo_direction_robust_support_fragile_branches",
  "loo_direction_fragile_branches",
  "loo_not_evaluable_branches",
  "loo_min_direction_fraction",
  "loo_min_support_retention",
  "loo_min_model_studies",
  "loo_support_robust_branch_names",
  "loo_fragile_branch_names",
  "loo_worst_omitted_studies",
  "transport_agreement_class",
  "transport_consensus_class",
  "transport_stability_score",
  "transport_review_priority",
  "transport_stability_component",
  "transport_fragility_component",
  "transport_relevant_supported_branches",
  "transport_evaluable_supported_branches",
  "transport_evaluable_axis_tests",
  "transport_possible_axis_tests",
  "transport_min_evaluable_axes",
  "transport_support_robust_branches",
  "transport_direction_robust_support_fragile_branches",
  "transport_direction_fragile_branches",
  "transport_not_evaluable_branches",
  "transport_min_direction_fraction",
  "transport_min_support_retention",
  "transport_min_remaining_studies",
  "transport_coverage_fraction",
  "transport_axis_coverage_fraction",
  "transport_evaluable_axis_names",
  "transport_support_robust_axis_names",
  "transport_fragile_axis_names",
  "transport_worst_omitted_groups",
  "transport_worst_omitted_studies",
  "context_interaction_gene_class",
  "context_interaction_top_class",
  "context_interaction_top_score",
  "context_interaction_review_priority",
  "context_interaction_component",
  "context_interaction_top_axis",
  "context_interaction_top_value",
  "context_interaction_top_branch_class",
  "context_interaction_context_delta",
  "context_interaction_context_delta_lower",
  "context_interaction_context_delta_upper",
  "context_interaction_complement_delta",
  "context_interaction_complement_delta_lower",
  "context_interaction_complement_delta_upper",
  "context_interaction_delta",
  "context_interaction_delta_q",
	  "context_interaction_context_studies",
	  "context_interaction_complement_studies",
	  "context_interaction_context_study_names",
	  "context_interaction_complement_study_names",
	  "context_interaction_summary",
	  "partial_pooling_gene_class",
	  "partial_pooling_top_class",
	  "partial_pooling_top_survival_class",
	  "partial_pooling_review_priority",
	  "partial_pooling_component",
	  "partial_pooling_top_axis",
	  "partial_pooling_top_value",
	  "partial_pooling_top_branch_class",
	  "partial_pooling_top_delta",
	  "partial_pooling_top_lower",
	  "partial_pooling_top_upper",
	  "partial_pooling_top_q",
	  "partial_pooling_top_raw_delta",
	  "partial_pooling_top_weight",
	  "partial_pooling_top_shrinkage_fraction",
	  "partial_pooling_rows",
	  "partial_pooling_pooled_supported_rows",
	  "partial_pooling_pooled_review_rows",
	  "partial_pooling_raw_supported_rows",
	  "partial_pooling_raw_survived_rows",
	  "partial_pooling_raw_shrunk_rows",
	  "partial_pooling_single_context_raw_signal_rows",
	  "partial_pooling_pooled_emergent_rows",
	  "partial_pooling_top_context_summary",
	  "partial_pooling_context_values",
	  "hierarchical_evidence_class", "hierarchical_priority_score",
  "hierarchical_uncalibrated_evidence_class",
  "hierarchical_uncalibrated_priority_score",
  "hierarchical_direction",
  "hierarchical_clean_cds_pct_change",
  "hierarchical_clean_cds_pct_lower",
  "hierarchical_clean_cds_pct_upper",
  "hierarchical_ci_excludes_zero",
  "hierarchical_uncalibrated_ci_excludes_zero",
  "hierarchical_calibrated_clean_cds_pct_lower",
  "hierarchical_calibrated_clean_cds_pct_upper",
  "hierarchical_calibrated_ci_excludes_zero",
  "hierarchical_calibrated_half_width_log2fc",
  "hierarchical_normal_half_width_log2fc",
  "hierarchical_conformal_q95_log2fc",
  "hierarchical_calibration_scope",
  "hierarchical_calibration_predictions",
  "hierarchical_direction_fraction",
  "hierarchical_direction_stable",
  "hierarchical_evidence_weight",
  "hierarchical_prior_source",
  "hierarchical_component",
  "hierarchical_loo_rows",
  "hierarchical_loo_mae_log2fc",
  "hierarchical_loo_rmse_log2fc",
  "hierarchical_loo_interval_coverage",
  "hierarchical_loo_calibrated_interval_coverage",
  "hierarchical_loo_direction_accuracy",
  "residual_risk_class",
  "residual_heterogeneity_score",
  "residual_predictions",
  "residual_mae_log2fc",
  "residual_rmse_log2fc",
  "residual_fraction_over_design_q90",
  "residual_direction_accuracy",
  "residual_worst_abs_error_log2fc",
  "residual_worst_context_key",
  "residual_top_studies",
  "residual_top_cell_lines",
  "residual_top_tissues",
  "residual_top_conditions",
  "context_calibration_status",
  "context_calibration_predictions",
  "context_calibration_peer_corrected_fraction",
  "context_calibration_original_rmse_log2fc",
  "context_calibration_context_rmse_log2fc",
  "context_calibration_rmse_delta_log2fc",
  "context_calibration_original_q95_log2fc",
  "context_calibration_context_q95_log2fc",
  "context_calibration_q95_delta_log2fc",
  "context_calibration_fraction_improved",
  "context_calibration_top_offset_level",
  "context_calibration_top_worse_context",
  "context_calibration_top_improved_contexts",
  "context_model_status",
  "context_model_predictions",
  "context_model_original_rmse_log2fc",
  "context_model_rmse_log2fc",
  "context_model_rmse_delta_log2fc",
  "context_model_original_q95_log2fc",
  "context_model_q95_log2fc",
  "context_model_q95_delta_log2fc",
  "context_model_fraction_improved",
  "context_model_mean_abs_offset_log2fc",
  "context_model_worst_context",
  "context_model_top_improved_contexts",
  "postviral_atlas_condition_evidence_score",
  "postviral_gene_design_evidence_class",
  "postviral_loo_effect_direction_stable",
  "postviral_study_loo_effect_direction_stable",
  "candidate_rank_v4", "novelty_class", "known_uorf_regulation",
  "manual_m_source_gene", "hidden_uorf_candidate", "has_overlapping_uorf",
  "n_uorfs", "n_overlapping_uorfs", "n_uorfs_with_m_support",
  "primary_clean_cds_usage", "n_rdg_features", "n_measured_uorf_ouorf",
  "n_measured_internal_orf", "n_measured_nte_ntt",
  "transcript_qc_class", "cds_exon_qc_status", "cds_exon_qc_min_ratio",
  "browser_review_rank", "browser_review_feature", "browser_review_reason",
  "browser_review_question", "rdg_png_file", "state_rdg_png_file",
  "relative_usage_matrix_file", "matched_control_table_file",
  "downstream_tis_feature", "downstream_tis_review_class",
  "downstream_tis_evidence_score", "nearest_designs", "warnings"
)
ordered_cols <- c(intersect(preferred_cols, names(cards)),
                  setdiff(names(cards), preferred_cols))
cards <- cards[, ..ordered_cols]

cards_file <- file.path(output_dir, "dominant_rdg_atlas_cards.csv")
fwrite(cards, cards_file)

branch_context_file <- file.path(output_dir,
                                 "dominant_rdg_atlas_branch_contexts.csv")
if (nrow(joint_branch_rows)) {
  branch_keep <- intersect(c(
    "gene_symbol", "tx_id", "design_family", "study", "AUTHOR",
    "case_condition", "TISSUE", "CELL_LINE", "GENE", "INHIBITOR",
    "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
    "control_conditions", "match_scope", "sample_set_signature",
    "n_case_groups", "n_control_groups", "n_case_runs_merged",
    "n_control_runs_merged", "case_branch_total_counts",
    "control_branch_total_counts", "joint_review_class",
    "joint_review_priority", "joint_top_branch_class",
    "joint_top_branch_delta", "joint_top_branch_q",
    "joint_l1_allocation_shift", "joint_max_abs_branch_delta",
    "joint_js_divergence", "joint_chisq_q", "joint_delta_clean_CDS",
    "joint_delta_leader_uORF", "joint_delta_overlapping_uORF",
    "joint_delta_other", "joint_case_prop_clean_CDS",
    "joint_control_prop_clean_CDS", "joint_case_prop_leader_uORF",
    "joint_control_prop_leader_uORF",
    "joint_case_prop_overlapping_uORF",
    "joint_control_prop_overlapping_uORF", "joint_case_prop_other",
    "joint_control_prop_other", "joint_evaluable", "joint_strict",
    "joint_top_motif_alignment", "joint_motif_alignment_score",
    "joint_top_motif_prior_signed", "joint_top_motif_prior_strength"
  ), names(joint_branch_rows))
  branch_contexts <- joint_branch_rows[, ..branch_keep]
  card_context <- cards[, .(
    gene_symbol, tx_id, design_family,
    atlas_evidence_class,
    atlas_confidence_score,
    atlas_browser_review_priority_score,
    atlas_clean_cds_pct_change,
    pooling_proof_tier,
    pooling_claim_language,
    pooling_proof_component,
    pooling_hazard_class,
    pooling_hazard_score,
    pooling_hazard_sources,
    pooling_exact_context_fraction,
    pooling_joint_strict_fraction,
    pooling_top_protocol_classes,
    joint_branch_allocation_component,
    context_interaction_gene_class,
    context_interaction_top_axis,
	    context_interaction_top_value,
	    context_interaction_top_branch_class,
	    context_interaction_delta,
	    context_interaction_summary,
	    partial_pooling_gene_class,
	    partial_pooling_top_class,
	    partial_pooling_top_axis,
	    partial_pooling_top_value,
	    partial_pooling_top_branch_class,
	    partial_pooling_top_delta,
	    partial_pooling_top_q,
	    warnings
	  )]
  branch_contexts <- merge(
    branch_contexts,
    card_context,
    by = c("gene_symbol", "tx_id", "design_family"),
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  branch_contexts[, branch_context_label := paste(
    study, CELL_LINE, TISSUE, case_condition,
    paste("vs", control_conditions),
    sep = " | "
  )]
  setorder(branch_contexts, gene_symbol, tx_id, design_family,
           -joint_review_priority, -joint_l1_allocation_shift)
  branch_contexts[, branch_context_rank := seq_len(.N),
                  by = .(gene_symbol, tx_id, design_family)]
  setorder(branch_contexts, gene_symbol, tx_id, design_family,
           branch_context_rank)
  context_preferred <- c(
    "gene_symbol", "tx_id", "design_family", "branch_context_rank",
    "branch_context_label", "study", "AUTHOR", "case_condition",
    "TISSUE", "CELL_LINE", "GENE", "INHIBITOR", "FRACTION",
    "Cancer_type", "Sex", "TIMEPOINT", "control_conditions",
    "match_scope", "sample_set_signature", "n_case_groups",
    "n_control_groups", "n_case_runs_merged", "n_control_runs_merged",
    "case_branch_total_counts", "control_branch_total_counts",
    "joint_review_class", "joint_review_priority",
    "joint_top_branch_class", "joint_top_branch_delta",
    "joint_top_branch_q", "joint_l1_allocation_shift",
    "joint_max_abs_branch_delta", "joint_js_divergence",
    "joint_chisq_q", "joint_delta_clean_CDS",
    "joint_delta_leader_uORF", "joint_delta_overlapping_uORF",
    "joint_delta_other", "joint_case_prop_clean_CDS",
    "joint_control_prop_clean_CDS", "joint_case_prop_leader_uORF",
    "joint_control_prop_leader_uORF",
    "joint_case_prop_overlapping_uORF",
    "joint_control_prop_overlapping_uORF", "joint_case_prop_other",
    "joint_control_prop_other", "joint_evaluable", "joint_strict",
    "joint_top_motif_alignment", "joint_motif_alignment_score",
    "joint_top_motif_prior_signed", "joint_top_motif_prior_strength",
    "atlas_evidence_class", "atlas_confidence_score",
	    "atlas_browser_review_priority_score", "atlas_clean_cds_pct_change",
	    "pooling_proof_tier", "pooling_claim_language",
	    "pooling_proof_component", "pooling_hazard_class",
	    "pooling_hazard_score", "pooling_hazard_sources",
	    "pooling_exact_context_fraction", "pooling_joint_strict_fraction",
	    "pooling_top_protocol_classes",
	    "joint_branch_allocation_component", "context_interaction_gene_class",
	    "context_interaction_top_axis", "context_interaction_top_value",
	    "context_interaction_top_branch_class", "context_interaction_delta",
	    "context_interaction_summary", "partial_pooling_gene_class",
	    "partial_pooling_top_class", "partial_pooling_top_axis",
	    "partial_pooling_top_value", "partial_pooling_top_branch_class",
	    "partial_pooling_top_delta", "partial_pooling_top_q", "warnings"
	  )
  branch_contexts <- branch_contexts[
    , c(intersect(context_preferred, names(branch_contexts)),
        setdiff(names(branch_contexts), context_preferred)),
    with = FALSE
  ]
} else {
  branch_contexts <- data.table(
    gene_symbol = character(),
    tx_id = character(),
    design_family = character(),
    branch_context_rank = integer(),
    branch_context_label = character()
  )
}
fwrite(branch_contexts, branch_context_file)

top_cards <- cards[
  atlas_evidence_class != "exploratory_low_support" |
    atlas_confidence_score >= 0.45 |
    (!is.na(sparse_onoff_score) & sparse_onoff_score >= 0.5)
]
fwrite(top_cards, file.path(output_dir, "dominant_rdg_atlas_cards_review.csv"))

validation_queue <- cards[
  atlas_evidence_class != "exploratory_low_support" |
    atlas_browser_review_priority_score >= 0.55 |
    is.finite(browser_validation_score)
]
setorder(validation_queue, -atlas_browser_review_priority_score,
         -browser_validation_component, -atlas_confidence_score,
         -abs_atlas_clean_cds_pct_change)
fwrite(
  validation_queue,
  file.path(output_dir, "dominant_rdg_atlas_browser_validation_queue.csv")
)

viral_review <- cards[
  design_family == "viral_infection" &
    (atlas_evidence_class != "exploratory_low_support" |
       atlas_browser_review_priority_score >= 0.45 |
       is.finite(browser_validation_score))
]
setorder(viral_review, -atlas_browser_review_priority_score,
         -browser_validation_component, -atlas_confidence_score,
         -abs_atlas_clean_cds_pct_change)
fwrite(
  viral_review,
  file.path(output_dir, "dominant_rdg_atlas_viral_infection_review.csv")
)

summary_metrics <- data.table(
  metric = c(
    "atlas_rows",
    "review_rows",
    "genes",
    "design_families",
    "atlas_ready_replicated_interval_rows",
    "replicated_strict_review_rows",
    "single_design_strict_review_rows",
    "sparse_onoff_review_rows",
    "browser_supported_rows",
    "branch_usage_supported_rows",
    "grouped_branch_usage_supported_rows",
    "grouped_branch_usage_replicated_exact_rows",
    "count_aware_branch_allocation_supported_rows",
    "count_aware_branch_allocation_replicated_rows",
    "viral_count_aware_branch_allocation_supported_rows",
    "viral_count_aware_branch_allocation_replicated_rows",
    "joint_branch_allocation_supported_rows",
    "joint_branch_allocation_strict_rows",
    "viral_joint_branch_allocation_supported_rows",
    "viral_joint_branch_allocation_strict_rows",
    "hierarchical_joint_allocation_supported_rows",
    "hierarchical_joint_allocation_replicated_rows",
    "viral_hierarchical_joint_allocation_supported_rows",
    "viral_hierarchical_joint_allocation_replicated_rows",
    "dirichlet_multinomial_supported_rows",
    "dirichlet_multinomial_replicated_rows",
    "viral_dirichlet_multinomial_supported_rows",
    "viral_dirichlet_multinomial_replicated_rows",
    "model_agreement_consensus_rows",
    "model_agreement_disagreement_rows",
    "viral_model_agreement_consensus_rows",
    "viral_model_agreement_disagreement_rows",
    "loo_consensus_support_robust_rows",
    "loo_consensus_direction_robust_rows",
    "loo_consensus_not_evaluable_rows",
    "loo_consensus_fragile_rows",
    "viral_loo_consensus_support_robust_rows",
    "viral_loo_consensus_direction_robust_rows",
    "viral_loo_consensus_not_evaluable_rows",
    "viral_loo_consensus_fragile_rows",
    "transport_consensus_support_robust_rows",
    "transport_consensus_direction_robust_rows",
    "transport_consensus_not_evaluable_rows",
    "transport_consensus_fragile_rows",
    "viral_transport_consensus_support_robust_rows",
    "viral_transport_consensus_direction_robust_rows",
    "viral_transport_consensus_not_evaluable_rows",
    "viral_transport_consensus_fragile_rows",
    "context_interaction_review_rows",
    "context_direction_reversal_rows",
    "context_specific_supported_rows",
    "context_magnitude_modulated_rows",
    "viral_context_interaction_review_rows",
	    "viral_context_direction_reversal_rows",
	    "viral_context_specific_supported_rows",
	    "viral_context_magnitude_modulated_rows",
	    "partial_pooling_supported_rows",
	    "partial_pooling_review_rows",
	    "partial_pooling_single_context_rows",
	    "partial_pooling_fragile_rows",
	    "viral_partial_pooling_supported_rows",
	    "viral_partial_pooling_review_rows",
	    "viral_partial_pooling_single_context_rows",
	    "viral_partial_pooling_fragile_rows",
	    "hierarchical_atlas_ready_rows",
    "hierarchical_interval_rows",
    "hierarchical_uncalibrated_atlas_ready_rows",
    "hierarchical_uncalibrated_interval_rows",
    "high_hierarchical_residual_rows",
    "viral_high_hierarchical_residual_rows",
    "context_calibration_improved_rows",
    "viral_context_calibration_improved_rows",
    "context_model_improved_rows",
    "viral_context_model_improved_rows",
    "viral_infection_review_rows",
    "median_confidence_score",
    "median_browser_review_priority_score",
    "max_abs_clean_cds_pct_change_review"
  ),
  value = c(
    nrow(cards),
    nrow(top_cards),
    uniqueN(cards$gene_symbol),
    uniqueN(cards$design_family),
    sum(cards$atlas_evidence_class == "atlas_ready_replicated_interval",
        na.rm = TRUE),
    sum(cards$atlas_evidence_class == "replicated_strict_review",
        na.rm = TRUE),
    sum(cards$atlas_evidence_class == "single_design_strict_review",
        na.rm = TRUE),
    sum(cards$atlas_evidence_class == "sparse_onoff_review",
        na.rm = TRUE),
    sum(is.finite(cards$browser_validation_score) &
          cards$browser_validation_score > 0, na.rm = TRUE),
    sum(grepl("branch_usage_shift", cards$branch_usage_top_support_class,
              fixed = TRUE), na.rm = TRUE),
    sum(grepl("group_branch_shift",
              cards$grouped_branch_usage_top_support_class,
              fixed = TRUE), na.rm = TRUE),
    sum(cards$grouped_branch_usage_top_support_class ==
          "replicated_exact_group_branch_shift", na.rm = TRUE),
    sum(grepl("count_aware_branch_shift",
              cards$count_aware_top_support_class,
              fixed = TRUE), na.rm = TRUE),
    sum(cards$count_aware_top_support_class ==
          "replicated_count_aware_branch_shift", na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          grepl("count_aware_branch_shift",
                cards$count_aware_top_support_class,
                fixed = TRUE), na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$count_aware_top_support_class ==
          "replicated_count_aware_branch_shift", na.rm = TRUE),
    sum(cards$joint_top_review_class %chin%
          c("strict_joint_branch_shift",
            "joint_branch_shift_review",
            "motif_joint_discordant_review"),
        na.rm = TRUE),
    sum(cards$joint_top_review_class == "strict_joint_branch_shift",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$joint_top_review_class %chin%
          c("strict_joint_branch_shift",
            "joint_branch_shift_review",
            "motif_joint_discordant_review"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$joint_top_review_class == "strict_joint_branch_shift",
        na.rm = TRUE),
    sum(cards$hierarchical_joint_top_support_class %chin%
          c("hierarchical_joint_replicated_shift",
            "hierarchical_joint_replicated_review",
            "hierarchical_joint_single_study_strong"),
        na.rm = TRUE),
    sum(cards$hierarchical_joint_top_support_class %chin%
          c("hierarchical_joint_replicated_shift",
            "hierarchical_joint_replicated_review"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$hierarchical_joint_top_support_class %chin%
          c("hierarchical_joint_replicated_shift",
            "hierarchical_joint_replicated_review",
            "hierarchical_joint_single_study_strong"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$hierarchical_joint_top_support_class %chin%
          c("hierarchical_joint_replicated_shift",
            "hierarchical_joint_replicated_review"),
        na.rm = TRUE),
    sum(cards$dm_top_support_class %chin%
          c("dm_replicated_shift",
            "dm_replicated_review",
            "dm_single_study_strong"),
        na.rm = TRUE),
    sum(cards$dm_top_support_class %chin%
          c("dm_replicated_shift",
            "dm_replicated_review"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$dm_top_support_class %chin%
          c("dm_replicated_shift",
            "dm_replicated_review",
            "dm_single_study_strong"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$dm_top_support_class %chin%
          c("dm_replicated_shift",
            "dm_replicated_review"),
        na.rm = TRUE),
    sum(cards$model_agreement_tier %chin%
          c("A_browser_validated_consensus",
            "B_replicated_branch_consensus",
            "C_multi_model_consensus"),
        na.rm = TRUE),
    sum(cards$model_agreement_tier == "D_model_disagreement_review",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$model_agreement_tier %chin%
          c("A_browser_validated_consensus",
            "B_replicated_branch_consensus",
            "C_multi_model_consensus"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$model_agreement_tier == "D_model_disagreement_review",
        na.rm = TRUE),
    sum(cards$loo_consensus_class == "consensus_loo_support_robust",
        na.rm = TRUE),
    sum(cards$loo_consensus_class == "consensus_loo_direction_robust",
        na.rm = TRUE),
    sum(cards$loo_consensus_class == "consensus_loo_not_evaluable",
        na.rm = TRUE),
    sum(cards$loo_consensus_class == "consensus_loo_fragile",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$loo_consensus_class == "consensus_loo_support_robust",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$loo_consensus_class == "consensus_loo_direction_robust",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$loo_consensus_class == "consensus_loo_not_evaluable",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$loo_consensus_class == "consensus_loo_fragile",
        na.rm = TRUE),
    sum(cards$transport_consensus_class ==
          "consensus_transport_support_robust", na.rm = TRUE),
    sum(cards$transport_consensus_class ==
          "consensus_transport_direction_robust", na.rm = TRUE),
    sum(cards$transport_consensus_class ==
          "consensus_transport_not_evaluable", na.rm = TRUE),
    sum(cards$transport_consensus_class ==
          "consensus_transport_fragile", na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$transport_consensus_class ==
          "consensus_transport_support_robust", na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$transport_consensus_class ==
          "consensus_transport_direction_robust", na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$transport_consensus_class ==
          "consensus_transport_not_evaluable", na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$transport_consensus_class ==
          "consensus_transport_fragile", na.rm = TRUE),
    sum(cards$context_interaction_gene_class %chin%
          c("context_direction_reversal",
            "context_specific_supported_shift",
            "context_magnitude_modulated"),
        na.rm = TRUE),
    sum(cards$context_interaction_gene_class == "context_direction_reversal",
        na.rm = TRUE),
    sum(cards$context_interaction_gene_class ==
          "context_specific_supported_shift", na.rm = TRUE),
    sum(cards$context_interaction_gene_class ==
          "context_magnitude_modulated", na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$context_interaction_gene_class %chin%
          c("context_direction_reversal",
            "context_specific_supported_shift",
            "context_magnitude_modulated"),
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$context_interaction_gene_class == "context_direction_reversal",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$context_interaction_gene_class ==
          "context_specific_supported_shift", na.rm = TRUE),
	    sum(cards$design_family == "viral_infection" &
	          cards$context_interaction_gene_class ==
	          "context_magnitude_modulated", na.rm = TRUE),
	    sum(cards$partial_pooling_gene_class ==
	          "partial_pooled_context_supported", na.rm = TRUE),
	    sum(cards$partial_pooling_gene_class ==
	          "partial_pooled_context_review", na.rm = TRUE),
	    sum(cards$partial_pooling_gene_class ==
	          "single_context_raw_signal_not_pooled", na.rm = TRUE),
	    sum(cards$partial_pooling_gene_class ==
	          "raw_context_signal_fragile_after_pooling", na.rm = TRUE),
	    sum(cards$design_family == "viral_infection" &
	          cards$partial_pooling_gene_class ==
	          "partial_pooled_context_supported", na.rm = TRUE),
	    sum(cards$design_family == "viral_infection" &
	          cards$partial_pooling_gene_class ==
	          "partial_pooled_context_review", na.rm = TRUE),
	    sum(cards$design_family == "viral_infection" &
	          cards$partial_pooling_gene_class ==
	          "single_context_raw_signal_not_pooled", na.rm = TRUE),
	    sum(cards$design_family == "viral_infection" &
	          cards$partial_pooling_gene_class ==
	          "raw_context_signal_fragile_after_pooling", na.rm = TRUE),
	    sum(cards$hierarchical_evidence_class ==
          "hierarchical_calibrated_atlas_ready",
        na.rm = TRUE),
    sum(cards$hierarchical_calibrated_ci_excludes_zero == TRUE, na.rm = TRUE),
    sum(cards$hierarchical_uncalibrated_evidence_class ==
          "hierarchical_atlas_ready", na.rm = TRUE),
    sum(cards$hierarchical_uncalibrated_ci_excludes_zero == TRUE,
        na.rm = TRUE),
    sum(cards$residual_risk_class == "high_residual_heterogeneity",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$residual_risk_class == "high_residual_heterogeneity",
        na.rm = TRUE),
    sum(cards$context_calibration_status == "context_reduces_residual_error",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$context_calibration_status == "context_reduces_residual_error",
        na.rm = TRUE),
    sum(cards$context_model_status == "context_model_reduces_residual_error",
        na.rm = TRUE),
    sum(cards$design_family == "viral_infection" &
          cards$context_model_status == "context_model_reduces_residual_error",
        na.rm = TRUE),
    nrow(viral_review),
    median_safe(cards$atlas_confidence_score),
    median_safe(cards$atlas_browser_review_priority_score),
    max_safe(abs(top_cards$atlas_clean_cds_pct_change))
  )
)
summary_metrics <- rbind(
  summary_metrics,
  data.table(
    metric = c(
      "pooling_exact_project_matched_control_rows",
      "pooling_exact_or_replicated_branch_model_rows",
      "pooling_pooled_recurrence_or_transport_rows",
      "pooling_context_specific_or_partial_rows",
      "pooling_high_confounding_hazard_rows",
      "pooling_moderate_confounding_hazard_rows",
      "pooling_start_site_protocol_rows",
      "median_pooling_proof_component",
      "median_pooling_hazard_score"
    ),
    value = c(
      sum(cards$pooling_proof_tier ==
            "A_exact_strict_project_matched_control", na.rm = TRUE),
      sum(cards$pooling_proof_tier ==
            "B_exact_or_replicated_branch_model", na.rm = TRUE),
      sum(cards$pooling_proof_tier ==
            "C_pooled_recurrence_or_transport", na.rm = TRUE),
      sum(cards$pooling_proof_tier ==
            "D_context_specific_or_partial_pooled", na.rm = TRUE),
      sum(cards$pooling_hazard_class ==
            "high_pooling_confounding_hazard", na.rm = TRUE),
      sum(cards$pooling_hazard_class ==
            "moderate_pooling_confounding_hazard", na.rm = TRUE),
      sum(cards$pooling_any_start_site_enriched_protocol == TRUE,
          na.rm = TRUE),
      median_safe(cards$pooling_proof_component),
      median_safe(cards$pooling_hazard_score)
    )
  ),
  fill = TRUE
)
fwrite(summary_metrics,
       file.path(output_dir, "dominant_rdg_atlas_summary_metrics.csv"))

if (nrow(cards)) {
  heat_dt <- cards[
    atlas_evidence_class != "exploratory_low_support" |
      atlas_confidence_score >= 0.35
  ]
  if (nrow(heat_dt) > 0) {
    keep_genes <- heat_dt[
      , .(score = max(atlas_confidence_score *
                        pmin(abs(atlas_clean_cds_pct_change) / 100, 2),
                      na.rm = TRUE)),
      by = gene_symbol
    ][order(-score)][seq_len(min(.N, 35)), gene_symbol]
    heat_dt <- heat_dt[gene_symbol %chin% keep_genes]
    heat_dt[, gene_symbol := factor(gene_symbol, levels = rev(keep_genes))]
    family_order <- heat_dt[
      , .(score = max(atlas_confidence_score, na.rm = TRUE)),
      by = design_family
    ][order(-score), design_family]
    heat_dt[, design_family := factor(design_family, levels = family_order)]
    heat_dt[, capped_pct := pmax(pmin(atlas_clean_cds_pct_change, 300), -100)]

    p_heat <- ggplot(
      heat_dt,
      aes(x = design_family, y = gene_symbol, fill = capped_pct,
          text = paste0(
            "Gene: ", gene_symbol,
            "<br>Design: ", design_family,
            "<br>CDS buffering: ", round(atlas_clean_cds_pct_change, 1), "%",
            "<br>Interval: [", round(atlas_clean_cds_pct_lower, 1), ", ",
            round(atlas_clean_cds_pct_upper, 1), "]",
            "<br>Evidence: ", atlas_evidence_class,
            "<br>Confidence: ", round(atlas_confidence_score, 2),
            "<br>Top feature: ", top_shifted_feature_id
          ))
    ) +
      geom_tile(color = "grey88", linewidth = 0.15) +
      geom_point(aes(size = atlas_confidence_score),
                 shape = 21, fill = "white", color = "grey20",
                 stroke = 0.2, alpha = 0.65) +
      scale_fill_gradient2(
        low = "#2b6cb0", mid = "white", high = "#8b0000",
        midpoint = 0, name = "Clean-CDS\nchange (%)",
        limits = c(-100, 300), oob = scales::squish
      ) +
      scale_size_continuous(range = c(0.4, 2.2), name = "Confidence") +
      labs(
        title = "RDG Atlas Cards",
        subtitle = "Clean-CDS buffering by gene and design family; size is atlas confidence",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )

    ggsave(file.path(figure_dir, "dominant_rdg_atlas_buffering_heatmap.png"),
           p_heat, width = 12.5, height = 8.5, dpi = 180)
    ggsave(file.path(figure_dir, "dominant_rdg_atlas_buffering_heatmap.pdf"),
           p_heat, width = 12.5, height = 8.5)
    if (exists("dominant_save_ggplotly") &&
        requireNamespace("plotly", quietly = TRUE)) {
      dominant_save_ggplotly(
        p_heat,
        file.path(figure_dir, "dominant_rdg_atlas_buffering_heatmap.html"),
        title = "RDG atlas buffering heatmap"
      )
    }
  }

  scatter_dt <- cards[is.finite(atlas_clean_cds_pct_change) &
                        is.finite(atlas_confidence_score)]
  if (nrow(scatter_dt) > 0) {
    p_scatter <- ggplot(
      scatter_dt,
      aes(x = atlas_clean_cds_pct_change, y = atlas_confidence_score,
          color = atlas_evidence_class,
          text = paste0(
            "Gene: ", gene_symbol,
            "<br>Design: ", design_family,
            "<br>CDS buffering: ", round(atlas_clean_cds_pct_change, 1), "%",
            "<br>Evidence: ", atlas_evidence_class,
            "<br>Warnings: ", warnings
          ))
    ) +
      geom_vline(xintercept = 0, color = "grey78", linewidth = 0.3) +
      geom_point(alpha = 0.75, size = 2) +
      scale_x_continuous(labels = function(x) paste0(x, "%")) +
      coord_cartesian(xlim = quantile(scatter_dt$atlas_clean_cds_pct_change,
                                      c(0.02, 0.98), na.rm = TRUE)) +
      labs(
        title = "RDG Atlas Confidence Versus Clean-CDS Effect",
        subtitle = "Large effects without replication remain review prompts, not atlas claims",
        x = "Clean-CDS buffering effect",
        y = "Atlas confidence"
      ) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")
    ggsave(file.path(figure_dir, "dominant_rdg_atlas_confidence_vs_effect.png"),
           p_scatter, width = 9.5, height = 6.4, dpi = 180)
    ggsave(file.path(figure_dir, "dominant_rdg_atlas_confidence_vs_effect.pdf"),
           p_scatter, width = 9.5, height = 6.4)
    if (exists("dominant_save_ggplotly") &&
        requireNamespace("plotly", quietly = TRUE)) {
      dominant_save_ggplotly(
        p_scatter,
        file.path(figure_dir, "dominant_rdg_atlas_confidence_vs_effect.html"),
        title = "RDG atlas confidence versus effect"
      )
    }
  }

  viral_plot_dt <- viral_review[
    is.finite(atlas_browser_review_priority_score)
  ][seq_len(min(.N, 35))]
  if (nrow(viral_plot_dt) > 0) {
    viral_plot_dt[, plot_label := paste0(gene_symbol, " / ",
                                         fifelse(
                                           !is.na(top_shifted_feature_id) &
                                             nzchar(top_shifted_feature_id),
                                           top_shifted_feature_id,
                                           "feature?"
                                         ))]
    viral_plot_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
    viral_plot_dt[, browser_supported := is.finite(browser_validation_score) &
                    browser_validation_score > 0]
    viral_plot_dt[, count_aware_supported :=
                    grepl("count_aware_branch_shift",
                          count_aware_top_support_class,
                          fixed = TRUE)]
    p_viral <- ggplot(
      viral_plot_dt,
      aes(x = plot_label, y = atlas_browser_review_priority_score,
          fill = pmax(pmin(atlas_clean_cds_pct_change, 300), -100),
          text = paste0(
            "Gene: ", gene_symbol,
            "<br>Feature: ", top_shifted_feature_id,
            "<br>CDS buffering: ", round(atlas_clean_cds_pct_change, 1), "%",
            "<br>Review priority: ",
            round(atlas_browser_review_priority_score, 2),
            "<br>Evidence: ", atlas_evidence_class,
            "<br>Count-aware support: ", count_aware_top_support_class,
            "<br>Count-aware feature: ", count_aware_top_feature_id,
            "<br>Browser validation: ", browser_validation_status,
            "<br>Warnings: ", warnings
          ))
    ) +
      geom_col(width = 0.75, color = "grey35", linewidth = 0.15) +
      geom_point(
        data = viral_plot_dt[browser_supported == TRUE],
        y = 0.04,
        shape = 21,
        size = 2.2,
        color = "grey10",
        fill = "white"
      ) +
      geom_point(
        data = viral_plot_dt[count_aware_supported == TRUE],
        y = 0.10,
        shape = 23,
        size = 2.4,
        color = "#5b2a7a",
        fill = "#ead7f6"
      ) +
      coord_flip() +
      scale_fill_gradient2(
        low = "#2b6cb0", mid = "white", high = "#8b0000",
        midpoint = 0, limits = c(-100, 300), oob = scales::squish,
        name = "Clean-CDS\nchange (%)"
      ) +
      labs(
        title = "Viral-Infection RDG Atlas Review Priority",
        subtitle = "Circle = browser supported; purple diamond = count-aware branch-allocation support",
        x = NULL,
        y = "Review priority score"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )
    ggsave(file.path(figure_dir, "dominant_rdg_atlas_viral_review_priority.png"),
           p_viral, width = 10.2, height = 8.4, dpi = 180)
    ggsave(file.path(figure_dir, "dominant_rdg_atlas_viral_review_priority.pdf"),
           p_viral, width = 10.2, height = 8.4)
    if (exists("dominant_save_ggplotly") &&
        requireNamespace("plotly", quietly = TRUE)) {
      dominant_save_ggplotly(
        p_viral,
        file.path(figure_dir, "dominant_rdg_atlas_viral_review_priority.html"),
        title = "RDG atlas viral review priority"
      )
    }
  }
}

message("Saved RDG atlas cards: ", cards_file)
message("Saved RDG atlas review cards: ",
        file.path(output_dir, "dominant_rdg_atlas_cards_review.csv"))
message("Saved RDG atlas browser-validation queue: ",
        file.path(output_dir, "dominant_rdg_atlas_browser_validation_queue.csv"))
message("Saved RDG atlas viral-infection review: ",
        file.path(output_dir, "dominant_rdg_atlas_viral_infection_review.csv"))
message("Saved RDG atlas summary metrics: ",
        file.path(output_dir, "dominant_rdg_atlas_summary_metrics.csv"))
