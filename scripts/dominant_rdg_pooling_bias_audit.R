#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) && dir.exists(file.path(candidate, "scripts"))) {
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
output_dir <- file.path(analysis_dir, "dominant_rdg_pooling_bias_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

condition_helper <- file.path(analysis_dir, "scripts",
                              "dominant_condition_families.R")
if (file.exists(condition_helper)) source(condition_helper)

message("RDG pooling and protocol-bias audit")
message("  1. Formalize evidence tiers for pooled branch claims")
message("  2. Quantify coverage gained by merged replicate groups")
message("  3. Flag what pooling loses: relaxed controls, context shrinkage, single-run dominance")
message("  4. Summarize inhibitor/protocol and condition-family bias hazards")

read_dt <- function(path, required = FALSE, select = NULL) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing required input: ", path, call. = FALSE)
    return(data.table())
  }
  if (is.null(select)) return(fread(path, showProgress = FALSE))
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  keep <- intersect(select, header)
  if (!length(keep)) return(data.table())
  fread(path, select = keep, showProgress = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_median <- function(x, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else median(x)
}

safe_mean <- function(x, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else mean(x)
}

safe_quantile <- function(x, prob = 0.90, default = NA_real_) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE,
                             type = 8))
}

safe_fraction <- function(x) {
  if (!length(x)) return(NA_real_)
  mean(x, na.rm = TRUE)
}

clip01 <- function(x) pmin(pmax(safe_num(x), 0), 1)

norm_label <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[is.na(y)] <- ""
  y <- gsub("[^a-z0-9]+", "_", y)
  gsub("^_+|_+$", "", y)
}

collapse_top <- function(x, n = 5L) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & !x %chin% c("MISSING", "NA")]
  if (!length(x)) return("")
  tab <- sort(table(x), decreasing = TRUE)
  paste(names(tab)[seq_len(min(length(tab), n))], collapse = ";")
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

num_col <- function(dt, column, default = NA_real_) {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- safe_num(dt[[column]])
  x[!is.finite(x)] <- default
  x
}

chr_col <- function(dt, column, default = "") {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- as.character(dt[[column]])
  x[is.na(x)] <- default
  x
}

metric_lookup <- function(dt, name, default = NA_real_) {
  if (!nrow(dt) || !all(c("metric", "value") %in% names(dt))) return(default)
  x <- safe_num(dt[metric == name, value][1])
  if (!length(x) || !is.finite(x)) default else x
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

condition_bias_family <- function(condition, fraction = "", inhibitor = "",
                                  cell_line = "", tissue = "", gene = "",
                                  cancer_type = "") {
  if (exists("dominant_classify_design_family", mode = "function")) {
    return(dominant_classify_design_family(
      condition = condition,
      fraction = fraction,
      inhibitor = inhibitor,
      cell_line = cell_line,
      tissue = tissue,
      gene = gene,
      cancer_type = cancer_type
    ))
  }
  text <- norm_label(paste(condition, fraction, inhibitor, cell_line, tissue,
                           gene, cancer_type, sep = " "))
  fcase(
    control_label_like(condition), "control",
    grepl("infect|virus|viral|sars|covid|influenza|hiv|hcv|rsv", text),
    "viral_infection",
    grepl("starv|depriv|amino_acid|no_b?caa|nobcaa|nutrient|serum_free",
          text),
    "nutrient_starvation",
    grepl("lactimidomycin|harringtonine|chx|cycloheximide|puro|arsen|er_stress|upr",
          text),
    "isr_er_translation_stress",
    default = "other_non_control"
  )
}

claim_language <- data.table(
  claim_tier = c(
    "A_project_specific_matched_contrast",
    "B_merged_replicate_measurement_support",
    "C_partial_pooled_context_effect",
    "D_design_family_recurrence_screen",
    "E_single_context_evidence",
    "F_all_merged_visual_background"
  ),
  allowed_language = c(
    "condition-specific RDG effect candidate",
    "merged-replicate branch measurement support",
    "context-modulated RDG effect candidate",
    "pooled design-family recurrence screen",
    "single-context branch evidence",
    "visual background only"
  ),
  minimum_evidence = c(
    "Exact project/study/cell/tissue/protocol/time matched case-control contrast with merged replicates and branch count gates.",
    "Same-condition replicate group improves branch counts; not itself a biological effect.",
    "Context-interaction row survives empirical-Bayes shrinkage across sibling contexts.",
    "Signal recurs after study-aware aggregation, model agreement, LOO/grouped omission, or design-family pooling.",
    "Strong local row lacks sibling contexts or shrinks poorly; inspect but do not generalize.",
    "All-run or all-condition profile used to see transcript shape; never used as a condition-specific proof."
  ),
  primary_loss = c(
    "Still not independent validation; can retain local protocol artifacts.",
    "Replicate heterogeneity and single-run dominance are compressed.",
    "Sibling contexts may still share protocol, lab, tissue, or virus biases.",
    "Condition, protocol, tissue, cell line, time point, and publication effects can be mixed.",
    "No transportability estimate; sensitive to one context or sample family.",
    "Condition labels and matched denominators are absent."
  ),
  recommended_action = c(
    "Prioritize browser validation and experimental replication.",
    "Use for branch count estimation; audit run dominance before claiming replication.",
    "Inspect matching context and complement profiles in the browser.",
    "Use as a candidate generator; seek exact WT/control evidence next.",
    "Treat as a targeted review lead or acquisition need.",
    "Use only as a coverage/shape reference."
  )
)
fwrite(claim_language,
       file.path(output_dir, "pooling_claim_language.csv"))

normalization_limits <- data.table(
  quantity = c(
    "branch_allocation_delta",
    "clean_cds_log2fc_or_pct_change",
    "merged_replicate_counts",
    "partial_pooled_context_delta",
    "all_merged_profile"
  ),
  normalized_or_controlled_by = c(
    "Feature counts divided by same gene/transcript branch denominator.",
    "Clean-CDS expression normalization and matched/control-aware log2FC summaries.",
    "Summed counts inside a same-condition metadata stratum.",
    "Empirical-Bayes shrinkage across context values within gene/design/branch/axis.",
    "No condition-specific normalization; all selected runs are pooled for shape."
  ),
  controls_well = c(
    "Library depth and gene-level total load within the allocation denominator.",
    "Depth and broad clean-CDS expression scale.",
    "Short-uORF sparsity and low per-run counts.",
    "Some single-context noise and weak sibling-context outliers.",
    "Transcript coordinate visibility and qualitative branch shape."
  ),
  does_not_control = c(
    "Inhibitor start peaks, footprint length, p-shift, isoform errors, and condition-wide biology.",
    "Protocol/inhibitor bias, host shutoff, stalling, or total ribosome-load compensation.",
    "Run heterogeneity, single-run dominance, or mixed biological states inside the stratum.",
    "Causal tissue/protocol separation when sibling contexts share confounders.",
    "Any case/control or condition-specific effect."
  ),
  audit_output = c(
    "pooling_contrast_tier_summary.csv",
    "design_family_pooling_bias_summary.csv",
    "replicate_pooling_loss_summary.csv",
    "partial_pooling_tradeoff_summary.csv",
    "pooling_claim_language.csv"
  )
)
fwrite(normalization_limits,
       file.path(output_dir, "normalization_limits.csv"))

coverage_metrics <- read_dt(file.path(
  analysis_dir, "dominant_rdg_branch_coverage_qc",
  "branch_coverage_summary_metrics.csv"
), required = TRUE)
coverage_gate <- read_dt(file.path(
  analysis_dir, "dominant_rdg_branch_coverage_qc",
  "branch_coverage_gate_summary.csv"
))

single_small_ge10 <- metric_lookup(coverage_metrics,
                                   "single_small_uorf_fraction_ge_10")
single_small_ge30 <- metric_lookup(coverage_metrics,
                                   "single_small_uorf_fraction_ge_30")
merged_small_ge10 <- metric_lookup(coverage_metrics,
                                   "merged_n2_small_uorf_fraction_ge_10")
merged_small_ge30 <- metric_lookup(coverage_metrics,
                                   "merged_n2_small_uorf_fraction_ge_30")
single_min_gate <- metric_lookup(
  coverage_metrics,
  "single_allocation_feature_minimal_gate_fraction"
)
merged_min_gate <- metric_lookup(
  coverage_metrics,
  "merged_n2_allocation_feature_minimal_gate_fraction"
)

pooling_gain_loss <- data.table(
  domain = c(
    "small_uorf_coverage", "small_uorf_coverage",
    "small_uorf_coverage", "small_uorf_coverage",
    "allocation_gate", "allocation_gate", "allocation_gate",
    "replicate_group_scale", "replicate_group_scale"
  ),
  metric = c(
    "single_run_fraction_ge_10",
    "merged_n2_fraction_ge_10",
    "single_run_fraction_ge_30",
    "merged_n2_fraction_ge_30",
    "single_run_minimal_gate_fraction",
    "merged_n2_minimal_gate_fraction",
    "merged_minus_single_minimal_gate",
    "merged_condition_groups",
    "merged_condition_groups_n2plus"
  ),
  value_numeric = c(
    single_small_ge10,
    merged_small_ge10,
    single_small_ge30,
    merged_small_ge30,
    single_min_gate,
    merged_min_gate,
    merged_min_gate - single_min_gate,
    metric_lookup(coverage_metrics, "merged_condition_groups"),
    metric_lookup(coverage_metrics, "merged_condition_groups_n2plus")
  ),
  interpretation = c(
    "Single runs with short-uORF support at >=10 reads.",
    "Merged same-condition groups with short-uORF support at >=10 reads.",
    "Single runs with short-uORF support at the stricter 30-read gate.",
    "Merged same-condition groups with short-uORF support at the stricter 30-read gate.",
    "Single-run allocation rows passing the current minimal branch model gate.",
    "Merged same-condition n>=2 rows passing the current minimal branch model gate.",
    "Coverage gain from replicate pooling; positive values support merged-replicate-first modeling.",
    "All same-condition metadata groups available for allocation features.",
    "Same-condition metadata groups with at least two runs."
  )
)
fwrite(pooling_gain_loss,
       file.path(output_dir, "pooling_gain_loss_summary.csv"))

grouped_cols <- c(
  "stratum_id", "case_condition", "perturbation_class", "design_family",
  "n_case_groups", "n_case_runs_merged", "case_branch_total_counts",
  "case_sum_raw_counts", "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
  "n_control_groups", "n_control_runs_merged", "control_branch_total_counts",
  "control_sum_raw_counts", "control_conditions", "grouped_usage_delta",
  "grouped_branch_usage_log_or", "grouped_branch_usage_log_or_se",
  "grouped_branch_usage_p", "match_scope", "gene_symbol", "tx_id",
  "feature_id", "branch_feature_class", "grouped_branch_usage_q",
  "grouped_min_runs_merged", "grouped_min_branch_total_counts",
  "grouped_max_feature_counts", "grouped_min_feature_counts",
  "grouped_branch_usage_evaluable", "grouped_branch_usage_strict",
  "grouped_branch_usage_relaxed"
)
grouped <- read_dt(file.path(
  analysis_dir, "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_contrast_table.csv"
), required = TRUE, select = grouped_cols)

if (nrow(grouped)) {
  grouped[, grouped_branch_usage_evaluable :=
            bool_col(.SD, "grouped_branch_usage_evaluable")]
  grouped[, grouped_branch_usage_strict :=
            bool_col(.SD, "grouped_branch_usage_strict")]
  grouped[, grouped_branch_usage_relaxed :=
            bool_col(.SD, "grouped_branch_usage_relaxed")]
  grouped[, grouped_usage_delta := num_col(.SD, "grouped_usage_delta")]
  grouped[, grouped_branch_usage_q := num_col(.SD, "grouped_branch_usage_q")]
  grouped[, control_is_wt_like := control_label_like(control_conditions)]
  grouped[, proof_tier := fcase(
    match_scope == "exact_context" &
      grouped_branch_usage_strict == TRUE &
      control_is_wt_like == TRUE,
    "A_exact_strict_project_matched_control",
    match_scope == "exact_context" &
      grouped_branch_usage_evaluable == TRUE,
    "B_exact_project_evaluable",
    match_scope != "exact_context" &
      grouped_branch_usage_evaluable == TRUE,
    "C_relaxed_study_cell_gene_pool",
    default = "D_low_count_or_unmatched"
  )]
  grouped[, strong_evaluable_signal :=
            grouped_branch_usage_evaluable == TRUE &
            is.finite(grouped_branch_usage_q) &
            grouped_branch_usage_q <= 0.10 &
            abs(grouped_usage_delta) >= 0.05]
  contrast_tier_summary <- grouped[, .(
    contrast_rows = .N,
    genes = uniqueN(gene_symbol),
    transcripts = uniqueN(tx_id),
    studies = uniqueN(study),
    exact_context_rows = sum(match_scope == "exact_context", na.rm = TRUE),
    relaxed_context_rows = sum(match_scope != "exact_context", na.rm = TRUE),
    evaluable_rows = sum(grouped_branch_usage_evaluable == TRUE, na.rm = TRUE),
    strict_rows = sum(grouped_branch_usage_strict == TRUE, na.rm = TRUE),
    strong_signal_rows = sum(strong_evaluable_signal == TRUE, na.rm = TRUE),
    wt_like_control_rows = sum(control_is_wt_like == TRUE, na.rm = TRUE),
    median_abs_usage_delta = safe_median(abs(grouped_usage_delta)),
    q90_abs_usage_delta = safe_quantile(abs(grouped_usage_delta), 0.90),
    top_case_conditions = collapse_top(case_condition),
    top_inhibitors = collapse_top(INHIBITOR),
    top_tissues = collapse_top(TISSUE),
    top_cell_lines = collapse_top(CELL_LINE)
  ), by = .(proof_tier, design_family, branch_feature_class)]
  setorder(contrast_tier_summary, proof_tier, design_family,
           branch_feature_class)
} else {
  contrast_tier_summary <- data.table()
}
fwrite(contrast_tier_summary,
       file.path(output_dir, "pooling_contrast_tier_summary.csv"))

grouped_review <- read_dt(file.path(
  analysis_dir, "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_review.csv"
))
if (nrow(grouped_review)) {
  grouped_review_summary <- grouped_review[, .(
    review_rows = .N,
    genes = uniqueN(gene_symbol),
    median_score = safe_median(grouped_branch_usage_top_score),
    max_score = max(num_col(.SD, "grouped_branch_usage_top_score"),
                    na.rm = TRUE),
    median_abs_top_delta =
      safe_median(abs(grouped_branch_usage_top_weighted_usage_delta)),
    exact_strict_top_rows =
      sum(grouped_branch_usage_top_strict_rows > 0, na.rm = TRUE),
    top_genes = collapse_top(gene_symbol, 8L)
  ), by = .(design_family, grouped_branch_usage_top_support_class,
            grouped_branch_usage_top_feature_class)]
  setorder(grouped_review_summary, -review_rows, design_family)
} else {
  grouped_review_summary <- data.table()
}
fwrite(grouped_review_summary,
       file.path(output_dir, "grouped_review_support_summary.csv"))

branch_context_cols <- c(
  "gene_symbol", "tx_id", "design_family", "branch_context_rank", "study",
  "case_condition", "control_conditions", "TISSUE", "CELL_LINE", "GENE",
  "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
  "match_scope", "n_case_runs_merged", "n_control_runs_merged",
  "case_branch_total_counts", "control_branch_total_counts",
  "joint_review_class", "joint_review_priority", "joint_top_branch_class",
  "joint_top_branch_delta", "joint_l1_allocation_shift", "joint_strict"
)
branch_contexts <- read_dt(file.path(
  analysis_dir, "dominant_rdg_atlas",
  "dominant_rdg_atlas_branch_contexts.csv"
), select = branch_context_cols)

if (nrow(branch_contexts)) {
  branch_contexts[, joint_strict := bool_col(.SD, "joint_strict")]
  branch_contexts[, protocol_bias_class :=
                    protocol_bias_class(INHIBITOR, FRACTION, case_condition)]
  branch_contexts[, joint_top_branch_delta :=
                    num_col(.SD, "joint_top_branch_delta")]
  branch_contexts[, joint_l1_allocation_shift :=
                    num_col(.SD, "joint_l1_allocation_shift")]
  design_family_pooling_bias <- branch_contexts[, {
    exact_fraction <- mean(match_scope == "exact_context", na.rm = TRUE)
    strict_fraction <- mean(joint_strict == TRUE, na.rm = TRUE)
    n_inhibitors <- uniqueN(INHIBITOR[!is.na(INHIBITOR) &
                                       INHIBITOR != "MISSING"])
    n_tissues <- uniqueN(TISSUE[!is.na(TISSUE) & TISSUE != "MISSING"])
    n_cell_lines <- uniqueN(CELL_LINE[!is.na(CELL_LINE) &
                                       CELL_LINE != "MISSING"])
    n_conditions <- uniqueN(case_condition[!is.na(case_condition) &
                                             case_condition != "MISSING"])
    hazard <- clip01(
      0.25 * (1 - exact_fraction) +
        0.15 * clip01((n_inhibitors - 1) / 3) +
        0.15 * clip01((n_tissues - 1) / 6) +
        0.15 * clip01((n_cell_lines - 1) / 6) +
        0.15 * clip01((n_conditions - 1) / 8) +
        0.15 * (1 - strict_fraction)
    )
    .(
      branch_context_rows = .N,
      gene_designs = uniqueN(paste(gene_symbol, tx_id, sep = "|")),
      studies = uniqueN(study),
      exact_context_rows = sum(match_scope == "exact_context", na.rm = TRUE),
      relaxed_context_rows = sum(match_scope != "exact_context", na.rm = TRUE),
      exact_context_fraction = exact_fraction,
      joint_strict_fraction = strict_fraction,
      inhibitors = n_inhibitors,
      tissues = n_tissues,
      cell_lines = n_cell_lines,
      case_conditions = n_conditions,
      protocol_classes = uniqueN(protocol_bias_class),
      median_case_runs = safe_median(n_case_runs_merged),
      median_control_runs = safe_median(n_control_runs_merged),
      median_abs_branch_delta = safe_median(abs(joint_top_branch_delta)),
      median_l1_allocation_shift = safe_median(joint_l1_allocation_shift),
      pooling_hazard_score = hazard,
      pooling_hazard_class = fcase(
        hazard >= 0.65, "high_pooling_confounding_hazard",
        hazard >= 0.40, "moderate_pooling_confounding_hazard",
        default = "lower_pooling_confounding_hazard"
      ),
      top_case_conditions = collapse_top(case_condition),
      top_inhibitors = collapse_top(INHIBITOR),
      top_protocol_classes = collapse_top(protocol_bias_class),
      top_tissues = collapse_top(TISSUE),
      top_cell_lines = collapse_top(CELL_LINE)
    )
  }, by = design_family]
  setorder(design_family_pooling_bias, -pooling_hazard_score,
           -branch_context_rows)

  protocol_context_summary <- branch_contexts[, .(
    branch_context_rows = .N,
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family, sep = "|")),
    studies = uniqueN(study),
    exact_context_fraction = mean(match_scope == "exact_context",
                                  na.rm = TRUE),
    joint_strict_fraction = mean(joint_strict == TRUE, na.rm = TRUE),
    median_abs_branch_delta = safe_median(abs(joint_top_branch_delta)),
    median_l1_allocation_shift = safe_median(joint_l1_allocation_shift),
    top_genes = collapse_top(gene_symbol, 8L),
    top_case_conditions = collapse_top(case_condition),
    top_inhibitors = collapse_top(INHIBITOR),
    top_tissues = collapse_top(TISSUE)
  ), by = .(design_family, protocol_bias_class)]
  setorder(protocol_context_summary, design_family, -branch_context_rows)
} else {
  design_family_pooling_bias <- data.table()
  protocol_context_summary <- data.table()
}
fwrite(design_family_pooling_bias,
       file.path(output_dir, "design_family_pooling_bias_summary.csv"))
fwrite(protocol_context_summary,
       file.path(output_dir, "protocol_context_bias_summary.csv"))

partial_cols <- c(
  "gene_symbol", "tx_id", "design_family", "branch_class", "context_axis",
  "context_value", "context_studies", "complement_studies",
  "interaction_delta", "interaction_delta_q", "raw_context_supported",
  "partial_pooling_weight", "pooled_interaction_delta",
  "pooled_interaction_q", "pooled_abs_delta", "shrinkage_fraction",
  "partial_pooled_context_class", "partial_pooling_survival_class",
  "partial_pooling_review_priority", "context_top_labels",
  "complement_top_labels"
)
partial_rows <- read_dt(file.path(
  analysis_dir, "dominant_rdg_context_partial_pooling",
  "dominant_rdg_context_partial_pooling_rows.csv"
), select = partial_cols)

if (nrow(partial_rows)) {
  partial_rows[, raw_context_supported :=
                 bool_col(.SD, "raw_context_supported")]
  for (column in c("interaction_delta", "interaction_delta_q",
                   "partial_pooling_weight", "pooled_interaction_delta",
                   "pooled_interaction_q", "pooled_abs_delta",
                   "shrinkage_fraction", "partial_pooling_review_priority")) {
    partial_rows[, (column) := num_col(.SD, column)]
  }
  partial_pooling_tradeoff <- partial_rows[, .(
    context_rows = .N,
    genes = uniqueN(gene_symbol),
    context_values = uniqueN(context_value),
    raw_supported_rows = sum(raw_context_supported == TRUE, na.rm = TRUE),
    pooled_supported_rows = sum(partial_pooled_context_class %chin% c(
      "partial_pooled_direction_reversal",
      "partial_pooled_supported_context_shift"
    ), na.rm = TRUE),
    pooled_review_rows = sum(partial_pooled_context_class ==
                               "partial_pooled_review", na.rm = TRUE),
    single_context_rows = sum(partial_pooling_survival_class ==
                                "single_context_raw_signal_not_pooled",
                              na.rm = TRUE),
    raw_shrunk_rows = sum(partial_pooling_survival_class ==
                            "raw_signal_shrunk_by_partial_pooling",
                          na.rm = TRUE),
    pooled_emergent_rows = sum(partial_pooling_survival_class ==
                                 "pooled_context_emerges", na.rm = TRUE),
    median_raw_abs_delta = safe_median(abs(interaction_delta)),
    median_pooled_abs_delta = safe_median(abs(pooled_interaction_delta)),
    median_shrinkage_fraction = safe_median(shrinkage_fraction),
    median_context_studies = safe_median(context_studies),
    median_complement_studies = safe_median(complement_studies),
    top_genes = collapse_top(gene_symbol, 8L),
    top_context_values = collapse_top(context_value)
  ), by = .(design_family, context_axis, branch_class)]
  setorder(partial_pooling_tradeoff, -raw_supported_rows,
           -pooled_supported_rows)

  partial_pooling_review <- partial_rows[
    raw_context_supported == TRUE |
      partial_pooled_context_class != "partial_pooled_neutral"
  ][
    order(-partial_pooling_review_priority)
  ][seq_len(min(.N, 300L))]
} else {
  partial_pooling_tradeoff <- data.table()
  partial_pooling_review <- data.table()
}
fwrite(partial_pooling_tradeoff,
       file.path(output_dir, "partial_pooling_tradeoff_summary.csv"))
fwrite(partial_pooling_review,
       file.path(output_dir, "partial_pooling_pool_loss_review.csv"))

replicate_summary <- read_dt(file.path(
  analysis_dir, "dominant_rdg_qualified_inspection",
  "qualified_replicate_summary.csv"
))
if (nrow(replicate_summary)) {
  replicate_summary[, exact_max_run_fraction :=
                      num_col(.SD, "exact_max_run_fraction")]
  replicate_summary[, exact_signal_run_fraction :=
                      num_col(.SD, "exact_signal_run_fraction")]
  replicate_pooling_loss <- replicate_summary[, .(
    profile_groups = .N,
    median_exact_run_count = safe_median(exact_run_count),
    median_max_run_fraction = safe_median(exact_max_run_fraction),
    q90_max_run_fraction = safe_quantile(exact_max_run_fraction, 0.90),
    fraction_single_run_dominant_ge_0_80 =
      mean(exact_max_run_fraction >= 0.80, na.rm = TRUE),
    fraction_single_run_dominant_ge_0_95 =
      mean(exact_max_run_fraction >= 0.95, na.rm = TRUE),
    median_signal_run_fraction = safe_median(exact_signal_run_fraction),
    fraction_half_or_fewer_runs_with_signal =
      mean(exact_signal_run_fraction <= 0.50, na.rm = TRUE),
    median_total_top_feature_counts =
      safe_median(exact_total_top_feature_counts),
    top_group_labels = collapse_top(group_label, 6L)
  ), by = .(group_role)]
} else {
  replicate_pooling_loss <- data.table()
}
fwrite(replicate_pooling_loss,
       file.path(output_dir, "replicate_pooling_loss_summary.csv"))

frame_file <- file.path(
  analysis_dir, "dominant_cds_frame_bumpiness_qc",
  "cds_frame_bumpiness_run_metrics.csv"
)
frame_cols <- c(
  "Run", "study", "CONDITION", "INHIBITOR", "FRACTION", "TISSUE",
  "CELL_LINE", "GENE", "Cancer_type", "tis_peak_strength",
  "global_canonical_frame_fraction", "global_possible_misshift",
  "metadata_possible_misshift", "top_readlength"
)
frame_rows <- read_dt(frame_file, select = frame_cols)
if (nrow(frame_rows)) {
  frame_runs <- unique(frame_rows, by = "Run")
  frame_runs[, protocol_bias_class :=
               protocol_bias_class(INHIBITOR, FRACTION, CONDITION)]
  frame_runs[, condition_bias_family := condition_bias_family(
    CONDITION, FRACTION, INHIBITOR, CELL_LINE, TISSUE, GENE, Cancer_type
  )]
  frame_runs[, global_possible_misshift :=
               bool_col(.SD, "global_possible_misshift")]
  frame_runs[, metadata_possible_misshift :=
               bool_col(.SD, "metadata_possible_misshift")]
  frame_runs[, tis_peak_strength := num_col(.SD, "tis_peak_strength")]
  frame_runs[, global_canonical_frame_fraction :=
               num_col(.SD, "global_canonical_frame_fraction")]
  protocol_bias_run_summary <- frame_runs[, .(
    runs = .N,
    studies = uniqueN(study),
    median_tis_peak_strength = safe_median(tis_peak_strength),
    q90_tis_peak_strength = safe_quantile(tis_peak_strength, 0.90),
    median_global_canonical_frame_fraction =
      safe_median(global_canonical_frame_fraction),
    q10_global_canonical_frame_fraction =
      safe_quantile(global_canonical_frame_fraction, 0.10),
    global_misshift_run_fraction =
      mean(global_possible_misshift == TRUE, na.rm = TRUE),
    metadata_misshift_run_fraction =
      mean(metadata_possible_misshift == TRUE, na.rm = TRUE),
    top_readlengths = collapse_top(top_readlength),
    top_conditions = collapse_top(CONDITION),
    top_inhibitors = collapse_top(INHIBITOR),
    top_fractions = collapse_top(FRACTION),
    top_tissues = collapse_top(TISSUE),
    top_cell_lines = collapse_top(CELL_LINE)
  ), by = .(protocol_bias_class, condition_bias_family)]
  setorder(protocol_bias_run_summary, -runs)
} else {
  frame_runs <- data.table()
  protocol_bias_run_summary <- data.table()
}
fwrite(protocol_bias_run_summary,
       file.path(output_dir, "protocol_bias_run_summary.csv"))

condition_protocol_hypotheses <- data.table(
  condition_or_protocol = c(
    "lactimidomycin_or_harringtonine",
    "cycloheximide",
    "viral_infection",
    "amino_acid_or_nutrient_starvation",
    "fractionated_or_polysome_protocol",
    "all_merged_background"
  ),
  expected_bias_or_biology = c(
    "Start-site enriched footprints can exaggerate initiation peaks and upstream branch starts.",
    "Elongation freezing can preserve CDS/translon footprints differently from start-enriched protocols.",
    "Host shutoff, immune induction, tissue composition, virus type, and time point can shift the whole translation distribution.",
    "Stalling, ribosome collision, GCN2/ISR activation, and altered elongation can create branch-like peaks.",
    "Different ribosome pools are sampled, so branch allocation may reflect fraction selection.",
    "Useful for transcript shape, not for condition-specific effects."
  ),
  audit_signal = c(
    "High TIS peak strength; protocol class differs between case and control; start-branch signal lacks matched protocol control.",
    "Protocol class dominates the evidence family or frame/QC metrics differ between groups.",
    "High context-interaction or grouped-omission sensitivity; lung/cell-line/tissue effects differ from complements.",
    "Nutrient-starvation branch effects with high bumpiness/frame/stalling risk or allocation-output disagreement.",
    "FRACTION terms differ between case and control or dominate a design family.",
    "No matched case/control denominator."
  ),
  claim_constraint = c(
    "Require same-inhibitor WT/control and browser profile before calling an RDG effect.",
    "Do not compare directly to start-enriched protocols without a protocol covariate or exact match.",
    "Call recurrent viral or context-modulated unless exact local controls and transport tests agree.",
    "Separate initiation-control claims from stalling/collision hypotheses.",
    "Treat as protocol-conditioned evidence, not a universal biological state.",
    "Never use as proof of condition direction."
  )
)
fwrite(condition_protocol_hypotheses,
       file.path(output_dir, "condition_protocol_hypotheses.csv"))

summary_metrics <- rbindlist(list(
  data.table(metric = "single_small_uorf_fraction_ge_30",
             value = single_small_ge30),
  data.table(metric = "merged_n2_small_uorf_fraction_ge_30",
             value = merged_small_ge30),
  data.table(metric = "single_allocation_minimal_gate_fraction",
             value = single_min_gate),
  data.table(metric = "merged_n2_allocation_minimal_gate_fraction",
             value = merged_min_gate),
  data.table(metric = "grouped_contrast_rows",
             value = nrow(grouped)),
  data.table(metric = "exact_strict_project_matched_rows",
             value = if (nrow(grouped)) {
               sum(grouped$proof_tier ==
                     "A_exact_strict_project_matched_control", na.rm = TRUE)
             } else 0),
  data.table(metric = "relaxed_evaluable_pool_rows",
             value = if (nrow(grouped)) {
               sum(grouped$proof_tier ==
                     "C_relaxed_study_cell_gene_pool", na.rm = TRUE)
             } else 0),
  data.table(metric = "partial_pooling_rows",
             value = nrow(partial_rows)),
  data.table(metric = "partial_pooling_pooled_supported_rows",
             value = if (nrow(partial_rows)) {
               sum(partial_rows$partial_pooled_context_class %chin% c(
                 "partial_pooled_direction_reversal",
                 "partial_pooled_supported_context_shift"
               ), na.rm = TRUE)
             } else 0),
  data.table(metric = "partial_pooling_single_context_rows",
             value = if (nrow(partial_rows)) {
               sum(partial_rows$partial_pooling_survival_class ==
                     "single_context_raw_signal_not_pooled", na.rm = TRUE)
             } else 0),
  data.table(metric = "partial_pooling_raw_shrunk_rows",
             value = if (nrow(partial_rows)) {
               sum(partial_rows$partial_pooling_survival_class ==
                     "raw_signal_shrunk_by_partial_pooling", na.rm = TRUE)
             } else 0),
  data.table(metric = "protocol_run_rows",
             value = nrow(frame_runs)),
  data.table(metric = "start_site_enriched_protocol_runs",
             value = if (nrow(frame_runs)) {
               sum(frame_runs$protocol_bias_class %chin% c(
                 "start_site_enriched_ltm",
                 "start_site_enriched_harringtonine"
               ), na.rm = TRUE)
             } else 0),
  data.table(metric = "high_pooling_hazard_design_families",
             value = if (nrow(design_family_pooling_bias)) {
               sum(design_family_pooling_bias$pooling_hazard_class ==
                     "high_pooling_confounding_hazard", na.rm = TRUE)
             } else 0)
), fill = TRUE)
summary_metrics[, value := safe_num(value)]
fwrite(summary_metrics,
       file.path(output_dir, "pooling_bias_summary_metrics.csv"))

fmt_pct <- function(x) {
  if (!is.finite(x)) return("NA")
  sprintf("%.1f%%", 100 * x)
}

metric_value <- function(name) {
  summary_metrics[metric == name, value][1]
}

top_hazard <- if (nrow(design_family_pooling_bias)) {
  head(design_family_pooling_bias[
    order(-pooling_hazard_score),
    .(design_family, pooling_hazard_class, pooling_hazard_score,
      branch_context_rows, top_protocol_classes, top_tissues,
      top_case_conditions)
  ], 8L)
} else {
  data.table()
}

top_partial <- if (nrow(partial_pooling_review)) {
  head(partial_pooling_review[, .(
    gene_symbol, design_family, branch_class, context_axis, context_value,
    interaction_delta, pooled_interaction_delta,
    partial_pooled_context_class, partial_pooling_survival_class
  )], 8L)
} else {
  data.table()
}

report <- c(
  "# RDG Pooling And Protocol-Bias Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Core Interpretation",
  "",
  paste0(
    "- Short-uORF >=30-read support improves from ",
    fmt_pct(metric_value("single_small_uorf_fraction_ge_30")),
    " in single runs to ",
    fmt_pct(metric_value("merged_n2_small_uorf_fraction_ge_30")),
    " in merged n>=2 replicate groups."
  ),
  paste0(
    "- The allocation minimal gate improves from ",
    fmt_pct(metric_value("single_allocation_minimal_gate_fraction")),
    " to ",
    fmt_pct(metric_value("merged_n2_allocation_minimal_gate_fraction")),
    " after merged-replicate pooling."
  ),
  paste0(
    "- Grouped branch contrasts contain ",
    format(metric_value("grouped_contrast_rows"), big.mark = ","),
    " rows; ",
    format(metric_value("exact_strict_project_matched_rows"), big.mark = ","),
    " are exact strict project-matched rows and ",
    format(metric_value("relaxed_evaluable_pool_rows"), big.mark = ","),
    " are relaxed evaluable pooled rows."
  ),
  paste0(
    "- Partial pooling has ",
    format(metric_value("partial_pooling_pooled_supported_rows"),
           big.mark = ","),
    " pooled-supported rows, ",
    format(metric_value("partial_pooling_single_context_rows"),
           big.mark = ","),
    " single-context rows, and ",
    format(metric_value("partial_pooling_raw_shrunk_rows"), big.mark = ","),
    " raw-supported rows that shrink to fragile."
  ),
  paste0(
    "- Protocol QC spans ",
    format(metric_value("protocol_run_rows"), big.mark = ","),
    " unique runs, including ",
    format(metric_value("start_site_enriched_protocol_runs"), big.mark = ","),
    " LTM/harringtonine-like start-site-enriched runs."
  ),
  "",
  "## Claim Rules",
  "",
  "- Use `A_project_specific_matched_contrast` language only for exact local case/control evidence with adequate counts and merged replicates.",
  "- Treat replicate pooling as a measurement aid. It raises counts but can hide single-run dominance.",
  "- Treat design-family pooling as recurrence evidence. It can mix tissue, protocol, virus, lab, time point, and sample composition.",
  "- Use partial pooling to separate shrinkage-supported context effects from single-context or fragile rows.",
  "- Use all-merged coverage only as a transcript-shape background.",
  "",
  "## Top Pooling-Hazard Design Families",
  "",
  if (nrow(top_hazard)) {
    paste(capture.output(print(top_hazard)), collapse = "\n")
  } else {
    "No branch-context pooling hazard table was available."
  },
  "",
  "## Top Partial-Pooling Review Rows",
  "",
  if (nrow(top_partial)) {
    paste(capture.output(print(top_partial)), collapse = "\n")
  } else {
    "No partial-pooling review rows were available."
  },
  "",
  "## Output Files",
  "",
  "- `pooling_claim_language.csv`: claim-tier vocabulary.",
  "- `normalization_limits.csv`: what each normalized quantity does and does not control.",
  "- `pooling_gain_loss_summary.csv`: coverage gains from merged replicates.",
  "- `pooling_contrast_tier_summary.csv`: exact versus relaxed evidence tier counts.",
  "- `design_family_pooling_bias_summary.csv`: design-family confounding hazards.",
  "- `protocol_bias_run_summary.csv`: inhibitor/protocol and frame-QC summaries.",
  "- `partial_pooling_tradeoff_summary.csv`: pooled support, single-context support, and shrinkage loss.",
  "- `replicate_pooling_loss_summary.csv`: single-run dominance diagnostics from static inspection packets."
)
writeLines(report, file.path(output_dir, "pooling_bias_report.md"))

message("Saved pooling/protocol-bias audit outputs to: ", output_dir)
