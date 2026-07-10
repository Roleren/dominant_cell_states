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
output_dir <- file.path(analysis_dir, "dominant_rdg_joint_branch_allocation")
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

max_safe <- function(x, default = 0) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

mean_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

median_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

first_nonempty <- function(x, default = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else default
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

q_component <- function(q) {
  q <- safe_num(q)
  out <- rep(0, length(q))
  finite <- is.finite(q) & q > 0
  out[finite] <- pmin(-log10(q[finite]) / 3, 1)
  out[is.finite(q) & q == 0] <- 1
  out
}

cor_safe <- function(x, y) {
  x <- safe_num(x)
  y <- safe_num(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5 || uniqueN(x[keep]) < 2 || uniqueN(y[keep]) < 2) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

safe_js_divergence <- function(p, q) {
  p <- safe_num(p)
  q <- safe_num(q)
  if (!all(is.finite(p)) || !all(is.finite(q)) || sum(p) <= 0 || sum(q) <= 0) {
    return(NA_real_)
  }
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl <- function(a, b) {
    keep <- a > 0 & b > 0
    if (!any(keep)) return(0)
    sum(a[keep] * log2(a[keep] / b[keep]))
  }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

dirichlet_var <- function(a, total_a) {
  a * (total_a - a) / (total_a^2 * (total_a + 1))
}

message("Joint RDG branch-allocation model:")
message("  1. Aggregate grouped branch counts into clean CDS / leader uORF / overlapping uORF / other")
message("  2. Fit row-level Dirichlet posterior compositions for case and control")
message("  3. Add motif-prior sensitivity and atlas-ready review summaries")

grouped_file <- file.path(
  analysis_dir,
  "dominant_rdg_grouped_branch_usage",
  "dominant_rdg_grouped_branch_usage_contrast_table.csv"
)
motif_file <- file.path(
  analysis_dir,
  "dominant_rdg_mechanistic_motif_model",
  "dominant_rdg_mechanistic_motif_scores.csv"
)
atlas_file <- file.path(
  analysis_dir,
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_cards.csv"
)

dt <- read_dt(grouped_file, required = TRUE)
motif <- read_dt(motif_file)
atlas <- read_dt(atlas_file)

required <- c(
  "gene_symbol", "tx_id", "design_family", "study", "case_condition",
  "control_conditions", "match_scope", "sample_set_signature",
  "branch_feature_class", "feature_id", "case_sum_raw_counts",
  "control_sum_raw_counts", "case_branch_total_counts",
  "control_branch_total_counts", "n_case_runs_merged",
  "n_control_runs_merged"
)
missing <- setdiff(required, names(dt))
if (length(missing)) {
  stop("Grouped branch table is missing columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

allocation_classes <- c("clean_CDS", "leader_uORF", "overlapping_uORF")
joint_classes <- c(allocation_classes, "other")
class_suffix <- c(
  clean_CDS = "clean_CDS",
  leader_uORF = "leader_uORF",
  overlapping_uORF = "overlapping_uORF",
  other = "other"
)

numeric_cols <- c(
  "case_sum_raw_counts", "control_sum_raw_counts",
  "case_branch_total_counts", "control_branch_total_counts",
  "n_case_runs_merged", "n_control_runs_merged", "n_case_groups",
  "n_control_groups", "grouped_row_weight", "grouped_min_branch_total_counts",
  "grouped_max_feature_counts", "grouped_min_feature_counts"
)
for (column in intersect(numeric_cols, names(dt))) {
  dt[, (column) := finite_or_zero(get(column))]
}

dt <- dt[branch_feature_class %chin% allocation_classes]
if (!nrow(dt)) {
  warning("No allocation feature classes found in grouped branch table.")
  empty_metrics <- data.table(metric = "joint_allocation_rows", value = 0)
  fwrite(empty_metrics, file.path(
    output_dir,
    "dominant_rdg_joint_branch_allocation_summary_metrics.csv"
  ))
  quit(save = "no", status = 0)
}

meta_cols <- intersect(
  c("stratum_id", "case_condition", "perturbation_class", "design_family",
    "study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE", "INHIBITOR",
    "FRACTION", "Cancer_type", "Sex", "TIMEPOINT", "control_conditions",
    "match_scope", "sample_set_signature", "gene_symbol", "tx_id"),
  names(dt)
)
key_cols <- meta_cols

denom <- dt[, .(
  n_case_groups = max_safe(n_case_groups),
  n_control_groups = max_safe(n_control_groups),
  n_case_runs_merged = max_safe(n_case_runs_merged),
  n_control_runs_merged = max_safe(n_control_runs_merged),
  case_branch_total_counts = max_safe(case_branch_total_counts),
  control_branch_total_counts = max_safe(control_branch_total_counts),
  case_total_denominator_values = uniqueN(case_branch_total_counts),
  control_total_denominator_values = uniqueN(control_branch_total_counts),
  max_grouped_row_weight = max_safe(grouped_row_weight, default = 1),
  n_allocation_features = uniqueN(feature_id),
  allocation_feature_ids = collapse_unique(feature_id, 10L),
  allocation_feature_classes = collapse_unique(branch_feature_class, 4L)
), by = key_cols]

class_counts <- dt[, .(
  n_features_in_class = uniqueN(feature_id),
  feature_ids = collapse_unique(feature_id, 10L),
  case_class_count_sum = sum(pmax(finite_or_zero(case_sum_raw_counts), 0)),
  control_class_count_sum = sum(pmax(finite_or_zero(control_sum_raw_counts), 0)),
  case_class_count_max = max_safe(case_sum_raw_counts),
  control_class_count_max = max_safe(control_sum_raw_counts)
), by = c(key_cols, "branch_feature_class")]

joint <- copy(denom)
for (cls in allocation_classes) {
  suffix <- class_suffix[[cls]]
  sub <- class_counts[branch_feature_class == cls]
  keep <- setdiff(names(sub), "branch_feature_class")
  sub <- sub[, ..keep]
  rename_cols <- setdiff(names(sub), key_cols)
  setnames(sub, rename_cols, paste0(rename_cols, "_", suffix))
  joint <- merge(joint, sub, by = key_cols, all.x = TRUE, sort = FALSE)
}

for (cls in allocation_classes) {
  suffix <- class_suffix[[cls]]
  for (prefix in c("n_features_in_class", "case_class_count_sum",
                   "control_class_count_sum", "case_class_count_max",
                   "control_class_count_max")) {
    column <- paste0(prefix, "_", suffix)
    if (!column %in% names(joint)) joint[, (column) := 0]
    joint[is.na(get(column)), (column) := 0]
  }
  feature_col <- paste0("feature_ids_", suffix)
  if (!feature_col %in% names(joint)) joint[, (feature_col) := ""]
  joint[is.na(get(feature_col)), (feature_col) := ""]
}

case_raw_cols <- paste0("case_class_count_sum_", class_suffix[allocation_classes])
control_raw_cols <- paste0("control_class_count_sum_", class_suffix[allocation_classes])
case_raw <- as.matrix(joint[, ..case_raw_cols])
control_raw <- as.matrix(joint[, ..control_raw_cols])
case_raw[!is.finite(case_raw)] <- 0
control_raw[!is.finite(control_raw)] <- 0
case_raw <- pmax(case_raw, 0)
control_raw <- pmax(control_raw, 0)

case_total <- pmax(finite_or_zero(joint$case_branch_total_counts), 0)
control_total <- pmax(finite_or_zero(joint$control_branch_total_counts), 0)
case_sum_raw <- rowSums(case_raw)
control_sum_raw <- rowSums(control_raw)
case_scale <- ifelse(case_sum_raw > case_total & case_sum_raw > 0,
                     case_total / case_sum_raw, 1)
control_scale <- ifelse(control_sum_raw > control_total & control_sum_raw > 0,
                        control_total / control_sum_raw, 1)
case_scaled <- case_raw * case_scale
control_scaled <- control_raw * control_scale

for (i in seq_along(allocation_classes)) {
  cls <- allocation_classes[[i]]
  suffix <- class_suffix[[cls]]
  joint[, (paste0("joint_case_count_", suffix)) := case_scaled[, i]]
  joint[, (paste0("joint_control_count_", suffix)) := control_scaled[, i]]
}
joint[, joint_case_raw_class_count_sum := case_sum_raw]
joint[, joint_control_raw_class_count_sum := control_sum_raw]
joint[, joint_case_class_count_scale := case_scale]
joint[, joint_control_class_count_scale := control_scale]
joint[, joint_case_counts_scaled := case_scale < 0.999999]
joint[, joint_control_counts_scaled := control_scale < 0.999999]

case_scaled_sum <- rowSums(case_scaled)
control_scaled_sum <- rowSums(control_scaled)
joint[, joint_case_count_other := pmax(case_total - case_scaled_sum, 0)]
joint[, joint_control_count_other := pmax(control_total - control_scaled_sum, 0)]

motif_required <- c(
  "gene_symbol", "tx_id", "inhibitory_overlap_burden",
  "reinitiation_distance_opportunity", "dense_uorf_collision_burden",
  "strong_start_leakage_burden", "downstream_rescue_opportunity",
  "overlap_without_rescue", "leakage_with_reinitiation",
  "dense_short_reinitiation", "rescue_after_overlap",
  "motif_complexity_index", "n_uorfs", "n_overlapping_uorfs",
  "uorf_count_bin", "overlap_class", "spacing_class",
  "last_leader_gap_class"
)
if (nrow(motif) && all(c("gene_symbol", "tx_id") %in% names(motif))) {
  motif_keep <- intersect(motif_required, names(motif))
  motif_small <- unique(motif[, ..motif_keep])
  joint <- merge(joint, motif_small, by = c("gene_symbol", "tx_id"),
                 all.x = TRUE, sort = FALSE)
}
for (column in setdiff(motif_required, c("gene_symbol", "tx_id",
                                         "uorf_count_bin", "overlap_class",
                                         "spacing_class",
                                         "last_leader_gap_class"))) {
  if (!column %in% names(joint)) joint[, (column) := 0]
  joint[is.na(get(column)), (column) := 0]
}
for (column in c("uorf_count_bin", "overlap_class", "spacing_class",
                 "last_leader_gap_class")) {
  if (!column %in% names(joint)) joint[, (column) := ""]
  joint[is.na(get(column)), (column) := ""]
}

joint[, motif_signed_clean_CDS :=
        clip01(0.45 * reinitiation_distance_opportunity +
                 0.35 * downstream_rescue_opportunity +
                 0.20 * leakage_with_reinitiation) -
        clip01(0.45 * overlap_without_rescue +
                 0.35 * dense_short_reinitiation +
                 0.20 * strong_start_leakage_burden)]
joint[, motif_signed_leader_uORF :=
        clip01(0.35 * dense_uorf_collision_burden +
                 0.30 * strong_start_leakage_burden +
                 0.20 * leakage_with_reinitiation +
                 0.15 * reinitiation_distance_opportunity)]
joint[, motif_signed_overlapping_uORF :=
        clip01(0.45 * inhibitory_overlap_burden +
                 0.35 * overlap_without_rescue +
                 0.20 * rescue_after_overlap)]
joint[, motif_signed_other := 0]

alpha_uniform <- rep(0.5, length(joint_classes))
names(alpha_uniform) <- joint_classes
joint[, motif_alpha_clean_CDS :=
        0.5 + 1.5 * clip01(motif_signed_clean_CDS +
                             pmax(downstream_rescue_opportunity, 0))]
joint[, motif_alpha_leader_uORF :=
        0.5 + 1.5 * clip01(motif_signed_leader_uORF)]
joint[, motif_alpha_overlapping_uORF :=
        0.5 + 1.5 * clip01(motif_signed_overlapping_uORF)]
joint[, motif_alpha_other := 0.5]

case_count_cols <- paste0("joint_case_count_", class_suffix[joint_classes])
control_count_cols <- paste0("joint_control_count_", class_suffix[joint_classes])
case_counts <- as.matrix(joint[, ..case_count_cols])
control_counts <- as.matrix(joint[, ..control_count_cols])
case_counts[!is.finite(case_counts)] <- 0
control_counts[!is.finite(control_counts)] <- 0

case_alpha_uniform <- sweep(case_counts, 2, alpha_uniform, "+")
control_alpha_uniform <- sweep(control_counts, 2, alpha_uniform, "+")
case_A <- rowSums(case_alpha_uniform)
control_A <- rowSums(control_alpha_uniform)
case_p <- case_alpha_uniform / case_A
control_p <- control_alpha_uniform / control_A
case_var <- case_alpha_uniform * (case_A - case_alpha_uniform) /
  (case_A^2 * (case_A + 1))
control_var <- control_alpha_uniform * (control_A - control_alpha_uniform) /
  (control_A^2 * (control_A + 1))
delta <- case_p - control_p
delta_se <- sqrt(case_var + control_var)
delta_z <- delta / delta_se
delta_p <- 2 * pnorm(-abs(delta_z))
delta_p[!is.finite(delta_p)] <- NA_real_

for (i in seq_along(joint_classes)) {
  cls <- joint_classes[[i]]
  suffix <- class_suffix[[cls]]
  joint[, (paste0("joint_case_prop_", suffix)) := case_p[, i]]
  joint[, (paste0("joint_control_prop_", suffix)) := control_p[, i]]
  joint[, (paste0("joint_delta_", suffix)) := delta[, i]]
  joint[, (paste0("joint_delta_se_", suffix)) := delta_se[, i]]
  joint[, (paste0("joint_delta_p_", suffix)) := delta_p[, i]]
  joint[, (paste0("joint_delta_q_", suffix)) := p.adjust(delta_p[, i],
                                                        method = "BH")]
}

motif_alpha_cols <- paste0("motif_alpha_", class_suffix[joint_classes])
motif_alpha <- as.matrix(joint[, ..motif_alpha_cols])
motif_alpha[!is.finite(motif_alpha)] <- 0.5
case_alpha_motif <- case_counts + motif_alpha
control_alpha_motif <- control_counts + motif_alpha
case_p_motif <- case_alpha_motif / rowSums(case_alpha_motif)
control_p_motif <- control_alpha_motif / rowSums(control_alpha_motif)
delta_motif <- case_p_motif - control_p_motif
for (i in seq_along(joint_classes)) {
  cls <- joint_classes[[i]]
  suffix <- class_suffix[[cls]]
  joint[, (paste0("motif_prior_joint_delta_", suffix)) := delta_motif[, i]]
}

delta_abs <- abs(delta)
dominant_idx <- max.col(delta_abs, ties.method = "first")
row_index <- seq_len(nrow(joint))
dominant_class <- joint_classes[dominant_idx]
joint[, joint_top_branch_class := dominant_class]
joint[, joint_top_branch_delta := delta[cbind(row_index, dominant_idx)]]
joint[, joint_top_branch_p := delta_p[cbind(row_index, dominant_idx)]]
joint[, joint_top_branch_se := delta_se[cbind(row_index, dominant_idx)]]
joint[, joint_top_branch_q := {
  q_cols <- paste0("joint_delta_q_", class_suffix[joint_classes])
  q_mat <- as.matrix(.SD)
  q_mat[cbind(row_index, dominant_idx)]
}, .SDcols = paste0("joint_delta_q_", class_suffix[joint_classes])]

joint[, joint_l1_allocation_shift := rowSums(delta_abs) / 2]
joint[, joint_max_abs_branch_delta := apply(delta_abs, 1, max)]
joint[, joint_js_divergence := vapply(
  seq_len(.N),
  function(i) safe_js_divergence(case_p[i, ], control_p[i, ]),
  numeric(1)
)]
joint[, joint_chisq_stat := {
  var_delta <- delta_se^2
  rowSums(ifelse(is.finite(delta) & is.finite(var_delta) & var_delta > 0,
                 delta^2 / var_delta, 0))
}]
joint[, joint_chisq_p := pchisq(joint_chisq_stat, df = length(joint_classes) - 1,
                                lower.tail = FALSE)]
joint[!is.finite(joint_chisq_p), joint_chisq_p := NA_real_]
joint[, joint_chisq_q := p.adjust(joint_chisq_p, method = "BH")]

joint[, joint_min_runs_merged :=
        pmin(n_case_runs_merged, n_control_runs_merged)]
joint[, joint_min_branch_total_counts :=
        pmin(case_branch_total_counts, control_branch_total_counts)]
max_class_count <- pmax(
  joint$joint_case_count_clean_CDS, joint$joint_control_count_clean_CDS,
  joint$joint_case_count_leader_uORF, joint$joint_control_count_leader_uORF,
  joint$joint_case_count_overlapping_uORF,
  joint$joint_control_count_overlapping_uORF
)
max_both_class_count <- pmax(
  pmin(joint$joint_case_count_clean_CDS, joint$joint_control_count_clean_CDS),
  pmin(joint$joint_case_count_leader_uORF, joint$joint_control_count_leader_uORF),
  pmin(joint$joint_case_count_overlapping_uORF,
       joint$joint_control_count_overlapping_uORF)
)
joint[, joint_max_class_count := max_class_count]
joint[, joint_max_both_class_count := max_both_class_count]
joint[, joint_passes_replicate_gate := joint_min_runs_merged >= 2]
joint[, joint_passes_total_gate := joint_min_branch_total_counts >= 100]
joint[, joint_passes_class_gate_any := joint_max_class_count >= 30]
joint[, joint_passes_class_gate_both := joint_max_both_class_count >= 30]
joint[, joint_evaluable :=
        joint_passes_replicate_gate &
        joint_passes_total_gate &
        joint_passes_class_gate_any &
        is.finite(joint_l1_allocation_shift) &
        is.finite(joint_chisq_p)]
joint[, joint_strict :=
        joint_evaluable &
        joint_passes_class_gate_both &
        match_scope == "exact_context"]

dominant_prior_cols <- paste0("motif_signed_", class_suffix[joint_classes])
prior_mat <- as.matrix(joint[, ..dominant_prior_cols])
prior_mat[!is.finite(prior_mat)] <- 0
dominant_prior <- prior_mat[cbind(row_index, dominant_idx)]
joint[, joint_top_motif_prior_signed := dominant_prior]
joint[, joint_top_motif_prior_strength := abs(dominant_prior)]
joint[, joint_top_motif_alignment := fcase(
  joint_top_branch_class == "other" |
    joint_top_motif_prior_strength < 0.20,
  "weak_or_no_architecture_prior",
  abs(joint_top_branch_delta) < 0.02,
  "weak_observed_shift",
  sign(joint_top_branch_delta) == sign(joint_top_motif_prior_signed),
  "motif_joint_concordant",
  sign(joint_top_branch_delta) != sign(joint_top_motif_prior_signed),
  "motif_joint_discordant",
  default = "motif_joint_ambiguous"
)]
joint[, joint_motif_alignment_score := clip01(
  joint_top_motif_prior_strength * pmin(abs(joint_top_branch_delta) / 0.15, 1)
)]

joint[, joint_row_weight := pmax(
  1,
  log1p(joint_max_class_count) *
    sqrt(pmax(joint_min_branch_total_counts, 1)) *
    fifelse(match_scope == "exact_context", 1, 0.55)
)]
joint[!is.finite(joint_row_weight), joint_row_weight := 1]
joint[, joint_review_priority := clip01(
  0.22 * pmin(joint_l1_allocation_shift / 0.20, 1) +
    0.16 * pmin(joint_max_abs_branch_delta / 0.12, 1) +
    0.15 * q_component(joint_chisq_q) +
    0.12 * pmin(joint_min_runs_merged / 3, 1) +
    0.12 * pmin(log1p(joint_max_class_count) / log1p(300), 1) +
    0.10 * as.numeric(joint_strict == TRUE) +
    0.08 * joint_motif_alignment_score +
    0.05 * as.numeric(design_family == "viral_infection")
)]
joint[, joint_review_class := fcase(
  !joint_evaluable,
  "low_count_or_no_grouped_control",
  joint_strict == TRUE &
    is.finite(joint_chisq_q) &
    joint_chisq_q <= 0.10 &
    joint_l1_allocation_shift >= 0.08 &
    joint_max_abs_branch_delta >= 0.04,
  "strict_joint_branch_shift",
  joint_evaluable == TRUE &
    is.finite(joint_chisq_p) &
    joint_chisq_p <= 0.01 &
    joint_l1_allocation_shift >= 0.10 &
    joint_max_abs_branch_delta >= 0.05,
  "joint_branch_shift_review",
  joint_evaluable == TRUE &
    joint_top_motif_alignment == "motif_joint_discordant" &
    joint_l1_allocation_shift >= 0.08,
  "motif_joint_discordant_review",
  joint_evaluable == TRUE &
    joint_top_motif_prior_strength >= 0.45 &
    joint_l1_allocation_shift < 0.04,
  "high_motif_prior_weak_joint_evidence",
  default = "joint_evaluable_no_strong_shift"
)]
joint[joint_evaluable != TRUE,
      joint_review_priority := pmin(joint_review_priority, 0.20)]

row_file <- file.path(
  output_dir,
  "dominant_rdg_joint_branch_allocation_rows.csv"
)
fwrite(joint, row_file)

class_long <- melt(
  joint,
  id.vars = key_cols,
  measure.vars = list(
    joint_case_prop = paste0("joint_case_prop_", class_suffix[joint_classes]),
    joint_control_prop = paste0("joint_control_prop_", class_suffix[joint_classes]),
    joint_delta = paste0("joint_delta_", class_suffix[joint_classes]),
    joint_delta_se = paste0("joint_delta_se_", class_suffix[joint_classes]),
    joint_delta_p = paste0("joint_delta_p_", class_suffix[joint_classes]),
    joint_delta_q = paste0("joint_delta_q_", class_suffix[joint_classes]),
    joint_case_count = paste0("joint_case_count_", class_suffix[joint_classes]),
    joint_control_count = paste0("joint_control_count_", class_suffix[joint_classes])
  ),
  variable.name = "branch_class_index"
)
class_long[, branch_class := joint_classes[branch_class_index]]
class_long[, branch_class_index := NULL]
class_long <- merge(
  class_long,
  joint[, c(key_cols, "joint_evaluable", "joint_strict",
            "joint_l1_allocation_shift", "joint_js_divergence",
            "joint_chisq_p", "joint_chisq_q", "joint_review_class",
            "joint_review_priority"), with = FALSE],
  by = key_cols,
  all.x = TRUE,
  sort = FALSE
)
class_file <- file.path(
  output_dir,
  "dominant_rdg_joint_branch_allocation_class_rows.csv"
)
fwrite(class_long, class_file)

gene_design <- joint[, {
  top_idx <- which.max(joint_review_priority)
  if (!length(top_idx) || is.na(top_idx)) top_idx <- 1L
  evaluable <- joint_evaluable == TRUE
  strict <- joint_strict == TRUE
  data.table(
    joint_rows = .N,
    joint_evaluable_rows = sum(evaluable, na.rm = TRUE),
    joint_strict_rows = sum(strict, na.rm = TRUE),
    joint_strict_studies = uniqueN(study[strict]),
    joint_evaluable_studies = uniqueN(study[evaluable]),
    joint_strict_cell_lines = uniqueN(CELL_LINE[strict]),
    joint_weighted_l1_shift = weighted_mean_safe(
      joint_l1_allocation_shift[evaluable],
      joint_row_weight[evaluable]
    ),
    joint_median_l1_shift = median_safe(joint_l1_allocation_shift[evaluable]),
    joint_max_l1_shift = max_safe(joint_l1_allocation_shift[evaluable]),
    joint_weighted_clean_CDS_delta = weighted_mean_safe(
      joint_delta_clean_CDS[evaluable],
      joint_row_weight[evaluable]
    ),
    joint_weighted_leader_uORF_delta = weighted_mean_safe(
      joint_delta_leader_uORF[evaluable],
      joint_row_weight[evaluable]
    ),
    joint_weighted_overlapping_uORF_delta = weighted_mean_safe(
      joint_delta_overlapping_uORF[evaluable],
      joint_row_weight[evaluable]
    ),
    joint_weighted_other_delta = weighted_mean_safe(
      joint_delta_other[evaluable],
      joint_row_weight[evaluable]
    ),
    joint_max_abs_clean_CDS_delta = max_safe(abs(joint_delta_clean_CDS[evaluable])),
    joint_max_abs_leader_uORF_delta =
      max_safe(abs(joint_delta_leader_uORF[evaluable])),
    joint_max_abs_overlapping_uORF_delta =
      max_safe(abs(joint_delta_overlapping_uORF[evaluable])),
    joint_max_abs_other_delta = max_safe(abs(joint_delta_other[evaluable])),
    joint_min_q = max_safe(-log10(pmax(joint_chisq_q[evaluable], 1e-300))),
    joint_scaled_row_fraction = mean(
      joint_case_counts_scaled == TRUE | joint_control_counts_scaled == TRUE,
      na.rm = TRUE
    ),
    joint_top_branch_class = joint_top_branch_class[top_idx],
    joint_top_branch_delta = joint_top_branch_delta[top_idx],
    joint_top_branch_q = joint_top_branch_q[top_idx],
    joint_top_l1_shift = joint_l1_allocation_shift[top_idx],
    joint_top_js_divergence = joint_js_divergence[top_idx],
    joint_top_review_class = joint_review_class[top_idx],
    joint_top_review_priority = joint_review_priority[top_idx],
    joint_top_motif_alignment = joint_top_motif_alignment[top_idx],
    joint_top_motif_prior_signed = joint_top_motif_prior_signed[top_idx],
    joint_top_motif_prior_strength = joint_top_motif_prior_strength[top_idx],
    joint_top_context = paste(study[top_idx], case_condition[top_idx],
                              "vs", control_conditions[top_idx], sep = " | "),
    joint_top_studies = collapse_unique(study, 8L),
    joint_top_conditions = collapse_unique(case_condition, 8L)
  )
}, by = .(gene_symbol, tx_id, design_family)]

if (nrow(atlas) && all(c("gene_symbol", "tx_id", "design_family") %in% names(atlas))) {
  atlas_keep <- intersect(names(atlas), c(
    "gene_symbol", "tx_id", "design_family",
    "atlas_evidence_class", "atlas_browser_review_priority_score",
    "atlas_clean_cds_pct_change", "browser_validation_status",
    "browser_validation_label", "count_aware_top_support_class",
    "count_aware_top_score", "hierarchical_evidence_class",
    "hierarchical_priority_score"
  ))
  gene_design <- merge(gene_design, unique(atlas[, ..atlas_keep]),
                       by = c("gene_symbol", "tx_id", "design_family"),
                       all.x = TRUE, sort = FALSE)
}

setorder(gene_design, -joint_top_review_priority, -joint_max_l1_shift,
         gene_symbol, design_family)
gene_file <- file.path(
  output_dir,
  "dominant_rdg_joint_branch_allocation_gene_design_summary.csv"
)
fwrite(gene_design, gene_file)

review <- gene_design[
  joint_top_review_class != "low_count_or_no_grouped_control" &
    (joint_top_review_class != "joint_evaluable_no_strong_shift" |
       joint_top_review_priority >= 0.55)
]
setorder(review, -joint_top_review_priority, -joint_max_l1_shift,
         gene_symbol, design_family)
review_file <- file.path(
  output_dir,
  "dominant_rdg_joint_branch_allocation_review.csv"
)
fwrite(review, review_file)

viral_review <- gene_design[design_family == "viral_infection"]
setorder(viral_review, -joint_top_review_priority, -joint_max_l1_shift,
         gene_symbol)
viral_file <- file.path(
  output_dir,
  "dominant_rdg_joint_branch_allocation_viral_review.csv"
)
fwrite(viral_review, viral_file)

summary <- rbindlist(list(
  joint[, .(
    joint_allocation_rows = .N,
    joint_evaluable_rows = sum(joint_evaluable, na.rm = TRUE),
    joint_strict_rows = sum(joint_strict, na.rm = TRUE),
    joint_review_rows = sum(joint_review_class !=
                              "low_count_or_no_grouped_control" &
                              (joint_review_class !=
                                 "joint_evaluable_no_strong_shift" |
                                 joint_review_priority >= 0.55),
                            na.rm = TRUE),
    joint_strict_shift_rows =
      sum(joint_review_class == "strict_joint_branch_shift", na.rm = TRUE),
    joint_shift_review_rows =
      sum(joint_review_class == "joint_branch_shift_review", na.rm = TRUE),
    motif_discordant_review_rows =
      sum(joint_review_class == "motif_joint_discordant_review", na.rm = TRUE),
    high_prior_weak_joint_rows =
      sum(joint_review_class == "high_motif_prior_weak_joint_evidence",
          na.rm = TRUE),
    scaled_rows = sum(joint_case_counts_scaled | joint_control_counts_scaled,
                      na.rm = TRUE),
    median_l1_shift = median_safe(joint_l1_allocation_shift[joint_evaluable]),
    max_l1_shift = max_safe(joint_l1_allocation_shift[joint_evaluable]),
    spearman_motif_strength_vs_l1 =
      cor_safe(joint_top_motif_prior_strength[joint_evaluable],
               joint_l1_allocation_shift[joint_evaluable]),
    genes = uniqueN(gene_symbol),
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    design_families = uniqueN(design_family)
  )][, scope := "all"],
  joint[, .(
    joint_allocation_rows = .N,
    joint_evaluable_rows = sum(joint_evaluable, na.rm = TRUE),
    joint_strict_rows = sum(joint_strict, na.rm = TRUE),
    joint_review_rows = sum(joint_review_class !=
                              "low_count_or_no_grouped_control" &
                              (joint_review_class !=
                                 "joint_evaluable_no_strong_shift" |
                                 joint_review_priority >= 0.55),
                            na.rm = TRUE),
    joint_strict_shift_rows =
      sum(joint_review_class == "strict_joint_branch_shift", na.rm = TRUE),
    joint_shift_review_rows =
      sum(joint_review_class == "joint_branch_shift_review", na.rm = TRUE),
    motif_discordant_review_rows =
      sum(joint_review_class == "motif_joint_discordant_review", na.rm = TRUE),
    high_prior_weak_joint_rows =
      sum(joint_review_class == "high_motif_prior_weak_joint_evidence",
          na.rm = TRUE),
    scaled_rows = sum(joint_case_counts_scaled | joint_control_counts_scaled,
                      na.rm = TRUE),
    median_l1_shift = median_safe(joint_l1_allocation_shift[joint_evaluable]),
    max_l1_shift = max_safe(joint_l1_allocation_shift[joint_evaluable]),
    spearman_motif_strength_vs_l1 =
      cor_safe(joint_top_motif_prior_strength[joint_evaluable],
               joint_l1_allocation_shift[joint_evaluable]),
    genes = uniqueN(gene_symbol),
    gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    design_families = uniqueN(design_family)
  ), by = .(scope = paste0("design_family:", design_family))]
), fill = TRUE)
setcolorder(summary, c("scope", setdiff(names(summary), "scope")))
summary_file <- file.path(
  output_dir,
  "dominant_rdg_joint_branch_allocation_summary_metrics.csv"
)
fwrite(summary, summary_file)

plot_dt <- review[joint_top_review_priority >= 0.35]
if (nrow(plot_dt)) {
  plot_dt <- plot_dt[order(-joint_top_review_priority)][seq_len(min(.N, 45))]
  plot_dt[, plot_label := paste0(gene_symbol, " / ", design_family)]
  plot_dt[, plot_label := factor(plot_label, levels = rev(plot_label))]
  p_top <- ggplot(
    plot_dt,
    aes(
      x = joint_top_l1_shift,
      y = plot_label,
      fill = joint_top_branch_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Design: ", design_family,
        "<br>Top class: ", joint_top_branch_class,
        "<br>Top delta: ", round(joint_top_branch_delta, 3),
        "<br>L1 shift: ", round(joint_top_l1_shift, 3),
        "<br>q: ", signif(joint_top_branch_q, 3),
        "<br>Review: ", joint_top_review_class,
        "<br>Motif alignment: ", joint_top_motif_alignment,
        "<br>Context: ", joint_top_context
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
      title = "Joint RDG Branch-Allocation Shifts",
      subtitle = "Dirichlet posterior composition over clean CDS, leader uORF, overlapping uORF, and residual branch mass",
      x = "Total allocation shift (L1 / 2)",
      y = NULL,
      fill = "Dominant branch"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "joint_branch_allocation_top_shifts.png"),
         p_top, width = 10.8, height = 8.4, dpi = 180)
  ggsave(file.path(figure_dir, "joint_branch_allocation_top_shifts.pdf"),
         p_top, width = 10.8, height = 8.4)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_top,
      file.path(figure_dir, "joint_branch_allocation_top_shifts.html"),
      title = "Joint branch allocation shifts"
    )
  }
}

viral_plot <- viral_review[joint_top_review_priority >= 0.25]
if (nrow(viral_plot)) {
  viral_plot <- viral_plot[order(-joint_top_review_priority)][seq_len(min(.N, 35))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p_viral <- ggplot(
    viral_plot,
    aes(
      x = joint_top_review_priority,
      y = gene_symbol,
      fill = joint_top_branch_class,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Review: ", joint_top_review_class,
        "<br>Priority: ", round(joint_top_review_priority, 3),
        "<br>Top class: ", joint_top_branch_class,
        "<br>Top delta: ", round(joint_top_branch_delta, 3),
        "<br>L1 shift: ", round(joint_top_l1_shift, 3),
        "<br>Motif alignment: ", joint_top_motif_alignment,
        "<br>Context: ", joint_top_context
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
      title = "Viral-Infection Joint RDG Allocation Review",
      subtitle = "Higher priority means stronger replicated allocation shift, count support, and motif/context support",
      x = "Joint allocation review priority",
      y = NULL,
      fill = "Dominant branch"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "joint_branch_allocation_viral_review.png"),
         p_viral, width = 9.8, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "joint_branch_allocation_viral_review.pdf"),
         p_viral, width = 9.8, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_viral,
      file.path(figure_dir, "joint_branch_allocation_viral_review.html"),
      title = "Viral joint branch allocation review"
    )
  }
}

heat_dt <- copy(gene_design[order(-joint_top_review_priority)][seq_len(min(.N, 40))])
if (nrow(heat_dt)) {
  heat_dt[, plot_label := paste0(gene_symbol, " / ", design_family)]
  heat_long <- melt(
    heat_dt,
    id.vars = c("plot_label", "gene_symbol", "design_family",
                "joint_top_review_priority", "joint_top_review_class"),
    measure.vars = c("joint_weighted_clean_CDS_delta",
                     "joint_weighted_leader_uORF_delta",
                     "joint_weighted_overlapping_uORF_delta"),
    variable.name = "branch",
    value.name = "weighted_delta"
  )
  heat_long[, branch := factor(
    branch,
    levels = c("joint_weighted_leader_uORF_delta",
               "joint_weighted_overlapping_uORF_delta",
               "joint_weighted_clean_CDS_delta"),
    labels = c("leader uORF", "overlapping uORF", "clean CDS")
  )]
  heat_long[, plot_label := factor(plot_label, levels = rev(unique(heat_dt$plot_label)))]
  p_heat <- ggplot(
    heat_long,
    aes(
      x = branch,
      y = plot_label,
      fill = weighted_delta,
      text = paste0(
        "Gene/design: ", plot_label,
        "<br>Branch: ", branch,
        "<br>Weighted delta: ", round(weighted_delta, 3),
        "<br>Review: ", joint_top_review_class,
        "<br>Priority: ", round(joint_top_review_priority, 3)
      )
    )
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#1f5aa6", mid = "white", high = "#8b1a1a",
      midpoint = 0, limits = c(-0.35, 0.35), oob = scales::squish
    ) +
    labs(
      title = "Joint RDG Allocation Delta Matrix",
      subtitle = "Positive values mean higher case allocation than matched control",
      x = NULL,
      y = NULL,
      fill = "Delta"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(figure_dir, "joint_branch_allocation_delta_heatmap.png"),
         p_heat, width = 9.6, height = 8.8, dpi = 180)
  ggsave(file.path(figure_dir, "joint_branch_allocation_delta_heatmap.pdf"),
         p_heat, width = 9.6, height = 8.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_heat,
      file.path(figure_dir, "joint_branch_allocation_delta_heatmap.html"),
      title = "Joint branch allocation delta matrix"
    )
  }
}

message("Saved joint branch-allocation outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
