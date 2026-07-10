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
output_dir <- file.path(analysis_dir,
                        "dominant_rdg_allocation_output_disagreements")
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

first_finite <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (!length(values)) NA_real_ else values[[1]]
}

safe_max <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

safe_mean <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else mean(x)
}

collapse_unique <- function(x, n = 5L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

collapse_contexts <- function(dt, n = 4L) {
  if (!nrow(dt)) return("")
  label_col <- if ("branch_context_label" %in% names(dt)) {
    "branch_context_label"
  } else if ("study" %in% names(dt)) {
    "study"
  } else {
    names(dt)[[1]]
  }
  labels <- unique(as.character(dt[[label_col]]))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  paste(head(labels, n), collapse = " | ")
}

wrap_text <- function(x, width = 92L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

atlas <- read_dt(file.path(analysis_dir, "dominant_rdg_atlas",
                           "dominant_rdg_atlas_cards.csv"),
                 required = TRUE)
agreement <- read_dt(file.path(analysis_dir, "dominant_rdg_model_agreement",
                               "dominant_rdg_model_agreement_gene_design.csv"),
                     required = FALSE)
branch_contexts <- read_dt(file.path(analysis_dir, "dominant_rdg_atlas",
                                     "dominant_rdg_atlas_branch_contexts.csv"),
                           required = FALSE)
study_evidence <- read_dt(file.path(
  analysis_dir,
  "dominant_rdg_model_agreement_loo",
  "dominant_rdg_model_agreement_study_evidence.csv"
), required = FALSE)

if (!nrow(atlas)) stop("Atlas card table is empty.", call. = FALSE)
atlas[, gene_symbol := as.character(gene_symbol)]
atlas[, design_family := as.character(design_family)]

key_cols <- c("gene_symbol", "tx_id", "design_family")
agreement_keep <- intersect(
  c(
    key_cols,
    "clean_cds_pct_change",
    "clean_cds_ci_excludes_zero",
    "clean_cds_direction_stable",
    "joint_weighted_clean_CDS_delta",
    "joint_weighted_leader_uORF_delta",
    "joint_weighted_overlapping_uORF_delta",
    "hierarchical_joint_delta_clean_CDS",
    "hierarchical_joint_delta_leader_uORF",
    "hierarchical_joint_delta_overlapping_uORF",
    "dm_delta_clean_CDS",
    "dm_delta_leader_uORF",
    "dm_delta_overlapping_uORF",
    "allocation_output_disagreement",
    "high_branch_heterogeneity",
    "model_disagreement_flags",
    "clean_cds_output_allocation_alignment",
    "model_agreement_score",
    "model_agreement_tier",
    "model_agreement_top_context"
  ),
  names(agreement)
)
if (nrow(agreement) && length(agreement_keep) > length(key_cols)) {
  setnames(
    agreement[, ..agreement_keep],
    setdiff(agreement_keep, key_cols),
    paste0("agreement_", setdiff(agreement_keep, key_cols)),
    skip_absent = TRUE
  )
  agreement_prefixed <- agreement[, ..agreement_keep]
  setnames(
    agreement_prefixed,
    setdiff(agreement_keep, key_cols),
    paste0("agreement_", setdiff(agreement_keep, key_cols)),
    skip_absent = TRUE
  )
  atlas <- merge(
    atlas,
    agreement_prefixed,
    by = key_cols,
    all.x = TRUE,
    sort = FALSE
  )
}

atlas[, allocation_output_disagreement_flag := bool_col(
  atlas, "allocation_output_disagreement"
) | bool_col(atlas, "agreement_allocation_output_disagreement") |
  grepl("allocation_vs_output", chr_col(atlas, "model_disagreement_flags"),
        ignore.case = TRUE) |
  grepl("allocation_vs_output",
        chr_col(atlas, "agreement_model_disagreement_flags"),
        ignore.case = TRUE) |
  chr_col(atlas, "clean_cds_output_allocation_alignment") ==
  "opposite_direction" |
  chr_col(atlas, "agreement_clean_cds_output_allocation_alignment") ==
  "opposite_direction"]

audit <- atlas[allocation_output_disagreement_flag == TRUE]
if (!nrow(audit)) {
  stop("No allocation-versus-output disagreement rows were found.",
       call. = FALSE)
}

audit[, clean_cds_output_pct := mapply(
  first_finite,
  num_col(audit, "atlas_clean_cds_pct_change"),
  num_col(audit, "agreement_clean_cds_pct_change"),
  SIMPLIFY = TRUE
)]
audit[, clean_cds_output_pct_lower := mapply(
  first_finite,
  num_col(audit, "atlas_clean_cds_pct_lower"),
  num_col(audit, "strict_random_effects_pct_change_lower"),
  SIMPLIFY = TRUE
)]
audit[, clean_cds_output_pct_upper := mapply(
  first_finite,
  num_col(audit, "atlas_clean_cds_pct_upper"),
  num_col(audit, "strict_random_effects_pct_change_upper"),
  SIMPLIFY = TRUE
)]
audit[, clean_cds_output_log2fc := fifelse(
  is.finite(clean_cds_output_pct) & clean_cds_output_pct > -99,
  log2(1 + clean_cds_output_pct / 100),
  NA_real_
)]

audit[, clean_cds_allocation_delta := mapply(
  first_finite,
  num_col(audit, "hierarchical_joint_delta_clean_CDS"),
  num_col(audit, "dm_delta_clean_CDS"),
  num_col(audit, "joint_weighted_clean_CDS_delta"),
  num_col(audit, "count_aware_weighted_usage_delta_clean_CDS"),
  num_col(audit, "agreement_hierarchical_joint_delta_clean_CDS"),
  num_col(audit, "agreement_dm_delta_clean_CDS"),
  num_col(audit, "agreement_joint_weighted_clean_CDS_delta"),
  SIMPLIFY = TRUE
)]
audit[, leader_uorf_allocation_delta := mapply(
  first_finite,
  num_col(audit, "hierarchical_joint_delta_leader_uORF"),
  num_col(audit, "dm_delta_leader_uORF"),
  num_col(audit, "joint_weighted_leader_uORF_delta"),
  num_col(audit, "count_aware_weighted_usage_delta_leader_uORF"),
  num_col(audit, "agreement_hierarchical_joint_delta_leader_uORF"),
  num_col(audit, "agreement_dm_delta_leader_uORF"),
  num_col(audit, "agreement_joint_weighted_leader_uORF_delta"),
  SIMPLIFY = TRUE
)]
audit[, overlapping_uorf_allocation_delta := mapply(
  first_finite,
  num_col(audit, "hierarchical_joint_delta_overlapping_uORF"),
  num_col(audit, "dm_delta_overlapping_uORF"),
  num_col(audit, "joint_weighted_overlapping_uORF_delta"),
  num_col(audit, "count_aware_weighted_usage_delta_overlapping_uORF"),
  num_col(audit, "agreement_hierarchical_joint_delta_overlapping_uORF"),
  num_col(audit, "agreement_dm_delta_overlapping_uORF"),
  num_col(audit, "agreement_joint_weighted_overlapping_uORF_delta"),
  SIMPLIFY = TRUE
)]
audit[, other_allocation_delta := mapply(
  first_finite,
  num_col(audit, "hierarchical_joint_delta_other"),
  num_col(audit, "joint_weighted_other_delta"),
  SIMPLIFY = TRUE
)]

audit[, upstream_allocation_delta := rowSums(
  cbind(
    fifelse(is.finite(leader_uorf_allocation_delta),
            leader_uorf_allocation_delta, 0),
    fifelse(is.finite(overlapping_uorf_allocation_delta),
            overlapping_uorf_allocation_delta, 0)
  )
)]
audit[, output_allocation_direction_pair := paste0(
  fifelse(clean_cds_output_pct > 0, "output_up",
          fifelse(clean_cds_output_pct < 0, "output_down", "output_neutral")),
  "__",
  fifelse(clean_cds_allocation_delta > 0, "allocation_up",
          fifelse(clean_cds_allocation_delta < 0, "allocation_down",
                  "allocation_neutral"))
)]
audit[, total_load_compensation_proxy := fifelse(
  is.finite(clean_cds_output_log2fc) &
    is.finite(clean_cds_allocation_delta),
  clean_cds_output_log2fc - clean_cds_allocation_delta,
  NA_real_
)]
audit[, upstream_shift_proxy := fifelse(
  is.finite(upstream_allocation_delta),
  upstream_allocation_delta,
  NA_real_
)]
audit[, disagreement_magnitude_score := clip01(
  0.45 * clip01(abs(clean_cds_output_log2fc) / 1) +
    0.35 * clip01(abs(clean_cds_allocation_delta) / 0.15) +
    0.20 * clip01(abs(upstream_shift_proxy) / 0.15)
)]

audit[, heterogeneity_risk_score := clip01(
  0.25 * as.numeric(bool_col(audit, "high_branch_heterogeneity") |
                      bool_col(audit, "agreement_high_branch_heterogeneity")) +
    0.20 * clip01(num_col(audit, "atlas_i2", 0)) +
    0.15 * clip01(num_col(audit, "random_effects_i2", 0)) +
    0.15 * clip01(num_col(audit, "hierarchical_joint_top_i2", 0)) +
    0.10 * as.numeric(!bool_col(audit, "atlas_direction_stable")) +
    0.10 * as.numeric(!bool_col(audit, "atlas_interval_supported")) +
    0.05 * as.numeric(num_col(audit, "strict_supported_studies", 0) < 2)
)]
audit[, compensation_evidence_score := clip01(
  0.35 * clip01(clean_cds_output_pct / 100) +
    0.25 * clip01(abs(clean_cds_allocation_delta) / 0.12) +
    0.20 * clip01(upstream_shift_proxy / 0.12) +
    0.10 * as.numeric(bool_col(audit, "atlas_interval_supported")) +
    0.10 * clip01(num_col(audit, "model_agreement_score",
                          num_col(audit, "agreement_model_agreement_score", 0)))
)]
audit[, disagreement_audit_class := fcase(
  clean_cds_output_pct >= 30 &
    clean_cds_allocation_delta <= -0.03 &
    compensation_evidence_score >= 0.55 &
    heterogeneity_risk_score < 0.45,
  "strong_buffering_or_total_load_compensation_candidate",
  clean_cds_output_pct >= 20 &
    clean_cds_allocation_delta <= -0.02 &
    heterogeneity_risk_score < 0.65,
  "possible_buffering_or_total_load_compensation",
  heterogeneity_risk_score >= 0.45,
  "aggregation_or_heterogeneity_risk",
  abs(clean_cds_output_pct) < 10 &
    abs(clean_cds_allocation_delta) >= 0.04,
  "weak_output_strong_allocation_review",
  default = "low_magnitude_disagreement_review"
)]
audit[, disagreement_review_priority := clip01(
  0.45 * compensation_evidence_score +
    0.25 * disagreement_magnitude_score +
    0.15 * clip01(num_col(audit, "atlas_browser_review_priority_score") / 1.5) +
    0.10 * clip01(num_col(audit, "context_interaction_component", 0)) +
    0.05 * as.numeric(bool_col(audit, "browser_supported")) -
    0.20 * heterogeneity_risk_score
)]
setorder(audit, -disagreement_review_priority, gene_symbol, design_family)
audit[, disagreement_review_rank := seq_len(.N)]

context_audit <- data.table()
if (nrow(branch_contexts)) {
  branch_contexts[, gene_symbol := as.character(gene_symbol)]
  branch_contexts[, design_family := as.character(design_family)]
  context_audit <- merge(
    branch_contexts,
    audit[, .(gene_symbol, tx_id, design_family,
              disagreement_review_rank,
              disagreement_audit_class)],
    by = key_cols,
    all.x = FALSE,
    sort = FALSE
  )
  if (nrow(context_audit)) {
    context_audit[, context_abs_disagreement := abs(
      num_col(context_audit, "atlas_clean_cds_pct_change", 0) / 100
    ) + abs(num_col(context_audit, "joint_delta_clean_CDS", 0))]
    setorder(context_audit, disagreement_review_rank,
             -joint_review_priority, -context_abs_disagreement)
    context_audit <- context_audit[
      ,
      head(.SD, 6L),
      by = .(gene_symbol, tx_id, design_family)
    ]
  }
}

study_audit <- data.table()
if (nrow(study_evidence)) {
  study_evidence[, gene_symbol := as.character(gene_symbol)]
  study_evidence[, design_family := as.character(design_family)]
  study_audit <- merge(
    study_evidence,
    audit[, .(gene_symbol, tx_id, design_family,
              disagreement_review_rank,
              disagreement_audit_class)],
    by = key_cols,
    all.x = FALSE,
    sort = FALSE
  )
  if (nrow(study_audit)) {
    setorder(study_audit, disagreement_review_rank, -study_top_priority)
    study_audit <- study_audit[
      ,
      head(.SD, 5L),
      by = .(gene_symbol, tx_id, design_family)
    ]
  }
}

if (nrow(context_audit)) {
  context_summary <- context_audit[
    ,
    .(
      n_context_rows_exported = .N,
      top_context_labels = collapse_contexts(.SD, 4L),
      top_context_clean_cds_delta = safe_max(abs(joint_delta_clean_CDS), 0),
      top_context_l1_shift = safe_max(joint_l1_allocation_shift, 0)
    ),
    by = key_cols
  ]
  audit <- merge(audit, context_summary, by = key_cols, all.x = TRUE,
                 sort = FALSE)
} else {
  audit[, `:=`(
    n_context_rows_exported = 0L,
    top_context_labels = "",
    top_context_clean_cds_delta = NA_real_,
    top_context_l1_shift = NA_real_
  )]
}

audit[is.na(n_context_rows_exported), n_context_rows_exported := 0L]
audit[is.na(top_context_labels), top_context_labels := ""]

review_cols <- intersect(
  c(
    "disagreement_review_rank", "gene_symbol", "tx_id", "design_family",
    "disagreement_audit_class", "disagreement_review_priority",
    "compensation_evidence_score", "heterogeneity_risk_score",
    "disagreement_magnitude_score", "clean_cds_output_pct",
    "clean_cds_output_pct_lower", "clean_cds_output_pct_upper",
    "clean_cds_output_log2fc", "clean_cds_allocation_delta",
    "leader_uorf_allocation_delta", "overlapping_uorf_allocation_delta",
    "other_allocation_delta", "upstream_allocation_delta",
    "total_load_compensation_proxy", "output_allocation_direction_pair",
    "top_shifted_feature_id", "top_shifted_feature_type",
    "hierarchical_joint_top_branch_class", "dm_top_branch_class",
    "model_agreement_tier", "model_agreement_class",
    "model_disagreement_flags", "clean_cds_output_allocation_alignment",
    "atlas_i2", "random_effects_i2", "hierarchical_joint_top_i2",
    "strict_supported_studies", "supported_studies", "n_studies",
    "best_case_condition", "best_control_conditions", "best_study",
    "top_context_labels", "rdg_png_file", "relative_usage_matrix_file",
    "warnings"
  ),
  names(audit)
)
review <- audit[, ..review_cols]
setorder(review, disagreement_review_rank)

summary_metrics <- data.table(
  metric = c(
    "n_disagreement_rows",
    "n_strong_compensation_candidates",
    "n_possible_compensation_candidates",
    "n_heterogeneity_risk_rows",
    "top_review_gene",
    "top_review_design",
    "top_review_class",
    "mean_compensation_evidence_score",
    "mean_heterogeneity_risk_score"
  ),
  value = c(
    as.character(nrow(audit)),
    as.character(sum(
      audit$disagreement_audit_class ==
        "strong_buffering_or_total_load_compensation_candidate"
    )),
    as.character(sum(
      audit$disagreement_audit_class ==
        "possible_buffering_or_total_load_compensation"
    )),
    as.character(sum(
      audit$disagreement_audit_class ==
        "aggregation_or_heterogeneity_risk"
    )),
    audit$gene_symbol[[1]],
    audit$design_family[[1]],
    audit$disagreement_audit_class[[1]],
    sprintf("%.3f", safe_mean(audit$compensation_evidence_score, 0)),
    sprintf("%.3f", safe_mean(audit$heterogeneity_risk_score, 0))
  )
)

fwrite(
  audit,
  file.path(output_dir,
            "dominant_rdg_allocation_output_disagreement_audit.csv")
)
fwrite(
  review,
  file.path(output_dir,
            "dominant_rdg_allocation_output_disagreement_review.csv")
)
fwrite(
  context_audit,
  file.path(output_dir,
            "dominant_rdg_allocation_output_disagreement_contexts.csv")
)
fwrite(
  study_audit,
  file.path(output_dir,
            "dominant_rdg_allocation_output_disagreement_study_evidence.csv")
)
fwrite(
  summary_metrics,
  file.path(output_dir,
            "dominant_rdg_allocation_output_disagreement_summary_metrics.csv")
)

plot_dt <- copy(review)
plot_dt[, gene_design := paste(gene_symbol, design_family, sep = " | ")]
plot_dt[, gene_design := factor(gene_design, levels = rev(gene_design))]
class_colors <- c(
  strong_buffering_or_total_load_compensation_candidate = "#1f9ed6",
  possible_buffering_or_total_load_compensation = "#1b9e77",
  weak_output_strong_allocation_review = "#cc4ccf",
  aggregation_or_heterogeneity_risk = "#f8766d",
  low_magnitude_disagreement_review = "#8c8c00"
)
class_labels <- c(
  strong_buffering_or_total_load_compensation_candidate = "strong compensation",
  possible_buffering_or_total_load_compensation = "possible compensation",
  weak_output_strong_allocation_review = "weak output / strong allocation",
  aggregation_or_heterogeneity_risk = "heterogeneity risk",
  low_magnitude_disagreement_review = "low magnitude"
)
p_scatter <- ggplot(
  plot_dt,
  aes(
    x = clean_cds_allocation_delta,
    y = clean_cds_output_pct,
    color = disagreement_audit_class,
    size = disagreement_review_priority
  )
  ) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_point(alpha = 0.88) +
  geom_text(aes(label = gene_symbol), nudge_y = 4, size = 3,
            show.legend = FALSE) +
  scale_color_manual(values = class_colors, labels = class_labels,
                     drop = FALSE) +
  scale_size_continuous(limits = c(0, 1), range = c(2.5, 6)) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    size = "none"
  ) +
  labs(
    title = "Allocation-Versus-Output Disagreement Audit",
    subtitle = wrap_text(
      "Rows in the upper-left quadrant have higher clean-CDS output despite lower clean-CDS branch allocation."
    ),
    x = "clean-CDS branch allocation delta",
    y = "clean-CDS output change (%)",
    color = "audit class"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9.5)
  )
ggsave(
  file.path(figure_dir,
            "allocation_output_disagreement_scatter.png"),
  p_scatter,
  width = 10.5,
  height = 6.8,
  dpi = 180
)
if (exists("dominant_save_ggplotly") &&
    requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    p_scatter,
    file.path(figure_dir,
              "allocation_output_disagreement_scatter.html"),
    title = "RDG allocation-output disagreement audit"
  )
}

heat_dt <- melt(
  review[, .(
    gene_symbol,
    design_family,
    clean_cds_allocation_delta,
    leader_uorf_allocation_delta,
    overlapping_uorf_allocation_delta,
    other_allocation_delta
  )],
  id.vars = c("gene_symbol", "design_family"),
  variable.name = "branch_metric",
  value.name = "allocation_delta"
)
heat_dt[, gene_design := paste(gene_symbol, design_family, sep = " | ")]
heat_dt[, gene_design := factor(gene_design, levels = rev(unique(
  plot_dt$gene_design
)))]
heat_dt[, branch_metric := factor(
  branch_metric,
  levels = c(
    "clean_cds_allocation_delta",
    "leader_uorf_allocation_delta",
    "overlapping_uorf_allocation_delta",
    "other_allocation_delta"
  ),
  labels = c("clean CDS", "leader uORF", "overlapping uORF", "other")
)]
p_heat <- ggplot(
  heat_dt,
  aes(x = branch_metric, y = gene_design, fill = allocation_delta)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(
    is.finite(allocation_delta),
    sprintf("%+.3f", allocation_delta),
    ""
  )), size = 2.7) +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#8b0000",
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    title = "Branch Allocation Shifts Behind Output Disagreements",
    subtitle = wrap_text(
      "Negative clean-CDS allocation paired with positive output suggests total-load compensation, buffering, or aggregation artifacts."
    ),
    x = NULL,
    y = NULL,
    fill = "allocation delta"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(figure_dir,
            "allocation_output_disagreement_branch_heatmap.png"),
  p_heat,
  width = 10.5,
  height = 5.8,
  dpi = 180
)
if (exists("dominant_save_ggplotly") &&
    requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    p_heat,
    file.path(figure_dir,
              "allocation_output_disagreement_branch_heatmap.html"),
    title = "RDG allocation-output branch heatmap"
  )
}

message("Saved allocation-output disagreement audit to: ", output_dir)
