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
output_dir <- file.path(analysis_dir, "dominant_rdg_motif_prior_calibration")
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

max_safe <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) 0 else max(x)
}

cor_safe <- function(x, y) {
  x <- safe_num(x)
  y <- safe_num(y)
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5 || length(unique(x[keep])) < 2 ||
      length(unique(y[keep])) < 2) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

collapse_unique <- function(x, n = 6L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

q_to_component <- function(q) {
  q <- safe_num(q)
  out <- rep(0.15, length(q))
  finite <- is.finite(q) & q > 0
  out[finite] <- pmin(1, -log10(q[finite]) / 3)
  out[is.finite(q) & q == 0] <- 1
  out
}

message("RDG motif-prior branch calibration:")
message("  1. Join bounded motif scores to count-aware branch rows")
message("  2. Compare motif prior strength to observed branch-allocation evidence")
message("  3. Rank concordant and discordant gene-designs for review")

motif_file <- file.path(
  analysis_dir, "dominant_rdg_mechanistic_motif_model",
  "dominant_rdg_mechanistic_motif_scores.csv"
)
branch_file <- file.path(
  analysis_dir, "dominant_rdg_branch_allocation",
  "dominant_rdg_branch_allocation_posterior_contrasts.csv"
)
atlas_file <- file.path(
  analysis_dir, "dominant_rdg_atlas",
  "dominant_rdg_atlas_cards.csv"
)

motif <- read_dt(motif_file, required = TRUE)
branch <- read_dt(branch_file, required = TRUE)
atlas <- read_dt(atlas_file)

motif_required <- c(
  "gene_symbol", "tx_id", "inhibitory_overlap_burden",
  "reinitiation_distance_opportunity", "dense_uorf_collision_burden",
  "strong_start_leakage_burden", "downstream_rescue_opportunity",
  "overlap_without_rescue", "leakage_with_reinitiation",
  "dense_short_reinitiation", "rescue_after_overlap",
  "motif_complexity_index"
)
missing <- setdiff(motif_required, names(motif))
if (length(missing)) {
  stop("Motif score table is missing required columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

branch_required <- c(
  "gene_symbol", "tx_id", "design_family", "study", "case_condition",
  "control_conditions", "branch_feature_class", "feature_id",
  "feature_type", "count_aware_usage_delta",
  "count_aware_usage_delta_q", "count_aware_evaluable",
  "count_aware_strict", "count_aware_row_weight",
  "case_sum_raw_counts", "control_sum_raw_counts"
)
missing <- setdiff(branch_required, names(branch))
if (length(missing)) {
  stop("Branch-allocation table is missing required columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

dt <- branch[
  branch_feature_class %chin% c("clean_CDS", "leader_uORF",
                                "overlapping_uORF") &
    count_aware_evaluable == TRUE
]
dt <- merge(dt, motif[, ..motif_required],
            by = c("gene_symbol", "tx_id"), all.x = TRUE, sort = FALSE)
for (column in setdiff(motif_required, c("gene_symbol", "tx_id"))) {
  dt[is.na(get(column)) | !is.finite(finite_or_zero(get(column))),
     (column) := 0]
}

dt[, `:=`(
  delta = finite_or_zero(count_aware_usage_delta),
  q = safe_num(count_aware_usage_delta_q),
  row_weight = pmax(finite_or_zero(count_aware_row_weight), 1),
  case_feature_counts = finite_or_zero(case_sum_raw_counts),
  control_feature_counts = finite_or_zero(control_sum_raw_counts),
  strict_row = count_aware_strict == TRUE
)]
dt[, q_component := q_to_component(q)]
dt[, observed_abs_delta := abs(delta)]
dt[, observed_branch_evidence_score := clip01(
  observed_abs_delta / 0.15
) * clip01(log1p(row_weight) / log1p(80)) *
  pmax(q_component, 0.15)]

dt[, motif_prior_signed := fcase(
  branch_feature_class == "clean_CDS",
  clip01(0.45 * reinitiation_distance_opportunity +
           0.35 * downstream_rescue_opportunity +
           0.20 * leakage_with_reinitiation) -
    clip01(0.45 * overlap_without_rescue +
             0.35 * dense_short_reinitiation +
             0.20 * strong_start_leakage_burden),
  branch_feature_class == "overlapping_uORF",
  clip01(0.45 * inhibitory_overlap_burden +
           0.35 * overlap_without_rescue +
           0.20 * rescue_after_overlap),
  branch_feature_class == "leader_uORF",
  clip01(0.35 * dense_uorf_collision_burden +
           0.30 * strong_start_leakage_burden +
           0.20 * leakage_with_reinitiation +
           0.15 * reinitiation_distance_opportunity),
  default = 0
)]
dt[, motif_prior_strength := abs(motif_prior_signed)]
dt[, motif_prior_direction := fifelse(
  motif_prior_signed > 0.10,
  1L,
  fifelse(motif_prior_signed < -0.10, -1L, 0L)
)]
dt[, observed_direction := fifelse(
  delta > 0.02,
  1L,
  fifelse(delta < -0.02, -1L, 0L)
)]
dt[, motif_prior_concordance_class := fcase(
  motif_prior_strength < 0.20,
  "low_prior",
  observed_branch_evidence_score < 0.12,
  "high_prior_weak_branch_evidence",
  motif_prior_direction != 0L & observed_direction != 0L &
    motif_prior_direction == observed_direction,
  "motif_branch_concordant",
  motif_prior_direction != 0L & observed_direction != 0L &
    motif_prior_direction != observed_direction,
  "motif_branch_discordant",
  default = "ambiguous_direction"
)]
dt[, `:=`(
  motif_branch_concordance_score = fifelse(
    motif_prior_concordance_class == "motif_branch_concordant",
    motif_prior_strength * observed_branch_evidence_score,
    0
  ),
  motif_branch_discordance_score = fifelse(
    motif_prior_concordance_class == "motif_branch_discordant",
    motif_prior_strength * observed_branch_evidence_score,
    0
  ),
  motif_unseen_prior_score = fifelse(
    motif_prior_concordance_class == "high_prior_weak_branch_evidence",
    motif_prior_strength * (1 - observed_branch_evidence_score),
    0
  )
)]

branch_rows_file <- file.path(
  output_dir,
  "dominant_rdg_motif_prior_branch_rows.csv"
)
fwrite(dt, branch_rows_file)

summary <- rbindlist(list(
  dt[, .(
    n_rows = .N,
    n_gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    mean_prior_strength = mean_safe(motif_prior_strength),
    median_prior_strength = median_safe(motif_prior_strength),
    mean_observed_abs_delta = mean_safe(observed_abs_delta),
    median_observed_abs_delta = median_safe(observed_abs_delta),
    spearman_prior_vs_abs_delta =
      cor_safe(motif_prior_strength, observed_abs_delta),
    spearman_prior_vs_evidence_score =
      cor_safe(motif_prior_strength, observed_branch_evidence_score),
    concordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_concordant"),
    discordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_discordant"),
    high_prior_weak_evidence_rows =
      sum(motif_prior_concordance_class ==
            "high_prior_weak_branch_evidence"),
    strict_rows = sum(strict_row, na.rm = TRUE)
  )][, scope := "all"],
  dt[, .(
    n_rows = .N,
    n_gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    mean_prior_strength = mean_safe(motif_prior_strength),
    median_prior_strength = median_safe(motif_prior_strength),
    mean_observed_abs_delta = mean_safe(observed_abs_delta),
    median_observed_abs_delta = median_safe(observed_abs_delta),
    spearman_prior_vs_abs_delta =
      cor_safe(motif_prior_strength, observed_abs_delta),
    spearman_prior_vs_evidence_score =
      cor_safe(motif_prior_strength, observed_branch_evidence_score),
    concordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_concordant"),
    discordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_discordant"),
    high_prior_weak_evidence_rows =
      sum(motif_prior_concordance_class ==
            "high_prior_weak_branch_evidence"),
    strict_rows = sum(strict_row, na.rm = TRUE)
  ), by = .(scope = paste0("branch_feature_class:", branch_feature_class))],
  dt[, .(
    n_rows = .N,
    n_gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
    mean_prior_strength = mean_safe(motif_prior_strength),
    median_prior_strength = median_safe(motif_prior_strength),
    mean_observed_abs_delta = mean_safe(observed_abs_delta),
    median_observed_abs_delta = median_safe(observed_abs_delta),
    spearman_prior_vs_abs_delta =
      cor_safe(motif_prior_strength, observed_abs_delta),
    spearman_prior_vs_evidence_score =
      cor_safe(motif_prior_strength, observed_branch_evidence_score),
    concordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_concordant"),
    discordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_discordant"),
    high_prior_weak_evidence_rows =
      sum(motif_prior_concordance_class ==
            "high_prior_weak_branch_evidence"),
    strict_rows = sum(strict_row, na.rm = TRUE)
  ), by = .(scope = paste0("design_family:", design_family))]
), fill = TRUE)
setcolorder(summary, c("scope", setdiff(names(summary), "scope")))
fwrite(summary, file.path(output_dir,
                          "dominant_rdg_motif_prior_calibration_summary.csv"))

gene_design <- dt[, {
  top_idx <- which.max(pmax(motif_branch_concordance_score,
                           motif_branch_discordance_score,
                           motif_unseen_prior_score))
  if (!length(top_idx) || is.na(top_idx)) top_idx <- 1L
  data.table(
    motif_prior_rows = .N,
    motif_prior_strict_rows = sum(strict_row, na.rm = TRUE),
    motif_prior_max_strength = max_safe(motif_prior_strength),
    motif_prior_mean_strength = mean_safe(motif_prior_strength),
    motif_observed_max_abs_delta = max_safe(observed_abs_delta),
    motif_observed_mean_abs_delta = mean_safe(observed_abs_delta),
    motif_observed_max_evidence_score =
      max_safe(observed_branch_evidence_score),
    motif_concordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_concordant"),
    motif_discordant_rows =
      sum(motif_prior_concordance_class == "motif_branch_discordant"),
    motif_high_prior_weak_evidence_rows =
      sum(motif_prior_concordance_class ==
            "high_prior_weak_branch_evidence"),
    motif_concordance_score = max_safe(motif_branch_concordance_score),
    motif_discordance_score = max_safe(motif_branch_discordance_score),
    motif_unseen_prior_score = max_safe(motif_unseen_prior_score),
    motif_net_alignment_score =
      max_safe(motif_branch_concordance_score) -
      max_safe(motif_branch_discordance_score),
    motif_top_feature_id = feature_id[top_idx],
    motif_top_feature_class = branch_feature_class[top_idx],
    motif_top_concordance_class =
      motif_prior_concordance_class[top_idx],
    motif_top_delta = delta[top_idx],
    motif_top_prior_signed = motif_prior_signed[top_idx],
    motif_top_prior_strength = motif_prior_strength[top_idx],
    motif_top_observed_evidence_score =
      observed_branch_evidence_score[top_idx],
    motif_top_context = paste(study[top_idx], case_condition[top_idx],
                              "vs", control_conditions[top_idx],
                              sep = " | "),
    motif_top_studies = collapse_unique(study, 6L),
    motif_top_conditions = collapse_unique(case_condition, 6L)
  )
}, by = .(gene_symbol, tx_id, design_family)]

if (nrow(atlas)) {
  atlas_keep <- intersect(names(atlas), c(
    "gene_symbol", "tx_id", "design_family",
    "atlas_evidence_class", "atlas_browser_review_priority_score",
    "atlas_clean_cds_pct_change", "browser_validation_status",
    "browser_validation_label", "browser_validation_notes",
    "count_aware_top_support_class", "count_aware_top_score",
    "hierarchical_evidence_class", "hierarchical_priority_score",
    "browser_review_rank", "browser_review_question"
  ))
  atlas_small <- unique(atlas[, ..atlas_keep])
  gene_design <- merge(gene_design, atlas_small,
                       by = c("gene_symbol", "tx_id", "design_family"),
                       all.x = TRUE, sort = FALSE)
}

gene_design[, motif_prior_review_class := fcase(
  motif_concordance_score >= 0.12 &
    motif_top_concordance_class == "motif_branch_concordant",
  "motif_prior_concordant_branch_evidence",
  motif_discordance_score >= 0.12 &
    motif_top_concordance_class == "motif_branch_discordant",
  "motif_prior_discordant_branch_evidence",
  motif_unseen_prior_score >= 0.45,
  "high_motif_prior_weak_branch_evidence",
  motif_prior_max_strength >= 0.45 &
    motif_observed_max_evidence_score < 0.12,
  "motif_prior_unresolved",
  default = "motif_prior_low_or_ambiguous"
)]
gene_design[, motif_prior_review_priority := clip01(
  0.35 * motif_concordance_score +
    0.35 * motif_discordance_score +
    0.20 * motif_unseen_prior_score +
    0.10 * clip01(finite_or_zero(hierarchical_priority_score))
)]
gene_design[, abs_motif_net_alignment_score := abs(motif_net_alignment_score)]
setorder(gene_design, -motif_prior_review_priority,
         -abs_motif_net_alignment_score, gene_symbol, design_family)

gene_design_file <- file.path(
  output_dir,
  "dominant_rdg_motif_prior_gene_design_calibration.csv"
)
fwrite(gene_design, gene_design_file)
fwrite(
  gene_design[design_family == "viral_infection"],
  file.path(output_dir,
            "dominant_rdg_motif_prior_viral_review.csv")
)
fwrite(
  gene_design[motif_prior_review_class != "motif_prior_low_or_ambiguous"],
  file.path(output_dir,
            "dominant_rdg_motif_prior_review_queue.csv")
)

plot_dt <- dt[
  motif_prior_strength >= 0.05 &
    is.finite(observed_branch_evidence_score)
]
if (nrow(plot_dt)) {
  plot_dt[, plot_class := fifelse(
    motif_prior_concordance_class == "motif_branch_concordant",
    "concordant",
    fifelse(motif_prior_concordance_class == "motif_branch_discordant",
            "discordant",
            fifelse(motif_prior_concordance_class ==
                      "high_prior_weak_branch_evidence",
                    "high prior, weak evidence",
                    "other"))
  )]
  p <- ggplot(
    plot_dt,
    aes(x = motif_prior_strength,
        y = observed_branch_evidence_score,
        color = plot_class,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Design: ", design_family,
          "<br>Feature: ", feature_id,
          "<br>Class: ", branch_feature_class,
          "<br>Prior: ", round(motif_prior_strength, 3),
          "<br>Delta: ", round(delta, 3),
          "<br>Evidence score: ",
          round(observed_branch_evidence_score, 3),
          "<br>Concordance: ", motif_prior_concordance_class
        ))
  ) +
    geom_point(alpha = 0.45, size = 1.6) +
    facet_wrap(~branch_feature_class) +
    scale_color_manual(values = c(
      "concordant" = "#0f766e",
      "discordant" = "#b91c1c",
      "high prior, weak evidence" = "#b45309",
      "other" = "#6b7280"
    )) +
    labs(
      title = "RDG Motif Priors Versus Count-Aware Branch Evidence",
      subtitle = "Rows are branch-allocation contrasts; priors are gene-level motif scores",
      x = "Motif prior strength",
      y = "Observed branch-evidence score",
      color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "motif_prior_vs_branch_evidence.png"),
         p, width = 10.4, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "motif_prior_vs_branch_evidence.pdf"),
         p, width = 10.4, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "motif_prior_vs_branch_evidence.html"),
      title = "Motif priors versus branch evidence"
    )
  }
}

viral_plot <- gene_design[design_family == "viral_infection"]
if (nrow(viral_plot)) {
  actionable_viral <- viral_plot[
    motif_prior_review_class != "motif_prior_low_or_ambiguous"
  ]
  if (nrow(actionable_viral)) {
    viral_plot <- actionable_viral
  }
  viral_plot <- viral_plot[
    order(-motif_prior_review_priority)
  ][seq_len(min(.N, 30))]
  viral_plot[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  viral_plot[, review_plot_class := fifelse(
    grepl("concordant", motif_prior_review_class),
    "concordant",
    fifelse(grepl("discordant", motif_prior_review_class),
            "discordant",
            fifelse(grepl("weak", motif_prior_review_class),
                    "high prior, weak evidence",
                    "other"))
  )]
  p2 <- ggplot(
    viral_plot,
    aes(x = motif_prior_review_priority, y = gene_symbol,
        fill = review_plot_class,
        text = paste0(
          "Gene: ", gene_symbol,
          "<br>Review class: ", motif_prior_review_class,
          "<br>Review priority: ",
          round(motif_prior_review_priority, 3),
          "<br>Net alignment: ", round(motif_net_alignment_score, 3),
          "<br>Top feature: ", motif_top_feature_id,
          "<br>Top class: ", motif_top_feature_class,
          "<br>Prior strength: ", round(motif_top_prior_strength, 3),
          "<br>Delta: ", round(motif_top_delta, 3),
          "<br>Evidence score: ",
          round(motif_top_observed_evidence_score, 3),
          "<br>Context: ", motif_top_context
        ))
  ) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(
      "concordant" = "#0f766e",
      "discordant" = "#b91c1c",
      "high prior, weak evidence" = "#b45309",
      "other" = "#6b7280"
    )) +
    labs(
      title = "Viral-Infection Motif Prior Review Queue",
      subtitle = "High-prior weak evidence marks missing branch evidence; discordant marks a model-failure review target",
      x = "Motif-prior review priority",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_dir, "motif_prior_viral_alignment.png"),
         p2, width = 9.6, height = 7.2, dpi = 180)
  ggsave(file.path(figure_dir, "motif_prior_viral_alignment.pdf"),
         p2, width = 9.6, height = 7.2)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p2,
      file.path(figure_dir, "motif_prior_viral_alignment.html"),
      title = "Viral motif prior calibration"
    )
  }
}

message("Saved motif-prior calibration outputs in: ", output_dir)
print(summary[scope %chin% c("all", "design_family:viral_infection")])
