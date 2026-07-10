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
output_dir <- file.path(
  analysis_dir, "dominant_rdg_therapeutic_perturbation_predictions"
)
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

num_col <- function(dt, column, default = 0) {
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

split_genes <- function(x) {
  out <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  toupper(out[nzchar(out)])
}

collapse_unique <- function(x, n = 4L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

safe_mean <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else mean(x)
}

safe_median <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else median(x)
}

safe_max <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

weighted_mean_safe <- function(x, w = NULL, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (is.null(w)) w <- rep(1, length(x))
  w <- suppressWarnings(as.numeric(w))
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(default)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

first_finite <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (!length(values)) 0 else values[[1]]
}

invert_interval <- function(lower, upper, sign = -1, scale = 1) {
  lower <- suppressWarnings(as.numeric(lower))
  upper <- suppressWarnings(as.numeric(upper))
  sign <- suppressWarnings(as.numeric(sign))
  if (!is.finite(lower) || !is.finite(upper) || !is.finite(sign)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  values <- c(sign * scale * lower, sign * scale * upper)
  c(lower = min(values), upper = max(values))
}

wrap_text <- function(x, width = 92L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

therapeutic_dir <- file.path(analysis_dir, "dominant_rdg_therapeutic_hypotheses")
scores <- read_dt(
  file.path(therapeutic_dir,
            "dominant_rdg_therapeutic_hypothesis_scores.csv"),
  required = TRUE
)
atlas <- read_dt(
  file.path(analysis_dir, "dominant_rdg_atlas",
            "dominant_rdg_atlas_cards.csv"),
  required = TRUE
)
if (!nrow(scores)) stop("Therapeutic-hypothesis score table is empty.",
                        call. = FALSE)
if (!nrow(atlas)) stop("Atlas card table is empty.", call. = FALSE)

scores[, therapeutic_hypothesis := as.character(therapeutic_hypothesis)]
if ("gene_symbol" %in% names(scores)) scores[, gene_symbol := NULL]
atlas[, gene_symbol := toupper(as.character(gene_symbol))]
atlas[, design_family := as.character(design_family)]

prediction_rules <- data.table(
  therapeutic_hypothesis = c(
    "ISR_eIF2B_ATF4_attenuation",
    "JAK_STAT_IFN_attenuation",
    "IFIH1_MDA5_axis_modulation",
    "HIF_hypoxia_angiogenic_modulation",
    "mTOR_translation_capacity_modulation",
    "mitochondrial_OXPHOS_support",
    "immune_antigen_complement_axis"
  ),
  disease_design_families = c(
    "isr_er_translation_stress; viral_infection; nutrient_starvation",
    "viral_infection; acute_interferon_cytokine",
    "viral_infection; acute_interferon_cytokine",
    "hypoxia_oxidative_stress; tumor_cancer_context; viral_infection",
    "nutrient_starvation; drug_or_stimulus; viral_infection; gain_of_function",
    "viral_infection; hypoxia_oxidative_stress; nutrient_starvation; tumor_cancer_context",
    "viral_infection; acute_interferon_cytokine; tumor_cancer_context"
  ),
  perturbation_assumption = c(
    "attenuate chronic ISR/eIF2B-ATF4 disease-like effect toward baseline",
    "attenuate persistent JAK/STAT IFN disease-like effect toward baseline",
    "stratify or dampen MDA5/IFIH1 disease-like effect toward baseline",
    "normalize hypoxia/HIF disease-like branch effects toward baseline",
    "context-dependent translation-capacity modulation toward baseline",
    "restore/normalize mitochondrial disease-like branch effects toward baseline",
    "stratify or attenuate immune-antigen/complement disease-like effects"
  ),
  perturbation_sign = c(-1, -1, -1, -1, -1, -1, -1),
  half_normalization_fraction = c(0.50, 0.50, 0.40, 0.40, 0.30, 0.35, 0.35),
  validation_assay = c(
    "ISR perturbation in infected/stressed cells; measure DDIT3/ATF4/PPP1R15A RDG allocation and protein output.",
    "JAK/STAT-axis perturbation in post-acute IFN-high cells; measure IFIH1/STAT1/IRF7 RDG allocation and CDS output.",
    "IFIH1/MDA5 stratified perturbation; compare lung-like and non-lung infected contexts before considering intervention.",
    "Hypoxia or HIF perturbation; separate oxygen-state biology from stress-linked translation rewiring.",
    "mTOR/translation-capacity perturbation; first test whether branch effects are causal or broad translation passengers.",
    "Respiration-state perturbation; measure OXPHOS branch allocation and downstream CDS buffering.",
    "Immune-enrichment audit plus perturbation; separate cell-composition effects from direct RDG rewiring."
  )
)

hypothesis_genes <- scores[
  ,
  .(gene_symbol = split_genes(primary_genes)),
  by = .(
    therapeutic_hypothesis,
    hypothesis_label,
    therapeutic_hypothesis_score,
    trial_readiness_class,
    intervention_stage,
    example_interventions,
    primary_genes,
    caveat
  )
]
hypothesis_genes <- merge(
  hypothesis_genes,
  prediction_rules,
  by = "therapeutic_hypothesis",
  all.x = TRUE,
  sort = FALSE
)
hypothesis_designs <- hypothesis_genes[
  ,
  .(design_family = trimws(unlist(strsplit(
    disease_design_families, ";", fixed = TRUE
  )))),
  by = names(hypothesis_genes)
]
hypothesis_designs <- hypothesis_designs[nzchar(design_family)]

warnings_chr <- chr_col(atlas, "warnings")
atlas[, interval_supported := bool_col(atlas, "atlas_interval_supported")]
atlas[, direction_stable := bool_col(atlas, "atlas_direction_stable")]
atlas[, interval_excludes_zero := (
  num_col(atlas, "atlas_clean_cds_pct_lower") > 0 |
    num_col(atlas, "atlas_clean_cds_pct_upper") < 0
)]
atlas[, browser_supported := num_col(atlas, "browser_validation_score") > 0 |
        grepl("browser_supported", warnings_chr, ignore.case = TRUE)]
atlas[, sample_component := clip01(
  log10(num_col(atlas, "total_case_samples") +
          num_col(atlas, "total_control_samples") + 1) / 3
)]
atlas[, row_prediction_evidence_score := clip01(
  0.14 * clip01(num_col(atlas, "atlas_confidence_score")) +
    0.12 * clip01(num_col(atlas, "atlas_browser_review_priority_score") / 1.5) +
    0.14 * clip01(num_col(atlas, "model_agreement_component")) +
    0.10 * clip01(num_col(atlas, "loo_stability_component")) +
    0.10 * clip01(num_col(atlas, "transport_stability_component")) +
    0.10 * clip01(num_col(atlas, "context_interaction_component")) +
    0.08 * clip01(num_col(atlas, "hierarchical_joint_allocation_component")) +
    0.08 * clip01(num_col(atlas, "dirichlet_multinomial_component")) +
    0.06 * as.numeric(interval_supported) +
    0.04 * as.numeric(browser_supported) +
    0.04 * sample_component
)]

atlas[, branch_delta_clean_CDS := mapply(
  first_finite,
  num_col(atlas, "hierarchical_joint_delta_clean_CDS"),
  num_col(atlas, "dm_delta_clean_CDS"),
  num_col(atlas, "joint_weighted_clean_CDS_delta"),
  num_col(atlas, "count_aware_weighted_usage_delta_clean_CDS"),
  SIMPLIFY = TRUE
)]
atlas[, branch_delta_leader_uORF := mapply(
  first_finite,
  num_col(atlas, "hierarchical_joint_delta_leader_uORF"),
  num_col(atlas, "dm_delta_leader_uORF"),
  num_col(atlas, "joint_weighted_leader_uORF_delta"),
  num_col(atlas, "count_aware_weighted_usage_delta_leader_uORF"),
  SIMPLIFY = TRUE
)]
atlas[, branch_delta_overlapping_uORF := mapply(
  first_finite,
  num_col(atlas, "hierarchical_joint_delta_overlapping_uORF"),
  num_col(atlas, "dm_delta_overlapping_uORF"),
  num_col(atlas, "joint_weighted_overlapping_uORF_delta"),
  num_col(atlas, "count_aware_weighted_usage_delta_overlapping_uORF"),
  SIMPLIFY = TRUE
)]

atlas_keep <- intersect(
  c(
    "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
    "atlas_clean_cds_pct_change", "atlas_clean_cds_pct_lower",
    "atlas_clean_cds_pct_upper", "abs_atlas_clean_cds_pct_change",
    "atlas_direction", "n_studies", "strict_supported_studies",
    "total_case_samples", "total_control_samples", "best_case_condition",
    "best_control_conditions", "best_study", "top_shifted_feature_id",
    "top_shifted_feature_type", "top_shifted_feature_family",
    "hierarchical_joint_top_branch_class", "dm_top_branch_class",
    "joint_top_branch_class", "branch_delta_clean_CDS",
    "branch_delta_leader_uORF", "branch_delta_overlapping_uORF",
    "model_agreement_tier", "loo_consensus_class",
    "transport_consensus_class", "context_interaction_gene_class",
    "context_interaction_top_axis", "context_interaction_top_value",
    "context_interaction_top_branch_class", "context_interaction_summary",
    "browser_validation_label", "warnings", "interval_supported",
    "interval_excludes_zero", "browser_supported",
    "row_prediction_evidence_score"
  ),
  names(atlas)
)

candidate_rows <- merge(
  hypothesis_designs,
  atlas[, ..atlas_keep],
  by = c("gene_symbol", "design_family"),
  all.x = TRUE,
  sort = FALSE,
  allow.cartesian = TRUE
)
candidate_rows[, missing_prediction_row := is.na(tx_id)]
candidate_rows[missing_prediction_row == TRUE,
               row_prediction_evidence_score := 0]
candidate_rows[missing_prediction_row == TRUE,
               atlas_clean_cds_pct_change := NA_real_]
candidate_rows[missing_prediction_row == TRUE,
               atlas_clean_cds_pct_lower := NA_real_]
candidate_rows[missing_prediction_row == TRUE,
               atlas_clean_cds_pct_upper := NA_real_]
candidate_rows[is.na(interval_supported), interval_supported := FALSE]
candidate_rows[is.na(interval_excludes_zero), interval_excludes_zero := FALSE]
candidate_rows[is.na(browser_supported), browser_supported := FALSE]

setorder(
  candidate_rows,
  therapeutic_hypothesis,
  gene_symbol,
  -row_prediction_evidence_score,
  -abs_atlas_clean_cds_pct_change
)
gene_predictions <- candidate_rows[
  ,
  .SD[1],
  by = .(therapeutic_hypothesis, gene_symbol)
]
gene_predictions[, disease_like_cds_pct_change :=
                   suppressWarnings(as.numeric(atlas_clean_cds_pct_change))]
gene_predictions[, disease_like_cds_pct_lower :=
                   suppressWarnings(as.numeric(atlas_clean_cds_pct_lower))]
gene_predictions[, disease_like_cds_pct_upper :=
                   suppressWarnings(as.numeric(atlas_clean_cds_pct_upper))]
gene_predictions[!is.finite(disease_like_cds_pct_change),
                 disease_like_cds_pct_change := NA_real_]
gene_predictions[, full_normalization_cds_pct_change :=
                   perturbation_sign * disease_like_cds_pct_change]
gene_predictions[, half_normalization_cds_pct_change :=
                   half_normalization_fraction *
                   full_normalization_cds_pct_change]

full_intervals <- mapply(
  invert_interval,
  gene_predictions$disease_like_cds_pct_lower,
  gene_predictions$disease_like_cds_pct_upper,
  gene_predictions$perturbation_sign,
  MoreArgs = list(scale = 1),
  SIMPLIFY = TRUE
)
half_intervals <- mapply(
  invert_interval,
  gene_predictions$disease_like_cds_pct_lower,
  gene_predictions$disease_like_cds_pct_upper,
  gene_predictions$perturbation_sign,
  gene_predictions$half_normalization_fraction,
  SIMPLIFY = TRUE
)
gene_predictions[, full_normalization_cds_pct_lower := full_intervals["lower", ]]
gene_predictions[, full_normalization_cds_pct_upper := full_intervals["upper", ]]
gene_predictions[, half_normalization_cds_pct_lower := half_intervals["lower", ]]
gene_predictions[, half_normalization_cds_pct_upper := half_intervals["upper", ]]

for (branch in c("clean_CDS", "leader_uORF", "overlapping_uORF")) {
  src <- paste0("branch_delta_", branch)
  full <- paste0("full_normalization_branch_delta_", branch)
  half <- paste0("half_normalization_branch_delta_", branch)
  gene_predictions[, (full) := perturbation_sign * suppressWarnings(
    as.numeric(get(src))
  )]
  gene_predictions[, (half) := half_normalization_fraction * get(full)]
}

gene_predictions[, top_branch_class := fifelse(
  !is.na(hierarchical_joint_top_branch_class) &
    nzchar(hierarchical_joint_top_branch_class),
  hierarchical_joint_top_branch_class,
  fifelse(
    !is.na(dm_top_branch_class) & nzchar(dm_top_branch_class),
    dm_top_branch_class,
    fifelse(
      !is.na(joint_top_branch_class) & nzchar(joint_top_branch_class),
      joint_top_branch_class,
      as.character(top_shifted_feature_type)
    )
  )
)]
gene_predictions[, top_branch_class := as.character(top_branch_class)]
gene_predictions[is.na(top_branch_class), top_branch_class := ""]

gene_predictions[, quantitative_prediction_confidence := clip01(
  0.45 * row_prediction_evidence_score +
    0.20 * clip01(suppressWarnings(as.numeric(therapeutic_hypothesis_score))) +
    0.13 * as.numeric(interval_excludes_zero) +
    0.09 * as.numeric(browser_supported) +
    0.08 * clip01(suppressWarnings(as.numeric(n_studies)) / 5) +
    0.05 * clip01(abs(disease_like_cds_pct_change) / 100)
)]
gene_predictions[, prediction_support_class := fcase(
  missing_prediction_row,
  "missing_target_design",
  quantitative_prediction_confidence >= 0.70 & interval_excludes_zero,
  "quantitative_review_ready",
  quantitative_prediction_confidence >= 0.55,
  "directional_prediction",
  quantitative_prediction_confidence >= 0.40,
  "exploratory_prediction",
  default = "low_support_prediction"
)]
gene_predictions[, abs_half_normalization_cds_pct_change :=
                   abs(half_normalization_cds_pct_change)]
setorder(
  gene_predictions,
  -quantitative_prediction_confidence,
  -abs_half_normalization_cds_pct_change,
  therapeutic_hypothesis,
  gene_symbol
)
gene_predictions[, perturbation_prediction_rank := seq_len(.N)]

hypothesis_summary <- gene_predictions[
  missing_prediction_row == FALSE,
  .(
    n_prediction_genes = .N,
    n_quantitative_review_ready = sum(
      prediction_support_class == "quantitative_review_ready"
    ),
    n_directional_or_better = sum(
      prediction_support_class %chin% c(
        "quantitative_review_ready", "directional_prediction"
      )
    ),
    mean_prediction_confidence = safe_mean(
      quantitative_prediction_confidence
    ),
    median_abs_full_normalization_cds_pct_change = safe_median(
      abs(full_normalization_cds_pct_change)
    ),
    median_abs_half_normalization_cds_pct_change = safe_median(
      abs(half_normalization_cds_pct_change)
    ),
    top_gene_predictions = {
      x <- .SD[order(-quantitative_prediction_confidence)]
      paste(
        head(
          paste0(
            x$gene_symbol,
            " ",
            sprintf("%+.1f%%", x$half_normalization_cds_pct_change),
            " (conf ",
            sprintf("%.2f", x$quantitative_prediction_confidence),
            ")"
          ),
          6L
        ),
        collapse = "; "
      )
    },
    top_validation_assays = collapse_unique(validation_assay, 2L)
  ),
  by = .(
    therapeutic_hypothesis,
    hypothesis_label,
    therapeutic_hypothesis_score,
    trial_readiness_class,
    perturbation_assumption
  )
]
hypothesis_summary[, hypothesis_prediction_score := clip01(
  0.35 * clip01(suppressWarnings(as.numeric(therapeutic_hypothesis_score))) +
    0.35 * mean_prediction_confidence +
    0.15 * clip01(n_directional_or_better / pmax(1, n_prediction_genes)) +
    0.15 * clip01(median_abs_half_normalization_cds_pct_change / 50)
)]
setorder(hypothesis_summary, -hypothesis_prediction_score)
hypothesis_summary[, perturbation_hypothesis_rank := seq_len(.N)]
setcolorder(
  hypothesis_summary,
  c(
    "perturbation_hypothesis_rank",
    setdiff(names(hypothesis_summary), "perturbation_hypothesis_rank")
  )
)

review <- gene_predictions[
  prediction_support_class %chin% c(
    "quantitative_review_ready", "directional_prediction",
    "exploratory_prediction"
  )
][
  order(-quantitative_prediction_confidence,
        -abs(half_normalization_cds_pct_change))
]
review <- review[
  ,
  .(
    perturbation_prediction_rank,
    hypothesis_label,
    gene_symbol,
    design_family,
    prediction_support_class,
    quantitative_prediction_confidence,
    disease_like_cds_pct_change,
    disease_like_cds_pct_lower,
    disease_like_cds_pct_upper,
    full_normalization_cds_pct_change,
    full_normalization_cds_pct_lower,
    full_normalization_cds_pct_upper,
    half_normalization_cds_pct_change,
    half_normalization_cds_pct_lower,
    half_normalization_cds_pct_upper,
    top_branch_class,
    half_normalization_branch_delta_clean_CDS,
    half_normalization_branch_delta_leader_uORF,
    half_normalization_branch_delta_overlapping_uORF,
    best_case_condition,
    best_control_conditions,
    best_study,
    context_interaction_summary,
    perturbation_assumption,
    validation_assay,
    caveat
  )
]

summary_metrics <- data.table(
  metric = c(
    "n_gene_predictions",
    "n_non_missing_gene_predictions",
    "n_quantitative_review_ready",
    "n_directional_or_better",
    "top_hypothesis",
    "top_hypothesis_prediction_score",
    "top_gene_prediction"
  ),
  value = c(
    as.character(nrow(gene_predictions)),
    as.character(sum(!gene_predictions$missing_prediction_row)),
    as.character(sum(
      gene_predictions$prediction_support_class == "quantitative_review_ready"
    )),
    as.character(sum(
      gene_predictions$prediction_support_class %chin% c(
        "quantitative_review_ready", "directional_prediction"
      )
    )),
    hypothesis_summary$hypothesis_label[[1]],
    sprintf("%.3f", hypothesis_summary$hypothesis_prediction_score[[1]]),
    paste0(
      review$gene_symbol[[1]], " / ", review$hypothesis_label[[1]], " / ",
      sprintf("%+.1f%% half-normalization CDS",
              review$half_normalization_cds_pct_change[[1]])
    )
  )
)

gene_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_perturbation_gene_predictions.csv"
)
summary_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_perturbation_hypothesis_summary.csv"
)
review_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_perturbation_review.csv"
)
metrics_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_perturbation_summary_metrics.csv"
)
fwrite(gene_predictions, gene_file)
fwrite(hypothesis_summary, summary_file)
fwrite(review, review_file)
fwrite(summary_metrics, metrics_file)

plot_review <- head(copy(review), 24L)
plot_review[, gene_label := paste(gene_symbol, hypothesis_label, sep = " | ")]
plot_review[, gene_label := factor(gene_label, levels = rev(gene_label))]
p_effects <- ggplot(
  plot_review,
  aes(
    x = half_normalization_cds_pct_change,
    y = gene_label,
    color = prediction_support_class
  )
) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_errorbar(
    aes(
      xmin = half_normalization_cds_pct_lower,
      xmax = half_normalization_cds_pct_upper
    ),
    width = 0.18,
    alpha = 0.65,
    orientation = "y"
  ) +
  geom_point(aes(size = quantitative_prediction_confidence), alpha = 0.9) +
  scale_color_manual(
    values = c(
      quantitative_review_ready = "#1b9e77",
      directional_prediction = "#d95f02",
      exploratory_prediction = "#7570b3",
      low_support_prediction = "#bdbdbd"
    ),
    drop = FALSE
  ) +
  scale_size_continuous(limits = c(0, 1), range = c(1.8, 5.2)) +
  guides(
    size = "none",
    color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    title = "Therapeutic Perturbation Prediction Targets",
    subtitle = wrap_text(
      "Effects are half-normalization targets: the predicted CDS change if a perturbation reverses half of the disease-like atlas effect."
    ),
    x = "predicted CDS buffering change under half-normalization",
    y = NULL,
    color = "support"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(figure_dir, "therapeutic_perturbation_top_gene_effects.png"),
  p_effects,
  width = 11,
  height = 7.5,
  dpi = 180
)
if (exists("dominant_save_ggplotly") &&
    requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    p_effects,
    file.path(figure_dir, "therapeutic_perturbation_top_gene_effects.html"),
    title = "RDG therapeutic perturbation target effects"
  )
}

plot_summary <- copy(hypothesis_summary)
plot_summary[, hypothesis_label := factor(
  hypothesis_label,
  levels = rev(hypothesis_summary$hypothesis_label)
)]
p_summary <- ggplot(
  plot_summary,
  aes(
    x = hypothesis_prediction_score,
    y = hypothesis_label,
    fill = trial_readiness_class
  )
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = sprintf("%.2f", hypothesis_prediction_score)),
    hjust = -0.12,
    size = 3.2
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(
    values = c(
      clinical_biomarker_trial_candidate = "#1b9e77",
      preclinical_mechanistic_perturbation_candidate = "#d95f02",
      target_discovery_or_stratification_candidate = "#7570b3",
      exploratory_preclinical_candidate = "#4d4d4d",
      background_or_marker_only = "#bdbdbd"
    ),
    labels = c(
      clinical_biomarker_trial_candidate = "clinical biomarker",
      preclinical_mechanistic_perturbation_candidate = "preclinical mechanism",
      target_discovery_or_stratification_candidate = "target discovery",
      exploratory_preclinical_candidate = "exploratory",
      background_or_marker_only = "marker/background"
    ),
    drop = FALSE
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "Perturbation-Prediction Priority",
    subtitle = wrap_text(
      "Hypothesis score combines therapeutic ranking, gene-level quantitative confidence, number of directional predictions, and expected half-normalization effect size."
    ),
    x = "prediction priority score",
    y = NULL,
    fill = "review class"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(figure_dir, "therapeutic_perturbation_hypothesis_summary.png"),
  p_summary,
  width = 10.5,
  height = 5.8,
  dpi = 180
)
if (exists("dominant_save_ggplotly") &&
    requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    p_summary,
    file.path(figure_dir, "therapeutic_perturbation_hypothesis_summary.html"),
    title = "RDG therapeutic perturbation hypothesis summary"
  )
}

message("Saved therapeutic perturbation predictions to: ", output_dir)
