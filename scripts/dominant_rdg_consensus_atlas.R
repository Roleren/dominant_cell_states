#!/usr/bin/env Rscript

# Consolidate the RDG atlas evidence layers into a compact consensus queue.
# This is not a new statistical model. It is a deterministic review layer that
# makes agreement, fragility, context dependence, and therapeutic relevance
# visible in one table.

dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
output_dir <- file.path(analysis_dir, "dominant_rdg_consensus")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

atlas_file <- file.path(analysis_dir, "dominant_rdg_atlas",
                        "dominant_rdg_atlas_cards.csv")
aux_file <- file.path(analysis_dir, "dominant_rdg_auxiliary_state_audit",
                      "dominant_rdg_auxiliary_state_atlas_annotations.csv")
therapeutic_gene_file <- file.path(
  analysis_dir,
  "dominant_rdg_therapeutic_hypotheses",
  "dominant_rdg_therapeutic_hypothesis_gene_evidence.csv"
)
perturb_file <- file.path(
  analysis_dir,
  "dominant_rdg_therapeutic_perturbation_predictions",
  "dominant_rdg_therapeutic_perturbation_gene_predictions.csv"
)

candidate_file <- file.path(output_dir, "dominant_rdg_consensus_candidates.csv")
review_file <- file.path(output_dir, "dominant_rdg_consensus_review_queue.csv")
evidence_file <- file.path(output_dir, "dominant_rdg_consensus_evidence_matrix.csv")
gene_file <- file.path(output_dir, "dominant_rdg_consensus_gene_summary.csv")
metrics_file <- file.path(output_dir, "dominant_rdg_consensus_summary_metrics.csv")

if (!file.exists(atlas_file)) {
  stop("Missing atlas cards file: ", atlas_file, call. = FALSE)
}

read_optional <- function(path) {
  if (!file.exists(path)) return(data.table())
  fread(path, showProgress = FALSE)
}

bounded01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- 0
  pmin(pmax(x, 0), 1)
}

scale_cap <- function(x, cap) {
  bounded01(suppressWarnings(as.numeric(x)) / cap)
}

num_col <- function(dt, col, default = 0) {
  if (!col %chin% names(dt)) return(rep(default, nrow(dt)))
  x <- suppressWarnings(as.numeric(dt[[col]]))
  x[!is.finite(x)] <- default
  x
}

chr_col <- function(dt, col, default = "") {
  if (!col %chin% names(dt)) return(rep(default, nrow(dt)))
  x <- as.character(dt[[col]])
  x[is.na(x)] <- default
  x
}

logical_col <- function(dt, col, default = FALSE) {
  if (!col %chin% names(dt)) return(rep(default, nrow(dt)))
  x <- dt[[col]]
  if (is.logical(x)) {
    x[is.na(x)] <- default
    return(x)
  }
  x <- tolower(as.character(x))
  fifelse(x %chin% c("true", "t", "1", "yes"), TRUE,
          fifelse(x %chin% c("false", "f", "0", "no"), FALSE, default))
}

pmax0 <- function(...) {
  values <- list(...)
  if (length(values) == 0) return(numeric())
  mat <- do.call(cbind, lapply(values, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[!is.finite(x)] <- 0
    x
  }))
  apply(mat, 1, max, na.rm = TRUE)
}

tag_string <- function(...) {
  tags <- list(...)
  n <- length(tags[[1]])
  out <- character(n)
  for (i in seq_len(n)) {
    present <- names(tags)[vapply(tags, function(x) isTRUE(x[[i]]), logical(1))]
    out[[i]] <- paste(present, collapse = "; ")
  }
  out
}

first_nonempty <- function(...) {
  values <- list(...)
  n <- length(values[[1]])
  out <- rep("", n)
  for (x in values) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    take <- !nzchar(out) & nzchar(x)
    out[take] <- x[take]
  }
  out
}

save_plot <- function(plot, name, width, height) {
  png_file <- file.path(figure_dir, paste0(name, ".png"))
  pdf_file <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220,
         bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  if (exists("dominant_save_ggplotly", mode = "function") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      plot,
      file.path(figure_dir, paste0(name, ".html")),
      title = name
    )
  }
}

atlas <- fread(atlas_file, showProgress = FALSE)
required_atlas <- c("gene_symbol", "tx_id", "design_family")
if (!all(required_atlas %chin% names(atlas))) {
  stop("Atlas cards are missing columns: ",
       paste(setdiff(required_atlas, names(atlas)), collapse = ", "),
       call. = FALSE)
}

keys <- c("gene_symbol", "tx_id", "design_family")

aux <- read_optional(aux_file)
if (nrow(aux)) {
  keep <- intersect(
    c(keys, "aux_top_state", "aux_top_overlap_class",
      "aux_top_abs_delta_z", "aux_top_overlap_score",
      "aux_top_weighted_overlap_score", "aux_top_readiness_weight",
      "aux_top_context_class", "aux_top_context",
      "aux_top_recommendation_label", "aux_audit_warning",
      "aux_audit_priority"),
    names(aux)
  )
  aux <- unique(aux[, ..keep], by = keys)
  atlas <- merge(atlas, aux, by = keys, all.x = TRUE, sort = FALSE)
}

ther_gene <- read_optional(therapeutic_gene_file)
if (nrow(ther_gene) && all(c("gene_symbol", "therapeutic_hypothesis") %chin%
                           names(ther_gene))) {
  ther_gene[, therapeutic_gene_score := pmax0(
    num_col(ther_gene, "max_gene_evidence_score"),
    num_col(ther_gene, "max_viral_evidence_score"),
    num_col(ther_gene, "mean_top_gene_design_score")
  )]
  setorder(ther_gene, gene_symbol, -therapeutic_gene_score)
  ther_gene <- ther_gene[
    ,
    .SD[1],
    by = gene_symbol
  ][
    ,
    .(
      gene_symbol,
      therapeutic_gene_hypothesis = therapeutic_hypothesis,
      therapeutic_gene_label = chr_col(.SD, "hypothesis_label"),
      therapeutic_gene_score
    )
  ]
  atlas <- merge(atlas, ther_gene, by = "gene_symbol", all.x = TRUE,
                 sort = FALSE)
}

perturb <- read_optional(perturb_file)
if (nrow(perturb) && all(keys %chin% names(perturb))) {
  perturb[, perturbation_component := pmax0(
    num_col(perturb, "row_prediction_evidence_score"),
    num_col(perturb, "quantitative_prediction_confidence"),
    scale_cap(abs(num_col(perturb, "half_normalization_cds_pct_change")), 100)
  )]
  setorder(perturb, gene_symbol, tx_id, design_family,
           -perturbation_component,
           perturbation_prediction_rank)
  perturb <- perturb[
    ,
    .SD[1],
    by = keys
  ][
    ,
    .(
      gene_symbol, tx_id, design_family,
      perturbation_hypothesis = therapeutic_hypothesis,
      perturbation_label = chr_col(.SD, "hypothesis_label"),
      perturbation_component,
      perturbation_support_class = chr_col(.SD, "prediction_support_class"),
      half_normalization_cds_pct_change =
        num_col(.SD, "half_normalization_cds_pct_change", NA_real_),
      abs_half_normalization_cds_pct_change =
        num_col(.SD, "abs_half_normalization_cds_pct_change", NA_real_),
      validation_assay = chr_col(.SD, "validation_assay")
    )
  ]
  atlas <- merge(atlas, perturb, by = keys, all.x = TRUE, sort = FALSE)
}

atlas[, atlas_component := bounded01(num_col(.SD, "atlas_confidence_score"))]
atlas[, output_effect_component :=
        scale_cap(abs(num_col(.SD, "atlas_clean_cds_pct_change")), 120)]
atlas[, browser_component := pmax0(
  num_col(.SD, "browser_validation_component"),
  num_col(.SD, "browser_validation_score")
)]
atlas[, branch_component := pmax0(
  num_col(.SD, "branch_usage_component"),
  num_col(.SD, "grouped_branch_usage_component"),
  num_col(.SD, "count_aware_branch_allocation_component"),
  num_col(.SD, "joint_branch_allocation_component"),
  num_col(.SD, "hierarchical_joint_allocation_component"),
  num_col(.SD, "dirichlet_multinomial_component")
)]
atlas[, model_component := pmax0(
  num_col(.SD, "model_agreement_component"),
  fifelse(grepl("browser_validated|replicated|multi_model",
                chr_col(.SD, "model_agreement_tier"),
                ignore.case = TRUE), 1, 0),
  fifelse(grepl("consensus", chr_col(.SD, "model_agreement_class"),
                ignore.case = TRUE), 0.75, 0)
)]
atlas[, robustness_component := pmax0(
  num_col(.SD, "loo_stability_component"),
  num_col(.SD, "transport_stability_component"),
  fifelse(grepl("support_robust", chr_col(.SD, "loo_consensus_class"),
                ignore.case = TRUE), 1, 0),
  fifelse(grepl("direction_robust", chr_col(.SD, "loo_consensus_class"),
                ignore.case = TRUE), 0.65, 0),
  fifelse(grepl("support_robust", chr_col(.SD, "transport_consensus_class"),
                ignore.case = TRUE), 1, 0),
  fifelse(grepl("direction_robust", chr_col(.SD, "transport_consensus_class"),
                ignore.case = TRUE), 0.65, 0)
)]
atlas[, context_component := pmax0(
  num_col(.SD, "context_interaction_component"),
  num_col(.SD, "partial_pooling_component"),
  fifelse(grepl("partial_pooled_context_supported",
                chr_col(.SD, "partial_pooling_gene_class"),
                ignore.case = TRUE), 1, 0),
  fifelse(grepl("context_direction_reversal|context_specific|magnitude_modulated",
                chr_col(.SD, "context_interaction_gene_class"),
                ignore.case = TRUE), 0.75, 0),
  fifelse(grepl("single_context_raw_signal",
                chr_col(.SD, "partial_pooling_gene_class"),
                ignore.case = TRUE), 0.45, 0)
)]
atlas[, pooling_component :=
        bounded01(num_col(.SD, "pooling_proof_component"))]
atlas[, pooling_hazard_component :=
        bounded01(num_col(.SD, "pooling_hazard_score"))]
atlas[, aux_context_component := bounded01(num_col(.SD, "aux_audit_priority"))]
atlas[, therapeutic_component := pmax0(
  num_col(.SD, "perturbation_component"),
  num_col(.SD, "therapeutic_gene_score")
)]
atlas[, sparse_component2 := bounded01(num_col(.SD, "sparse_onoff_score"))]
atlas[, novelty_component := fifelse(
  !grepl("known", chr_col(.SD, "novelty_class"), ignore.case = TRUE) &
    (num_col(.SD, "hidden_uorf_candidate") > 0 |
       logical_col(.SD, "has_overlapping_uorf") |
       num_col(.SD, "n_overlapping_uorfs") > 0),
  0.5,
  0
)]

atlas[, fragility_component := pmax0(
  num_col(.SD, "loo_fragility_component"),
  num_col(.SD, "transport_fragility_component"),
  num_col(.SD, "model_disagreement_component"),
  fifelse(grepl("fragile|not_evaluable|non_consensus",
                paste(chr_col(.SD, "loo_consensus_class"),
                      chr_col(.SD, "transport_consensus_class"),
                      chr_col(.SD, "partial_pooling_gene_class")),
                ignore.case = TRUE), 0.45, 0),
  fifelse(logical_col(.SD, "high_branch_heterogeneity"), 0.35, 0),
  fifelse(logical_col(.SD, "allocation_output_disagreement"), 0.35, 0),
  0.85 * pooling_hazard_component,
  scale_cap(num_col(.SD, "residual_heterogeneity_score"), 1)
)]
atlas[, warning_penalty := pmin(0.30, 0.16 * fragility_component +
                                 0.05 * aux_context_component +
                                 0.07 * pooling_hazard_component)]

atlas[, consensus_score := pmin(pmax(
  0.15 * atlas_component +
    0.22 * branch_component +
    0.18 * model_component +
    0.15 * robustness_component +
    0.10 * browser_component +
    0.07 * output_effect_component +
    0.05 * context_component +
    0.04 * pooling_component +
    0.04 * therapeutic_component +
    0.02 * sparse_component2 +
    0.02 * novelty_component -
    warning_penalty,
  0), 1)]

browser_supported <- atlas$browser_component > 0
multi_model <- atlas$model_component >= 0.75
branch_supported <- atlas$branch_component >= 0.55
robust_supported <- atlas$robustness_component >= 0.65
context_specific <- atlas$context_component >= 0.65
therapeutic_ready <- atlas$therapeutic_component >= 0.55
exact_project_matched <- chr_col(atlas, "pooling_proof_tier") ==
  "A_exact_strict_project_matched_control"
pooling_confounded <- atlas$pooling_hazard_component >= 0.40
fragile <- atlas$fragility_component >= 0.5 |
  atlas$pooling_hazard_component >= 0.65
buffering_discordant <- logical_col(atlas, "allocation_output_disagreement")

atlas[, consensus_class := fcase(
  browser_supported & multi_model & branch_supported & robust_supported,
  "A_browser_validated_consensus",
  multi_model & branch_supported & robust_supported,
  "B_replicated_robust_consensus",
  multi_model & branch_supported,
  "C_multi_model_branch_consensus",
  context_specific & !robust_supported,
  "D_context_specific_review",
  buffering_discordant,
  "E_buffering_compensation_review",
  therapeutic_ready & (branch_supported | multi_model),
  "F_therapeutic_translation_review",
  atlas_component >= 0.65 | branch_supported | consensus_score >= 0.45,
  "G_general_atlas_review",
  default = "H_low_support_background"
)]
atlas[fragile & consensus_class %chin% c("G_general_atlas_review",
                                         "H_low_support_background"),
      consensus_class := "I_fragile_or_confounded_review"]

atlas[, evidence_tags := tag_string(
  browser_supported = browser_component > 0,
  branch_supported = branch_component >= 0.55,
  multi_model = model_component >= 0.75,
  loo_or_transport_robust = robustness_component >= 0.65,
  context_specific = context_component >= 0.65,
  partial_pooled = grepl("partial_pooled_context_supported",
                         chr_col(.SD, "partial_pooling_gene_class"),
                         ignore.case = TRUE),
  single_context = grepl("single_context_raw_signal",
                         chr_col(.SD, "partial_pooling_gene_class"),
                         ignore.case = TRUE),
  exact_project_matched = chr_col(.SD, "pooling_proof_tier") ==
    "A_exact_strict_project_matched_control",
  pooling_hazard = pooling_hazard_component >= 0.40,
  therapeutic_candidate = therapeutic_component >= 0.55,
  auxiliary_state_warning = aux_context_component >= 0.45,
  large_cds_effect = output_effect_component >= 0.50,
  sparse_onoff = sparse_component2 >= 0.50,
  allocation_output_disagreement =
    logical_col(.SD, "allocation_output_disagreement"),
  fragile = fragility_component >= 0.50
)]

atlas[, consensus_warning_summary := first_nonempty(
  chr_col(.SD, "aux_audit_warning"),
  chr_col(.SD, "warnings")
)]

atlas[, consensus_review_question := fifelse(
  consensus_class == "A_browser_validated_consensus",
  "Use as positive-control/anchor; verify branch and CDS effect direction in matched browser views.",
  fifelse(
    consensus_class == "B_replicated_robust_consensus",
    "High-priority browser validation: replicated branch and robustness layers agree.",
    fifelse(
      consensus_class == "C_multi_model_branch_consensus",
      "Inspect branch redistribution and study/context heterogeneity before biological interpretation.",
      fifelse(
        consensus_class == "D_context_specific_review",
        "Inspect the named tissue/cell-line context against its complement; do not generalize yet.",
        fifelse(
          consensus_class == "E_buffering_compensation_review",
          "Check whether branch allocation and net CDS output disagree because of total-load buffering or alternate routing.",
          fifelse(
            consensus_class == "F_therapeutic_translation_review",
            "Review as perturbation-assay target; confirm the predicted CDS-buffering direction first.",
            "Review only if the gene/design is biologically important or needed as a negative-control context."
          )
        )
      )
    )
  )
)]

atlas[, consensus_rank := frank(-consensus_score, ties.method = "first")]
setorder(atlas, consensus_rank)

output_cols <- intersect(c(
  "consensus_rank", "gene_symbol", "tx_id", "design_family",
  "consensus_class", "consensus_score", "evidence_tags",
  "atlas_evidence_class", "atlas_confidence_score",
  "atlas_browser_review_priority_score", "atlas_clean_cds_pct_change",
  "atlas_clean_cds_pct_lower", "atlas_clean_cds_pct_upper",
  "atlas_direction", "n_studies", "strict_supported_studies",
  "total_case_samples", "total_control_samples",
  "best_case_condition", "best_control_conditions", "best_study",
  "top_shifted_feature_id", "top_shifted_feature_type",
  "top_shifted_feature_family",
  "joint_top_branch_class", "joint_top_branch_delta",
  "hierarchical_joint_top_branch_class", "hierarchical_joint_top_delta",
  "dm_top_branch_class", "dm_top_delta",
  "model_agreement_tier", "model_agreement_class",
  "loo_consensus_class", "transport_consensus_class",
  "pooling_proof_tier", "pooling_claim_language",
  "pooling_proof_component", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources",
  "pooling_top_protocol_classes", "pooling_exact_context_fraction",
  "pooling_joint_strict_fraction",
  "context_interaction_gene_class", "context_interaction_top_axis",
  "context_interaction_top_value", "context_interaction_top_branch_class",
  "context_interaction_summary", "partial_pooling_gene_class",
  "partial_pooling_top_axis", "partial_pooling_top_value",
  "partial_pooling_top_branch_class",
  "browser_validation_label", "browser_validation_context",
  "aux_top_state", "aux_top_overlap_class", "aux_audit_warning",
  "therapeutic_gene_hypothesis", "therapeutic_gene_label",
  "perturbation_hypothesis", "perturbation_label",
  "perturbation_support_class", "half_normalization_cds_pct_change",
  "validation_assay",
  "atlas_component", "branch_component", "model_component",
  "robustness_component", "browser_component", "output_effect_component",
  "context_component", "aux_context_component", "therapeutic_component",
  "pooling_component", "pooling_hazard_component",
  "sparse_component2", "novelty_component", "fragility_component",
  "warning_penalty", "consensus_warning_summary",
  "consensus_review_question", "warnings"
), names(atlas))

candidate <- copy(atlas[, ..output_cols])
score_cols <- intersect(c(
  "consensus_score", "atlas_confidence_score",
  "atlas_browser_review_priority_score", "atlas_clean_cds_pct_change",
  "atlas_clean_cds_pct_lower", "atlas_clean_cds_pct_upper",
  "joint_top_branch_delta", "hierarchical_joint_top_delta", "dm_top_delta",
  "half_normalization_cds_pct_change", "pooling_proof_component",
  "pooling_hazard_score", "pooling_exact_context_fraction",
  "pooling_joint_strict_fraction", "atlas_component",
  "branch_component", "model_component", "robustness_component",
  "browser_component", "output_effect_component", "context_component",
  "aux_context_component", "therapeutic_component",
  "pooling_component", "pooling_hazard_component", "sparse_component2",
  "novelty_component", "fragility_component", "warning_penalty"
), names(candidate))
candidate[, (score_cols) := lapply(.SD, function(x) round(as.numeric(x), 4)),
          .SDcols = score_cols]

review <- candidate[
  consensus_class != "H_low_support_background" |
    consensus_score >= 0.35 |
    grepl("viral", design_family, ignore.case = TRUE)
]
review[, review_rank := seq_len(.N)]
setcolorder(review, c("review_rank", setdiff(names(review), "review_rank")))

component_cols <- c(
  atlas_component = "Atlas clean-CDS",
  branch_component = "Branch redistribution",
  model_component = "Model agreement",
  robustness_component = "LOO/transport robustness",
  browser_component = "Browser support",
  output_effect_component = "CDS effect size",
  context_component = "Context specificity",
  pooling_component = "Pooling proof",
  pooling_hazard_component = "Pooling hazard",
  aux_context_component = "Auxiliary-state warning",
  therapeutic_component = "Therapeutic relevance",
  fragility_component = "Fragility/disagreement",
  warning_penalty = "Penalty"
)
component_cols <- component_cols[names(component_cols) %chin% names(atlas)]
evidence_matrix <- melt(
  atlas[
    ,
    c(keys, "consensus_rank", "consensus_class", "consensus_score",
      names(component_cols)),
    with = FALSE
  ],
  id.vars = c(keys, "consensus_rank", "consensus_class", "consensus_score"),
  variable.name = "component",
  value.name = "component_score",
  variable.factor = FALSE
)
evidence_matrix[, component_label := unname(component_cols[component])]
evidence_matrix[, component_score := round(as.numeric(component_score), 4)]
evidence_matrix[, component_order := match(component, names(component_cols))]
setorder(evidence_matrix, consensus_rank, component_order)
evidence_matrix[, component_order := NULL]

gene_summary <- atlas[
  ,
  .(
    best_consensus_rank = min(consensus_rank),
    best_consensus_score = max(consensus_score, na.rm = TRUE),
    best_consensus_class = consensus_class[which.max(consensus_score)],
    n_design_families = uniqueN(design_family),
    n_consensus_rows = sum(consensus_score >= 0.50, na.rm = TRUE),
    n_browser_supported_rows = sum(browser_component > 0, na.rm = TRUE),
    n_viral_rows = sum(grepl("viral", design_family, ignore.case = TRUE)),
    best_viral_score = suppressWarnings(max(
      fifelse(grepl("viral", design_family, ignore.case = TRUE),
              consensus_score, NA_real_),
      na.rm = TRUE
    )),
    best_design_family = design_family[which.max(consensus_score)],
    best_review_question = consensus_review_question[which.max(consensus_score)]
  ),
  by = gene_symbol
]
gene_summary[!is.finite(best_viral_score), best_viral_score := NA_real_]
gene_summary[, best_consensus_score := round(best_consensus_score, 4)]
gene_summary[, best_viral_score := round(best_viral_score, 4)]
setorder(gene_summary, best_consensus_rank)

metrics <- rbindlist(list(
  data.table(metric = "atlas_rows", value = nrow(atlas)),
  data.table(metric = "review_queue_rows", value = nrow(review)),
  data.table(metric = "genes", value = uniqueN(atlas$gene_symbol)),
  data.table(metric = "consensus_score_ge_0_50",
             value = sum(atlas$consensus_score >= 0.50, na.rm = TRUE)),
  data.table(metric = "consensus_score_ge_0_70",
             value = sum(atlas$consensus_score >= 0.70, na.rm = TRUE)),
  data.table(metric = "viral_review_rows",
             value = sum(grepl("viral", review$design_family,
                               ignore.case = TRUE))),
  data.table(metric = "browser_supported_rows",
             value = sum(atlas$browser_component > 0, na.rm = TRUE)),
  data.table(metric = "class_A_browser_validated_consensus",
             value = sum(atlas$consensus_class ==
                           "A_browser_validated_consensus")),
  data.table(metric = "class_B_replicated_robust_consensus",
             value = sum(atlas$consensus_class ==
                           "B_replicated_robust_consensus")),
  data.table(metric = "class_C_multi_model_branch_consensus",
             value = sum(atlas$consensus_class ==
                           "C_multi_model_branch_consensus")),
  data.table(metric = "class_D_context_specific_review",
             value = sum(atlas$consensus_class ==
                           "D_context_specific_review")),
  data.table(metric = "class_E_buffering_compensation_review",
             value = sum(atlas$consensus_class ==
                           "E_buffering_compensation_review")),
  data.table(metric = "class_F_therapeutic_translation_review",
             value = sum(atlas$consensus_class ==
                           "F_therapeutic_translation_review")),
  data.table(metric = "class_G_general_atlas_review",
             value = sum(atlas$consensus_class ==
                           "G_general_atlas_review")),
  data.table(metric = "class_H_low_support_background",
             value = sum(atlas$consensus_class ==
                           "H_low_support_background")),
  data.table(metric = "class_I_fragile_or_confounded_review",
             value = sum(atlas$consensus_class ==
                           "I_fragile_or_confounded_review")),
  data.table(metric = "pooling_exact_project_matched_rows",
             value = sum(exact_project_matched, na.rm = TRUE)),
  data.table(metric = "pooling_hazard_ge_0_40_rows",
             value = sum(pooling_confounded, na.rm = TRUE)),
  data.table(metric = "pooling_hazard_ge_0_65_rows",
             value = sum(atlas$pooling_hazard_component >= 0.65,
                         na.rm = TRUE))
), fill = TRUE)
metrics[, value := as.numeric(value)]

fwrite(candidate, candidate_file)
fwrite(review, review_file)
fwrite(evidence_matrix, evidence_file)
fwrite(gene_summary, gene_file)
fwrite(metrics, metrics_file)

top_plot <- atlas[consensus_rank <= 35]
top_plot[, plot_label := paste(gene_symbol, design_family, sep = " | ")]
top_plot[, plot_label := factor(plot_label, levels = rev(plot_label))]
p_bar <- ggplot(top_plot, aes(consensus_score, plot_label,
                              fill = consensus_class)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = sprintf("%.2f", consensus_score)),
    hjust = -0.08,
    size = 2.8
  ) +
  scale_x_continuous(limits = c(0, min(1, max(top_plot$consensus_score) + 0.12)),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Top RDG consensus candidates",
    subtitle = "Consensus combines atlas output, branch redistribution, model agreement, robustness, browser support, context evidence, and penalties.",
    x = "Consensus score",
    y = NULL,
    fill = "Class"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  )
save_plot(p_bar, "dominant_rdg_consensus_top_candidates", 11.5, 8.6)

heat_dt <- evidence_matrix[consensus_rank <= 35]
heat_dt[, row_label := paste(gene_symbol, design_family, sep = " | ")]
heat_dt[, row_label := factor(row_label, levels = rev(unique(row_label)))]
heat_dt[, component_label := factor(component_label,
                                    levels = unname(component_cols))]
p_heat <- ggplot(heat_dt, aes(component_label, row_label,
                              fill = component_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2F6C9F",
    mid = "#F7F7F7",
    high = "#9E2F2F",
    midpoint = 0.5,
    limits = c(0, 1),
    oob = scales::squish
  ) +
  labs(
    title = "Consensus evidence components",
    subtitle = "Blue = weak, white = intermediate, red = strong. Penalty and fragility are shown as evidence against over-interpretation.",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_minimal(base_size = 10.2) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title.position = "plot"
  )
save_plot(p_heat, "dominant_rdg_consensus_evidence_heatmap", 11.8, 8.8)

class_counts <- atlas[
  ,
  .(N = .N),
  by = consensus_class
][order(-N)]
class_counts[, consensus_class := factor(consensus_class,
                                         levels = rev(consensus_class))]
p_class <- ggplot(class_counts, aes(N, consensus_class,
                                    fill = consensus_class)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = N), hjust = -0.08, size = 3.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "RDG consensus class distribution",
    x = "Atlas rows",
    y = NULL
  ) +
  theme_minimal(base_size = 10.8) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )
save_plot(p_class, "dominant_rdg_consensus_class_counts", 9.5, 5.6)

stopifnot(nrow(candidate) == nrow(atlas))
stopifnot(candidate$consensus_rank[1] == 1)
stopifnot(all(diff(candidate$consensus_rank) == 1))
stopifnot(!any(is.na(candidate$consensus_score)))
stopifnot(all(candidate$consensus_score >= -1e-12 &
                candidate$consensus_score <= 1 + 1e-12))
stopifnot(!any(is.na(candidate$consensus_class) |
                !nzchar(candidate$consensus_class)))

message("RDG consensus atlas:")
message("  atlas rows: ", nrow(candidate))
message("  review rows: ", nrow(review))
message("  genes: ", nrow(gene_summary))
message("  top candidate: ",
        paste(candidate[1, .(gene_symbol, design_family, consensus_score)],
              collapse = " | "))
message("  outputs: ", output_dir)
