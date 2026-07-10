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
output_dir <- file.path(analysis_dir, "dominant_rdg_therapeutic_hypotheses")
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

safe_max <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

safe_mean <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) default else mean(x)
}

collapse_unique <- function(x, n = 5L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

collapse_top_genes <- function(gene, score, n = 6L) {
  keep <- !is.na(gene) & nzchar(as.character(gene))
  if (!any(keep)) return("")
  dt <- data.table(
    gene_symbol = as.character(gene[keep]),
    score = suppressWarnings(as.numeric(score[keep]))
  )
  dt[!is.finite(score), score := 0]
  dt <- dt[order(-score, gene_symbol)]
  dt <- dt[!duplicated(gene_symbol)]
  paste(
    head(paste0(dt$gene_symbol, "=", sprintf("%.2f", dt$score)), n),
    collapse = "; "
  )
}

wrap_text <- function(x, width = 88L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

top_n_mean <- function(x, n = 5L) {
  x <- suppressWarnings(as.numeric(x))
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  if (!length(x)) return(0)
  mean(head(x, n))
}

split_genes <- function(x) {
  out <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  out[nzchar(out)]
}

hypothesis_definitions <- data.table(
  therapeutic_hypothesis = c(
    "JAK_STAT_IFN_attenuation",
    "ISR_eIF2B_ATF4_attenuation",
    "IFIH1_MDA5_axis_modulation",
    "HIF_hypoxia_angiogenic_modulation",
    "mTOR_translation_capacity_modulation",
    "mitochondrial_OXPHOS_support",
    "immune_antigen_complement_axis"
  ),
  hypothesis_label = c(
    "JAK/STAT IFN attenuation",
    "ISR/eIF2B-ATF4 attenuation",
    "IFIH1/MDA5 axis modulation",
    "HIF/hypoxia angiogenic modulation",
    "mTOR/translation-capacity modulation",
    "Mitochondrial/OXPHOS support",
    "Antigen/complement inflammation axis"
  ),
  intervention_stage = c(
    "clinical_adjacent_biomarker_hypothesis",
    "mechanistic_preclinical_hypothesis",
    "target_discovery_not_trial_ready",
    "exploratory_context_specific",
    "exploratory_context_specific",
    "exploratory_context_specific",
    "target_discovery_not_trial_ready"
  ),
  example_interventions = c(
    "baricitinib; ruxolitinib; tofacitinib",
    "eIF2B activators/ISRIB-like compounds; DNL343-like; fosigotifator-like",
    "direct MDA5/MAVS modulation; downstream IFN dampening",
    "HIF/VEGF pathway modulation",
    "mTORC1 or translation-initiation modulation",
    "mitochondrial support or respiration-state perturbation",
    "immune-axis stratification; complement/antigen-presentation follow-up"
  ),
  primary_genes = c(
    "IFIH1; STAT1; IRF7; IFIT1; IFIT3; ISG15; MX1; OAS1; OAS2; CCL2; TNFAIP3",
    "ATF4; DDIT3; PPP1R15A; ATF3; DDIT4; TRIB3; ASNS; HSPA5; HERPUD1; CEBPB",
    "IFIH1; STAT1; IRF7; IFIT1; IFIT3; ISG15; MX1; OAS1; OAS2",
    "VEGFA; EGLN3; NDRG1; SLC2A1; SLC2A3; LDHA; PDK1; HK2; ENO1",
    "RPS6; RPLP0; RPL32; RPL13A; RPS3; EIF3E; EIF4B; EEF2; PABPC1; MYC; NPM1; NCL; FBL",
    "NDUFA1; NDUFS1; NDUFV1; SDHB; UQCRC1; COX5A; COX4I1; COX8A; ATP5F1A; ATP5MC1; ATP5MC3; TOMM20; TFAM; MRPL12; MRPS12",
    "HLA-DRB1; CFB; C1QB; C1QC; MS4A1; VCAM1; VWF; CCL2; TNFAIP3; B2M"
  ),
  desired_biological_direction = c(
    "Lower persistent antiviral/IFN pressure after acute infection has cleared.",
    "Lower chronic ATF4/DDIT3 ISR output or restore productive CDS buffering.",
    "Stratify or dampen persistent dsRNA-sensing/MDA5 signaling with high safety caution.",
    "Separate true hypoxic angiogenic adaptation from stress-linked translation rewiring.",
    "Test whether global translation capacity is a cause or passenger of RDG shifts.",
    "Test whether respiration state modifies CDS buffering and postviral phenotype output.",
    "Stratify immune-cell or inflammatory contamination versus real RDG reprogramming."
  ),
  clinical_readiness_score = c(1.00, 0.55, 0.25, 0.35, 0.35, 0.25, 0.30),
  mechanistic_alignment_score = c(0.75, 1.00, 0.85, 0.50, 0.35, 0.35, 0.55),
  risk_penalty = c(0.20, 0.25, 0.35, 0.25, 0.20, 0.15, 0.25),
  preferred_context = c(
    "post-viral / Long-COVID-like inflammatory samples",
    "ISR, ER-stress, arsenite, viral-infection, and postviral fatigue contexts",
    "lung-biased viral-infection contexts and IFIH1-high samples",
    "hypoxia, cancer, and viral samples with HIF-like state shifts",
    "high translation-capacity or MYC/mTOR-linked contexts",
    "OXPHOS-enriched contexts and fatigue-like metabolic signatures",
    "immune-enriched viral/fatigue samples and tissue/cell-mixture audits"
  ),
  caveat = c(
    "Biomarker/trial-adjacent only; immune suppression is unsafe during uncontrolled acute infection.",
    "Preclinical perturbation candidate; current evidence is RDG/mechanistic, not patient outcome evidence.",
    "Strong stratification target but not a direct treatment claim; antiviral defense safety risk is high.",
    "Exploratory; pathway is broad and disease-context dependent.",
    "Exploratory; broad translation perturbation is likely nonspecific.",
    "Exploratory; support-like interventions must be separated from causal disease mechanisms.",
    "Likely mixed biology and cell-composition signal; needs metadata and browser validation."
  )
)

atlas_file <- file.path(analysis_dir, "dominant_rdg_atlas",
                        "dominant_rdg_atlas_cards.csv")
atlas <- read_dt(atlas_file, required = TRUE)
if (!nrow(atlas)) stop("Atlas card table is empty: ", atlas_file,
                       call. = FALSE)
atlas[, gene_symbol := toupper(as.character(gene_symbol))]
atlas[, design_family := as.character(design_family)]

warnings_chr <- chr_col(atlas, "warnings")
model_class <- paste(
  chr_col(atlas, "model_agreement_class"),
  chr_col(atlas, "model_agreement_tier"),
  chr_col(atlas, "model_support_signature")
)
transport_class <- paste(
  chr_col(atlas, "transport_consensus_class"),
  chr_col(atlas, "transport_agreement_class")
)
loo_class <- paste(
  chr_col(atlas, "loo_consensus_class"),
  chr_col(atlas, "loo_agreement_class")
)
context_class <- paste(
  chr_col(atlas, "context_interaction_gene_class"),
  chr_col(atlas, "context_interaction_top_class")
)

atlas[, `:=`(
  is_viral_design = grepl("viral|fatigue", design_family, ignore.case = TRUE),
  replicated_consensus = grepl("replicated|multi_model|consensus",
                               model_class, ignore.case = TRUE),
  transport_robust = grepl("support_robust|direction_robust",
                           transport_class, ignore.case = TRUE),
  loo_robust = grepl("support_robust|direction_robust",
                     loo_class, ignore.case = TRUE),
  context_interaction_supported = grepl(
    "direction_reversal|specific_supported|magnitude_modulated",
    context_class,
    ignore.case = TRUE
  ),
  browser_supported = num_col(atlas, "browser_validation_score") > 0 |
    grepl("browser_supported", warnings_chr, ignore.case = TRUE),
  branch_supported = grepl(
    "branch_allocation_supported|joint_branch_allocation_supported|dirichlet_multinomial_supported|hierarchical_joint_allocation_supported",
    warnings_chr,
    ignore.case = TRUE
  ),
  interval_supported = bool_col(atlas, "atlas_interval_supported"),
  direction_stable = bool_col(atlas, "atlas_direction_stable")
)]

atlas[, atlas_priority_component := clip01(
  num_col(atlas, "atlas_browser_review_priority_score") / 1.5
)]
atlas[, sample_component := clip01(
  log10(num_col(atlas, "total_case_samples") +
          num_col(atlas, "total_control_samples") + 1) / 3
)]
atlas[, gene_design_evidence_score := clip01(
  0.10 * clip01(num_col(atlas, "atlas_confidence_score")) +
    0.08 * atlas_priority_component +
    0.12 * clip01(num_col(atlas, "model_agreement_component")) +
    0.07 * clip01(num_col(atlas, "branch_usage_component")) +
    0.10 * clip01(num_col(atlas, "joint_branch_allocation_component")) +
    0.10 * clip01(num_col(atlas, "hierarchical_joint_allocation_component")) +
    0.06 * clip01(num_col(atlas, "dirichlet_multinomial_component")) +
    0.08 * clip01(num_col(atlas, "transport_stability_component")) +
    0.08 * clip01(num_col(atlas, "context_interaction_component")) +
    0.06 * clip01(num_col(atlas, "browser_validation_component")) +
    0.08 * clip01(num_col(atlas, "postviral_atlas_condition_evidence_score")) +
    0.07 * sample_component
)]

hypothesis_genes <- hypothesis_definitions[
  ,
  .(gene_symbol = toupper(split_genes(primary_genes))),
  by = .(therapeutic_hypothesis, hypothesis_label)
]

atlas_keep_columns <- intersect(
  c(
    "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
    "atlas_confidence_score", "atlas_browser_review_priority_score",
    "atlas_direction", "atlas_clean_cds_pct_change",
    "n_designs", "n_studies", "strict_supported_studies",
    "total_case_samples", "total_control_samples",
    "best_case_condition", "best_control_conditions", "best_study",
    "top_shifted_feature_id", "top_shifted_feature_type",
    "model_agreement_tier", "model_agreement_component",
    "model_agreement_review_priority", "transport_consensus_class",
    "transport_stability_component", "loo_consensus_class",
    "loo_stability_component", "context_interaction_gene_class",
    "context_interaction_component", "context_interaction_top_axis",
    "context_interaction_top_value", "context_interaction_top_branch_class",
    "context_interaction_delta", "context_interaction_delta_q",
    "context_interaction_summary", "browser_validation_score",
    "browser_validation_label", "warnings", "is_viral_design",
    "replicated_consensus", "transport_robust", "loo_robust",
    "context_interaction_supported", "browser_supported",
    "branch_supported", "interval_supported", "direction_stable",
    "gene_design_evidence_score"
  ),
  names(atlas)
)

gene_design_evidence <- merge(
  hypothesis_genes,
  atlas[, ..atlas_keep_columns],
  by = "gene_symbol",
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)
gene_design_evidence[, missing_atlas_gene := is.na(design_family)]
for (column in c(
  "is_viral_design", "replicated_consensus", "transport_robust", "loo_robust",
  "context_interaction_supported", "browser_supported", "branch_supported",
  "interval_supported", "direction_stable"
)) {
  if (column %in% names(gene_design_evidence)) {
    idx <- is.na(gene_design_evidence[[column]])
    if (any(idx)) gene_design_evidence[idx, (column) := FALSE]
  }
}
for (column in c(
  "gene_design_evidence_score", "atlas_browser_review_priority_score",
  "atlas_confidence_score", "atlas_clean_cds_pct_change",
  "model_agreement_component", "transport_stability_component",
  "loo_stability_component", "context_interaction_component",
  "browser_validation_score", "total_case_samples", "total_control_samples"
)) {
  if (column %in% names(gene_design_evidence)) {
    value <- suppressWarnings(as.numeric(gene_design_evidence[[column]]))
    idx <- !is.finite(value)
    if (any(idx)) gene_design_evidence[idx, (column) := 0]
  }
}

gene_summary <- gene_design_evidence[
  ,
  {
    score <- suppressWarnings(as.numeric(gene_design_evidence_score))
    score[!is.finite(score)] <- 0
    viral_score <- score
    viral_score[!is_viral_design] <- 0
    best_idx <- which.max(score)
    if (!length(best_idx) || !is.finite(score[best_idx])) best_idx <- 1L
    .(
      n_design_rows = sum(!missing_atlas_gene),
      n_viral_design_rows = sum(is_viral_design & !missing_atlas_gene),
      max_gene_evidence_score = safe_max(score),
      max_viral_evidence_score = safe_max(viral_score),
      mean_top_gene_design_score = top_n_mean(score, 3L),
      replicated_designs = sum(replicated_consensus & !missing_atlas_gene),
      viral_replicated_designs = sum(
        replicated_consensus & is_viral_design & !missing_atlas_gene
      ),
      transport_robust_designs = sum(transport_robust & !missing_atlas_gene),
      loo_robust_designs = sum(loo_robust & !missing_atlas_gene),
      context_interaction_designs = sum(
        context_interaction_supported & !missing_atlas_gene
      ),
      viral_context_interaction_designs = sum(
        context_interaction_supported & is_viral_design & !missing_atlas_gene
      ),
      browser_supported_designs = sum(browser_supported & !missing_atlas_gene),
      branch_supported_designs = sum(branch_supported & !missing_atlas_gene),
      interval_supported_designs = sum(interval_supported & !missing_atlas_gene),
      direction_stable_designs = sum(direction_stable & !missing_atlas_gene),
      top_design_family = if ("design_family" %in% names(.SD)) {
        as.character(design_family[best_idx])
      } else "",
      top_case_condition = if ("best_case_condition" %in% names(.SD)) {
        as.character(best_case_condition[best_idx])
      } else "",
      top_study = if ("best_study" %in% names(.SD)) {
        as.character(best_study[best_idx])
      } else "",
      top_feature = if ("top_shifted_feature_id" %in% names(.SD)) {
        as.character(top_shifted_feature_id[best_idx])
      } else "",
      top_feature_type = if ("top_shifted_feature_type" %in% names(.SD)) {
        as.character(top_shifted_feature_type[best_idx])
      } else "",
      top_clean_cds_pct_change = if (
        "atlas_clean_cds_pct_change" %in% names(.SD)
      ) {
        suppressWarnings(as.numeric(atlas_clean_cds_pct_change[best_idx]))
      } else 0,
      top_context_summary = collapse_unique(context_interaction_summary, 2L),
      top_warnings = collapse_unique(warnings, 2L)
    )
  },
  by = .(therapeutic_hypothesis, hypothesis_label, gene_symbol)
]
gene_summary[is.na(top_design_family), top_design_family := ""]
gene_summary[is.na(top_case_condition), top_case_condition := ""]
gene_summary[is.na(top_study), top_study := ""]
gene_summary[is.na(top_feature), top_feature := ""]
gene_summary[is.na(top_feature_type), top_feature_type := ""]
gene_summary[!is.finite(top_clean_cds_pct_change),
             top_clean_cds_pct_change := 0]

score_summary <- gene_summary[
  ,
  .(
    n_primary_genes = .N,
    n_genes_in_atlas = sum(n_design_rows > 0),
    n_genes_with_any_evidence = sum(max_gene_evidence_score >= 0.35),
    n_viral_genes = sum(max_viral_evidence_score >= 0.35),
    n_consensus_genes = sum(replicated_designs > 0),
    n_viral_consensus_genes = sum(viral_replicated_designs > 0),
    n_transport_robust_genes = sum(transport_robust_designs > 0),
    n_loo_robust_genes = sum(loo_robust_designs > 0),
    n_context_interaction_genes = sum(context_interaction_designs > 0),
    n_viral_context_interaction_genes = sum(
      viral_context_interaction_designs > 0
    ),
    n_browser_supported_genes = sum(browser_supported_designs > 0),
    n_branch_supported_genes = sum(branch_supported_designs > 0),
    max_gene_evidence_score = safe_max(max_gene_evidence_score),
    mean_top5_gene_evidence_score = top_n_mean(max_gene_evidence_score, 5L),
    mean_top5_viral_evidence_score = top_n_mean(
      max_viral_evidence_score, 5L
    ),
    top_evidence_genes = collapse_top_genes(
      gene_symbol, max_gene_evidence_score, 7L
    ),
    top_viral_genes = collapse_top_genes(
      gene_symbol, max_viral_evidence_score, 7L
    ),
    top_contexts = collapse_unique(top_context_summary, 4L),
    representative_studies = collapse_unique(top_study, 5L)
  ),
  by = .(therapeutic_hypothesis, hypothesis_label)
]

scores <- merge(
  hypothesis_definitions,
  score_summary,
  by = c("therapeutic_hypothesis", "hypothesis_label"),
  all.x = TRUE,
  sort = FALSE
)
for (column in c(
  "n_primary_genes", "n_genes_in_atlas", "n_genes_with_any_evidence",
  "n_viral_genes", "n_consensus_genes", "n_viral_consensus_genes",
  "n_transport_robust_genes", "n_loo_robust_genes",
  "n_context_interaction_genes", "n_viral_context_interaction_genes",
  "n_browser_supported_genes", "n_branch_supported_genes",
  "max_gene_evidence_score", "mean_top5_gene_evidence_score",
  "mean_top5_viral_evidence_score"
)) {
  idx <- is.na(scores[[column]])
  if (any(idx)) scores[idx, (column) := 0]
}
for (column in c("top_evidence_genes", "top_viral_genes", "top_contexts",
                 "representative_studies")) {
  idx <- is.na(scores[[column]])
  if (any(idx)) scores[idx, (column) := ""]
}

scores[, gene_coverage_component := clip01(
  n_genes_in_atlas / pmax(3, pmin(n_primary_genes, 8))
)]
scores[, gene_strength_component := clip01(mean_top5_gene_evidence_score)]
scores[, consensus_component := clip01(
  n_consensus_genes / pmax(1, pmin(3, n_genes_in_atlas))
)]
scores[, robust_component := clip01(
  (n_transport_robust_genes + n_loo_robust_genes) /
    pmax(1, 2 * pmin(3, n_genes_in_atlas))
)]
scores[, browser_component := clip01(
  n_browser_supported_genes / pmax(1, pmin(2, n_genes_in_atlas))
)]
scores[, context_component := clip01(
  n_context_interaction_genes / pmax(1, pmin(2, n_genes_in_atlas))
)]
scores[, viral_gene_component := clip01(
  n_viral_genes / pmax(1, pmin(4, n_genes_in_atlas))
)]
scores[, viral_consensus_component := clip01(
  n_viral_consensus_genes / pmax(1, pmin(2, n_viral_genes))
)]
scores[, viral_context_component := clip01(
  n_viral_context_interaction_genes / pmax(1, pmin(2, n_viral_genes))
)]
scores[, transport_context_score := clip01(
  0.45 * robust_component +
    0.35 * context_component +
    0.20 * browser_component
)]
scores[, data_evidence_score := clip01(
  0.18 * gene_coverage_component +
    0.27 * gene_strength_component +
    0.20 * consensus_component +
    0.15 * robust_component +
    0.10 * context_component +
    0.10 * browser_component
)]
scores[, viral_specificity_score := clip01(
  0.40 * clip01(mean_top5_viral_evidence_score) +
    0.25 * viral_gene_component +
    0.20 * viral_consensus_component +
    0.15 * viral_context_component
)]
scores[, therapeutic_hypothesis_score := clip01(
  0.30 * data_evidence_score +
    0.25 * viral_specificity_score +
    0.18 * clinical_readiness_score +
    0.17 * mechanistic_alignment_score +
    0.10 * transport_context_score -
    0.15 * risk_penalty
)]
scores[, trial_readiness_class := fcase(
  intervention_stage == "target_discovery_not_trial_ready" &
    therapeutic_hypothesis_score >= 0.45,
  "target_discovery_or_stratification_candidate",
  therapeutic_hypothesis_score >= 0.62 &
    clinical_readiness_score >= 0.75,
  "clinical_biomarker_trial_candidate",
  therapeutic_hypothesis_score >= 0.58 &
    mechanistic_alignment_score >= 0.75,
  "preclinical_mechanistic_perturbation_candidate",
  therapeutic_hypothesis_score >= 0.48 &
    viral_specificity_score >= 0.60,
  "target_discovery_or_stratification_candidate",
  therapeutic_hypothesis_score >= 0.40,
  "exploratory_preclinical_candidate",
  default = "background_or_marker_only"
)]
scores[, therapeutic_review_rank := frank(
  -therapeutic_hypothesis_score,
  ties.method = "first"
)]
setorder(scores, therapeutic_review_rank)

summary_metrics <- data.table(
  metric = c(
    "n_hypotheses",
    "top_hypothesis",
    "top_hypothesis_score",
    "top_trial_readiness_class",
    "n_clinical_biomarker_trial_candidates",
    "n_preclinical_mechanistic_candidates",
    "n_target_discovery_candidates"
  ),
  value = c(
    as.character(nrow(scores)),
    scores$hypothesis_label[[1]],
    sprintf("%.3f", scores$therapeutic_hypothesis_score[[1]]),
    scores$trial_readiness_class[[1]],
    as.character(sum(
      scores$trial_readiness_class == "clinical_biomarker_trial_candidate"
    )),
    as.character(sum(
      scores$trial_readiness_class ==
        "preclinical_mechanistic_perturbation_candidate"
    )),
    as.character(sum(
      scores$trial_readiness_class ==
        "target_discovery_or_stratification_candidate"
    ))
  )
)

gene_evidence_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_hypothesis_gene_evidence.csv"
)
score_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_hypothesis_scores.csv"
)
review_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_hypothesis_review.csv"
)
summary_file <- file.path(
  output_dir,
  "dominant_rdg_therapeutic_hypothesis_summary_metrics.csv"
)
fwrite(gene_summary, gene_evidence_file)
fwrite(scores, score_file)
fwrite(
  scores[, .(
    therapeutic_review_rank,
    hypothesis_label,
    therapeutic_hypothesis_score,
    trial_readiness_class,
    intervention_stage,
    example_interventions,
    desired_biological_direction,
    n_genes_in_atlas,
    n_viral_genes,
    n_consensus_genes,
    n_context_interaction_genes,
    n_browser_supported_genes,
    top_evidence_genes,
    top_viral_genes,
    top_contexts,
    caveat
  )],
  review_file
)
fwrite(summary_metrics, summary_file)

score_plot <- copy(scores)
score_plot[, hypothesis_label := factor(
  hypothesis_label,
  levels = rev(scores$hypothesis_label)
)]
score_colors <- c(
  clinical_biomarker_trial_candidate = "#1b9e77",
  preclinical_mechanistic_perturbation_candidate = "#d95f02",
  target_discovery_or_stratification_candidate = "#7570b3",
  exploratory_preclinical_candidate = "#4d4d4d",
  background_or_marker_only = "#bdbdbd"
)
score_labels <- c(
  clinical_biomarker_trial_candidate = "clinical biomarker",
  preclinical_mechanistic_perturbation_candidate = "preclinical mechanism",
  target_discovery_or_stratification_candidate = "target discovery",
  exploratory_preclinical_candidate = "exploratory",
  background_or_marker_only = "marker/background"
)
p_scores <- ggplot(
  score_plot,
  aes(
    x = therapeutic_hypothesis_score,
    y = hypothesis_label,
    fill = trial_readiness_class
  )
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = sprintf("%.2f", therapeutic_hypothesis_score)),
    hjust = -0.12,
    size = 3.2
  ) +
  scale_fill_manual(values = score_colors, labels = score_labels,
                    drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.07))) +
  labs(
    title = "RDG Atlas Therapeutic-Hypothesis Prioritization",
    subtitle = wrap_text(paste(
      "Scores are preclinical/biomarker hypotheses from atlas evidence,",
      "viral specificity, robustness, readiness, and safety penalty."
    )),
    x = "hypothesis score",
    y = NULL,
    fill = "review class"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9.5),
    plot.margin = margin(8, 18, 8, 8)
  )
ggsave(
  file.path(figure_dir, "therapeutic_hypothesis_scores.png"),
  p_scores,
  width = 10.5,
  height = 5.6,
  dpi = 180
)
if (exists("dominant_save_ggplotly") &&
    requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    p_scores,
    file.path(figure_dir, "therapeutic_hypothesis_scores.html"),
    title = "RDG therapeutic hypothesis scores"
  )
}

matrix_cols <- c(
  data_evidence_score = "data evidence",
  viral_specificity_score = "viral specificity",
  clinical_readiness_score = "clinical readiness",
  mechanistic_alignment_score = "mechanistic alignment",
  transport_context_score = "transport/context",
  risk_penalty = "risk penalty",
  therapeutic_hypothesis_score = "final score"
)
matrix_dt <- melt(
  scores[, c("hypothesis_label", names(matrix_cols)), with = FALSE],
  id.vars = "hypothesis_label",
  variable.name = "score_component",
  value.name = "score"
)
matrix_dt[, score_component := factor(
  as.character(matrix_cols[as.character(score_component)]),
  levels = matrix_cols
)]
matrix_dt[, hypothesis_label := factor(
  hypothesis_label,
  levels = rev(scores$hypothesis_label)
)]
p_matrix <- ggplot(
  matrix_dt,
  aes(x = score_component, y = hypothesis_label, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.7) +
  scale_fill_gradientn(
    colors = c("#2c7bb6", "white", "#8b0000"),
    values = c(0, 0.5, 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Therapeutic-Hypothesis Evidence Matrix",
    subtitle = "Risk penalty is shown as a component; it is subtracted in the final score.",
    x = NULL,
    y = NULL,
    fill = "score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(figure_dir, "therapeutic_hypothesis_evidence_matrix.png"),
  p_matrix,
  width = 9.5,
  height = 5.5,
  dpi = 180
)
if (exists("dominant_save_ggplotly") &&
    requireNamespace("plotly", quietly = TRUE)) {
  dominant_save_ggplotly(
    p_matrix,
    file.path(figure_dir, "therapeutic_hypothesis_evidence_matrix.html"),
    title = "RDG therapeutic hypothesis evidence matrix"
  )
}

message("Saved therapeutic-hypothesis outputs to: ", output_dir)
