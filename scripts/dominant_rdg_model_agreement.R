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
output_dir <- file.path(analysis_dir, "dominant_rdg_model_agreement")
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

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) paths[1] else hit
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

finite_or_zero <- function(x) {
  x <- safe_num(x)
  x[!is.finite(x)] <- 0
  x
}

clip01 <- function(x) pmin(pmax(finite_or_zero(x), 0), 1)

logical_safe <- function(x) {
  if (is.logical(x)) return(fifelse(is.na(x), FALSE, x))
  tolower(as.character(x)) %chin% c("true", "t", "1", "yes")
}

collapse_unique <- function(x, n = 8L) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (!length(x)) return("")
  paste(head(x, n), collapse = "; ")
}

existing_columns <- function(dt, columns) {
  columns[columns %in% names(dt)]
}

select_columns <- function(dt, columns) {
  dt[, ..columns]
}

rename_columns <- function(dt, mapping) {
  present <- names(mapping)[names(mapping) %in% names(dt)]
  if (length(present)) setnames(dt, present, unname(mapping[present]))
  dt
}

column_or <- function(dt, column, default = NA) {
  if (column %in% names(dt)) dt[[column]] else rep(default, nrow(dt))
}

message("RDG model-agreement analysis:")
message("  1. Compare clean-CDS, joint, hierarchical-joint, and posterior-predictive support")
message("  2. Compare hierarchical-joint and Dirichlet-multinomial branch directions")
message("  3. Create compact consensus and disagreement review tiers")

hierarchical <- read_dt(
  file.path(analysis_dir, "dominant_rdg_hierarchical",
            "dominant_rdg_hierarchical_gene_design_effects.csv"),
  required = TRUE
)
joint <- read_dt(
  file.path(analysis_dir, "dominant_rdg_joint_branch_allocation",
            "dominant_rdg_joint_branch_allocation_gene_design_summary.csv"),
  required = TRUE
)
hierarchical_joint <- read_dt(
  file.path(analysis_dir, "dominant_rdg_hierarchical_joint_allocation",
            "dominant_rdg_hierarchical_joint_allocation_gene_design_summary.csv"),
  required = TRUE
)
dm <- read_dt(
  file.path(analysis_dir, "dominant_rdg_dirichlet_multinomial",
            "dominant_rdg_dirichlet_multinomial_gene_design_summary.csv"),
  required = TRUE
)
motif <- read_dt(
  file.path(analysis_dir, "dominant_rdg_motif_prior_calibration",
            "dominant_rdg_motif_prior_gene_design_calibration.csv")
)
browser <- read_dt(
  first_existing_file(c(
    file.path(analysis_dir, "dominant_rdg_atlas", "browser_validation_notes.csv"),
    file.path(results_dir, "curated_inputs", "browser_validation_notes.csv")
  ))
)

key <- c("gene_symbol", "tx_id", "design_family")
for (input in list(hierarchical, joint, hierarchical_joint, dm)) {
  missing <- setdiff(key, names(input))
  if (length(missing)) {
    stop("Model-agreement input is missing keys: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

clean_keep <- existing_columns(hierarchical, c(
  key,
  "hierarchical_evidence_class",
  "hierarchical_priority_score",
  "hierarchical_direction",
  "hierarchical_clean_cds_pct_change",
  "hierarchical_direction_fraction",
  "hierarchical_direction_stable",
  "hierarchical_ci_excludes_zero",
  "hierarchical_calibrated_ci_excludes_zero",
  "n_studies",
  "n_strict_studies"
))
agreement <- select_columns(hierarchical, clean_keep)
agreement <- rename_columns(agreement, c(
  hierarchical_evidence_class = "clean_cds_evidence_class",
  hierarchical_priority_score = "clean_cds_priority",
  hierarchical_direction = "clean_cds_direction",
  hierarchical_clean_cds_pct_change = "clean_cds_pct_change",
  hierarchical_direction_fraction = "clean_cds_direction_fraction",
  hierarchical_direction_stable = "clean_cds_direction_stable",
  hierarchical_ci_excludes_zero = "clean_cds_ci_excludes_zero",
  hierarchical_calibrated_ci_excludes_zero =
    "clean_cds_calibrated_ci_excludes_zero",
  n_studies = "clean_cds_studies",
  n_strict_studies = "clean_cds_strict_studies"
))
agreement <- unique(agreement, by = key)

joint_keep <- existing_columns(joint, c(
  key,
  "joint_top_review_class",
  "joint_top_review_priority",
  "joint_top_branch_class",
  "joint_top_branch_delta",
  "joint_top_l1_shift",
  "joint_strict_studies",
  "joint_evaluable_studies",
  "joint_weighted_clean_CDS_delta",
  "joint_weighted_leader_uORF_delta",
  "joint_weighted_overlapping_uORF_delta",
  "joint_top_motif_alignment",
  "joint_top_context"
))
agreement <- merge(
  agreement, select_columns(joint, joint_keep),
  by = key, all.x = TRUE, sort = FALSE
)

hierarchical_joint_keep <- existing_columns(hierarchical_joint, c(
  key,
  "hierarchical_joint_top_support_class",
  "hierarchical_joint_top_priority",
  "hierarchical_joint_top_branch_class",
  "hierarchical_joint_top_delta",
  "hierarchical_joint_top_studies",
  "hierarchical_joint_top_direction_fraction",
  "hierarchical_joint_top_i2",
  "hierarchical_joint_top_context",
  paste0("hierarchical_joint_delta_",
         c("clean_CDS", "leader_uORF", "overlapping_uORF")),
  paste0("hierarchical_joint_priority_",
         c("clean_CDS", "leader_uORF", "overlapping_uORF")),
  paste0("hierarchical_joint_support_class_",
         c("clean_CDS", "leader_uORF", "overlapping_uORF"))
))
agreement <- merge(
  agreement, select_columns(hierarchical_joint, hierarchical_joint_keep),
  by = key, all.x = TRUE, sort = FALSE
)

dm_keep <- existing_columns(dm, c(
  key,
  "dm_top_support_class",
  "dm_top_priority",
  "dm_top_branch_class",
  "dm_top_delta",
  "dm_top_studies",
  "dm_top_direction_fraction",
  "dm_top_i2",
  "dm_top_prior_source",
  "dm_top_context",
  paste0("dm_delta_", c("clean_CDS", "leader_uORF", "overlapping_uORF")),
  paste0("dm_priority_", c("clean_CDS", "leader_uORF", "overlapping_uORF")),
  paste0("dm_support_class_",
         c("clean_CDS", "leader_uORF", "overlapping_uORF"))
))
agreement <- merge(
  agreement, select_columns(dm, dm_keep),
  by = key, all.x = TRUE, sort = FALSE
)

if (nrow(motif)) {
  motif_keep <- existing_columns(motif, c(
    key,
    "motif_prior_review_class",
    "motif_prior_review_priority",
    "motif_top_feature_class",
    "motif_top_concordance_class",
    "motif_top_delta",
    "motif_top_context"
  ))
  agreement <- merge(
    agreement, select_columns(motif, motif_keep),
    by = key, all.x = TRUE, sort = FALSE
  )
}

if (nrow(browser)) {
  browser_keep <- existing_columns(browser, c(
    key,
    "browser_validation_status",
    "browser_validation_label",
    "validation_context",
    "validation_notes",
    "validation_score"
  ))
  browser_small <- select_columns(browser, browser_keep)
  if (!"validation_score" %in% names(browser_small)) {
    browser_small[, validation_score := NA_real_]
  }
  browser_small[, validation_score := safe_num(validation_score)]
  browser_small[!is.finite(validation_score), validation_score := 0.5]
  setorder(browser_small, gene_symbol, tx_id, design_family, -validation_score)
  browser_small <- browser_small[, .SD[1], by = key]
  browser_rename <- c(
    validation_context = "browser_validation_context",
    validation_notes = "browser_validation_notes",
    validation_score = "browser_validation_score"
  )
  browser_present <- names(browser_rename)[
    names(browser_rename) %in% names(browser_small)
  ]
  setnames(browser_small, browser_present, browser_rename[browser_present])
  agreement <- merge(agreement, browser_small, by = key,
                     all.x = TRUE, sort = FALSE)
}

for (column in c(
  "clean_cds_evidence_class", "clean_cds_direction",
  "joint_top_review_class", "joint_top_branch_class",
  "hierarchical_joint_top_support_class",
  "hierarchical_joint_top_branch_class",
  "dm_top_support_class", "dm_top_branch_class", "dm_top_prior_source",
  "motif_prior_review_class", "motif_top_concordance_class",
  "browser_validation_status", "browser_validation_label",
  "browser_validation_context", "browser_validation_notes"
)) {
  if (!column %in% names(agreement)) agreement[, (column) := ""]
  agreement[, (column) := fifelse(is.na(get(column)), "", as.character(get(column)))]
}
for (column in c(
  "clean_cds_priority", "clean_cds_pct_change",
  "clean_cds_direction_fraction", "clean_cds_studies",
  "clean_cds_strict_studies", "joint_top_review_priority",
  "joint_top_branch_delta", "joint_top_l1_shift", "joint_strict_studies",
  "joint_evaluable_studies", "hierarchical_joint_top_priority",
  "hierarchical_joint_top_delta", "hierarchical_joint_top_studies",
  "hierarchical_joint_top_direction_fraction", "hierarchical_joint_top_i2",
  "dm_top_priority", "dm_top_delta", "dm_top_studies",
  "dm_top_direction_fraction", "dm_top_i2",
  "motif_prior_review_priority", "browser_validation_score"
)) {
  if (!column %in% names(agreement)) agreement[, (column) := NA_real_]
  agreement[, (column) := safe_num(get(column))]
}
for (column in c(
  "clean_cds_direction_stable", "clean_cds_ci_excludes_zero",
  "clean_cds_calibrated_ci_excludes_zero"
)) {
  if (!column %in% names(agreement)) agreement[, (column) := FALSE]
  agreement[, (column) := logical_safe(get(column))]
}

agreement[, clean_cds_supported :=
  clean_cds_evidence_class != "hierarchical_low_support" &
    clean_cds_direction_stable &
    is.finite(clean_cds_priority) & clean_cds_priority >= 0.45]
agreement[, clean_cds_strong :=
  clean_cds_supported & clean_cds_ci_excludes_zero &
    clean_cds_priority >= 0.60]
agreement[, joint_supported := joint_top_review_class %chin% c(
  "strict_joint_branch_shift",
  "joint_branch_shift_review",
  "motif_joint_discordant_review"
)]
agreement[, joint_strict := joint_top_review_class == "strict_joint_branch_shift"]
agreement[, hierarchical_joint_supported :=
  hierarchical_joint_top_support_class %chin% c(
    "hierarchical_joint_replicated_shift",
    "hierarchical_joint_replicated_review",
    "hierarchical_joint_single_study_strong"
  )]
agreement[, hierarchical_joint_replicated :=
  hierarchical_joint_top_support_class %chin% c(
    "hierarchical_joint_replicated_shift",
    "hierarchical_joint_replicated_review"
  )]
agreement[, dm_supported := dm_top_support_class %chin% c(
  "dm_replicated_shift", "dm_replicated_review", "dm_single_study_strong"
)]
agreement[, dm_replicated := dm_top_support_class %chin% c(
  "dm_replicated_shift", "dm_replicated_review"
)]
agreement[, browser_supported :=
  is.finite(browser_validation_score) & browser_validation_score > 0]
agreement[, motif_concordant :=
  motif_prior_review_class == "motif_prior_concordant_branch_evidence"]
agreement[, motif_discordant :=
  motif_prior_review_class == "motif_prior_discordant_branch_evidence" |
    motif_top_concordance_class == "motif_branch_discordant"]

branch_levels <- c("clean_CDS", "leader_uORF", "overlapping_uORF")
hierarchical_supported_classes <- c(
  "hierarchical_joint_replicated_shift",
  "hierarchical_joint_replicated_review",
  "hierarchical_joint_single_study_strong"
)
hierarchical_replicated_classes <- c(
  "hierarchical_joint_replicated_shift",
  "hierarchical_joint_replicated_review"
)
dm_supported_classes <- c(
  "dm_replicated_shift", "dm_replicated_review", "dm_single_study_strong"
)
dm_replicated_classes <- c("dm_replicated_shift", "dm_replicated_review")

branch_comparison <- rbindlist(lapply(branch_levels, function(branch) {
  h_delta <- safe_num(column_or(
    agreement, paste0("hierarchical_joint_delta_", branch)
  ))
  h_priority <- safe_num(column_or(
    agreement, paste0("hierarchical_joint_priority_", branch)
  ))
  h_class <- as.character(column_or(
    agreement, paste0("hierarchical_joint_support_class_", branch), ""
  ))
  dm_delta <- safe_num(column_or(agreement, paste0("dm_delta_", branch)))
  dm_priority <- safe_num(column_or(agreement, paste0("dm_priority_", branch)))
  dm_class <- as.character(column_or(
    agreement, paste0("dm_support_class_", branch), ""
  ))
  joint_delta <- safe_num(column_or(
    agreement, paste0("joint_weighted_", branch, "_delta")
  ))
  h_supported <- h_class %chin% hierarchical_supported_classes
  dm_supported <- dm_class %chin% dm_supported_classes
  both_supported <- h_supported & dm_supported
  direction_agree <- both_supported & is.finite(h_delta) & is.finite(dm_delta) &
    sign(h_delta) == sign(dm_delta) & sign(h_delta) != 0
  direction_opposite <- both_supported & is.finite(h_delta) & is.finite(dm_delta) &
    sign(h_delta) != sign(dm_delta) & sign(h_delta) != 0 & sign(dm_delta) != 0
  clean_output_sign <- sign(safe_num(agreement$clean_cds_pct_change))
  allocation_sign <- sign(rowMeans(
    cbind(
      fifelse(h_supported & is.finite(h_delta), h_delta, NA_real_),
      fifelse(dm_supported & is.finite(dm_delta), dm_delta, NA_real_)
    ),
    na.rm = TRUE
  ))
  allocation_sign[!is.finite(allocation_sign)] <- NA_real_
  data.table(
    agreement[, ..key],
    branch_class = branch,
    hierarchical_joint_delta = h_delta,
    hierarchical_joint_priority = h_priority,
    hierarchical_joint_support_class = h_class,
    hierarchical_joint_supported = h_supported,
    hierarchical_joint_replicated =
      h_class %chin% hierarchical_replicated_classes,
    dm_delta = dm_delta,
    dm_priority = dm_priority,
    dm_support_class = dm_class,
    dm_supported = dm_supported,
    dm_replicated = dm_class %chin% dm_replicated_classes,
    joint_weighted_delta = joint_delta,
    both_supported = both_supported,
    direction_agree = direction_agree,
    direction_opposite = direction_opposite,
    effect_abs_difference = abs(h_delta - dm_delta),
    branch_agreement_class = fcase(
      direction_agree, "both_supported_same_direction",
      direction_opposite, "both_supported_opposite_direction",
      h_supported & !dm_supported, "hierarchical_joint_only",
      dm_supported & !h_supported, "dirichlet_multinomial_only",
      both_supported, "both_supported_zero_or_undefined",
      default = "neither_supported"
    ),
    clean_cds_output_allocation_alignment = fifelse(
      branch == "clean_CDS" & agreement$clean_cds_supported &
        is.finite(allocation_sign) & clean_output_sign != 0,
      fifelse(clean_output_sign == allocation_sign,
              "same_direction", "opposite_direction"),
      "not_evaluable"
    )
  )
}), use.names = TRUE, fill = TRUE)

branch_summary <- branch_comparison[, .(
  common_supported_branches = sum(both_supported),
  agreeing_common_branches = sum(direction_agree),
  disagreeing_common_branches = sum(direction_opposite),
  hierarchical_only_branches =
    sum(branch_agreement_class == "hierarchical_joint_only"),
  dm_only_branches =
    sum(branch_agreement_class == "dirichlet_multinomial_only"),
  common_supported_branch_names =
    collapse_unique(branch_class[both_supported], n = 3L),
  agreeing_common_branch_names =
    collapse_unique(branch_class[direction_agree], n = 3L),
  disagreeing_common_branch_names =
    collapse_unique(branch_class[direction_opposite], n = 3L),
  max_branch_model_abs_difference =
    if (any(is.finite(effect_abs_difference))) {
      max(effect_abs_difference, na.rm = TRUE)
    } else {
      NA_real_
    },
  clean_cds_output_allocation_alignment =
    clean_cds_output_allocation_alignment[branch_class == "clean_CDS"][1]
), by = key]

top_disagreement <- branch_comparison[
  is.finite(effect_abs_difference)
][order(-effect_abs_difference), .SD[1], by = key][, .(
  gene_symbol, tx_id, design_family,
  top_model_disagreement_branch = branch_class,
  top_model_disagreement_abs_difference = effect_abs_difference
)]
agreement <- merge(agreement, branch_summary, by = key, all.x = TRUE,
                   sort = FALSE)
agreement <- merge(agreement, top_disagreement, by = key, all.x = TRUE,
                   sort = FALSE)

agreement[, model_support_count :=
  as.integer(clean_cds_supported) + as.integer(joint_supported) +
    as.integer(hierarchical_joint_supported) + as.integer(dm_supported)]
agreement[, replicated_model_count :=
  as.integer(clean_cds_strong) + as.integer(hierarchical_joint_replicated) +
    as.integer(dm_replicated)]
agreement[, branch_model_consensus :=
  agreeing_common_branches > 0 & disagreeing_common_branches == 0]
agreement[, branch_model_direction_disagreement :=
  disagreeing_common_branches > 0]
agreement[, different_branch_support :=
  hierarchical_joint_supported & dm_supported &
    common_supported_branches == 0]
agreement[, allocation_output_disagreement :=
  clean_cds_output_allocation_alignment == "opposite_direction"]
agreement[, high_branch_heterogeneity :=
  (is.finite(hierarchical_joint_top_i2) & hierarchical_joint_top_i2 > 0.75) |
    (is.finite(dm_top_i2) & dm_top_i2 > 0.75)]

agreement[, clean_cds_component :=
  fifelse(clean_cds_supported, clip01(clean_cds_priority), 0)]
agreement[, joint_component :=
  fifelse(joint_supported, clip01(joint_top_review_priority), 0)]
agreement[, hierarchical_joint_component :=
  fifelse(hierarchical_joint_supported,
          clip01(hierarchical_joint_top_priority), 0)]
agreement[, dm_component :=
  fifelse(dm_supported, clip01(dm_top_priority), 0)]
agreement[, browser_component :=
  fifelse(browser_supported, clip01(browser_validation_score), 0)]
agreement[, model_agreement_score := clip01(
  0.15 * clean_cds_component +
    0.15 * joint_component +
    0.25 * hierarchical_joint_component +
    0.25 * dm_component +
    0.20 * browser_component +
    0.08 * model_support_count / 4 +
    0.08 * as.numeric(hierarchical_joint_replicated & dm_replicated &
                        branch_model_consensus) +
    0.05 * as.numeric(branch_model_consensus) +
    0.04 * as.numeric(motif_concordant) -
    0.12 * as.numeric(branch_model_direction_disagreement) -
    0.05 * as.numeric(allocation_output_disagreement) -
    0.04 * as.numeric(motif_discordant)
)]

agreement[, model_agreement_class := fcase(
  browser_supported & model_support_count >= 2 &
    !branch_model_direction_disagreement,
  "browser_validated_consensus",
  hierarchical_joint_replicated & dm_replicated &
    branch_model_consensus & !allocation_output_disagreement,
  "replicated_branch_consensus",
  model_support_count >= 3 & !branch_model_direction_disagreement &
    !allocation_output_disagreement,
  "multi_model_consensus",
  branch_model_direction_disagreement,
  "branch_model_direction_disagreement",
  different_branch_support,
  "different_branch_support_review",
  allocation_output_disagreement,
  "allocation_output_disagreement_review",
  replicated_model_count >= 1,
  "replicated_single_model_review",
  model_support_count >= 1 | browser_supported,
  "single_model_review",
  default = "low_support"
)]
agreement[, model_agreement_tier := fcase(
  model_agreement_class == "browser_validated_consensus",
  "A_browser_validated_consensus",
  model_agreement_class == "replicated_branch_consensus",
  "B_replicated_branch_consensus",
  model_agreement_class == "multi_model_consensus",
  "C_multi_model_consensus",
  model_agreement_class %chin% c(
    "branch_model_direction_disagreement",
    "different_branch_support_review",
    "allocation_output_disagreement_review"
  ),
  "D_model_disagreement_review",
  model_agreement_class %chin% c(
    "replicated_single_model_review", "single_model_review"
  ),
  "E_single_model_review",
  default = "F_low_support"
)]
agreement[, model_disagreement_flag :=
  model_agreement_tier == "D_model_disagreement_review"]
agreement[, model_agreement_review_priority := pmin(
  1.5,
  model_agreement_score +
    0.12 * as.numeric(design_family == "viral_infection") +
    0.10 * as.numeric(model_disagreement_flag) +
    0.10 * as.numeric(browser_supported) +
    0.04 * as.numeric(high_branch_heterogeneity)
)]

agreement[, model_support_signature := vapply(seq_len(.N), function(i) {
  supported <- c(
    clean_CDS_output = clean_cds_supported[i],
    joint_allocation = joint_supported[i],
    hierarchical_joint = hierarchical_joint_supported[i],
    dirichlet_multinomial = dm_supported[i],
    browser = browser_supported[i],
    motif_concordant = motif_concordant[i]
  )
  names(supported)[supported] |>
    paste(collapse = "; ") |>
    (\(x) if (nzchar(x)) x else "none")()
}, character(1))]
agreement[, model_disagreement_flags := vapply(seq_len(.N), function(i) {
  flags <- c(
    branch_direction = branch_model_direction_disagreement[i],
    different_supported_branches = different_branch_support[i],
    allocation_vs_output = allocation_output_disagreement[i],
    motif_discordant = motif_discordant[i],
    high_branch_heterogeneity = high_branch_heterogeneity[i]
  )
  names(flags)[flags] |>
    paste(collapse = "; ") |>
    (\(x) if (nzchar(x)) x else "none")()
}, character(1))]
agreement[, model_agreement_top_context := fcase(
  browser_supported & nzchar(browser_validation_context),
  browser_validation_context,
  hierarchical_joint_supported & nzchar(hierarchical_joint_top_context),
  hierarchical_joint_top_context,
  dm_supported & nzchar(dm_top_context),
  dm_top_context,
  joint_supported & nzchar(joint_top_context),
  joint_top_context,
  default = ""
)]

setorder(agreement, -model_agreement_review_priority, -model_agreement_score,
         gene_symbol, design_family)
setorder(branch_comparison, gene_symbol, design_family, branch_class)

review <- agreement[
  model_agreement_tier != "F_low_support" |
    model_agreement_review_priority >= 0.45
]
setorder(review, -model_agreement_review_priority, -model_agreement_score)
viral_review <- agreement[design_family == "viral_infection"]
setorder(viral_review, -model_agreement_review_priority,
         -model_agreement_score)

summarize_scope <- function(dt, scope) {
  data.table(
    scope = scope,
    gene_designs = nrow(dt),
    genes = uniqueN(dt$gene_symbol),
    supported_gene_designs = sum(dt$model_support_count > 0),
    browser_validated_consensus =
      sum(dt$model_agreement_class == "browser_validated_consensus"),
    replicated_branch_consensus =
      sum(dt$model_agreement_class == "replicated_branch_consensus"),
    multi_model_consensus =
      sum(dt$model_agreement_class == "multi_model_consensus"),
    model_disagreement_reviews = sum(dt$model_disagreement_flag),
    single_model_reviews =
      sum(dt$model_agreement_tier == "E_single_model_review"),
    low_support = sum(dt$model_agreement_tier == "F_low_support"),
    hierarchical_dm_common_supported =
      sum(dt$common_supported_branches > 0, na.rm = TRUE),
    hierarchical_dm_direction_disagreements =
      sum(dt$branch_model_direction_disagreement, na.rm = TRUE),
    allocation_output_disagreements =
      sum(dt$allocation_output_disagreement, na.rm = TRUE),
    median_agreement_score =
      if (nrow(dt)) median(dt$model_agreement_score, na.rm = TRUE) else NA_real_,
    max_agreement_score =
      if (nrow(dt)) max(dt$model_agreement_score, na.rm = TRUE) else NA_real_
  )
}
summary_metrics <- rbindlist(list(
  summarize_scope(agreement, "all"),
  summarize_scope(agreement[design_family == "viral_infection"],
                  "design_family:viral_infection")
))

summarize_branches <- function(dt, scope) {
  dt[, .(
    comparison_rows = .N,
    hierarchical_supported = sum(hierarchical_joint_supported),
    dm_supported = sum(dm_supported),
    both_supported = sum(both_supported),
    both_supported_same_direction = sum(direction_agree),
    both_supported_opposite_direction = sum(direction_opposite),
    hierarchical_only =
      sum(branch_agreement_class == "hierarchical_joint_only"),
    dm_only =
      sum(branch_agreement_class == "dirichlet_multinomial_only"),
    pearson_delta_correlation = if (sum(both_supported) >= 3) {
      cor(hierarchical_joint_delta[both_supported],
          dm_delta[both_supported], method = "pearson")
    } else {
      NA_real_
    },
    spearman_delta_correlation = if (sum(both_supported) >= 3) {
      cor(hierarchical_joint_delta[both_supported],
          dm_delta[both_supported], method = "spearman")
    } else {
      NA_real_
    },
    median_both_supported_abs_difference = if (sum(both_supported) > 0) {
      median(effect_abs_difference[both_supported], na.rm = TRUE)
    } else {
      NA_real_
    },
    max_both_supported_abs_difference = if (sum(both_supported) > 0) {
      max(effect_abs_difference[both_supported], na.rm = TRUE)
    } else {
      NA_real_
    }
  ), by = branch_class][, scope := scope][]
}
branch_summary_metrics <- rbindlist(list(
  summarize_branches(branch_comparison, "all"),
  summarize_branches(
    branch_comparison[design_family == "viral_infection"],
    "design_family:viral_infection"
  )
), use.names = TRUE, fill = TRUE)
setcolorder(branch_summary_metrics, c(
  "scope", "branch_class",
  setdiff(names(branch_summary_metrics), c("scope", "branch_class"))
))

fwrite(branch_comparison, file.path(
  output_dir, "dominant_rdg_model_agreement_branch_comparison.csv"
))
fwrite(branch_summary_metrics, file.path(
  output_dir, "dominant_rdg_model_agreement_branch_summary.csv"
))
fwrite(agreement, file.path(
  output_dir, "dominant_rdg_model_agreement_gene_design.csv"
))
fwrite(review, file.path(
  output_dir, "dominant_rdg_model_agreement_review.csv"
))
fwrite(viral_review, file.path(
  output_dir, "dominant_rdg_model_agreement_viral_review.csv"
))
fwrite(summary_metrics, file.path(
  output_dir, "dominant_rdg_model_agreement_summary_metrics.csv"
))

scatter_dt <- branch_comparison[
  (hierarchical_joint_supported | dm_supported) &
    is.finite(hierarchical_joint_delta) & is.finite(dm_delta)
]
if (nrow(scatter_dt)) {
  scatter_dt[, plot_label := paste0(gene_symbol, " / ", design_family)]
  scatter_dt[, agreement_plot_class := fcase(
    branch_agreement_class == "both_supported_same_direction",
    "both supported",
    branch_agreement_class == "hierarchical_joint_only",
    "hierarchical only",
    branch_agreement_class == "dirichlet_multinomial_only",
    "posterior predictive only",
    branch_agreement_class == "both_supported_opposite_direction",
    "opposite direction",
    default = "other"
  )]
  label_dt <- scatter_dt[
    branch_agreement_class == "both_supported_opposite_direction" |
      effect_abs_difference >= quantile(effect_abs_difference, 0.985,
                                        na.rm = TRUE)
  ][order(-effect_abs_difference)][1:min(.N, 16)]
  p_scatter <- ggplot(
    scatter_dt,
    aes(
      x = hierarchical_joint_delta,
      y = dm_delta,
      color = branch_class,
      shape = agreement_plot_class,
      text = paste0(
        "Gene/design: ", plot_label,
        "<br>Branch: ", branch_class,
        "<br>Hierarchical delta: ", round(hierarchical_joint_delta, 3),
        "<br>Posterior-predictive delta: ", round(dm_delta, 3),
        "<br>Agreement: ", branch_agreement_class
      )
    )
  ) +
    geom_hline(yintercept = 0, color = "#d0d5dd", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "#d0d5dd", linewidth = 0.35) +
    geom_abline(slope = 1, intercept = 0, color = "#667085",
                linewidth = 0.45, linetype = "dashed") +
    geom_point(size = 2.6, alpha = 0.82) +
    geom_text(
      data = label_dt,
      aes(label = gene_symbol),
      size = 2.5, vjust = -0.7, check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      clean_CDS = "#1f5aa6",
      leader_uORF = "#c57b00",
      overlapping_uORF = "#9b1c31"
    )) +
    labs(
      title = "RDG Branch-Model Agreement",
      subtitle = "Study-aware hierarchical delta versus posterior-predictive delta",
      x = "Hierarchical joint allocation delta",
      y = "Dirichlet-multinomial posterior-predictive delta",
      color = "Branch",
      shape = "Agreement"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "model_agreement_branch_comparison.png"),
         p_scatter, width = 10.4, height = 7.6, dpi = 180)
  ggsave(file.path(figure_dir, "model_agreement_branch_comparison.pdf"),
         p_scatter, width = 10.4, height = 7.6)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_scatter,
      file.path(figure_dir, "model_agreement_branch_comparison.html"),
      title = "RDG branch-model agreement"
    )
  }
}

viral_matrix <- head(viral_review, 35)
if (nrow(viral_matrix)) {
  viral_matrix[, plot_label := paste0(gene_symbol, " / ", model_agreement_tier)]
  matrix_dt <- rbindlist(list(
    viral_matrix[, .(
      plot_label, gene_symbol, design_family,
      model = "Clean CDS output",
      state = fcase(clean_cds_strong, "replicated/strong",
                    clean_cds_supported, "supported", default = "none")
    )],
    viral_matrix[, .(
      plot_label, gene_symbol, design_family,
      model = "Joint allocation",
      state = fcase(joint_strict, "replicated/strong",
                    joint_supported, "supported", default = "none")
    )],
    viral_matrix[, .(
      plot_label, gene_symbol, design_family,
      model = "Hierarchical joint",
      state = fcase(hierarchical_joint_replicated, "replicated/strong",
                    hierarchical_joint_supported, "supported",
                    default = "none")
    )],
    viral_matrix[, .(
      plot_label, gene_symbol, design_family,
      model = "Posterior predictive",
      state = fcase(dm_replicated, "replicated/strong",
                    dm_supported, "supported", default = "none")
    )],
    viral_matrix[, .(
      plot_label, gene_symbol, design_family,
      model = "Browser",
      state = fifelse(browser_supported, "manual validation", "none")
    )],
    viral_matrix[, .(
      plot_label, gene_symbol, design_family,
      model = "Motif prior",
      state = fcase(motif_concordant, "motif concordant",
                    motif_discordant, "motif discordant", default = "none")
    )]
  ))
  matrix_dt[, model := factor(
    model,
    levels = c("Clean CDS output", "Joint allocation", "Hierarchical joint",
               "Posterior predictive", "Browser", "Motif prior")
  )]
  matrix_dt[, plot_label := factor(
    plot_label, levels = rev(unique(viral_matrix$plot_label))
  )]
  matrix_dt[, tile_label := fcase(
    state == "replicated/strong", "R",
    state == "supported", "S",
    state == "manual validation", "V",
    state == "motif concordant", "+",
    state == "motif discordant", "-",
    default = ""
  )]
  p_matrix <- ggplot(
    matrix_dt,
    aes(
      x = model,
      y = plot_label,
      fill = state,
      text = paste0(
        "Gene: ", gene_symbol,
        "<br>Model: ", model,
        "<br>State: ", state
      )
    )
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = tile_label), size = 3, color = "#202124") +
    scale_fill_manual(values = c(
      none = "#eef1f4",
      supported = "#8eb6d8",
      `replicated/strong` = "#245b9e",
      `manual validation` = "#6f2c91",
      `motif concordant` = "#4f8f62",
      `motif discordant` = "#c46b45"
    )) +
    labs(
      title = "Viral RDG Model-Agreement Matrix",
      subtitle = "R = replicated/strong, S = supported, V = browser validated",
      x = NULL,
      y = NULL,
      fill = "Evidence"
    ) +
    theme_minimal(base_size = 9.5) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_dir, "model_agreement_viral_matrix.png"),
         p_matrix, width = 10.8, height = 8.8, dpi = 180)
  ggsave(file.path(figure_dir, "model_agreement_viral_matrix.pdf"),
         p_matrix, width = 10.8, height = 8.8)
  if (exists("dominant_save_ggplotly") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      p_matrix,
      file.path(figure_dir, "model_agreement_viral_matrix.html"),
      title = "Viral RDG model-agreement matrix"
    )
  }
}

message("Saved RDG model-agreement outputs in: ", output_dir)
print(summary_metrics)
print(branch_summary_metrics)
