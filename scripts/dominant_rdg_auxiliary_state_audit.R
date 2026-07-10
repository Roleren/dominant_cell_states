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
output_dir <- file.path(analysis_dir, "dominant_rdg_auxiliary_state_audit")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
aux_score_file <- file.path(
  analysis_dir,
  "dominant_state_auxiliary_modules",
  "dominant_state_auxiliary_module_scores.csv"
)
aux_summary_file <- file.path(
  analysis_dir,
  "dominant_state_auxiliary_modules",
  "dominant_state_auxiliary_module_summary.csv"
)
joint_file <- file.path(
  analysis_dir,
  "dominant_rdg_joint_branch_allocation",
  "dominant_rdg_joint_branch_allocation_class_rows.csv"
)
atlas_file <- file.path(
  analysis_dir,
  "dominant_rdg_atlas",
  "dominant_rdg_atlas_cards.csv"
)

required_files <- c(metadata_file, aux_score_file, aux_summary_file,
                    joint_file, atlas_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required auxiliary-state audit input(s):\n",
       paste(missing_files, collapse = "\n"), call. = FALSE)
}

clip01 <- function(x) pmax(pmin(x, 1), 0)

mean_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

max_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

norm_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- "MISSING"
  x
}

split_conditions <- function(x) {
  x <- norm_text(x)
  parts <- unique(trimws(unlist(strsplit(x, "[;,|]", perl = TRUE))))
  parts <- parts[!is.na(parts) & nzchar(parts) & parts != "MISSING"]
  parts
}

make_key <- function(dt, cols) {
  values <- lapply(cols, function(col) {
    if (!col %chin% names(dt)) return(rep("MISSING", nrow(dt)))
    norm_text(dt[[col]])
  })
  do.call(paste, c(values, sep = " || "))
}

state_short_name <- function(x) {
  map <- c(
    "Antigen presentation / immune composition" = "Antigen",
    "DNA damage / p53 checkpoint" = "p53",
    "Mitochondrial UPR / mitonuclear stress" = "Mito UPR",
    "Host shutoff / antiviral translation repression" = "Host shutoff",
    "Proteostasis / Heat shock" = "Proteostasis",
    "Autophagy / Lysosome (TFEB)" = "Autophagy",
    "NRF2 / Oxidative-ferroptosis defense" = "NRF2",
    "Secretory / ER proteostasis" = "ER proteostasis"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 220, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, bg = "white")
  if (exists("dominant_save_ggplotly", mode = "function") &&
      requireNamespace("plotly", quietly = TRUE)) {
    dominant_save_ggplotly(
      plot, file.path(figure_dir, paste0(name, ".html")), title = name
    )
  }
}

message("Reading auxiliary scores, metadata, joint branch rows, and atlas cards")
aux_scores <- fread(aux_score_file, showProgress = FALSE)
aux_summary <- fread(aux_summary_file, showProgress = FALSE)
joint_rows <- fread(joint_file, showProgress = FALSE)
atlas_cards <- fread(atlas_file, showProgress = FALSE)

needed_aux <- c("Run", "candidate_state", "signed_score", "robust_z",
                "percentile")
if (!all(needed_aux %chin% names(aux_scores))) {
  stop("Auxiliary score file lacks columns: ",
       paste(setdiff(needed_aux, names(aux_scores)), collapse = ", "),
       call. = FALSE)
}

context_cols <- c("study", "AUTHOR", "CELL_LINE", "TISSUE", "GENE",
                  "INHIBITOR", "FRACTION", "Cancer_type", "Sex",
                  "TIMEPOINT")
metadata_cols <- unique(c("Run", "CONDITION", context_cols))
metadata_cols <- metadata_cols[metadata_cols %chin%
                                 names(fread(metadata_file, nrows = 0))]
metadata <- fread(metadata_file, select = metadata_cols, showProgress = FALSE)
if (!all(c("Run", "CONDITION", "study") %chin% names(metadata))) {
  stop("Metadata file must include Run, CONDITION, and study columns.",
       call. = FALSE)
}

for (col in setdiff(c("CONDITION", context_cols), names(metadata))) {
  metadata[, (col) := "MISSING"]
}
metadata[, Run := as.character(Run)]
metadata[, CONDITION := norm_text(CONDITION)]
metadata[, context_key := make_key(.SD, context_cols), .SDcols = context_cols]

aux_scores[, Run := as.character(Run)]
aux_scores <- aux_scores[, .(
  Run,
  candidate_state = as.character(candidate_state),
  aux_signed_score = as.numeric(signed_score),
  aux_robust_z = as.numeric(robust_z),
  aux_percentile = as.numeric(percentile)
)]

meta_aux <- merge(
  metadata[, c("Run", "CONDITION", "context_key", context_cols), with = FALSE],
  aux_scores,
  by = "Run",
  allow.cartesian = TRUE
)

stratum_cols <- unique(c(
  "stratum_id", "case_condition", "control_conditions", "match_scope",
  context_cols
))
stratum_cols <- stratum_cols[stratum_cols %chin% names(joint_rows)]
strata <- unique(joint_rows[, ..stratum_cols])
for (col in setdiff(context_cols, names(strata))) strata[, (col) := "MISSING"]
strata[, case_condition := norm_text(case_condition)]
strata[, control_conditions := norm_text(control_conditions)]
strata[, context_key := make_key(.SD, context_cols), .SDcols = context_cols]

case_map <- unique(strata[, .(
  stratum_id,
  context_key,
  CONDITION = case_condition,
  sample_role = "case"
)])
control_map <- unique(strata[, .(
  CONDITION = split_conditions(control_conditions)
), by = .(stratum_id, context_key)])
if (nrow(control_map)) control_map[, sample_role := "control"]
role_map <- unique(rbindlist(list(case_map, control_map), fill = TRUE))
role_map <- role_map[CONDITION != "MISSING"]

message("Reconstructing matched case/control sample sets for auxiliary states")
sample_roles <- merge(
  meta_aux,
  role_map,
  by = c("context_key", "CONDITION"),
  allow.cartesian = TRUE
)
sample_roles <- unique(
  sample_roles,
  by = c("Run", "candidate_state", "stratum_id", "sample_role")
)

stratum_aux <- sample_roles[, .(
  case_n = uniqueN(Run[sample_role == "case"]),
  control_n = uniqueN(Run[sample_role == "control"]),
  case_mean_z = mean_safe(aux_robust_z[sample_role == "case"]),
  control_mean_z = mean_safe(aux_robust_z[sample_role == "control"]),
  case_median_z = median_safe(aux_robust_z[sample_role == "case"]),
  control_median_z = median_safe(aux_robust_z[sample_role == "control"]),
  case_mean_score = mean_safe(aux_signed_score[sample_role == "case"]),
  control_mean_score = mean_safe(aux_signed_score[sample_role == "control"]),
  case_mean_percentile = mean_safe(aux_percentile[sample_role == "case"]),
  control_mean_percentile = mean_safe(aux_percentile[sample_role == "control"]),
  case_runs = paste(sort(unique(Run[sample_role == "case"])), collapse = "; "),
  control_runs = paste(sort(unique(Run[sample_role == "control"])), collapse = "; ")
), by = .(stratum_id, candidate_state)]

all_strata_states <- CJ(
  stratum_id = unique(strata$stratum_id),
  candidate_state = unique(aux_scores$candidate_state),
  unique = TRUE
)
stratum_aux <- merge(
  all_strata_states,
  stratum_aux,
  by = c("stratum_id", "candidate_state"),
  all.x = TRUE
)
stratum_aux[is.na(case_n), case_n := 0L]
stratum_aux[is.na(control_n), control_n := 0L]
stratum_aux[, `:=`(
  aux_delta_z = case_mean_z - control_mean_z,
  aux_delta_score = case_mean_score - control_mean_score,
  aux_delta_percentile = case_mean_percentile - control_mean_percentile,
  aux_evaluable = case_n >= 2L & control_n >= 2L
)]
stratum_aux[, abs_aux_delta_z := abs(aux_delta_z)]
stratum_aux[, aux_shift_strength := fcase(
  !aux_evaluable, "not_evaluable",
  abs_aux_delta_z >= 1.0, "strong_aux_shift",
  abs_aux_delta_z >= 0.5, "moderate_aux_shift",
  default = "low_aux_shift"
)]
strata_labels <- strata[, c(
  "stratum_id", "case_condition", "control_conditions", "match_scope",
  context_cols
), with = FALSE]
strata_labels <- unique(strata_labels, by = "stratum_id")
strata_labels[, stratum_label := paste0(
  study, " | ", CELL_LINE, " | ", TISSUE, " | ", case_condition,
  " vs ", control_conditions
)]
stratum_aux <- merge(
  stratum_aux,
  strata_labels,
  by = "stratum_id",
  all.x = TRUE
)

branch_keep <- c(
  "stratum_id", "gene_symbol", "tx_id", "design_family", "branch_class",
  "joint_delta", "joint_delta_se", "joint_delta_p", "joint_delta_q",
  "joint_case_count", "joint_control_count", "joint_evaluable",
  "joint_strict", "joint_l1_allocation_shift", "joint_js_divergence",
  "joint_chisq_q", "joint_review_class", "joint_review_priority"
)
branch_keep <- branch_keep[branch_keep %chin% names(joint_rows)]
branch_rows <- merge(
  joint_rows[, ..branch_keep],
  stratum_aux,
  by = "stratum_id",
  allow.cartesian = TRUE
)
branch_rows[, `:=`(
  gene_symbol = as.character(gene_symbol),
  tx_id = as.character(tx_id),
  branch_class = as.character(branch_class),
  joint_delta = as.numeric(joint_delta),
  joint_delta_q = as.numeric(joint_delta_q),
  joint_review_priority = as.numeric(joint_review_priority),
  joint_strict = as.logical(joint_strict)
)]
branch_rows[, branch_supported := joint_strict == TRUE &
              joint_review_class == "strict_joint_branch_shift"]
branch_rows[, branch_abs_delta := abs(joint_delta)]
branch_rows[, aux_branch_alignment := fcase(
  !aux_evaluable | !is.finite(joint_delta) | !is.finite(aux_delta_z),
  "not_evaluable",
  joint_delta == 0 | aux_delta_z == 0, "zero_or_tied",
  sign(joint_delta) == sign(aux_delta_z), "same_direction",
  default = "opposite_direction"
)]
branch_rows[, branch_aux_context_class := fcase(
  !aux_evaluable, "aux_not_evaluable",
  branch_supported & abs_aux_delta_z >= 1.0,
  "supported_branch_in_strong_aux_context",
  branch_supported & abs_aux_delta_z >= 0.5,
  "supported_branch_in_moderate_aux_context",
  branch_supported, "supported_branch_without_large_aux_context",
  !branch_supported & abs_aux_delta_z >= 1.0,
  "aux_context_shift_without_supported_branch",
  default = "low_aux_low_branch"
)]
branch_rows[, aux_branch_overlap_score :=
              clip01(abs_aux_delta_z / 2) *
              (0.4 + 0.6 * clip01(joint_review_priority)) *
              fifelse(branch_supported, 1.0, 0.55)]
branch_output_filter <- branch_rows$aux_evaluable == TRUE &
  (
    branch_rows$branch_supported == TRUE |
      branch_rows$abs_aux_delta_z >= 1.0 |
      branch_rows$branch_aux_context_class %chin% c(
        "supported_branch_in_strong_aux_context",
        "supported_branch_in_moderate_aux_context",
        "aux_context_shift_without_supported_branch"
      )
  )
branch_write_cols <- c(
  "stratum_id", "gene_symbol", "tx_id", "design_family", "branch_class",
  "joint_delta", "joint_delta_se", "joint_delta_q", "joint_case_count",
  "joint_control_count", "joint_evaluable", "joint_strict",
  "joint_review_class", "joint_review_priority", "candidate_state",
  "case_n", "control_n", "case_mean_z", "control_mean_z", "aux_delta_z",
  "aux_delta_score", "aux_delta_percentile", "abs_aux_delta_z",
  "aux_evaluable", "aux_shift_strength", "branch_supported",
  "branch_abs_delta", "aux_branch_alignment", "branch_aux_context_class",
  "aux_branch_overlap_score"
)
branch_write_cols <- intersect(branch_write_cols, names(branch_rows))
branch_rows_for_write <- branch_rows[branch_output_filter, ..branch_write_cols]

message("Summarizing auxiliary-state overlap per RDG gene/design")
gene_design_aux <- branch_rows[aux_evaluable == TRUE, {
  score <- aux_branch_overlap_score
  if (!any(is.finite(score))) score <- rep(NA_real_, .N)
  idx <- if (any(is.finite(score))) which.max(score) else 1L
  .(
    n_aux_evaluable_strata = uniqueN(stratum_id),
    n_branch_rows = .N,
    n_supported_branch_rows = sum(branch_supported, na.rm = TRUE),
    n_strong_aux_rows = sum(abs_aux_delta_z >= 1.0, na.rm = TRUE),
    n_moderate_aux_rows = sum(abs_aux_delta_z >= 0.5, na.rm = TRUE),
    max_abs_aux_delta_z = max_safe(abs_aux_delta_z),
    mean_abs_aux_delta_z = mean_safe(abs_aux_delta_z),
    max_aux_branch_overlap_score = max_safe(aux_branch_overlap_score),
    top_aux_delta_z = aux_delta_z[idx],
    top_aux_delta_score = aux_delta_score[idx],
    top_aux_delta_percentile = aux_delta_percentile[idx],
    top_branch_class = branch_class[idx],
    top_joint_delta = joint_delta[idx],
    top_joint_delta_q = joint_delta_q[idx],
    top_joint_review_class = joint_review_class[idx],
    top_joint_review_priority = joint_review_priority[idx],
    top_aux_context_class = branch_aux_context_class[idx],
    top_aux_branch_alignment = aux_branch_alignment[idx],
    top_context = stratum_label[idx],
    top_case_n = case_n[idx],
    top_control_n = control_n[idx],
    top_case_condition = case_condition[idx],
    top_control_conditions = control_conditions[idx],
    top_match_scope = match_scope[idx]
  )
}, by = .(gene_symbol, tx_id, design_family, candidate_state)]

gene_design_aux[, aux_overlap_class := fcase(
  n_supported_branch_rows > 0 & max_abs_aux_delta_z >= 1.0,
  "supported_rdg_in_strong_aux_context",
  n_supported_branch_rows > 0 & max_abs_aux_delta_z >= 0.5,
  "supported_rdg_in_moderate_aux_context",
  n_supported_branch_rows > 0,
  "supported_rdg_without_large_aux_context",
  max_abs_aux_delta_z >= 1.0,
  "aux_context_shift_only",
  default = "low_aux_overlap"
)]
gene_design_aux[, aux_state_short := state_short_name(candidate_state)]
gene_design_aux <- merge(
  gene_design_aux,
  aux_summary[, .(
    candidate_state,
    auxiliary_recommendation = recommendation,
    auxiliary_recommendation_label = recommendation_label,
    auxiliary_novelty_score = novelty_score,
    auxiliary_gene_fraction_present = gene_fraction_present
  )],
  by = "candidate_state",
  all.x = TRUE
)
gene_design_aux[, aux_readiness_weight := fcase(
  auxiliary_recommendation == "candidate_dominant_state_for_pipeline", 1.00,
  auxiliary_recommendation == "exploratory_screen_only", 0.75,
  auxiliary_recommendation == "needs_fst_marker_expansion", 0.45,
  auxiliary_recommendation == "mostly_existing_state_alias", 0.35,
  default = 0.50
)]
gene_design_aux[, max_weighted_aux_branch_overlap_score :=
                  max_aux_branch_overlap_score * aux_readiness_weight]
setorder(gene_design_aux, -max_weighted_aux_branch_overlap_score,
         -max_aux_branch_overlap_score,
         -max_abs_aux_delta_z, gene_symbol, design_family, candidate_state)

top_aux <- gene_design_aux[is.finite(max_aux_branch_overlap_score),
                           .SD[which.max(max_weighted_aux_branch_overlap_score)],
                           by = .(gene_symbol, tx_id, design_family)]
setnames(
  top_aux,
  old = c("candidate_state", "aux_state_short", "aux_overlap_class",
          "max_abs_aux_delta_z", "max_aux_branch_overlap_score",
          "max_weighted_aux_branch_overlap_score",
          "aux_readiness_weight",
          "top_aux_delta_z", "top_branch_class", "top_joint_delta",
          "top_joint_review_class", "top_joint_review_priority",
          "top_aux_context_class", "top_aux_branch_alignment", "top_context",
          "top_case_n", "top_control_n", "auxiliary_recommendation_label"),
  new = c("aux_top_state", "aux_top_state_short", "aux_top_overlap_class",
          "aux_top_abs_delta_z", "aux_top_overlap_score",
          "aux_top_weighted_overlap_score",
          "aux_top_readiness_weight",
          "aux_top_delta_z", "aux_top_branch_class", "aux_top_joint_delta",
          "aux_top_joint_review_class", "aux_top_joint_review_priority",
          "aux_top_context_class", "aux_top_branch_alignment",
          "aux_top_context", "aux_top_case_n", "aux_top_control_n",
          "aux_top_recommendation_label")
)

atlas_keep <- c(
  "gene_symbol", "tx_id", "design_family", "atlas_evidence_class",
  "atlas_confidence_score", "atlas_browser_review_priority_score",
  "atlas_direction", "atlas_clean_cds_pct_change",
  "model_agreement_class", "model_agreement_tier",
  "context_interaction_gene_class", "partial_pooling_gene_class",
  "browser_validation_status", "browser_validation_label"
)
atlas_keep <- atlas_keep[atlas_keep %chin% names(atlas_cards)]
atlas_annotations <- merge(
  atlas_cards[, ..atlas_keep],
  top_aux[, c(
    "gene_symbol", "tx_id", "design_family", "aux_top_state",
    "aux_top_state_short", "aux_top_overlap_class", "aux_top_abs_delta_z",
    "aux_top_overlap_score", "aux_top_weighted_overlap_score",
    "aux_top_readiness_weight", "aux_top_delta_z", "aux_top_branch_class",
    "aux_top_joint_delta", "aux_top_joint_review_class",
    "aux_top_joint_review_priority", "aux_top_context_class",
    "aux_top_branch_alignment", "aux_top_context", "aux_top_case_n",
    "aux_top_control_n", "aux_top_recommendation_label"
  ), with = FALSE],
  by = c("gene_symbol", "tx_id", "design_family"),
  all.x = TRUE
)
atlas_annotations[, aux_audit_warning := fcase(
  aux_top_overlap_class == "supported_rdg_in_strong_aux_context",
  "Strong auxiliary-state shift in matched RDG context; inspect as possible state/context driver.",
  aux_top_overlap_class == "supported_rdg_in_moderate_aux_context",
  "Moderate auxiliary-state shift in matched RDG context; treat as possible confounder.",
  aux_top_overlap_class == "aux_context_shift_only",
  "Auxiliary state shifts, but RDG branch support is weak in this audit.",
  default = ""
)]
atlas_annotations[, aux_audit_priority := fifelse(
  aux_top_overlap_class %chin% c(
    "supported_rdg_in_strong_aux_context",
    "supported_rdg_in_moderate_aux_context"
  ),
  clip01(aux_top_weighted_overlap_score) *
    (0.5 + 0.5 * clip01(atlas_browser_review_priority_score)),
  0
)]
setorder(atlas_annotations, -aux_audit_priority,
         -atlas_browser_review_priority_score, gene_symbol, design_family)

selected_states <- c(
  "Antigen presentation / immune composition",
  "DNA damage / p53 checkpoint",
  "Mitochondrial UPR / mitonuclear stress",
  "Host shutoff / antiviral translation repression"
)
selected_long <- gene_design_aux[candidate_state %chin% selected_states]
selected_wide <- dcast(
  selected_long,
  gene_symbol + tx_id + design_family ~ aux_state_short,
  value.var = "max_abs_aux_delta_z"
)

summary_by_state <- gene_design_aux[, .(
  n_gene_designs = uniqueN(paste(gene_symbol, tx_id, design_family)),
  n_supported_strong_aux = sum(
    aux_overlap_class == "supported_rdg_in_strong_aux_context", na.rm = TRUE
  ),
  n_supported_moderate_aux = sum(
    aux_overlap_class == "supported_rdg_in_moderate_aux_context", na.rm = TRUE
  ),
  n_aux_shift_only = sum(aux_overlap_class == "aux_context_shift_only",
                         na.rm = TRUE),
  median_max_abs_aux_delta_z = median_safe(max_abs_aux_delta_z),
  max_abs_aux_delta_z = max_safe(max_abs_aux_delta_z),
  max_aux_branch_overlap_score = max_safe(max_aux_branch_overlap_score),
  max_weighted_aux_branch_overlap_score =
    max_safe(max_weighted_aux_branch_overlap_score)
), by = .(candidate_state, aux_state_short, auxiliary_recommendation_label)]
setorder(summary_by_state, -n_supported_strong_aux, -n_supported_moderate_aux,
         -max_weighted_aux_branch_overlap_score)

summary_metrics <- data.table(
  metric = c(
    "strata_total",
    "strata_with_aux_case_control",
    "sample_role_rows",
    "branch_aux_rows",
    "branch_aux_rows_written",
    "gene_design_aux_rows",
    "atlas_rows",
    "atlas_rows_with_aux_annotation",
    "atlas_rows_strong_aux_overlap",
    "atlas_rows_moderate_aux_overlap",
    "atlas_rows_aux_shift_only",
    "top_raw_aux_overlap_state",
    "top_weighted_aux_overlap_state"
  ),
  value = c(
    as.character(nrow(strata)),
    as.character(uniqueN(stratum_aux[aux_evaluable == TRUE, stratum_id])),
    as.character(nrow(sample_roles)),
    as.character(nrow(branch_rows)),
    as.character(nrow(branch_rows_for_write)),
    as.character(nrow(gene_design_aux)),
    as.character(nrow(atlas_annotations)),
    as.character(sum(!is.na(atlas_annotations$aux_top_state))),
    as.character(sum(
      atlas_annotations$aux_top_overlap_class ==
        "supported_rdg_in_strong_aux_context",
      na.rm = TRUE
    )),
    as.character(sum(
      atlas_annotations$aux_top_overlap_class ==
        "supported_rdg_in_moderate_aux_context",
      na.rm = TRUE
    )),
    as.character(sum(
      atlas_annotations$aux_top_overlap_class == "aux_context_shift_only",
      na.rm = TRUE
    )),
    summary_by_state[
      which.max(max_aux_branch_overlap_score),
      aux_state_short
    ],
    summary_by_state[
      which.max(max_weighted_aux_branch_overlap_score),
      aux_state_short
    ]
  )
)

message("Writing auxiliary-state RDG audit tables")
fwrite(
  stratum_aux,
  file.path(output_dir, "dominant_rdg_auxiliary_state_stratum_scores.csv")
)
fwrite(
  branch_rows_for_write,
  file.path(output_dir, "dominant_rdg_auxiliary_state_branch_rows.csv")
)
fwrite(
  gene_design_aux,
  file.path(output_dir, "dominant_rdg_auxiliary_state_gene_design_summary.csv")
)
fwrite(
  atlas_annotations,
  file.path(output_dir, "dominant_rdg_auxiliary_state_atlas_annotations.csv")
)
fwrite(
  selected_wide,
  file.path(output_dir, "dominant_rdg_auxiliary_state_selected_state_matrix.csv")
)
fwrite(
  summary_by_state,
  file.path(output_dir, "dominant_rdg_auxiliary_state_summary_by_state.csv")
)
fwrite(
  summary_metrics,
  file.path(output_dir, "dominant_rdg_auxiliary_state_summary_metrics.csv")
)

if (nrow(gene_design_aux)) {
  plot_dt <- gene_design_aux[
    aux_overlap_class %chin% c(
      "supported_rdg_in_strong_aux_context",
      "supported_rdg_in_moderate_aux_context",
      "aux_context_shift_only"
    ) | max_aux_branch_overlap_score >= 0.25
  ]
  if (nrow(plot_dt)) {
    top_labels <- unique(
      plot_dt[order(-max_aux_branch_overlap_score),
              paste(gene_symbol, design_family, sep = " / ")][1:60]
    )
    plot_dt <- plot_dt[
      paste(gene_symbol, design_family, sep = " / ") %chin% top_labels
    ]
    plot_dt[, gene_design_label := paste(gene_symbol, design_family, sep = " / ")]
    plot_dt[, gene_design_label := factor(
      gene_design_label,
      levels = rev(unique(plot_dt[
        order(max_aux_branch_overlap_score),
        gene_design_label
      ]))
    )]
    plot_dt[, support_flag := ifelse(n_supported_branch_rows > 0,
                                     "supported branch", "aux shift only")]
    p_heat <- ggplot(
      plot_dt,
      aes(x = aux_state_short, y = gene_design_label, fill = top_aux_delta_z)
    ) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_point(
        aes(size = max_aux_branch_overlap_score, shape = support_flag),
        color = "#222222", alpha = 0.75, stroke = 0.25
      ) +
      scale_fill_gradient2(
        low = "#2B6CB0", mid = "white", high = "#8B1A1A",
        midpoint = 0, na.value = "#EFEFEF",
        name = "case-control\naux z"
      ) +
      scale_size_continuous(
        range = c(0.9, 4.2), limits = c(0, 1),
        name = "overlap\nscore"
      ) +
      scale_shape_manual(
        values = c("supported branch" = 21, "aux shift only" = 4),
        name = "RDG evidence"
      ) +
      labs(
        title = "Auxiliary State Overlap With RDG Branch Contexts",
        subtitle = "Color is case-control auxiliary-state z delta reconstructed from matched metadata; points mark overlap score with RDG branch evidence.",
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "right",
        plot.title = element_text(face = "bold")
      )
    save_plot(p_heat, "auxiliary_state_rdg_overlap_heatmap", 12, 10)
  }

  count_dt <- gene_design_aux[
    aux_overlap_class %chin% c(
      "supported_rdg_in_strong_aux_context",
      "supported_rdg_in_moderate_aux_context"
    ),
    .N,
    by = .(design_family, aux_state_short, aux_overlap_class)
  ]
  if (nrow(count_dt)) {
    count_dt[, design_family := factor(
      design_family,
      levels = unique(count_dt[, .(N = sum(N)), by = design_family][
        order(N), design_family
      ])
    )]
    p_counts <- ggplot(
      count_dt,
      aes(x = design_family, y = N, fill = aux_state_short)
    ) +
      geom_col(width = 0.75) +
      coord_flip() +
      facet_wrap(~ aux_overlap_class, ncol = 1, scales = "free_y") +
      scale_fill_brewer(type = "qual", palette = "Dark2", name = "Aux state") +
      labs(
        title = "RDG Candidates With Auxiliary-State Context Shifts",
        subtitle = "Counts are gene/design/state summaries where a supported RDG branch occurs in a moderate or strong auxiliary-state shift.",
        x = NULL,
        y = "Gene/design/state rows"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
      )
    save_plot(p_counts, "auxiliary_state_overlap_counts", 10, 7)
  }
}

scatter_dt <- atlas_annotations[
  is.finite(aux_top_abs_delta_z) &
    is.finite(atlas_clean_cds_pct_change) &
    is.finite(atlas_browser_review_priority_score)
]
if (nrow(scatter_dt)) {
  scatter_dt[, aux_top_state_short := factor(
    aux_top_state_short,
    levels = unique(scatter_dt[
      order(-aux_audit_priority),
      aux_top_state_short
    ])
  )]
  p_scatter <- ggplot(
    scatter_dt,
    aes(
      x = atlas_clean_cds_pct_change,
      y = aux_top_abs_delta_z,
      color = aux_top_state_short,
      size = atlas_browser_review_priority_score
    )
  ) +
    geom_hline(yintercept = c(0.5, 1.0), linetype = c("dotted", "dashed"),
               color = "#666666", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "#888888", linewidth = 0.3) +
    geom_point(alpha = 0.78) +
    scale_size_continuous(range = c(1.2, 5), name = "browser\npriority") +
    scale_color_brewer(type = "qual", palette = "Dark2", name = "Top aux state") +
    labs(
      title = "Auxiliary-State Shift Versus Atlas Clean-CDS Effect",
      subtitle = "Dashed lines mark moderate and strong auxiliary-state case-control deltas; points are atlas gene/design rows.",
      x = "Atlas clean-CDS effect (%)",
      y = "Top auxiliary-state absolute z delta"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  save_plot(p_scatter, "auxiliary_state_vs_atlas_effect", 10, 7)
}

message("Saved auxiliary-state RDG audit outputs to: ", output_dir)
