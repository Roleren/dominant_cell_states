#!/usr/bin/env Rscript

# Convert completed active-review rows into normalized labels and candidate
# calibration inputs. This script never edits the curated atlas validation file.

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
})

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
output_dir <- file.path(analysis_dir, "dominant_rdg_review_feedback")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

template_file <- file.path(
  analysis_dir,
  "dominant_rdg_review_links",
  "dominant_rdg_review_batch_template_with_links.csv"
)
draft_file <- file.path(
  output_dir,
  "dominant_rdg_review_label_drafts.csv"
)

labels_out <- file.path(
  output_dir,
  "dominant_rdg_review_feedback_labels.csv"
)
training_out <- file.path(
  output_dir,
  "dominant_rdg_review_feedback_training_set.csv"
)
browser_notes_out <- file.path(
  output_dir,
  "dominant_rdg_review_feedback_browser_validation_notes_candidate.csv"
)
summary_out <- file.path(
  output_dir,
  "dominant_rdg_review_feedback_summary.csv"
)
batch_summary_out <- file.path(
  output_dir,
  "dominant_rdg_review_feedback_batch_summary.csv"
)
report_out <- file.path(
  output_dir,
  "dominant_rdg_review_feedback_report.md"
)

safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(data.table())
  data.table::fread(path, showProgress = FALSE, nThread = 1, ...)
}

read_required <- function(path, label) {
  dt <- safe_fread(path)
  if (nrow(dt) == 0) {
    stop("Missing or empty ", label, ": ", path, call. = FALSE)
  }
  dt
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

has_text <- function(x) {
  x <- clean_text(x)
  nzchar(x) & !tolower(x) %chin% c("na", "nan", "null", "missing")
}

norm_text <- function(x) {
  x <- tolower(clean_text(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

contains_any <- function(x, pattern) {
  grepl(pattern, norm_text(x), perl = TRUE)
}

ensure_columns <- function(dt, columns, value = "") {
  for (column in columns) {
    if (!column %in% names(dt)) dt[, (column) := value]
  }
  invisible(dt)
}

direction_status <- function(x) {
  y <- norm_text(x)
  fifelse(
    !nzchar(y), "",
    fifelse(
      grepl("^(yes|y|true|verified|match|matched|same|expected|pass|passed)$",
            y) |
        grepl("verified|matches_expected|direction_match|confirmed", y),
      "verified",
      fifelse(
        grepl("^(no|n|false|opposite|mismatch|failed|fail)$", y) |
          grepl("opposite|mismatch|wrong_direction|direction_fail", y),
        "opposite",
        fifelse(
          grepl("unclear|ambiguous|mixed|partial|weak", y),
          "unclear",
          fifelse(grepl("not_applicable|not_relevant|na", y),
                  "not_applicable", "other")
        )
      )
    )
  )
}

direction_score <- function(x) {
  fifelse(
    x == "verified", 1,
    fifelse(x == "opposite", -1,
            fifelse(x %chin% c("", "unclear", "not_applicable"),
                    NA_real_, 0))
  )
}

label_from_fields <- function(manual_status, browser_grade, branch_dir,
                              clean_cds_dir, context_match, translon_status,
                              confounder_status, negative_control, batch_role) {
  manual <- norm_text(manual_status)
  browser <- norm_text(browser_grade)
  translon <- norm_text(translon_status)
  confounder <- norm_text(confounder_status)
  neg <- norm_text(negative_control)
  role <- norm_text(batch_role)
  direction <- norm_text(paste(branch_dir, clean_cds_dir, context_match))
  combined <- norm_text(paste(
    manual, browser, translon, confounder, neg, direction, role
  ))
  combined_without_confounder <- norm_text(paste(
    manual, browser, translon, neg, direction, role
  ))
  confounder_positive <- grepl(
    "^(possible_confounder|confounded|wrong_context)$|metadata|context_wrong",
    confounder
  )

  if (nzchar(neg) || grepl("negative_control", role)) {
    if (grepl("fail|failed|unexpected|positive|present|strong_signal|signal_present",
              neg)) {
      return("negative_control_fail")
    }
    if (grepl("pass|passed|absent|no_signal|none|clean|negative", neg)) {
      return("negative_control_pass")
    }
  }

  if (confounder_positive ||
      grepl("confound|confounded|metadata|context_wrong|wrong_context",
            combined_without_confounder)) {
    return("confounded")
  }
  if (grepl("measurement|translon|isoform|wrong_isoform|wrong_transcript|p_shift|frameshift|frame_problem|overlap_problem",
            combined)) {
    return("measurement_problem")
  }
  if (grepl("not_reviewable|unreviewable|unclear|ambiguous|too_sparse|low_count|insufficient|mixed",
            manual)) {
    return("not_reviewable")
  }
  if (grepl("negative|reject|rejected|refute|refuted|absent|no_signal|none_visible|not_supported|unsupported|no_support|false|(^|_)no($|_)",
            manual)) {
    return("negative")
  }
  if (grepl("negative|absent|no_signal|none|not_supported|unsupported|no_support",
            browser)) {
    return("negative")
  }
  if (grepl("weak|partial|suggestive|maybe", manual) ||
      grepl("weak|partial|suggestive", browser)) {
    return("weak_positive")
  }
  if (grepl("positive|support|supported|valid|validated|confirm|confirmed|yes|true|interesting",
            manual) ||
      grepl("strong|clear|support|supported|positive|interesting", browser)) {
    return("positive")
  }
  if (grepl("opposite", direction)) {
    return("negative")
  }

  ""
}

label_score <- function(label) {
  fcase(
    label == "positive", 1,
    label == "weak_positive", 0.5,
    label == "negative", -1,
    label == "confounded", -0.25,
    label == "measurement_problem", -0.5,
    label == "negative_control_pass", -1,
    label == "negative_control_fail", 1,
    default = NA_real_
  )
}

label_weight <- function(label) {
  fcase(
    label %chin% c("positive", "negative"), 1,
    label == "weak_positive", 0.5,
    label %chin% c("negative_control_pass", "negative_control_fail"), 0.8,
    label == "confounded", 0.25,
    label == "measurement_problem", 0.2,
    default = 0
  )
}

browser_status_from_label <- function(label) {
  fcase(
    label == "positive", "browser_supported_routing_change",
    label == "weak_positive", "browser_weak_support",
    label == "negative", "browser_no_visible_signal",
    label == "confounded", "browser_signal_confounded",
    label == "measurement_problem", "measurement_review_issue",
    label == "negative_control_pass", "negative_control_no_signal",
    label == "negative_control_fail", "negative_control_signal_detected",
    label == "not_reviewable", "not_reviewable",
    default = ""
  )
}

first_nonempty <- function(...) {
  values <- list(...)
  out <- rep("", length(values[[1]]))
  for (value in values) {
    value <- clean_text(value)
    take <- !has_text(out) & has_text(value)
    out[take] <- value[take]
  }
  out
}

template <- read_required(template_file, "linked active-review template")

manual_cols <- c(
  "manual_review_status",
  "browser_signal_grade",
  "branch_direction_verified",
  "clean_cds_direction_verified",
  "context_match_verified",
  "translon_structure_status",
  "metadata_confounder_status",
  "negative_control_result",
  "reviewer_notes",
  "reviewer",
  "review_date",
  "next_action"
)
review_context_cols <- c(
  "target_label_class", "label_gap_class", "label_acquisition_score",
  "protocol_pooling_need", "pooling_proof_tier", "pooling_hazard_class",
  "pooling_hazard_score", "pooling_hazard_sources"
)
ensure_columns(template, manual_cols)
for (column in manual_cols) {
  set(template, j = column, value = clean_text(template[[column]]))
}
ensure_columns(template, review_context_cols)

drafts <- safe_fread(draft_file)
if (nrow(drafts)) {
  ensure_columns(drafts, c("global_review_order", manual_cols,
                           review_context_cols))
  for (column in manual_cols) {
    set(drafts, j = column, value = clean_text(drafts[[column]]))
  }
  for (column in review_context_cols) {
    set(drafts, j = column, value = clean_text(drafts[[column]]))
  }
  drafts[, global_review_order :=
           suppressWarnings(as.integer(global_review_order))]
  drafts <- drafts[is.finite(global_review_order)]
  setorder(drafts, global_review_order)
  drafts <- drafts[!duplicated(global_review_order, fromLast = TRUE)]
  setkey(template, global_review_order)
  setkey(drafts, global_review_order)
  for (column in manual_cols) {
    draft_values <- drafts[template, get(column)]
    take <- has_text(draft_values)
    if (any(take)) {
      set(template, i = which(take), j = column, value = draft_values[take])
    }
  }
  for (column in review_context_cols) {
    draft_values <- drafts[template, get(column)]
    take <- has_text(draft_values)
    if (any(take)) {
      set(template, i = which(take), j = column, value = draft_values[take])
    }
  }
  setkey(template, NULL)
}

feedback_matrix <- as.data.frame(lapply(template[, ..manual_cols], has_text))
template[, any_manual_feedback := rowSums(feedback_matrix) > 0]
template[, branch_direction_norm :=
           direction_status(branch_direction_verified)]
template[, clean_cds_direction_norm :=
           direction_status(clean_cds_direction_verified)]
template[, context_match_norm :=
           direction_status(context_match_verified)]
template[, branch_direction_score :=
           direction_score(branch_direction_norm)]
template[, clean_cds_direction_score :=
           direction_score(clean_cds_direction_norm)]
template[, context_match_score :=
           direction_score(context_match_norm)]
template[, manual_status_norm := norm_text(manual_review_status)]
template[, browser_signal_norm := norm_text(browser_signal_grade)]
template[, translon_status_norm := norm_text(translon_structure_status)]
template[, metadata_confounder_norm := norm_text(metadata_confounder_status)]
template[, negative_control_norm := norm_text(negative_control_result)]

template[, review_label_class := mapply(
  label_from_fields,
  manual_review_status,
  browser_signal_grade,
  branch_direction_verified,
  clean_cds_direction_verified,
  context_match_verified,
  translon_structure_status,
  metadata_confounder_status,
  negative_control_result,
  batch_role,
  USE.NAMES = FALSE
)]
template[template[["any_manual_feedback"]] == FALSE, review_label_class := ""]
template[, review_label_score := label_score(review_label_class)]
template[, calibration_weight := label_weight(review_label_class)]

template[, usable_for_motif_calibration :=
           review_label_class %chin% c("positive", "weak_positive",
                                       "negative") &
           !translon_status_norm %chin% c("measurement_problem",
                                          "wrong_isoform",
                                          "wrong_transcript")]
template[, usable_for_context_filtering :=
           review_label_class == "confounded" |
           context_match_norm %chin% c("opposite", "unclear") |
           nzchar(metadata_confounder_norm)]
template[, usable_for_counterfactual_calibration :=
           review_label_class %chin% c("positive", "weak_positive",
                                       "negative") &
           batch_role %chin% c("candidate_positive",
                               "counterfactual_review")]
template[, usable_for_negative_control_calibration :=
           review_label_class %chin% c("negative_control_pass",
                                       "negative_control_fail") |
           (batch_role == "negative_control" &
              review_label_class %chin% c("positive", "negative"))]
template[, usable_for_training :=
           calibration_weight > 0 & is.finite(review_label_score) &
           review_label_class != "not_reviewable"]

id_cols <- c(
  "global_review_order", "batch_id", "batch_name", "batch_review_order",
  "batch_role", "gene_symbol", "tx_id", "design_family",
  "review_packet_id", "observatory_url", "observatory_link_mode",
  "observatory_link_note", "action_label", "claim_stage",
  "review_reason", "review_priority_tier", review_context_cols
)
norm_cols <- c(
  "manual_status_norm", "browser_signal_norm", "branch_direction_norm",
  "clean_cds_direction_norm", "context_match_norm",
  "translon_status_norm", "metadata_confounder_norm",
  "negative_control_norm", "branch_direction_score",
  "clean_cds_direction_score", "context_match_score",
  "review_label_class", "review_label_score", "calibration_weight",
  "usable_for_motif_calibration", "usable_for_context_filtering",
  "usable_for_counterfactual_calibration",
  "usable_for_negative_control_calibration", "usable_for_training"
)
label_cols <- unique(c(id_cols, manual_cols, norm_cols))
ensure_columns(template, setdiff(label_cols, names(template)))

labels <- template[template[["any_manual_feedback"]] == TRUE, ..label_cols]
data.table::fwrite(labels, labels_out)

training_covariates <- c(
  "consensus_rank", "consensus_class", "consensus_score",
  "atlas_evidence_class", "atlas_clean_cds_pct_change",
  "main_branch_class", "main_branch_delta",
  "top_shifted_feature_id", "top_shifted_feature_type",
  "top_shifted_feature_family", "model_agreement_tier",
  "loo_consensus_class", "transport_consensus_class",
  "context_interaction_gene_class", "partial_pooling_gene_class",
  "aux_top_state", "aux_top_overlap_class",
  "therapeutic_gene_label", "perturbation_label",
  "perturbation_support_class", "half_normalization_cds_pct_change",
  "n_studies", "strict_supported_studies",
  "total_case_samples", "total_control_samples",
  "best_case_condition", "best_control_conditions", "best_study",
  "branch_context_label", "study", "AUTHOR", "case_condition",
  "control_conditions", "TISSUE", "CELL_LINE", "GENE", "INHIBITOR",
  "FRACTION", "Cancer_type", "Sex", "TIMEPOINT", "match_scope",
  "n_case_runs_merged", "n_control_runs_merged",
  "case_branch_total_counts", "control_branch_total_counts",
  "context_count_gate", "replicate_gate", "joint_review_class",
  "joint_review_priority", "joint_top_branch_class",
  "joint_top_branch_delta", "joint_l1_allocation_shift",
  "joint_delta_clean_CDS", "joint_delta_leader_uORF",
  "joint_delta_overlapping_uORF", "transcript_qc_class",
  "cds_exon_qc_status", "n_measured_features",
  "n_signal_supported_features", "n_low_count_features",
  "n_manual_features", "n_overlapping_features", "n_internal_features",
  "warning_item_count", "warning_items", "warnings",
  "observatory_gene_label", "observatory_case_n",
  "observatory_control_n", "observatory_context_label",
  "observatory_case_condition", "observatory_control_conditions"
)
training_cols <- unique(c(id_cols, norm_cols, intersect(training_covariates,
                                                        names(template))))
training <- template[usable_for_training == TRUE, ..training_cols]
data.table::fwrite(training, training_out)

browser_candidate_cols <- c(
  "gene_symbol", "tx_id", "design_family", "validation_context",
  "browser_validation_status", "browser_validation_label",
  "validated_by", "validation_date", "validation_notes",
  "validation_score", "source_review_file", "global_review_order",
  "batch_id", "review_label_class", "calibration_weight"
)
if (nrow(labels)) {
  browser_notes <- copy(template[template[["any_manual_feedback"]] == TRUE])
  browser_notes[, validation_context := first_nonempty(
    observatory_context_label,
    branch_context_label,
    review_packet_id
  )]
  browser_notes[, browser_validation_status :=
                  browser_status_from_label(review_label_class)]
  browser_notes[, browser_validation_label := fifelse(
    has_text(browser_signal_grade),
    clean_text(browser_signal_grade),
    review_label_class
  )]
  browser_notes[, validated_by := fifelse(
    has_text(reviewer), clean_text(reviewer), "manual_review_template"
  )]
  browser_notes[, validation_date := clean_text(review_date)]
  browser_notes[, validation_notes := clean_text(paste(
    "Active review", global_review_order,
    paste0("(", batch_role, ")"),
    reviewer_notes,
    next_action,
    observatory_link_mode,
    sep = "; "
  ))]
  browser_notes[, validation_score := review_label_score]
  browser_notes[, source_review_file := template_file]
  browser_notes <- browser_notes[, ..browser_candidate_cols]
} else {
  browser_notes <- as.data.table(setNames(
    replicate(length(browser_candidate_cols), character(0),
              simplify = FALSE),
    browser_candidate_cols
  ))
}
data.table::fwrite(browser_notes, browser_notes_out)

summary_dt <- rbindlist(list(
  data.table(metric = "review_template_rows", value = nrow(template)),
  data.table(metric = "labeled_review_rows", value = nrow(labels)),
  data.table(metric = "training_rows", value = nrow(training)),
  data.table(metric = "browser_validation_note_candidate_rows",
             value = nrow(browser_notes)),
  data.table(metric = "positive_rows",
             value = sum(template$review_label_class == "positive")),
  data.table(metric = "weak_positive_rows",
             value = sum(template$review_label_class == "weak_positive")),
  data.table(metric = "negative_rows",
             value = sum(template$review_label_class == "negative")),
  data.table(metric = "confounded_rows",
             value = sum(template$review_label_class == "confounded")),
  data.table(metric = "measurement_problem_rows",
             value = sum(template$review_label_class == "measurement_problem")),
  data.table(metric = "negative_control_pass_rows",
             value = sum(template$review_label_class == "negative_control_pass")),
  data.table(metric = "negative_control_fail_rows",
             value = sum(template$review_label_class == "negative_control_fail")),
  data.table(metric = "motif_calibration_rows",
             value = sum(template$usable_for_motif_calibration)),
  data.table(metric = "context_filtering_rows",
             value = sum(template$usable_for_context_filtering)),
  data.table(metric = "counterfactual_calibration_rows",
             value = sum(template$usable_for_counterfactual_calibration)),
  data.table(metric = "negative_control_calibration_rows",
             value = sum(template$usable_for_negative_control_calibration))
), fill = TRUE)
data.table::fwrite(summary_dt, summary_out)

batch_id_cols <- c("batch_id", "batch_name", "batch_role")
batch_summary <- template[, .(
  review_rows = .N,
  labeled_rows = sum(any_manual_feedback),
  training_rows = sum(usable_for_training),
  positive_rows = sum(review_label_class == "positive"),
  weak_positive_rows = sum(review_label_class == "weak_positive"),
  negative_rows = sum(review_label_class == "negative"),
  confounded_rows = sum(review_label_class == "confounded"),
  measurement_problem_rows = sum(review_label_class == "measurement_problem"),
  negative_control_pass_rows =
    sum(review_label_class == "negative_control_pass"),
  negative_control_fail_rows =
    sum(review_label_class == "negative_control_fail")
), by = batch_id_cols]
setorder(batch_summary, batch_id)
data.table::fwrite(batch_summary, batch_summary_out)

report <- c(
  "# RDG Review Feedback Ingest",
  "",
  paste0("Input template: `", template_file, "`"),
  paste0("Draft overlay: `", draft_file, "`"),
  paste0("Draft rows loaded: ", nrow(drafts)),
  "",
  "## Summary",
  "",
  paste0("- Review rows: ", nrow(template)),
  paste0("- Rows with manual feedback: ", nrow(labels)),
  paste0("- Calibration/training rows: ", nrow(training)),
  paste0("- Candidate browser-validation note rows: ", nrow(browser_notes)),
  "",
  "## Outputs",
  "",
  paste0("- `", labels_out, "`"),
  paste0("- `", training_out, "`"),
  paste0("- `", browser_notes_out, "`"),
  paste0("- `", summary_out, "`"),
  paste0("- `", batch_summary_out, "`"),
  "",
  "## Use",
  "",
  paste0(
    "Fill manual columns in `",
    file.path("dominant_rdg_review_links",
              "dominant_rdg_review_batch_template_with_links.csv"),
    "` or save draft labels from the atlas browser, rerun ",
    "`rdg_review_feedback`, then inspect the candidate browser ",
    "notes before manually merging trusted rows into ",
    "`dominant_rdg_atlas/browser_validation_notes.csv`."
  )
)
if (!nrow(labels)) {
  report <- c(
    report,
    "",
    "No manual feedback rows were detected in the linked review template yet."
  )
}
writeLines(report, report_out)

message("RDG review feedback ingest complete: ", nrow(labels),
        " labeled rows, ", nrow(training), " training rows.")
