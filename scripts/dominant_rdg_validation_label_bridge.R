#!/usr/bin/env Rscript

# Combine curated browser validation notes with active-review feedback labels.
# This is a readback/triage layer; it does not overwrite curated inputs.

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
results_dir <- Sys.getenv(
  "DOMINANT_CELL_STATES_RESULTS_DIR",
  unset = file.path(analysis_dir, "results")
)
output_dir <- file.path(analysis_dir, "dominant_rdg_validation_label_bridge")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) paths[1] else hit
}

browser_notes_file <- first_existing_file(c(
  file.path(analysis_dir, "dominant_rdg_atlas", "browser_validation_notes.csv"),
  file.path(results_dir, "curated_inputs", "browser_validation_notes.csv")
))
feedback_labels_file <- file.path(
  analysis_dir, "dominant_rdg_review_feedback",
  "dominant_rdg_review_feedback_labels.csv"
)
motif_file <- file.path(
  analysis_dir, "dominant_rdg_motif_prior_calibration",
  "dominant_rdg_motif_prior_gene_design_calibration.csv"
)
review_template_file <- file.path(
  analysis_dir, "dominant_rdg_review_links",
  "dominant_rdg_review_batch_template_with_links.csv"
)

label_set_out <- file.path(
  output_dir, "dominant_rdg_validation_label_set.csv"
)
gene_design_out <- file.path(
  output_dir, "dominant_rdg_validation_label_gene_design.csv"
)
motif_overlay_out <- file.path(
  output_dir, "dominant_rdg_validation_label_motif_overlay.csv"
)
gap_queue_out <- file.path(
  output_dir, "dominant_rdg_validation_label_gap_queue.csv"
)
balance_out <- file.path(
  output_dir, "dominant_rdg_validation_label_balance.csv"
)
summary_out <- file.path(
  output_dir, "dominant_rdg_validation_label_bridge_summary.csv"
)
report_out <- file.path(
  output_dir, "dominant_rdg_validation_label_bridge_report.md"
)

safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(data.table())
  data.table::fread(path, showProgress = FALSE, nThread = 1, ...)
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

num_col <- function(dt, column, default = NA_real_) {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- suppressWarnings(as.numeric(dt[[column]]))
  x[!is.finite(x)] <- default
  x
}

chr_col <- function(dt, column, default = "") {
  if (!column %in% names(dt)) return(rep(default, nrow(dt)))
  x <- clean_text(dt[[column]])
  x[!has_text(x)] <- default
  x
}

ensure_columns <- function(dt, columns, value = "") {
  for (column in columns) {
    if (!column %in% names(dt)) dt[, (column) := value]
  }
  invisible(dt)
}

clip01 <- function(x) pmax(0, pmin(1, x))

scale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  out <- rep(0, length(x))
  if (!any(ok)) return(out)
  rng <- range(x[ok])
  if (diff(rng) < 1e-12) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (x[ok] - rng[1]) / diff(rng)
  }
  out
}

weighted_mean_safe <- function(x, w) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

collapse_unique <- function(x) {
  x <- unique(clean_text(x))
  x <- x[has_text(x)]
  if (!length(x)) "" else paste(x, collapse = "; ")
}

browser_label_class <- function(status, label, notes, score) {
  text <- norm_text(paste(status, label, notes))
  score <- suppressWarnings(as.numeric(score))
  fcase(
    grepl("negative|no_signal|absent|unsupported|not_supported|reject|refute",
          text) |
      (is.finite(score) & score < 0),
    "negative",
    grepl("weak|partial|suggestive|ambiguous", text) |
      (is.finite(score) & score > 0 & score < 0.75),
    "weak_positive",
    grepl("support|supported|confirm|confirmed|valid|validated|interesting|positive",
          text) |
      (is.finite(score) & score >= 0.75),
    "positive",
    default = ""
  )
}

score_from_label <- function(label, fallback = NA_real_) {
  fallback <- suppressWarnings(as.numeric(fallback))
  fcase(
    label == "positive", fifelse(is.finite(fallback),
                                  clip01(fallback), 1),
    label == "weak_positive", fifelse(is.finite(fallback),
                                       pmax(0.25, pmin(0.75, fallback)),
                                       0.5),
    label %chin% c("negative", "negative_control_pass"), -1,
    label == "negative_control_fail", 1,
    label == "confounded", -0.25,
    label == "measurement_problem", -0.5,
    default = NA_real_
  )
}

empty_labels <- function() {
  data.table(
    label_source = character(),
    label_id = character(),
    gene_symbol = character(),
    tx_id = character(),
    design_family = character(),
    validation_context = character(),
    label_class = character(),
    label_score = numeric(),
    label_weight = numeric(),
    label_status = character(),
    label_text = character(),
    reviewer = character(),
    review_date = character(),
    review_notes = character(),
    source_file = character(),
    global_review_order = integer(),
    batch_id = character(),
    batch_role = character(),
    observatory_url = character()
  )
}

browser_notes <- safe_fread(browser_notes_file)
feedback_labels <- safe_fread(feedback_labels_file)
motif <- safe_fread(motif_file)
review_template <- safe_fread(review_template_file)

curated_labels <- empty_labels()
if (nrow(browser_notes)) {
  ensure_columns(
    browser_notes,
    c("gene_symbol", "tx_id", "design_family", "validation_context",
      "browser_validation_status", "browser_validation_label",
      "validated_by", "validation_date", "validation_notes",
      "validation_score")
  )
  browser_notes[, validation_score_num := num_col(.SD, "validation_score")]
  browser_notes[, label_class := browser_label_class(
    browser_validation_status,
    browser_validation_label,
    validation_notes,
    validation_score_num
  )]
  curated_labels <- browser_notes[has_text(label_class), .(
    label_source = "curated_browser_validation",
    label_id = "",
    gene_symbol = clean_text(gene_symbol),
    tx_id = clean_text(tx_id),
    design_family = clean_text(design_family),
    validation_context = clean_text(validation_context),
    label_class = clean_text(label_class),
    label_score = score_from_label(label_class, validation_score_num),
    label_weight = pmax(0.25, pmin(1, abs(score_from_label(
      label_class, validation_score_num
    )))),
    label_status = clean_text(browser_validation_status),
    label_text = clean_text(browser_validation_label),
    reviewer = clean_text(validated_by),
    review_date = clean_text(validation_date),
    review_notes = clean_text(validation_notes),
    source_file = browser_notes_file,
    global_review_order = NA_integer_,
    batch_id = "",
    batch_role = "",
    observatory_url = ""
  )]
}

active_labels <- empty_labels()
if (nrow(feedback_labels)) {
  ensure_columns(
    feedback_labels,
    c("gene_symbol", "tx_id", "design_family", "review_packet_id",
      "review_label_class", "review_label_score", "calibration_weight",
      "manual_review_status", "browser_signal_grade", "reviewer",
      "review_date", "reviewer_notes", "global_review_order",
      "batch_id", "batch_role", "observatory_url")
  )
  feedback_labels[, review_label_score_num :=
                    num_col(.SD, "review_label_score")]
  feedback_labels[, calibration_weight_num :=
                    num_col(.SD, "calibration_weight", default = 0)]
  active_labels <- feedback_labels[has_text(review_label_class), .(
    label_source = "active_review_feedback",
    label_id = "",
    gene_symbol = clean_text(gene_symbol),
    tx_id = clean_text(tx_id),
    design_family = clean_text(design_family),
    validation_context = clean_text(review_packet_id),
    label_class = clean_text(review_label_class),
    label_score = fifelse(
      is.finite(review_label_score_num),
      review_label_score_num,
      score_from_label(review_label_class)
    ),
    label_weight = fifelse(
      is.finite(calibration_weight_num) & calibration_weight_num > 0,
      calibration_weight_num,
      pmax(0.25, pmin(1, abs(score_from_label(review_label_class))))
    ),
    label_status = clean_text(manual_review_status),
    label_text = clean_text(browser_signal_grade),
    reviewer = clean_text(reviewer),
    review_date = clean_text(review_date),
    review_notes = clean_text(reviewer_notes),
    source_file = feedback_labels_file,
    global_review_order = suppressWarnings(as.integer(global_review_order)),
    batch_id = clean_text(batch_id),
    batch_role = clean_text(batch_role),
    observatory_url = clean_text(observatory_url)
  )]
}

labels <- rbindlist(list(curated_labels, active_labels), fill = TRUE)
if (nrow(labels)) {
  labels[, label_id := sprintf("VL%04d", .I)]
  setorder(labels, gene_symbol, tx_id, design_family, label_source,
           global_review_order)
}
data.table::fwrite(labels, label_set_out)

key_cols <- c("gene_symbol", "tx_id", "design_family")
label_counts <- function(x, value) sum(x == value, na.rm = TRUE)

if (nrow(labels)) {
  gene_design <- labels[, .(
    label_rows = .N,
    curated_browser_label_count =
      sum(label_source == "curated_browser_validation"),
    feedback_label_count = sum(label_source == "active_review_feedback"),
    positive_label_count = label_counts(label_class, "positive"),
    weak_positive_label_count = label_counts(label_class, "weak_positive"),
    negative_label_count = label_counts(label_class, "negative"),
    confounded_label_count = label_counts(label_class, "confounded"),
    measurement_problem_label_count =
      label_counts(label_class, "measurement_problem"),
    negative_control_pass_label_count =
      label_counts(label_class, "negative_control_pass"),
    negative_control_fail_label_count =
      label_counts(label_class, "negative_control_fail"),
    weighted_label_score = weighted_mean_safe(label_score, label_weight),
    max_label_weight = max(label_weight, na.rm = TRUE),
    label_classes = collapse_unique(label_class),
    label_sources = collapse_unique(label_source),
    validation_contexts = collapse_unique(validation_context),
    reviewers = collapse_unique(reviewer),
    latest_review_date = max(clean_text(review_date))
  ), by = key_cols]
  gene_design[, validation_label_class := fcase(
    positive_label_count > 0 &
      negative_label_count + negative_control_pass_label_count == 0,
    "positive_seed",
    positive_label_count == 0 &
      negative_label_count + negative_control_pass_label_count > 0,
    "negative_seed",
    positive_label_count > 0 &
      negative_label_count + negative_control_pass_label_count > 0,
    "mixed_positive_negative",
    confounded_label_count > 0,
    "confounded_seed",
    measurement_problem_label_count > 0,
    "measurement_problem_seed",
    default = "labeled_unclear"
  )]
} else {
  gene_design <- data.table(
    gene_symbol = character(),
    tx_id = character(),
    design_family = character(),
    label_rows = integer(),
    validation_label_class = character()
  )
}

motif_small_cols <- intersect(c(
  key_cols,
  "motif_prior_review_class", "motif_prior_review_priority",
  "motif_top_concordance_class", "motif_top_prior_strength",
  "motif_top_delta", "motif_net_alignment_score",
  "atlas_evidence_class", "atlas_browser_review_priority_score",
  "atlas_clean_cds_pct_change", "browser_review_rank",
  "browser_review_question"
), names(motif))
if (length(motif_small_cols)) {
  gene_design <- merge(
    gene_design,
    unique(motif[, ..motif_small_cols], by = key_cols),
    by = key_cols,
    all.x = TRUE
  )
}
setorder(gene_design, -label_rows, gene_symbol, design_family)
data.table::fwrite(gene_design, gene_design_out)

if (nrow(motif)) {
  overlay <- copy(motif)
  overlay <- merge(overlay, gene_design, by = key_cols, all.x = TRUE,
                   suffixes = c("", "_label"))
  count_cols <- c(
    "label_rows", "curated_browser_label_count", "feedback_label_count",
    "positive_label_count", "weak_positive_label_count",
    "negative_label_count", "confounded_label_count",
    "measurement_problem_label_count",
    "negative_control_pass_label_count",
    "negative_control_fail_label_count"
  )
  for (column in intersect(count_cols, names(overlay))) {
    overlay[is.na(get(column)), (column) := 0L]
  }
  ensure_columns(overlay, c("validation_label_class"))
  overlay[!has_text(validation_label_class),
          validation_label_class := "unlabeled"]
  overlay[, label_calibration_role := fcase(
    validation_label_class == "positive_seed", "positive_calibrator",
    validation_label_class == "negative_seed", "negative_calibrator",
    validation_label_class == "mixed_positive_negative",
    "conflict_or_threshold_calibrator",
    validation_label_class == "confounded_seed", "context_filter_calibrator",
    validation_label_class == "measurement_problem_seed",
    "measurement_filter_calibrator",
    default = "unlabeled_candidate"
  )]
  overlay[, motif_label_agreement := fcase(
    validation_label_class == "positive_seed" &
      grepl("concordant", motif_prior_review_class),
    "positive_label_motif_concordant",
    validation_label_class == "positive_seed" &
      grepl("discordant", motif_prior_review_class),
    "positive_label_motif_discordant",
    validation_label_class == "negative_seed" &
      motif_prior_review_priority >= 0.5,
    "negative_label_high_motif_priority",
    validation_label_class == "negative_seed",
    "negative_label_low_motif_priority",
    validation_label_class == "unlabeled" &
      motif_prior_review_priority >= 0.5,
    "unlabeled_high_motif_priority",
    default = "not_assessed"
  )]
  setorder(overlay, -label_rows, -motif_prior_review_priority)
} else {
  overlay <- data.table()
}
data.table::fwrite(overlay, motif_overlay_out)

target_balance <- data.table(
  label_class = c(
    "positive", "weak_positive", "negative", "confounded",
    "measurement_problem", "negative_control_pass",
    "negative_control_fail"
  ),
  target_minimum = c(8L, 5L, 12L, 8L, 5L, 10L, 2L),
  acquisition_role = c(
    "confirm candidate-positive browser signals",
    "capture borderline/suggestive browser signals",
    "falsify top-looking candidates",
    "teach context and metadata filters",
    "teach translon/measurement filters",
    "confirm negative-control absences",
    "catch unexpected negative-control failures"
  )
)
current_counts <- labels[, .N, by = label_class]
target_balance <- merge(target_balance, current_counts, by = "label_class",
                        all.x = TRUE)
target_balance[is.na(N), N := 0L]
setnames(target_balance, "N", "current_count")
target_balance[, label_deficit := pmax(target_minimum - current_count, 0L)]
target_balance[, deficit_fraction := fifelse(
  target_minimum > 0,
  label_deficit / target_minimum,
  0
)]
target_balance[, acquisition_priority := fcase(
  label_deficit == 0, "covered",
  current_count == 0, "missing",
  deficit_fraction >= 0.5, "thin",
  default = "needs_more"
)]
setorder(target_balance, -deficit_fraction, label_class)
data.table::fwrite(target_balance, balance_out)

if (nrow(review_template)) {
  review <- copy(review_template)
  ensure_columns(review, c(key_cols, "batch_role", "batch_id",
                          "batch_name", "review_question",
                          "review_packet_id", "observatory_url"))
  motif_nonkey_cols <- setdiff(motif_small_cols, key_cols)
  label_review_cols <- c(
    key_cols,
    setdiff(names(gene_design), c(key_cols, motif_nonkey_cols))
  )
  review <- merge(review, gene_design[, ..label_review_cols],
                  by = key_cols, all.x = TRUE)
  if (length(motif_small_cols)) {
    review <- merge(
      review,
      unique(motif[, ..motif_small_cols], by = key_cols),
      by = key_cols,
      all.x = TRUE,
      suffixes = c("", "_motif")
    )
  }
  for (column in intersect(count_cols, names(review))) {
    review[is.na(get(column)), (column) := 0L]
  }
  ensure_columns(review, c("validation_label_class"))
  review[!has_text(validation_label_class),
         validation_label_class := "unlabeled"]
  review[, already_labeled := label_rows > 0]
  review[, label_gap_class := fcase(
    batch_role == "negative_control", "need_negative_control_labels",
    batch_role == "context_audit", "need_context_filter_labels",
    batch_role == "measurement_review", "need_measurement_filter_labels",
    batch_role == "counterfactual_review", "need_counterfactual_labels",
    batch_role == "replication_followup", "need_replication_labels",
    default = "need_candidate_browser_labels"
  )]
  review[, target_label_class := fcase(
    batch_role == "negative_control", "negative_control_pass_or_fail",
    batch_role == "context_audit", "confounded_or_context_match",
    batch_role == "measurement_review", "measurement_problem_or_clean",
    batch_role == "counterfactual_review", "positive_or_negative",
    batch_role == "replication_followup", "positive_or_weak_or_negative",
    default = "positive_or_negative"
  )]
  score_base <- 0.35 * scale01(num_col(review, "acquisition_priority_score", 0)) +
    0.25 * scale01(num_col(review, "counterfactual_readiness_score", 0)) +
    0.20 * scale01(num_col(review, "biological_value_score", 0)) +
    0.20 * scale01(num_col(review, "motif_prior_review_priority", 0))
  class_weight <- fcase(
    review$batch_role == "negative_control", 1,
    review$batch_role == "context_audit", 0.9,
    review$batch_role == "measurement_review", 0.8,
    review$batch_role == "counterfactual_review", 0.75,
    review$batch_role == "candidate_positive", 0.7,
    default = 0.55
  )
  review[, label_acquisition_score := clip01(0.65 * score_base +
                                               0.35 * class_weight)]
  review[already_labeled == TRUE,
         label_acquisition_score := label_acquisition_score * 0.35]
  setorder(review, -label_acquisition_score, already_labeled,
           batch_id, global_review_order)
  review[, label_gap_rank := seq_len(.N)]
  gap_cols <- unique(c(
    "label_gap_rank", key_cols, "batch_id", "batch_name", "batch_role",
    "label_gap_class", "target_label_class", "label_acquisition_score",
    "already_labeled", "validation_label_class", "label_rows",
    "positive_label_count", "negative_label_count",
    "confounded_label_count", "measurement_problem_label_count",
    "negative_control_pass_label_count", "negative_control_fail_label_count",
    "motif_prior_review_class", "motif_prior_review_priority",
    "motif_top_concordance_class", "motif_top_prior_strength",
    "atlas_evidence_class", "atlas_clean_cds_pct_change",
    "consensus_class", "consensus_score", "review_question",
    "review_reason", "observatory_url", "observatory_link_mode",
    "observatory_link_note"
  ))
  gap_cols <- intersect(gap_cols, names(review))
  gap_queue <- review[, ..gap_cols]
} else {
  gap_queue <- data.table()
}
data.table::fwrite(gap_queue, gap_queue_out)

summary_dt <- rbindlist(list(
  data.table(metric = "curated_browser_label_rows",
             value = nrow(curated_labels)),
  data.table(metric = "active_feedback_label_rows",
             value = nrow(active_labels)),
  data.table(metric = "validation_label_rows", value = nrow(labels)),
  data.table(metric = "labeled_gene_designs", value = nrow(gene_design)),
  data.table(metric = "positive_label_rows",
             value = sum(labels$label_class == "positive")),
  data.table(metric = "weak_positive_label_rows",
             value = sum(labels$label_class == "weak_positive")),
  data.table(metric = "negative_label_rows",
             value = sum(labels$label_class == "negative")),
  data.table(metric = "confounded_label_rows",
             value = sum(labels$label_class == "confounded")),
  data.table(metric = "measurement_problem_label_rows",
             value = sum(labels$label_class == "measurement_problem")),
  data.table(metric = "negative_control_pass_rows",
             value = sum(labels$label_class == "negative_control_pass")),
  data.table(metric = "negative_control_fail_rows",
             value = sum(labels$label_class == "negative_control_fail")),
  data.table(metric = "motif_overlay_rows", value = nrow(overlay)),
  data.table(metric = "label_gap_queue_rows", value = nrow(gap_queue)),
  data.table(metric = "unlabeled_gap_queue_rows",
             value = if (nrow(gap_queue)) {
               sum(gap_queue$already_labeled == FALSE, na.rm = TRUE)
             } else {
               0L
             })
), fill = TRUE)
data.table::fwrite(summary_dt, summary_out)

top_gaps <- if (nrow(gap_queue)) {
  gap_queue[already_labeled == FALSE][
    seq_len(min(.N, 10)),
    paste0(
      "- ", label_gap_rank, ". ", gene_symbol, " | ", design_family,
      " | ", batch_role, " | score ",
      sprintf("%.3f", label_acquisition_score)
    )
  ]
} else {
  character()
}
if (!length(top_gaps)) top_gaps <- "- No unlabeled review rows available."

missing_classes <- target_balance[label_deficit > 0][
  ,
  paste0("- ", label_class, ": current ", current_count,
         ", target ", target_minimum,
         ", deficit ", label_deficit)
]
if (!length(missing_classes)) missing_classes <- "- All target label classes covered."

report <- c(
  "# RDG Validation Label Bridge",
  "",
  "## Summary",
  "",
  paste0("- Curated browser labels: ", nrow(curated_labels)),
  paste0("- Active-review feedback labels: ", nrow(active_labels)),
  paste0("- Unified validation labels: ", nrow(labels)),
  paste0("- Labeled gene-designs: ", nrow(gene_design)),
  paste0("- Motif overlay rows: ", nrow(overlay)),
  paste0("- Label-gap queue rows: ", nrow(gap_queue)),
  "",
  "## Missing Label Classes",
  "",
  missing_classes,
  "",
  "## Top Unlabeled Gap Rows",
  "",
  top_gaps,
  "",
  "## Outputs",
  "",
  paste0("- `", label_set_out, "`"),
  paste0("- `", gene_design_out, "`"),
  paste0("- `", motif_overlay_out, "`"),
  paste0("- `", gap_queue_out, "`"),
  paste0("- `", balance_out, "`"),
  paste0("- `", summary_out, "`")
)
writeLines(report, report_out)

message("RDG validation label bridge complete: ", nrow(labels),
        " labels over ", nrow(gene_design), " gene-designs.")
