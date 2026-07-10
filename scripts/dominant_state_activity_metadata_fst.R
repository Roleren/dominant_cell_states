#!/usr/bin/env Rscript

# Browser-ready dominant-state metadata layer.
# The primary score is an empirical within-state activity rank in [0, 1]:
# 0 = lowest scoring run for that state, 1 = highest scoring run for that state.

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
score_dir <- file.path(analysis_dir, "dominant_state_score_audit")
score_file <- file.path(score_dir, "dominant_state_score_audit_sample_scores_long.csv")
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
module_definition_file <- file.path(
  analysis_dir,
  "scripts",
  "dominant_state_module_definitions.R"
)

out_fst <- file.path(analysis_dir, "human_dominant_state_activity_metadata.fst")
out_long <- file.path(analysis_dir, "human_dominant_state_activity_metadata_long.csv")
out_manifest <- file.path(analysis_dir, "human_dominant_state_activity_metadata_manifest.csv")

if (!requireNamespace("fst", quietly = TRUE)) {
  stop("Package 'fst' is required to write: ", out_fst, call. = FALSE)
}
if (!file.exists(score_file)) {
  stop("Missing dominant-state audit score file: ", score_file, call. = FALSE)
}
if (!file.exists(state_file)) {
  stop("Missing clean-CDS dominant-state file: ", state_file, call. = FALSE)
}
if (!file.exists(module_definition_file)) {
  stop("Missing dominant-state module definitions: ", module_definition_file,
       call. = FALSE)
}
source(module_definition_file)

modules <- load_dominant_state_modules(analysis_dir)
baseline_state <- "Baseline (control)"
state_order <- modules[dominant_state != baseline_state, dominant_state]

state_dt <- fread(
  state_file,
  select = c("Run", "sample", "dominant_state", "dominant_state_score"),
  showProgress = FALSE
)
required_state_cols <- c("Run", "sample", "dominant_state", "dominant_state_score")
if (!all(required_state_cols %chin% names(state_dt))) {
  stop("State table is missing columns: ",
       paste(setdiff(required_state_cols, names(state_dt)), collapse = ", "),
       call. = FALSE)
}
if (anyDuplicated(state_dt$Run)) {
  stop("State table has duplicate Run values.", call. = FALSE)
}

scores <- fread(score_file, showProgress = FALSE)
if ("sample" %chin% names(scores)) {
  scores[, sample := NULL]
}
required_score_cols <- c(
  "Run", "dominant_state", "signed_score", "robust_z",
  "up_genes_present", "up_genes_total",
  "down_genes_present", "down_genes_total"
)
if (!all(required_score_cols %chin% names(scores))) {
  stop("Score audit table is missing columns: ",
       paste(setdiff(required_score_cols, names(scores)), collapse = ", "),
       call. = FALSE)
}
scores <- scores[dominant_state %chin% state_order]

missing_runs <- setdiff(state_dt$Run, scores$Run)
if (length(missing_runs)) {
  stop("Score audit table is missing runs from the state table: ",
       paste(head(missing_runs, 8), collapse = ", "), call. = FALSE)
}

activity_rank_0_1 <- function(x) {
  out <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  n <- sum(valid)
  if (n == 0) return(out)
  if (n == 1) {
    out[valid] <- 1
    return(out)
  }
  out[valid] <- (frank(x[valid], ties.method = "average") - 1) / (n - 1)
  out
}

minmax_0_1 <- function(x) {
  out <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  if (!any(valid)) return(out)
  limits <- range(x[valid], na.rm = TRUE)
  span <- diff(limits)
  if (!is.finite(span) || span <= 0) {
    out[valid] <- 1
    return(out)
  }
  out[valid] <- (x[valid] - limits[[1]]) / span
  out
}

safe_state_column <- function(state) {
  x <- iconv(state, from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  paste0("DCS_", x)
}

column_map <- data.table(
  dominant_state = state_order,
  activity_column = make.unique(safe_state_column(state_order), sep = "_")
)
scores <- merge(scores, column_map, by = "dominant_state", all.x = TRUE,
                sort = FALSE)

scores[
  ,
  `:=`(
    activity_score = activity_rank_0_1(signed_score),
    minmax_score = minmax_0_1(signed_score)
  ),
  by = dominant_state
]

scores <- merge(
  state_dt[, .(Run, run_order = .I, sample)],
  scores,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
scores[, state_order_index := match(dominant_state, state_order)]
setorder(scores, run_order, state_order_index)

activity_wide <- dcast(
  scores[, .(Run, activity_column, activity_score)],
  Run ~ activity_column,
  value.var = "activity_score"
)
setcolorder(activity_wide, c("Run", column_map$activity_column))

top_state <- scores[
  is.finite(activity_score),
  .SD[which.max(activity_score)],
  by = Run
][
  ,
  .(
    Run,
    dominant_state_top = dominant_state,
    dominant_state_top_activity = activity_score,
    dominant_state_top_raw_score = signed_score
  )
]

metadata <- merge(
  state_dt[, .(Run, run_order = .I)],
  activity_wide,
  by = "Run",
  all.x = TRUE,
  sort = FALSE
)
metadata <- merge(metadata, top_state, by = "Run", all.x = TRUE, sort = FALSE)
setorder(metadata, run_order)
metadata[, run_order := NULL]
setcolorder(
  metadata,
  c(
    "Run",
    column_map$activity_column,
    "dominant_state_top",
    "dominant_state_top_activity",
    "dominant_state_top_raw_score"
  )
)

round_double_columns <- function(dt, digits = 2) {
  double_cols <- names(dt)[vapply(dt, is.double, logical(1))]
  if (length(double_cols)) {
    dt[, (double_cols) := lapply(.SD, round, digits = digits),
       .SDcols = double_cols]
  }
  invisible(dt)
}

manifest <- scores[
  ,
  .(
    n_runs = uniqueN(Run),
    n_finite_scores = sum(is.finite(signed_score)),
    n_missing_scores = sum(!is.finite(signed_score)),
    raw_score_min = suppressWarnings(min(signed_score, na.rm = TRUE)),
    raw_score_median = suppressWarnings(median(signed_score, na.rm = TRUE)),
    raw_score_max = suppressWarnings(max(signed_score, na.rm = TRUE)),
    robust_z_min = suppressWarnings(min(robust_z, na.rm = TRUE)),
    robust_z_median = suppressWarnings(median(robust_z, na.rm = TRUE)),
    robust_z_max = suppressWarnings(max(robust_z, na.rm = TRUE)),
    up_genes_present = suppressWarnings(max(up_genes_present, na.rm = TRUE)),
    up_genes_total = suppressWarnings(max(up_genes_total, na.rm = TRUE)),
    down_genes_present = suppressWarnings(max(down_genes_present, na.rm = TRUE)),
    down_genes_total = suppressWarnings(max(down_genes_total, na.rm = TRUE))
  ),
  by = .(dominant_state, activity_column)
]
manifest[!is.finite(raw_score_min), raw_score_min := NA_real_]
manifest[!is.finite(raw_score_median), raw_score_median := NA_real_]
manifest[!is.finite(raw_score_max), raw_score_max := NA_real_]
manifest[!is.finite(robust_z_min), robust_z_min := NA_real_]
manifest[!is.finite(robust_z_median), robust_z_median := NA_real_]
manifest[!is.finite(robust_z_max), robust_z_max := NA_real_]
manifest[
  ,
  `:=`(
    score_method = paste(
      "clean-CDS signed module score; activity column is empirical",
      "within-state rank scaled to [0,1]"
    ),
    score_direction = "higher means more active relative to other runs",
    baseline_handling = paste(
      "Baseline (control) genes are used upstream for normalization and are",
      "not exported as a dominant-state activity column"
    )
  )
]
manifest[, state_order_index := match(dominant_state, state_order)]
setorder(manifest, state_order_index)
manifest[, state_order_index := NULL]

long_cols <- c(
  "Run", "sample", "dominant_state", "activity_column", "activity_score",
  "minmax_score", "signed_score", "robust_z", "up_genes_present",
  "up_genes_total", "down_genes_present", "down_genes_total"
)
long <- scores[order(run_order, state_order_index), ..long_cols]

round_double_columns(metadata, digits = 2)
round_double_columns(long, digits = 2)
round_double_columns(manifest, digits = 2)

fst::write_fst(metadata, out_fst, compress = 50)
fwrite(long, out_long)
fwrite(manifest, out_manifest)

read_back <- fst::read_fst(out_fst, as.data.table = TRUE)
if (!identical(read_back$Run, state_dt$Run)) {
  stop("FST read-back Run order does not match the clean-CDS state table.",
       call. = FALSE)
}
if (!identical(names(read_back), names(metadata))) {
  stop("FST read-back columns do not match the written metadata table.",
       call. = FALSE)
}

activity_matrix <- as.matrix(
  read_back[, column_map$activity_column, with = FALSE]
)
if (any(activity_matrix < 0 | activity_matrix > 1, na.rm = TRUE)) {
  stop("Activity scores outside [0,1] after FST read-back.", call. = FALSE)
}

message("Dominant-state activity metadata:")
message("  runs: ", nrow(metadata))
message("  activity states: ", length(column_map$activity_column))
message("  fst: ", out_fst)
message("  manifest: ", out_manifest)
