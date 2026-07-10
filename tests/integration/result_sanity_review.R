find_analysis_dir_local <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Could not find dominant_cell_states from: ", start, call. = FALSE)
    }
    here <- parent
  }
}

suppressPackageStartupMessages(library(data.table))

analysis_dir <- find_analysis_dir_local()
repo_root <- analysis_dir
results_dir <- file.path(analysis_dir, "results")
generated_root_names <- c(
  "^dominant_",
  "^human_dominant_",
  "^required_fst_",
  "^postviral_fatigue_outputs$",
  "^htmlwidget_libs$",
  "^readlength_inputs$"
)
generated_manual_files <- c(
  "manual_translon_candidates.csv",
  "manual_translon_requested_gene_uoorf_candidates.csv",
  "manual_translons_manifest.csv",
  "manual_translons_ranges.rds"
)
is_generated_root_entry <- function(x) {
  is.character(x) &&
    length(x) >= 1L &&
    any(vapply(generated_root_names, grepl, logical(1), x = x[[1]]))
}
is_analysis_root_arg <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    identical(normalizePath(x, mustWork = FALSE), analysis_dir)
}
is_generated_manual_entry <- function(parts) {
  length(parts) >= 3L &&
    identical(parts[[2]], "manual_translons") &&
    is.character(parts[[3]]) &&
    parts[[3]][[1]] %in% generated_manual_files
}
file.path <- function(..., fsep = .Platform$file.sep) {
  parts <- list(...)
  if (length(parts) >= 2L &&
      is_analysis_root_arg(parts[[1]]) &&
      (is_generated_root_entry(parts[[2]]) ||
       is_generated_manual_entry(parts))) {
    parts[[1]] <- results_dir
  }
  do.call(base::file.path, c(parts, list(fsep = fsep)))
}
out_dir <- file.path(analysis_dir, "results", "test_reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_dt <- function(path) {
  tryCatch(
    fread(path, showProgress = FALSE, nThread = 1),
    error = function(e) e
  )
}

metric <- function(name, value, status = "info", note = "") {
  data.table(metric = name, value = as.character(value), status = status,
             note = note)
}

csv_files <- list.files(analysis_dir, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE)
csv_inventory <- rbindlist(lapply(csv_files, function(path) {
  x <- read_dt(path)
  if (inherits(x, "error")) {
    return(data.table(
      file = path, rows = NA_integer_, cols = NA_integer_,
      size_bytes = file.info(path)$size,
      status = "read_error", note = conditionMessage(x)
    ))
  }
  status <- if (nrow(x) == 0) "empty" else "ok"
  note <- if (status == "empty" &&
              basename(path) ==
              "dominant_rdg_hierarchical_counterfactual_predictions.csv") {
    "expected when every gene x design family cell is observed"
  } else if (status == "empty" &&
             basename(path) %chin% c(
               "dominant_rdg_review_feedback_labels.csv",
               "dominant_rdg_review_feedback_training_set.csv",
               "dominant_rdg_review_feedback_browser_validation_notes_candidate.csv"
             )) {
    "expected until the linked active-review template contains manual labels"
  } else {
    ""
  }
  data.table(
    file = path, rows = nrow(x), cols = ncol(x),
    size_bytes = file.info(path)$size,
    status = status, note = note
  )
}), fill = TRUE)

figure_files <- list.files(
  analysis_dir,
  pattern = "\\.(png|pdf|html)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
private_widget_dirs <- list.dirs(analysis_dir, recursive = TRUE,
                                 full.names = TRUE)
private_widget_dirs <- private_widget_dirs[
  grepl("_files$", basename(private_widget_dirs)) &
    basename(private_widget_dirs) != "htmlwidget_libs"
]
figure_inventory <- data.table(
  file = figure_files,
  ext = tolower(tools::file_ext(figure_files)),
  size_bytes = file.info(figure_files)$size
)
figure_inventory[, status := fcase(
  !file.exists(file), "missing",
  is.na(size_bytes) | size_bytes == 0, "empty",
  ext %chin% c("png", "pdf") & size_bytes < 2500, "suspiciously_small",
  default = "ok"
)]

key_checks <- list()
add_key <- function(name, value, status = "info", note = "") {
  key_checks[[length(key_checks) + 1L]] <<- metric(name, value, status, note)
}

cards_file <- file.path(analysis_dir, "dominant_rdg_atlas",
                        "dominant_rdg_atlas_cards.csv")
candidates_file <- file.path(analysis_dir, "dominant_rdg_validation_dossiers",
                             "dominant_rdg_validation_dossier_candidates.csv")
template_file <- file.path(analysis_dir, "dominant_rdg_validation_dossiers",
                           "dominant_rdg_validation_dossier_manual_review_template.csv")
pages_file <- file.path(analysis_dir, "required_fst_pages.csv")
clean_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
activity_file <- file.path(analysis_dir,
                           "human_dominant_state_activity_metadata_long.csv")

if (file.exists(cards_file)) {
  cards <- fread(cards_file, showProgress = FALSE)
  dup_cards <- cards[, .N, by = .(gene_symbol, tx_id, design_family)][N > 1, .N]
  add_key("atlas_cards_rows", nrow(cards), if (nrow(cards) > 0) "ok" else "fail")
  add_key("atlas_card_duplicate_gene_tx_design_rows", dup_cards,
          if (dup_cards == 0) "ok" else "fail")
  add_key("atlas_card_genes", uniqueN(cards$gene_symbol), "ok")
  add_key("atlas_card_design_families", uniqueN(cards$design_family), "ok")
}

if (file.exists(candidates_file)) {
  candidates <- fread(candidates_file, showProgress = FALSE)
  primary <- candidates[selected_for_primary_review == TRUE]
  add_key("validation_candidate_rows", nrow(candidates),
          if (nrow(candidates) > 0) "ok" else "fail")
  add_key("primary_review_rows", nrow(primary),
          if (nrow(primary) > 0) "ok" else "fail")
  add_key("primary_review_genes", uniqueN(primary$gene_symbol), "ok")
  top_candidates <- primary[
    order(primary_review_rank)
  ][1:min(.N, 20), .(
    primary_review_rank, gene_symbol, tx_id, design_family,
    consensus_class, atlas_clean_cds_pct_change,
    strict_supported_studies, total_case_samples, total_control_samples,
    warning_penalty
  )]
  fwrite(top_candidates, file.path(out_dir,
                                   "result_sanity_top_primary_candidates.csv"))
}

if (file.exists(template_file)) {
  template <- fread(template_file, showProgress = FALSE)
  review_cols <- c(
    "validation_status", "browser_signal_grade",
    "branch_direction_verified", "clean_cds_direction_verified",
    "context_match_verified", "likely_failure_mode", "reviewer_notes",
    "reviewer", "review_date", "next_action"
  )
  missing_review_cols <- setdiff(review_cols, names(template))
  add_key("manual_review_template_rows", nrow(template),
          if (nrow(template) > 0) "ok" else "fail")
  add_key("manual_review_template_missing_review_columns",
          length(missing_review_cols),
          if (!length(missing_review_cols)) "ok" else "fail",
          paste(missing_review_cols, collapse = "; "))
}

if (file.exists(pages_file)) {
  pages <- fread(pages_file, showProgress = FALSE)
  add_key("required_fst_pages", nrow(pages),
          if (nrow(pages) > 0) "ok" else "fail")
  if ("local_file_exists" %in% names(pages)) {
    add_key("required_fst_pages_local_present",
            sum(pages$local_file_exists == TRUE, na.rm = TRUE),
            if (any(pages$local_file_exists == FALSE, na.rm = TRUE)) {
              "warn"
            } else if (all(is.na(pages$local_file_exists))) {
              "warn"
            } else {
              "ok"
            },
            paste0("missing=", sum(pages$local_file_exists == FALSE,
                                   na.rm = TRUE),
                   "; unchecked=", sum(is.na(pages$local_file_exists))))
  }
}

if (file.exists(clean_file)) {
  clean <- fread(clean_file, showProgress = FALSE)
  add_key("clean_cds_runs", nrow(clean),
          if (nrow(clean) == 3857) "ok" else "warn",
          "Expected 3857 for current all_samples-Homo_sapiens experiment.")
  add_key("clean_cds_missing_run_ids",
          sum(is.na(clean$Run) | !nzchar(as.character(clean$Run))),
          if (sum(is.na(clean$Run) | !nzchar(as.character(clean$Run))) == 0) {
            "ok"
          } else {
            "fail"
          })
}

if (file.exists(activity_file)) {
  activity <- fread(activity_file, showProgress = FALSE)
  add_key("state_activity_long_rows", nrow(activity),
          if (nrow(activity) > 0) "ok" else "fail")
  add_key("state_activity_states", uniqueN(activity$dominant_state), "ok")
}

add_key("csv_files", nrow(csv_inventory), "ok")
add_key("csv_read_errors", csv_inventory[status == "read_error", .N],
        if (csv_inventory[status == "read_error", .N] == 0) "ok" else "fail")
add_key("csv_empty_files", csv_inventory[status == "empty", .N],
        if (csv_inventory[status == "empty" & note == "", .N] == 0) {
          "ok"
        } else {
          "warn"
        })
add_key("figure_files", nrow(figure_inventory), "ok")
add_key("figure_empty_files", figure_inventory[status == "empty", .N],
        if (figure_inventory[status == "empty", .N] == 0) "ok" else "fail")
add_key("figure_suspiciously_small_files",
        figure_inventory[status == "suspiciously_small", .N],
        if (figure_inventory[status == "suspiciously_small", .N] == 0) {
          "ok"
        } else {
          "warn"
        })
add_key("private_htmlwidget_lib_dirs", length(private_widget_dirs),
        if (!length(private_widget_dirs)) "ok" else "warn")

key_checks <- rbindlist(key_checks, fill = TRUE)

fwrite(csv_inventory, file.path(out_dir, "result_sanity_csv_inventory.csv"))
fwrite(figure_inventory, file.path(out_dir, "result_sanity_figure_inventory.csv"))
fwrite(key_checks, file.path(out_dir, "result_sanity_key_checks.csv"))

report_file <- file.path(out_dir, "result_sanity_report.md")
warn_fail <- key_checks[status %chin% c("warn", "fail")]
report <- c(
  "# Dominant-State Result Sanity Review",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Summary",
  "",
  paste0("- CSV files read: ", nrow(csv_inventory)),
  paste0("- CSV read errors: ", csv_inventory[status == "read_error", .N]),
  paste0("- Empty CSV files: ", csv_inventory[status == "empty", .N],
         " (expected empties are annotated in the inventory)"),
  paste0("- Figure/HTML files: ", nrow(figure_inventory)),
  paste0("- Empty figure/HTML files: ", figure_inventory[status == "empty", .N]),
  paste0("- Suspiciously small static figures: ",
         figure_inventory[status == "suspiciously_small", .N]),
  paste0("- Private htmlwidget dependency folders: ",
         length(private_widget_dirs)),
  "",
  "## Key Checks",
  "",
  paste0("- ", key_checks$metric, ": ", key_checks$value,
         " [", key_checks$status, "]",
         ifelse(nzchar(key_checks$note), paste0(" - ", key_checks$note), "")),
  "",
  "## Warnings / Failures",
  "",
  if (nrow(warn_fail)) {
    paste0("- ", warn_fail$metric, ": ", warn_fail$value,
           " [", warn_fail$status, "]",
           ifelse(nzchar(warn_fail$note), paste0(" - ", warn_fail$note), ""))
  } else {
    "- None."
  },
  "",
  "## Files Written",
  "",
  "- `result_sanity_csv_inventory.csv`",
  "- `result_sanity_figure_inventory.csv`",
  "- `result_sanity_key_checks.csv`",
  "- `result_sanity_top_primary_candidates.csv`"
)
writeLines(report, report_file, useBytes = TRUE)

message("Result sanity review written to: ", report_file)
