#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) source(dominant_loader)

find_repo_root <- function(start = getwd()) {
  if (exists("dominant_find_ribocrypt_repo")) {
    repo <- dominant_find_ribocrypt_repo(start)
    if (nzchar(repo)) return(repo)
  }
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        identical(dominant_description_package(here), "RiboCrypt")) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find RiboCrypt repo root from: ", start, call. = FALSE)
}

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

repo_root <- find_repo_root()
analysis_dir <- find_analysis_dir()

setup_runtime_env <- function(analysis_dir) {
  runtime_root <- file.path(analysis_dir, ".runtime")
  cache_dir <- file.path(runtime_root, "xdg-cache")
  config_dir <- file.path(runtime_root, "xdg-config")
  bfc_dir <- file.path(runtime_root, "biocfilecache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(bfc_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(
    XDG_CACHE_HOME = cache_dir,
    XDG_CONFIG_HOME = config_dir,
    BFC_CACHE = bfc_dir
  )
}

setup_runtime_env(analysis_dir)

if (!"RiboCrypt" %in% loadedNamespaces()) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load local RiboCrypt.",
         call. = FALSE)
  }
  devtools::load_all(repo_root, quiet = TRUE)
}

as_bool <- function(x, default = FALSE) {
  if (!length(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(x) %chin% c("1", "true", "t", "yes", "y")
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--force")) {
      out$force <- TRUE
    } else if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=", "", arg)
      out[[gsub("-", "_", key)]] <- value
    }
  }
  out
}

opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))

default_queue <- file.path(
  analysis_dir,
  "dominant_rdg_review_links",
  "dominant_rdg_review_batch_template_with_links.csv"
)
queue_file <- opts$queue %||% Sys.getenv("RDG_QUALIFIED_QUEUE",
                                         unset = default_queue)
output_dir <- opts$output_dir %||% Sys.getenv(
  "RDG_QUALIFIED_OUTPUT_DIR",
  unset = file.path(analysis_dir, "dominant_rdg_qualified_inspection")
)
limit_raw <- opts$limit %||% Sys.getenv("RDG_QUALIFIED_LIMIT", unset = "")
row_raw <- opts$row %||% Sys.getenv("RDG_QUALIFIED_ROW", unset = "")
force <- isTRUE(opts$force) ||
  as_bool(Sys.getenv("RDG_QUALIFIED_FORCE", unset = ""), FALSE)
stop_on_error <- as_bool(Sys.getenv("RDG_QUALIFIED_STOP_ON_ERROR", unset = ""),
                         FALSE)
max_profile_cells <- suppressWarnings(as.numeric(
  Sys.getenv("RDG_QUALIFIED_MAX_PROFILE_CELLS", unset = "20000000")
))
if (!is.finite(max_profile_cells) || max_profile_cells <= 0) {
  max_profile_cells <- 20000000
}
max_plot_positions <- suppressWarnings(as.integer(
  Sys.getenv("RDG_QUALIFIED_MAX_PLOT_POSITIONS", unset = "3500")
))
if (!is.finite(max_plot_positions) || max_plot_positions <= 0L) {
  max_plot_positions <- 3500L
}

limit_n <- suppressWarnings(as.integer(limit_raw))
if (length(limit_n) == 0 || is.na(limit_n) || limit_n <= 0L) {
  limit_n <- NA_integer_
}
row_n <- suppressWarnings(as.integer(row_raw))
if (length(row_n) == 0 || is.na(row_n) || row_n <= 0L) {
  row_n <- NA_integer_
}

profile_dir <- file.path(output_dir, "profiles")
feature_dir <- file.path(output_dir, "feature_allocation")
flow_dir <- file.path(output_dir, "rdg_flow")
replicate_dir <- file.path(output_dir, "replicate_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(flow_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
feature_expression_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
feature_annotation_file <- file.path(analysis_dir, "dominant_rdg_outputs",
                                     "rdg_feature_annotation.csv")
rdg_edges_file <- file.path(analysis_dir, "dominant_rdg_outputs",
                            "rdg_edges.csv")

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

scalar_chr <- function(x, default = "") {
  if (!length(x) || is.na(x[[1]])) return(default)
  as.character(x[[1]])
}

scalar_num <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x[[1]]))
  if (!length(value) || !is.finite(value)) default else value
}

max_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

mean_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

median_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

sd_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x)
}

coalesce_col <- function(row, names, default = "") {
  for (name in names) {
    if (name %in% names(row)) {
      value <- scalar_chr(row[[name]])
      if (nzchar(value)) return(value)
    }
  }
  default
}

first_finite_num <- function(row, names, default = NA_real_) {
  for (name in names) {
    if (name %in% names(row)) {
      value <- suppressWarnings(as.numeric(row[[name]][[1]]))
      if (is.finite(value)) return(value)
    }
  }
  default
}

real_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  !x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null",
              "None", "none", "MISSING", "missing", "unknown",
              "Unknown")
}

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !real_value(x)] <- ""
  x
}

split_values <- function(x) {
  x <- clean_value(x)
  if (!nzchar(x)) return(character())
  out <- unlist(strsplit(x, "\\s*[,;]\\s*"), use.names = FALSE)
  out <- clean_value(out)
  unique(out[nzchar(out)])
}

equal_filter <- function(dt, column, value) {
  keep <- dt[[column]] == value
  keep[is.na(keep)] <- FALSE
  keep
}

short_label <- function(x, max_chars = 82L) {
  x <- gsub("\\s+", " ", as.character(x))
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 3L), "...")
  x
}

extract_hash <- function(url) {
  if (!grepl("#", url, fixed = TRUE)) return("")
  sub("^[^#]*(#.*)$", "\\1", url)
}

parse_observatory_state <- function(url) {
  hash <- extract_hash(url)
  if (!nzchar(hash)) return(NULL)
  RiboCrypt:::parse_observatory_url_hash(hash)
}

selection_label <- function(labels, id) {
  value <- labels[[id]]
  if (is.null(value)) id else scalar_chr(value, id)
}

selection_role <- function(label) {
  lower <- tolower(label)
  if (grepl("control", lower, fixed = TRUE)) return("control")
  if (grepl("case", lower, fixed = TRUE)) return("case")
  if (grepl("all merged", lower, fixed = TRUE)) return("all_merged")
  "selection"
}

state_selection_groups <- function(state) {
  selections <- state$selections$plot_selections %||% list()
  labels <- state$selections$labels %||% list()
  ids <- as.character(state$selections$index %||% names(selections))
  if (!length(ids)) ids <- names(selections)
  ids <- ids[nzchar(ids)]
  groups <- lapply(ids, function(id) {
    runs <- unique(as.character(selections[[id]] %||% character()))
    runs <- runs[!is.na(runs) & nzchar(runs)]
    label <- selection_label(labels, id)
    list(
      id = id,
      label = label,
      role = selection_role(label),
      runs = runs
    )
  })
  groups <- groups[lengths(lapply(groups, `[[`, "runs")) > 0]
  groups
}

collection_available_runs <- function(collection_path) {
  if (!file.exists(collection_path)) return(character())
  if (!identical(basename(collection_path), "coverage_index.fst")) {
    return(colnames(fst::read_fst(collection_path, from = 1, to = 1)))
  }
  index_row <- fst::read_fst(collection_path, from = 1, to = 1)
  if (!"file_forward" %in% names(index_row)) return(character())
  file_forward <- file.path(dirname(collection_path), basename(index_row$file_forward[[1]]))
  if (!file.exists(file_forward)) return(character())
  fst::metadata_fst(file_forward)$columnNames
}

display_region_from_state <- function(tx_region, browser_state) {
  out <- tx_region
  leader <- suppressWarnings(as.numeric(browser_state$extendLeaders %||% 0))
  trailer <- suppressWarnings(as.numeric(browser_state$extendTrailers %||% 0))
  if (is.finite(leader) && leader != 0) out <- ORFik::extendLeaders(out, leader)
  if (is.finite(trailer) && trailer != 0) out <- ORFik::extendTrailers(out, trailer)
  out
}

filter_metadata_runs <- function(metadata, available_runs, study = "",
                                 tissue = "", cell_line = "",
                                 condition = character()) {
  dt <- metadata[Run %chin% available_runs]
  study <- clean_value(study)
  tissue <- clean_value(tissue)
  cell_line <- clean_value(cell_line)
  condition <- clean_value(condition)
  condition <- condition[nzchar(condition)]
  if (nzchar(study) && "study" %chin% names(dt)) {
    keep <- equal_filter(dt, "study", study)
    dt <- dt[keep]
  }
  if (nzchar(tissue) && "TISSUE" %chin% names(dt)) {
    keep <- equal_filter(dt, "TISSUE", tissue)
    dt <- dt[keep]
  }
  if (nzchar(cell_line) && "CELL_LINE" %chin% names(dt)) {
    keep <- equal_filter(dt, "CELL_LINE", cell_line)
    dt <- dt[keep]
  }
  if (length(condition) && "CONDITION" %chin% names(dt)) {
    keep <- dt[["CONDITION"]] %chin% condition
    keep[is.na(keep)] <- FALSE
    dt <- dt[keep]
  }
  unique(dt$Run)
}

make_group_row <- function(row_index, group_id, group_label, group_role,
                           group_scope, comparison_anchor, priority, runs,
                           available_runs) {
  runs <- unique(as.character(runs))
  runs <- runs[!is.na(runs) & nzchar(runs)]
  present <- intersect(runs, available_runs)
  data.table(
    source_row = row_index,
    group_id = group_id,
    group_label = group_label,
    group_role = group_role,
    group_scope = group_scope,
    comparison_anchor = comparison_anchor,
    profile_priority = priority,
    requested_run_count = length(runs),
    available_run_count = length(present),
    missing_run_count = length(setdiff(runs, available_runs)),
    runs = list(present),
    requested_runs = list(runs)
  )
}

build_context_groups <- function(row, row_index, state, metadata, available_runs) {
  exact <- state_selection_groups(state)
  exact <- exact[vapply(exact, function(x) x$role != "all_merged", logical(1))]
  pieces <- list()
  for (i in seq_along(exact)) {
    role <- exact[[i]]$role
    scope <- paste0("exact_", role)
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index = row_index,
      group_id = paste0("exact_", role),
      group_label = paste0("Exact ", role, ": ", exact[[i]]$label),
      group_role = role,
      group_scope = scope,
      comparison_anchor = "exact_case_control",
      priority = if (role %chin% c("case", "control")) 1L else 2L,
      runs = exact[[i]]$runs,
      available_runs = available_runs
    )
  }

  study <- coalesce_col(row, c("study", "best_study"), default = "")
  tissue <- coalesce_col(row, c("TISSUE"), default = "")
  cell_line <- coalesce_col(row, c("CELL_LINE"), default = "")
  case_condition <- coalesce_col(row, c("case_condition", "best_case_condition",
                                        "observatory_case_condition"), default = "")
  control_conditions <- split_values(coalesce_col(
    row,
    c("control_conditions", "best_control_conditions",
      "observatory_control_conditions"),
    default = ""
  ))

  study_case <- filter_metadata_runs(metadata, available_runs, study = study,
                                     tissue = tissue, cell_line = cell_line,
                                     condition = case_condition)
  study_control <- filter_metadata_runs(metadata, available_runs, study = study,
                                        tissue = tissue, cell_line = cell_line,
                                        condition = control_conditions)
  tissue_case <- filter_metadata_runs(metadata, available_runs, tissue = tissue,
                                      condition = case_condition)
  tissue_control <- filter_metadata_runs(metadata, available_runs, tissue = tissue,
                                         condition = control_conditions)
  tissue_all <- filter_metadata_runs(metadata, available_runs, tissue = tissue)
  condition_all <- filter_metadata_runs(metadata, available_runs,
                                        condition = case_condition)
  study_all <- filter_metadata_runs(metadata, available_runs, study = study)

  if (length(study_case)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "study_context_case",
      paste0("Study context case: ", study, " | ", tissue, " | ",
             cell_line, " | ", case_condition),
      "case", "study_context_case", "study_context_case_control", 2L,
      study_case, available_runs
    )
  }
  if (length(study_control)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "study_context_control",
      paste0("Study context control: ", study, " | ", tissue, " | ",
             cell_line, " | ", paste(control_conditions, collapse = "/")),
      "control", "study_context_control", "study_context_case_control", 2L,
      study_control, available_runs
    )
  }
  if (length(tissue_case)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "tissue_case_condition",
      paste0(tissue, " ", case_condition, " across studies"),
      "context_case", "tissue_case_condition", "tissue_case_vs_tissue_all",
      3L, tissue_case, available_runs
    )
  }
  if (length(tissue_control)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "tissue_control_condition",
      paste0(tissue, " ", paste(control_conditions, collapse = "/"),
             " across studies"),
      "context_control", "tissue_control_condition",
      "tissue_control_vs_tissue_all", 3L, tissue_control, available_runs
    )
  }
  if (length(tissue_all)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "tissue_general",
      paste0(tissue, " all available conditions"),
      "context_general", "tissue_general", "broad_context", 4L,
      tissue_all, available_runs
    )
  }
  if (length(condition_all)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "case_condition_global",
      paste0(case_condition, " across all tissues"),
      "context_case", "case_condition_global", "condition_vs_all_merged", 5L,
      condition_all, available_runs
    )
  }
  if (length(study_all)) {
    pieces[[length(pieces) + 1L]] <- make_group_row(
      row_index, "study_general",
      paste0(study, " all available conditions"),
      "context_general", "study_general", "broad_context", 5L,
      study_all, available_runs
    )
  }
  pieces[[length(pieces) + 1L]] <- make_group_row(
    row_index, "all_merged_available",
    "All merged available runs for transcript",
    "all_merged", "all_merged", "global_baseline", 6L,
    available_runs, available_runs
  )

  groups <- rbindlist(pieces, fill = TRUE)
  groups <- groups[available_run_count > 0]
  if (!nrow(groups)) return(groups)
  groups[, run_signature := vapply(runs, function(x) {
    paste(sort(unique(x)), collapse = ";")
  }, character(1))]
  groups[, duplicate_of := group_id[1L], by = run_signature]
  groups[, duplicate_run_set := duplicate_of != group_id]
  groups[, group_label_short := short_label(paste0(group_label, " (n=",
                                                   available_run_count, ")"))]
  setorder(groups, profile_priority, group_id)
  groups[]
}

choose_profile_groups <- function(groups, display_width, max_profile_cells) {
  if (!nrow(groups)) return(groups)
  groups[, profile_include := FALSE]
  groups[, profile_skip_reason := ""]
  groups[startsWith(group_id, "exact_"), profile_include := TRUE]
  groups[duplicate_run_set == TRUE & profile_include == FALSE,
         profile_skip_reason := "duplicate_run_set"]
  included_runs <- unique(unlist(groups[profile_include == TRUE, runs],
                                 use.names = FALSE))
  candidate_idx <- which(groups$profile_include == FALSE &
                           groups$duplicate_run_set == FALSE)
  candidate_idx <- candidate_idx[order(groups$profile_priority[candidate_idx],
                                       groups$group_id[candidate_idx])]
  for (i in candidate_idx) {
    proposed_runs <- unique(c(included_runs, groups$runs[[i]]))
    proposed_cells <- length(proposed_runs) * display_width
    if (proposed_cells <= max_profile_cells) {
      groups[i, profile_include := TRUE]
      included_runs <- proposed_runs
    } else {
      groups[i, profile_skip_reason := paste0(
        "profile_cells_limit_", max_profile_cells
      )]
    }
  }
  groups[, profile_union_run_count := length(included_runs)]
  groups[, profile_union_cell_count := length(included_runs) * display_width]
  groups[]
}

profile_groups_from_reads <- function(reads, groups) {
  mat <- as.matrix(reads)
  rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    cols <- intersect(groups$runs[[i]], colnames(mat))
    counts <- if (length(cols)) {
      matrixStats::rowSums2(mat, cols = match(cols, colnames(mat)),
                            na.rm = TRUE)
    } else {
      rep(0, nrow(mat))
    }
    total <- sum(counts, na.rm = TRUE)
    data.table(
      source_row = groups$source_row[[i]],
      group_id = groups$group_id[[i]],
      group_label = groups$group_label_short[[i]],
      group_role = groups$group_role[[i]],
      group_scope = groups$group_scope[[i]],
      run_count = length(cols),
      position = seq_len(nrow(mat)),
      frame = factor((seq_len(nrow(mat)) - 1L) %% 3L, levels = 0:2),
      count = as.numeric(counts),
      norm_per_1k = if (is.finite(total) && total > 0) {
        as.numeric(counts) / total * 1000
      } else {
        rep(0, length(counts))
      }
    )
  }), fill = TRUE)
}

profile_metrics <- function(profile_dt) {
  if (!nrow(profile_dt)) return(data.table())
  profile_dt[, {
    counts <- count
    total <- sum(counts, na.rm = TRUE)
    ord <- sort(counts, decreasing = TRUE)
    top3 <- sum(head(ord, 3L), na.rm = TRUE)
    .(
      run_count = unique(run_count)[1L],
      display_width = .N,
      display_total_counts = total,
      max_position_count = max(counts, na.rm = TRUE),
      max_position_fraction = if (total > 0) max(counts, na.rm = TRUE) / total else NA_real_,
      top3_position_fraction = if (total > 0) top3 / total else NA_real_,
      nonzero_positions = sum(counts > 0, na.rm = TRUE)
    )
  }, by = .(source_row, group_id, group_label, group_role, group_scope)]
}

feature_track <- function(feature_annotation, edges, tx_id, leader_extension,
                          display_width, profile_dt) {
  keep_fx <- equal_filter(feature_annotation, "tx_id", tx_id)
  fx <- copy(feature_annotation[keep_fx])
  if (!nrow(fx)) return(list(features = data.table(), edges = data.table()))
  fx[, `:=`(
    xmin = as.numeric(tx_start) + leader_extension,
    xmax = as.numeric(tx_end) + leader_extension
  )]
  fx <- fx[is.finite(xmin) & is.finite(xmax) & xmax >= 1 & xmin <= display_width]
  fx[, `:=`(
    xmin = pmax(1, xmin),
    xmax = pmin(display_width, xmax)
  )]
  if (nrow(profile_dt)) {
    groups <- unique(profile_dt[, .(group_id, group_label)])
    y <- profile_dt[, .(ymax_group = max(norm_per_1k, 0.01, na.rm = TRUE)),
                    by = group_id]
    strip <- fx[rep(seq_len(nrow(fx)), nrow(groups))]
    strip[, group_id := rep(groups$group_id, each = nrow(fx))]
    strip[, group_label := rep(groups$group_label, each = nrow(fx))]
    strip <- y[strip, on = "group_id"]
    strip[, `:=`(
      ymin = -0.20 * ymax_group,
      ymax = -0.05 * ymax_group,
      label_y = -0.125 * ymax_group
    )]
  } else {
    strip <- fx
  }

  keep_edges <- equal_filter(edges, "tx_id", tx_id)
  ed <- copy(edges[keep_edges])
  if (nrow(ed)) {
    ed[, `:=`(
      x = as.numeric(x) + leader_extension,
      xend = as.numeric(xend) + leader_extension
    )]
    ed <- ed[is.finite(x) & is.finite(xend) & pmax(x, xend) >= 1 &
               pmin(x, xend) <= display_width]
    ed[, `:=`(
      x = pmax(1, pmin(display_width, x)),
      xend = pmax(1, pmin(display_width, xend))
    )]
  }
  list(features = strip, edges = ed)
}

plot_profile <- function(profile_dt, strip_dt, row, top_feature, output_file,
                         max_plot_positions) {
  if (!nrow(profile_dt)) return(FALSE)
  display_width <- max(profile_dt$position, na.rm = TRUE)
  bin_size <- max(1L, ceiling(display_width / max_plot_positions))
  plot_dt <- copy(profile_dt)
  if (bin_size > 1L) {
    plot_dt[, position_bin := ((position - 1L) %/% bin_size) * bin_size + 1L]
    plot_dt <- plot_dt[, .(
      count = sum(count, na.rm = TRUE),
      norm_per_1k = sum(norm_per_1k, na.rm = TRUE),
      run_count = unique(run_count)[1L]
    ), by = .(source_row, group_id, group_label, group_role, group_scope,
              position = position_bin)]
  }
  plot_dt[, group_label := factor(group_label, levels = unique(group_label))]
  if (nrow(strip_dt)) {
    strip_dt[, group_label := factor(group_label, levels = levels(plot_dt$group_label))]
  }

  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  tx_id <- coalesce_col(row, c("tx_id"), default = "")
  design <- coalesce_col(row, c("design_family"), default = "")
  top_branch <- coalesce_col(row, c("joint_top_branch_class", "main_branch_class"),
                             default = "")
  expected <- coalesce_col(row, c("expected_browser_pattern"), default = "")
  title <- paste(gene, tx_id, design, sep = " | ")
  subtitle <- paste(c(
    paste0("Expected top branch: ", top_branch,
           if (nzchar(top_feature)) paste0("; top feature: ", top_feature) else ""),
    expected
  ), collapse = "\n")

  role_cols <- c(
    case = "#C44E52",
    control = "#4C78A8",
    context_case = "#E15759",
    context_control = "#72B7B2",
    context_general = "#59A14F",
    all_merged = "#79706E",
    selection = "#B07AA1"
  )
  feature_cols <- c(
    clean_CDS = "#4C78A8",
    leader_uorf = "#F58518",
    overlapping_uorf = "#E45756",
    internal_orf = "#54A24B",
    nte_extension = "#B279A2",
    ntt_truncation = "#9D755D",
    cds_translon = "#BAB0AC",
    other_translon = "#BAB0AC"
  )

  p <- ggplot(plot_dt, aes(position, norm_per_1k)) +
    geom_line(aes(color = group_role), linewidth = 0.35) +
    facet_grid(group_label ~ ., scales = "free_y", switch = "y") +
    scale_y_continuous(expand = expansion(mult = c(0.24, 0.08))) +
    scale_color_manual(values = role_cols, drop = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Displayed transcript coordinate",
      y = "Normalized P-site coverage per 1k displayed counts",
      color = "Group"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, lineheight = 0.95),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 7.5),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  if (nrow(strip_dt)) {
    p <- p +
      geom_rect(
        data = strip_dt,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
            fill = feature_type),
        inherit.aes = FALSE,
        alpha = 0.75,
        color = "grey25",
        linewidth = 0.18
      ) +
      geom_rect(
        data = strip_dt[feature_id == top_feature],
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = "black",
        linewidth = 0.55
      ) +
      scale_fill_manual(values = feature_cols, drop = FALSE)
  }

  ggsave(output_file, p,
         width = 12.5,
         height = max(5.2, 2.5 + uniqueN(plot_dt$group_label) * 1.15),
         dpi = 135,
         limitsize = FALSE)
  TRUE
}

summarize_feature_counts <- function(row, row_index, groups, feature_expression,
                                     feature_annotation) {
  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  tx_value <- coalesce_col(row, c("tx_id"), default = "")
  fe <- feature_expression[
    data.table(gene_symbol = gene, tx_id = tx_value),
    on = c("gene_symbol", "tx_id"),
    nomatch = 0
  ]
  keep_ann <- feature_annotation[["gene_symbol"]] == gene &
    feature_annotation[["tx_id"]] == tx_value
  keep_ann[is.na(keep_ann)] <- FALSE
  ann <- feature_annotation[keep_ann]
  if (!nrow(fe) && !nrow(ann)) return(data.table())

  out <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    runs <- groups$runs[[i]]
    dt <- fe[Run %chin% runs]
    if (!nrow(dt)) {
      if (!nrow(ann)) return(NULL)
      return(ann[, .(
        source_row = row_index,
        group_id = groups$group_id[[i]],
        group_label = groups$group_label_short[[i]],
        group_role = groups$group_role[[i]],
        group_scope = groups$group_scope[[i]],
        run_count = length(runs),
        feature_id,
        feature_type,
        category,
        raw_counts_sum = 0,
        fpkm_like_mean = NA_real_,
        fpkm_like_median = NA_real_,
        signal_run_count = 0L,
        max_run_raw_counts = 0,
        max_run_fraction = NA_real_
      )])
    }
    dt[, .(
      source_row = row_index,
      group_id = groups$group_id[[i]],
      group_label = groups$group_label_short[[i]],
      group_role = groups$group_role[[i]],
      group_scope = groups$group_scope[[i]],
      run_count = length(intersect(runs, unique(dt$Run))),
      raw_counts_sum = sum(raw_counts, na.rm = TRUE),
      fpkm_like_mean = mean(fpkm_like, na.rm = TRUE),
      fpkm_like_median = median(fpkm_like, na.rm = TRUE),
      signal_run_count = sum(raw_counts > 0, na.rm = TRUE),
      max_run_raw_counts = max(raw_counts, 0, na.rm = TRUE)
    ), by = .(feature_id, feature_type, category)]
  }), fill = TRUE)

  if (!nrow(out)) return(out)
  out[, max_run_fraction := fifelse(raw_counts_sum > 0,
                                    max_run_raw_counts / raw_counts_sum,
                                    NA_real_)]
  out[, total_all_feature_counts := sum(raw_counts_sum, na.rm = TRUE),
      by = .(source_row, group_id)]
  out[, total_allocation_counts := sum(
    raw_counts_sum[feature_type %chin% c("leader_uorf", "overlapping_uorf",
                                         "clean_CDS")],
    na.rm = TRUE
  ), by = .(source_row, group_id)]
  out[, prop_all_feature := fifelse(total_all_feature_counts > 0,
                                    raw_counts_sum / total_all_feature_counts,
                                    NA_real_)]
  out[, prop_allocation := fifelse(
    total_allocation_counts > 0 &
      feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS"),
    raw_counts_sum / total_allocation_counts,
    NA_real_
  )]
  out[]
}

contrast_pair <- function(feature_counts, row_index, contrast_id,
                          numerator_group, denominator_group) {
  num <- feature_counts[group_id == numerator_group]
  den <- feature_counts[group_id == denominator_group]
  if (!nrow(num) || !nrow(den)) return(data.table())
  merged <- merge(
    num[, .(feature_id, feature_type, category,
            numerator_raw = raw_counts_sum,
            numerator_prop = prop_all_feature,
            numerator_allocation_prop = prop_allocation,
            numerator_max_run_fraction = max_run_fraction,
            numerator_run_count = run_count)],
    den[, .(feature_id, denominator_raw = raw_counts_sum,
            denominator_prop = prop_all_feature,
            denominator_allocation_prop = prop_allocation,
            denominator_max_run_fraction = max_run_fraction,
            denominator_run_count = run_count)],
    by = "feature_id",
    all = TRUE
  )
  merged[is.na(numerator_raw), numerator_raw := 0]
  merged[is.na(denominator_raw), denominator_raw := 0]
  merged[, `:=`(
    source_row = row_index,
    contrast_id = contrast_id,
    numerator_group = numerator_group,
    denominator_group = denominator_group,
    raw_log2_ratio = log2((numerator_raw + 1) / (denominator_raw + 1)),
    prop_delta = numerator_prop - denominator_prop,
    allocation_prop_delta = numerator_allocation_prop - denominator_allocation_prop
  )]
  setcolorder(merged, c("source_row", "contrast_id", "numerator_group",
                        "denominator_group"))
  merged[]
}

feature_contrasts <- function(feature_counts, row_index) {
  pairs <- list(
    c("exact_case_vs_exact_control", "exact_case", "exact_control"),
    c("study_context_case_vs_control", "study_context_case",
      "study_context_control"),
    c("tissue_case_vs_tissue_general", "tissue_case_condition",
      "tissue_general"),
    c("tissue_control_vs_tissue_general", "tissue_control_condition",
      "tissue_general"),
    c("case_condition_global_vs_all_merged", "case_condition_global",
      "all_merged_available"),
    c("exact_case_vs_all_merged", "exact_case", "all_merged_available")
  )
  rbindlist(lapply(pairs, function(pair) {
    contrast_pair(feature_counts, row_index, pair[[1]], pair[[2]], pair[[3]])
  }), fill = TRUE)
}

plot_feature_allocation <- function(feature_counts, row, top_feature, output_file) {
  if (!nrow(feature_counts)) return(FALSE)
  keep_groups <- c("exact_control", "exact_case", "study_context_control",
                   "study_context_case", "tissue_case_condition",
                   "tissue_general", "case_condition_global",
                   "all_merged_available")
  dt <- copy(feature_counts[group_id %chin% keep_groups])
  if (!nrow(dt)) return(FALSE)
  dt[, plot_prop_all_feature := fifelse(is.finite(prop_all_feature),
                                        prop_all_feature, 0)]
  dt[, group_label := factor(group_label, levels = unique(group_label))]
  dt[, feature_id := factor(feature_id, levels = unique(feature_id[order(feature_type)]))]
  dt[, top_feature := feature_id == top_feature]
  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  tx_id <- coalesce_col(row, c("tx_id"), default = "")
  title <- paste(gene, tx_id, "feature allocation", sep = " | ")

  p <- ggplot(dt, aes(feature_id, plot_prop_all_feature, fill = group_label)) +
    geom_col(aes(color = top_feature),
             position = position_dodge2(width = 0.82, preserve = "single"),
             width = 0.75, linewidth = 0.18) +
    facet_grid(. ~ feature_type, scales = "free_x", space = "free_x") +
    scale_color_manual(
      values = c("FALSE" = "grey35", "TRUE" = "black"),
      guide = "none",
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = "Bars show each feature's share of all modeled feature counts in the group; black outline marks the expected top feature.",
      x = "RDG feature",
      y = "Share of modeled feature counts",
      fill = "Group"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8),
      axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 7)
    )

  ggsave(output_file, p, width = 12.5,
         height = max(4.8, 3.4 + uniqueN(dt$group_label) * 0.14),
         dpi = 135,
         limitsize = FALSE)
  TRUE
}

summarize_replicate_diagnostics <- function(row, row_index, groups,
                                            feature_expression, metadata,
                                            top_feature) {
  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  tx_value <- coalesce_col(row, c("tx_id"), default = "")
  fe <- feature_expression[
    data.table(gene_symbol = gene, tx_id = tx_value),
    on = c("gene_symbol", "tx_id"),
    nomatch = 0
  ]
  exact_groups <- groups[group_id %chin% c("exact_case", "exact_control")]
  if (!nrow(fe) || !nrow(exact_groups)) return(data.table())

  run_map <- rbindlist(lapply(seq_len(nrow(exact_groups)), function(i) {
    runs <- exact_groups$runs[[i]]
    if (!length(runs)) return(NULL)
    data.table(
      source_row = row_index,
      group_id = exact_groups$group_id[[i]],
      group_label = exact_groups$group_label_short[[i]],
      group_role = exact_groups$group_role[[i]],
      Run = runs,
      run_order = seq_along(runs)
    )
  }), fill = TRUE)
  if (!nrow(run_map)) return(data.table())

  dt <- fe[Run %chin% run_map$Run]
  if (!nrow(dt)) return(data.table())
  run_totals <- dt[, .(
    run_total_all_feature_counts = sum(raw_counts, na.rm = TRUE),
    run_total_allocation_counts = sum(
      raw_counts[feature_type %chin% c("leader_uorf", "overlapping_uorf",
                                       "clean_CDS")],
      na.rm = TRUE
    )
  ), by = Run]

  selected_features <- unique(c(top_feature, "clean_CDS"))
  selected_features <- selected_features[nzchar(selected_features)]
  out <- dt[feature_id %chin% selected_features]
  if (!nrow(out)) return(data.table())
  out <- run_map[out, on = "Run", nomatch = 0]
  out <- run_totals[out, on = "Run"]

  metadata_cols <- intersect(
    c("Run", "study", "AUTHOR", "CONDITION", "TISSUE", "CELL_LINE",
      "GENE", "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT",
      "Total_Cov_Signal", "Frame_usage_0", "Frame_usage_1",
      "Frame_usage_2", "top_readlength", "tis_peak_strength"),
    names(metadata)
  )
  if (length(metadata_cols) > 1L) {
    md <- unique(metadata[, ..metadata_cols], by = "Run")
    out <- md[out, on = "Run"]
  }

  out[, `:=`(
    feature_role = fifelse(feature_id == top_feature, "top_feature",
                           fifelse(feature_id == "clean_CDS", "clean_CDS",
                                   "other")),
    run_feature_prop_all = fifelse(
      run_total_all_feature_counts > 0,
      raw_counts / run_total_all_feature_counts,
      NA_real_
    ),
    run_feature_prop_allocation = fifelse(
      run_total_allocation_counts > 0 &
        feature_type %chin% c("leader_uorf", "overlapping_uorf", "clean_CDS"),
      raw_counts / run_total_allocation_counts,
      NA_real_
    )
  )]
  setorder(out, group_id, run_order, feature_role)
  out[]
}

summarize_replicate_groups <- function(replicate_dt, top_feature) {
  if (!nrow(replicate_dt)) return(data.table())
  dt <- replicate_dt[feature_id == top_feature]
  if (!nrow(dt)) return(data.table())
  dt[, {
    all_counts <- suppressWarnings(as.numeric(raw_counts))
    measured <- is.finite(all_counts)
    counts <- all_counts[measured]
    total <- if (length(counts)) sum(counts) else 0
    max_count <- if (length(counts)) max(counts) else NA_real_
    .(
      exact_run_count = as.integer(.N),
      exact_measured_run_count = as.integer(sum(measured)),
      exact_missing_run_count = as.integer(sum(!measured)),
      exact_signal_run_count = as.integer(sum(counts > 0)),
      exact_total_top_feature_counts = as.numeric(total),
      exact_median_top_feature_counts = if (length(counts)) {
        as.numeric(median(counts))
      } else {
        NA_real_
      },
      exact_max_top_feature_counts = as.numeric(max_count),
      exact_max_run_fraction = if (total > 0) max_count / total else NA_real_,
      exact_signal_run_fraction = as.numeric(sum(counts > 0)) / .N
    )
  }, by = .(source_row, group_id, group_label, group_role)]
}

summarize_replicate_effects <- function(replicate_diagnostics) {
  if (!nrow(replicate_diagnostics)) return(data.table())
  dt <- copy(replicate_diagnostics[
    feature_role == "top_feature" &
      group_id %chin% c("exact_case", "exact_control")
  ])
  if (!nrow(dt)) return(data.table())
  dt[, `:=`(
    raw_counts_num = suppressWarnings(as.numeric(raw_counts)),
    prop_all_num = suppressWarnings(as.numeric(run_feature_prop_all)),
    prop_allocation_num = suppressWarnings(as.numeric(run_feature_prop_allocation))
  )]

  group_summary <- dt[, {
    measured <- is.finite(raw_counts_num)
    .(
      replicate_run_count = as.integer(.N),
      replicate_measured_run_count = as.integer(sum(measured)),
      replicate_missing_run_count = as.integer(sum(!measured)),
      replicate_signal_run_count = as.integer(sum(raw_counts_num[measured] > 0)),
      replicate_mean_raw_counts = mean_or_na(raw_counts_num),
      replicate_median_raw_counts = median_or_na(raw_counts_num),
      replicate_sd_raw_counts = sd_or_na(raw_counts_num),
      replicate_mean_prop_all = mean_or_na(prop_all_num),
      replicate_median_prop_all = median_or_na(prop_all_num),
      replicate_sd_prop_all = sd_or_na(prop_all_num),
      replicate_mean_prop_allocation = mean_or_na(prop_allocation_num),
      replicate_median_prop_allocation = median_or_na(prop_allocation_num),
      replicate_sd_prop_allocation = sd_or_na(prop_allocation_num)
    )
  }, by = .(source_row, group_id)]

  out <- dcast(
    group_summary,
    source_row ~ group_id,
    value.var = c("replicate_run_count", "replicate_measured_run_count",
                  "replicate_missing_run_count", "replicate_signal_run_count",
                  "replicate_mean_raw_counts", "replicate_median_raw_counts",
                  "replicate_sd_raw_counts", "replicate_mean_prop_all",
                  "replicate_median_prop_all", "replicate_sd_prop_all",
                  "replicate_mean_prop_allocation",
                  "replicate_median_prop_allocation",
                  "replicate_sd_prop_allocation"),
    fill = NA
  )
  if (!nrow(out)) return(out)

  required <- c(
    "replicate_mean_raw_counts_exact_case",
    "replicate_mean_raw_counts_exact_control",
    "replicate_mean_prop_all_exact_case",
    "replicate_mean_prop_all_exact_control",
    "replicate_median_prop_all_exact_case",
    "replicate_median_prop_all_exact_control",
    "replicate_mean_prop_allocation_exact_case",
    "replicate_mean_prop_allocation_exact_control",
    "replicate_median_prop_allocation_exact_case",
    "replicate_median_prop_allocation_exact_control"
  )
  for (col in setdiff(required, names(out))) out[, (col) := NA_real_]

  out[, `:=`(
    replicate_mean_raw_log2_ratio = log2(
      (replicate_mean_raw_counts_exact_case + 1) /
        (replicate_mean_raw_counts_exact_control + 1)
    ),
    replicate_mean_prop_delta = replicate_mean_prop_all_exact_case -
      replicate_mean_prop_all_exact_control,
    replicate_median_prop_delta = replicate_median_prop_all_exact_case -
      replicate_median_prop_all_exact_control,
    replicate_mean_allocation_prop_delta =
      replicate_mean_prop_allocation_exact_case -
      replicate_mean_prop_allocation_exact_control,
    replicate_median_allocation_prop_delta =
      replicate_median_prop_allocation_exact_case -
      replicate_median_prop_allocation_exact_control
  )]
  out[]
}

replicate_measurement_flags <- function(replicate_summary) {
  if (!nrow(replicate_summary)) {
    return(data.table(source_row = integer(), replicate_measurement_flag = character()))
  }
  wide <- dcast(
    replicate_summary,
    source_row ~ group_id,
    value.var = c("exact_run_count", "exact_measured_run_count",
                  "exact_missing_run_count"),
    fill = NA
  )
  for (col in c("exact_run_count_exact_case", "exact_run_count_exact_control",
                "exact_measured_run_count_exact_case",
                "exact_measured_run_count_exact_control",
                "exact_missing_run_count_exact_case",
                "exact_missing_run_count_exact_control")) {
    if (!col %chin% names(wide)) wide[, (col) := NA_integer_]
  }
  missing_cols <- c("exact_missing_run_count_exact_case",
                    "exact_missing_run_count_exact_control")
  measured_cols <- c("exact_measured_run_count_exact_case",
                     "exact_measured_run_count_exact_control")
  run_cols <- c("exact_run_count_exact_case", "exact_run_count_exact_control")
  wide[, `:=`(
    exact_missing_run_count_total = rowSums(.SD, na.rm = TRUE)
  ), .SDcols = missing_cols]
  wide[, `:=`(
    exact_measured_run_count_total = rowSums(.SD, na.rm = TRUE)
  ), .SDcols = measured_cols]
  wide[, `:=`(
    exact_run_count_total = rowSums(.SD, na.rm = TRUE)
  ), .SDcols = run_cols]
  wide[, replicate_measurement_flag := fcase(
    exact_run_count_total > 0 &
      exact_measured_run_count_total == 0 &
      exact_missing_run_count_total > 0,
    "top_feature_exact_replicates_all_missing",
    (exact_run_count_exact_case > 0 &
       exact_measured_run_count_exact_case == 0 &
       exact_missing_run_count_exact_case > 0) |
      (exact_run_count_exact_control > 0 &
         exact_measured_run_count_exact_control == 0 &
         exact_missing_run_count_exact_control > 0),
    "top_feature_exact_replicate_group_missing",
    exact_missing_run_count_total > 0,
    "top_feature_exact_replicates_partly_missing",
    default = ""
  )]
  wide[nzchar(replicate_measurement_flag),
       .(source_row, replicate_measurement_flag)]
}

append_replicate_measurement_flags <- function(manifest, replicate_summary) {
  flag_dt <- replicate_measurement_flags(replicate_summary)
  if (!nrow(flag_dt) || !nrow(manifest)) return(manifest)
  for (i in seq_len(nrow(flag_dt))) {
    row_id <- flag_dt$source_row[[i]]
    flag <- flag_dt$replicate_measurement_flag[[i]]
    current <- manifest[source_row == row_id, flags]
    current <- if (length(current) && !is.na(current[[1]]) && nzchar(current[[1]])) {
      strsplit(current[[1]], ";", fixed = TRUE)[[1]]
    } else {
      character()
    }
    manifest[source_row == row_id, flags := paste(unique(c(current, flag)),
                                                  collapse = ";")]
  }
  manifest[]
}

append_replicate_effect_flags <- function(manifest, replicate_effects) {
  if (!nrow(replicate_effects) || !nrow(manifest)) return(manifest)
  effect_cols <- c("source_row", "replicate_mean_prop_delta")
  if (!all(effect_cols %chin% names(replicate_effects)) ||
      !"exact_top_feature_prop_delta" %chin% names(manifest)) {
    return(manifest)
  }
  effect_dt <- replicate_effects[, ..effect_cols]
  dt <- effect_dt[manifest, on = "source_row"]
  dt[, aggregate_delta := suppressWarnings(as.numeric(exact_top_feature_prop_delta))]
  dt[, replicate_delta := suppressWarnings(as.numeric(replicate_mean_prop_delta))]
  dt[, direction_discordant :=
       is.finite(aggregate_delta) & is.finite(replicate_delta) &
       abs(aggregate_delta) >= 1e-9 & abs(replicate_delta) >= 1e-9 &
       sign(aggregate_delta) != sign(replicate_delta)]
  discordant_rows <- dt[direction_discordant == TRUE, source_row]
  if (!length(discordant_rows)) return(manifest)
  for (row_id in discordant_rows) {
    current <- manifest[source_row == row_id, flags]
    current <- if (length(current) && !is.na(current[[1]]) && nzchar(current[[1]])) {
      strsplit(current[[1]], ";", fixed = TRUE)[[1]]
    } else {
      character()
    }
    manifest[source_row == row_id, flags := paste(
      unique(c(current, "aggregate_replicate_top_feature_direction_discordant")),
      collapse = ";"
    )]
  }
  manifest[]
}

plot_replicate_diagnostics <- function(replicate_dt, row, top_feature,
                                       output_file) {
  if (!nrow(replicate_dt)) return(FALSE)
  dt <- copy(replicate_dt[feature_role %chin% c("top_feature", "clean_CDS")])
  if (!nrow(dt)) return(FALSE)
  dt[, raw_counts_plot := as.numeric(raw_counts)]
  dt[, raw_counts_missing := !is.finite(raw_counts_plot)]
  dt[, run_axis := paste(group_role, run_order, Run, sep = "\n")]
  dt[, run_axis := factor(run_axis, levels = unique(run_axis))]
  dt[, feature_label := fifelse(
    feature_role == "top_feature",
    paste0("top feature: ", feature_id, " / ", feature_type),
    paste0(feature_id, " / ", feature_type)
  )]
  plot_dt <- melt(
    dt,
    id.vars = c("source_row", "group_id", "group_role", "Run", "run_order",
                "run_axis", "feature_id", "feature_type", "feature_label",
                "feature_role", "raw_counts_missing"),
    measure.vars = c("raw_counts_plot", "run_feature_prop_all"),
    variable.name = "metric",
    value.name = "value"
  )
  plot_dt[, value_missing := !is.finite(value)]
  plot_dt[!is.finite(value), value := 0]
  plot_dt[, metric := factor(
    metric,
    levels = c("raw_counts_plot", "run_feature_prop_all"),
    labels = c("Raw feature counts", "Share of modeled feature counts")
  )]
  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  tx_id <- coalesce_col(row, c("tx_id"), default = "")
  title <- paste(gene, tx_id, "exact replicate diagnostics", sep = " | ")
  subtitle <- paste0(
    "Exact case/control per-run counts and proportions; expected top feature: ",
    top_feature,
    "; x marks missing values"
  )

  role_cols <- c(case = "#C44E52", control = "#4C78A8")
  p <- ggplot(plot_dt, aes(run_axis, value, fill = group_role)) +
    geom_col(width = 0.82, color = "grey35", linewidth = 0.12) +
    geom_point(
      data = plot_dt[value_missing == TRUE],
      aes(run_axis, 0),
      inherit.aes = FALSE,
      shape = 4,
      color = "black",
      size = 1.7,
      stroke = 0.45
    ) +
    facet_grid(metric ~ feature_label, scales = "free_y") +
    scale_fill_manual(values = role_cols, drop = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Exact selected run",
      y = NULL,
      fill = "Group"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8),
      axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 6),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  ggsave(output_file, p,
         width = max(8.5, min(16, uniqueN(dt$Run) * 0.42 + 4.2)),
         height = 6.2,
         dpi = 135,
         limitsize = FALSE)
  TRUE
}

plot_rdg_flow <- function(feature_annotation, edges, row, tx_id,
                          leader_extension, display_width, top_feature,
                          output_file) {
  keep_fx <- equal_filter(feature_annotation, "tx_id", tx_id)
  fx <- copy(feature_annotation[keep_fx])
  if (!nrow(fx)) return(FALSE)
  fx[, `:=`(
    xmin = as.numeric(tx_start) + leader_extension,
    xmax = as.numeric(tx_end) + leader_extension,
    is_top_feature = feature_id == top_feature
  )]
  fx <- fx[is.finite(xmin) & is.finite(xmax) & xmax >= 1 & xmin <= display_width]
  fx[, `:=`(
    xmin = pmax(1, xmin),
    xmax = pmin(display_width, xmax)
  )]
  keep_edges <- equal_filter(edges, "tx_id", tx_id)
  ed <- copy(edges[keep_edges])
  if (nrow(ed)) {
    ed[, `:=`(
      x = as.numeric(x) + leader_extension,
      xend = as.numeric(xend) + leader_extension
    )]
    ed <- ed[is.finite(x) & is.finite(xend) & pmax(x, xend) >= 1 &
               pmin(x, xend) <= display_width]
    ed[, `:=`(
      x = pmax(1, pmin(display_width, x)),
      xend = pmax(1, pmin(display_width, xend))
    )]
  }
  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  design <- coalesce_col(row, c("design_family"), default = "")
  title <- paste(gene, tx_id, design, sep = " | ")
  branch <- coalesce_col(row, c("joint_top_branch_class", "main_branch_class"),
                         default = "")

  feature_cols <- c(
    clean_CDS = "#4C78A8",
    leader_uorf = "#F58518",
    overlapping_uorf = "#E45756",
    internal_orf = "#54A24B",
    nte_extension = "#B279A2",
    ntt_truncation = "#9D755D",
    cds_translon = "#BAB0AC",
    other_translon = "#BAB0AC"
  )
  edge_cols <- c(
    scanning = "#666666",
    translon_translation = "#E45756",
    cds_translation = "#4C78A8",
    posttermination_scanning = "#59A14F",
    competitive_exit = "#B279A2"
  )

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey78", linewidth = 0.3) +
    geom_segment(
      data = ed,
      aes(x = x, xend = xend, y = y, yend = yend, color = edge_type),
      linewidth = 0.55,
      alpha = 0.9,
      lineend = "round"
    ) +
    geom_rect(
      data = fx,
      aes(xmin = xmin, xmax = xmax, ymin = track_y - 0.07,
          ymax = track_y + 0.07, fill = feature_type),
      color = "grey20",
      linewidth = 0.22,
      alpha = 0.88
    ) +
    geom_rect(
      data = fx[is_top_feature == TRUE],
      aes(xmin = xmin, xmax = xmax, ymin = track_y - 0.10,
          ymax = track_y + 0.10),
      fill = NA,
      color = "black",
      linewidth = 0.7
    ) +
    geom_text(
      data = fx,
      aes(x = (xmin + xmax) / 2, y = track_y + 0.16,
          label = paste0(feature_id, "\n", start_codon, "/", stop_codon)),
      size = 2.5,
      lineheight = 0.82,
      check_overlap = TRUE
    ) +
    scale_fill_manual(values = feature_cols, drop = FALSE) +
    scale_color_manual(values = edge_cols, drop = FALSE) +
    labs(
      title = title,
      subtitle = paste0("RDG topology; expected top branch: ", branch,
                        if (nzchar(top_feature)) paste0("; top feature: ", top_feature) else ""),
      x = "Displayed transcript coordinate",
      y = "RDG track",
      fill = "Feature",
      color = "Edge"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  ggsave(output_file, p, width = 12.5, height = 4.8, dpi = 135,
         limitsize = FALSE)
  TRUE
}

sequence_feature_context <- function(df, tx_region, display_region, row_index,
                                     row, feature_annotation, tx_id) {
  gene <- coalesce_col(row, c("gene_symbol"), default = "")
  fx <- copy(feature_annotation[
    feature_annotation[["gene_symbol"]] == gene &
      feature_annotation[["tx_id"]] == tx_id
  ])
  if (!nrow(fx)) return(data.table())
  fa <- tryCatch(ORFik::findFa(df), error = function(e) NULL)
  mrna_seq <- tryCatch({
    if (is.null(fa)) stop("No FASTA found")
    as.character(RiboCrypt:::browser_tx_seqs_getSeq(fa, tx_region,
                                                    keep.names = TRUE))
  }, error = function(e) NA_character_)
  display_seq <- tryCatch({
    if (is.null(fa)) stop("No FASTA found")
    as.character(RiboCrypt:::browser_tx_seqs_getSeq(fa, display_region,
                                                    keep.names = TRUE))
  }, error = function(e) NA_character_)
  mrna_len <- if (!is.na(mrna_seq) && nzchar(mrna_seq)) nchar(mrna_seq) else NA_integer_
  display_len <- if (!is.na(display_seq) && nzchar(display_seq)) nchar(display_seq) else NA_integer_
  if (is.finite(mrna_len)) {
    fx[, start_context_sequence := vapply(tx_start, function(pos) {
      pos <- as.integer(pos)
      if (!is.finite(pos) || pos < 1L || pos > mrna_len) return(NA_character_)
      substr(mrna_seq, max(1L, pos - 6L), min(mrna_len, pos + 8L))
    }, character(1))]
  } else {
    fx[, start_context_sequence := NA_character_]
  }
  fx[, `:=`(
    source_row = row_index,
    mrna_sequence_available = is.finite(mrna_len),
    mrna_sequence_length = mrna_len,
    display_sequence_available = is.finite(display_len),
    display_sequence_length = display_len
  )]
  fx[, .(
    source_row,
    gene_symbol,
    tx_id,
    feature_id,
    feature_type,
    category,
    tx_start,
    tx_end,
    feature_bases,
    cds_overlap_bases,
    start_codon,
    stop_codon,
    start_codon_class,
    frame_relative_to_cds,
    peptide_length_aa,
    kozak_minus3,
    kozak_plus4,
    kozak_strength,
    coding_potential_proxy,
    start_context_sequence,
    mrna_sequence_available,
    mrna_sequence_length,
    display_sequence_available,
    display_sequence_length
  )]
}

make_review_flags <- function(row, row_index, feature_counts, contrasts,
                              profile_metric_dt, sequence_dt, top_feature) {
  branch_case <- first_finite_num(row, c("case_branch_total_counts"))
  branch_control <- first_finite_num(row, c("control_branch_total_counts"))
  branch_min <- suppressWarnings(min(branch_case, branch_control, na.rm = TRUE))
  if (!is.finite(branch_min)) branch_min <- NA_real_

  exact_profile <- profile_metric_dt[group_id %chin% c("exact_case", "exact_control")]
  display_min <- if (nrow(exact_profile)) {
    min(exact_profile$display_total_counts, na.rm = TRUE)
  } else {
    NA_real_
  }
  if (!is.finite(display_min)) display_min <- NA_real_
  max_spike <- if (nrow(exact_profile)) {
    max_or_na(exact_profile$max_position_fraction)
  } else {
    NA_real_
  }
  top3_spike <- if (nrow(exact_profile)) {
    max_or_na(exact_profile$top3_position_fraction)
  } else {
    NA_real_
  }

  top_fc <- feature_counts[feature_id == top_feature &
                             group_id %chin% c("exact_case", "exact_control")]
  top_min <- if (nrow(top_fc)) min(top_fc$raw_counts_sum, na.rm = TRUE) else NA_real_
  if (!is.finite(top_min)) top_min <- NA_real_
  top_max_run_fraction <- if (nrow(top_fc)) {
    max_or_na(top_fc$max_run_fraction)
  } else {
    NA_real_
  }

  top_contrast <- contrasts[contrast_id == "exact_case_vs_exact_control" &
                              feature_id == top_feature]
  exact_top_prop_delta <- if (nrow(top_contrast)) {
    top_contrast$prop_delta[[1]]
  } else {
    NA_real_
  }
  exact_top_log2 <- if (nrow(top_contrast)) {
    top_contrast$raw_log2_ratio[[1]]
  } else {
    NA_real_
  }

  seq_ok <- nrow(sequence_dt) > 0 &&
    any(sequence_dt$mrna_sequence_available == TRUE, na.rm = TRUE)
  top_feature_in_rdg <- nrow(sequence_dt[feature_id == top_feature]) > 0
  has_internal_or_nte <- nrow(sequence_dt[
    feature_type %chin% c("internal_orf", "nte_extension", "ntt_truncation")
  ]) > 0

  flags <- character()
  if (is.na(branch_min)) {
    flags <- c(flags, "branch_count_missing")
  } else if (branch_min < 30) {
    flags <- c(flags, "branch_count_below_30")
  } else if (branch_min < 100) {
    flags <- c(flags, "branch_count_30_to_99")
  } else if (branch_min < 300) {
    flags <- c(flags, "branch_count_100_to_299")
  }
  if (is.na(display_min)) {
    flags <- c(flags, "display_coverage_missing")
  } else if (display_min < 30) {
    flags <- c(flags, "display_coverage_below_30")
  } else if (display_min < 100) {
    flags <- c(flags, "display_coverage_30_to_99")
  } else if (display_min < 300) {
    flags <- c(flags, "display_coverage_100_to_299")
  }
  if (is.na(top_min)) {
    flags <- c(flags, "top_feature_count_missing")
  } else if (top_min < 30) {
    flags <- c(flags, "top_feature_count_below_30")
  } else if (top_min < 100) {
    flags <- c(flags, "top_feature_count_30_to_99")
  } else if (top_min < 300) {
    flags <- c(flags, "top_feature_count_100_to_299")
  }
  if (is.finite(max_spike) && max_spike >= 0.45) {
    flags <- c(flags, "single_position_spike_risk")
  }
  if (is.finite(top3_spike) && top3_spike >= 0.65) {
    flags <- c(flags, "top3_position_spike_risk")
  }
  if (is.finite(top_max_run_fraction) && top_max_run_fraction >= 0.70) {
    flags <- c(flags, "top_feature_single_run_dominated")
  }
  if (!seq_ok) flags <- c(flags, "mrna_sequence_unavailable")
  if (nzchar(top_feature) && !top_feature_in_rdg) {
    flags <- c(flags, "top_feature_not_in_rdg_annotation")
  }
  if (has_internal_or_nte) flags <- c(flags, "noncanonical_rdg_architecture_present")

  readiness <- if ("branch_count_missing" %chin% flags ||
                   "display_coverage_missing" %chin% flags ||
                   "top_feature_count_missing" %chin% flags ||
                   "mrna_sequence_unavailable" %chin% flags) {
    "incomplete_inputs"
  } else if (any(c("branch_count_below_30", "display_coverage_below_30",
                   "top_feature_count_below_30") %chin% flags)) {
    "not_claimable_from_coverage"
  } else if (any(c("single_position_spike_risk", "top3_position_spike_risk",
                   "top_feature_single_run_dominated") %chin% flags)) {
    "spike_risk_review_required"
  } else if (is.finite(branch_min) && branch_min >= 300 &&
             is.finite(display_min) && display_min >= 300 &&
             is.finite(top_min) && top_min >= 100) {
    "strong_visual_review_ready"
  } else if (is.finite(branch_min) && branch_min >= 100 &&
             is.finite(display_min) && display_min >= 100 &&
             is.finite(top_min) && top_min >= 30) {
    "qualified_visual_review_ready"
  } else {
    "weak_count_visual_review_only"
  }

  data.table(
    source_row = row_index,
    branch_min_case_control_counts = branch_min,
    display_min_exact_counts = display_min,
    top_feature_min_exact_counts = top_min,
    max_exact_position_fraction = max_spike,
    max_exact_top3_position_fraction = top3_spike,
    top_feature_max_run_fraction = top_max_run_fraction,
    exact_top_feature_prop_delta = exact_top_prop_delta,
    exact_top_feature_raw_log2_ratio = exact_top_log2,
    mrna_sequence_available = seq_ok,
    top_feature_in_rdg_annotation = top_feature_in_rdg,
    has_internal_or_nte_architecture = has_internal_or_nte,
    inspection_readiness = readiness,
    flags = paste(unique(flags), collapse = ";")
  )
}

render_qualified_row <- function(row, row_index, df_cache, metadata,
                                 feature_expression, feature_annotation,
                                 rdg_edges) {
  url <- coalesce_col(row, c("observatory_url"), default = "")
  state <- parse_observatory_state(url)
  if (is.null(state)) stop("Could not parse observatory_url")

  tx_id <- coalesce_col(row, c("tx_id"), default = state$browser$tx %||% "")
  gene <- coalesce_col(row, c("gene_symbol"), default = state$browser$gene %||% "")
  if (!nzchar(tx_id)) stop("Missing transcript id")

  exp_name <- scalar_chr(state$exp, "all_samples-Homo_sapiens")
  if (is.null(df_cache[[exp_name]])) {
    df_cache[[exp_name]] <- ORFik::read.experiment(exp_name, validate = FALSE)
  }
  df <- df_cache[[exp_name]]

  tx_region <- ORFik::loadRegion(df, "mrna", names.keep = tx_id)
  if (!length(tx_region)) stop("Transcript not found in experiment: ", tx_id)
  display_region <- display_region_from_state(tx_region, state$browser)
  display_width <- as.integer(ORFik::widthPerGroup(display_region)[[1]])
  leader_extension <- suppressWarnings(as.numeric(state$browser$extendLeaders %||% 0))
  if (!is.finite(leader_extension)) leader_extension <- 0

  collection_path <- RiboCrypt:::collection_path_from_exp(
    df, tx_id, grl_all = tx_region
  )
  available_runs <- collection_available_runs(collection_path)
  if (!length(available_runs)) stop("No coverage runs available for ", tx_id)

  groups <- build_context_groups(row, row_index, state, metadata, available_runs)
  if (!nrow(groups)) stop("No context groups have available runs")
  groups <- choose_profile_groups(groups, display_width, max_profile_cells)

  top_feature <- coalesce_col(
    row,
    c("top_shifted_feature_id", "branch_usage_top_feature_id",
      "grouped_branch_usage_top_feature_id", "count_aware_top_feature_id"),
    default = ""
  )

  profile_dt <- data.table()
  metric_dt <- data.table()
  profile_file <- ""
  profile_status <- "not_requested"
  plot_groups <- groups[profile_include == TRUE]
  if (nrow(plot_groups)) {
    profile_runs <- unique(unlist(plot_groups$runs, use.names = FALSE))
    reads <- RiboCrypt:::load_collection(
      collection_path,
      grl = display_region,
      columns = profile_runs
    )
    profile_dt <- profile_groups_from_reads(reads, plot_groups)
    metric_dt <- profile_metrics(profile_dt)
    track <- feature_track(feature_annotation, rdg_edges, tx_id,
                           leader_extension, display_width, profile_dt)
    row_rank <- coalesce_col(row, c("global_review_order", "label_gap_rank"),
                             default = row_index)
    file_stem <- paste(
      sprintf("%03d", as.integer(row_index)),
      safe_name(row_rank),
      safe_name(gene),
      safe_name(tx_id),
      sep = "_"
    )
    profile_file <- file.path(profile_dir, paste0(file_stem, "_profile.png"))
    if (force || !file.exists(profile_file)) {
      plot_profile(profile_dt, track$features, row, top_feature, profile_file,
                   max_plot_positions)
    }
    profile_status <- "ok"
  }

  feature_counts <- summarize_feature_counts(row, row_index, groups,
                                             feature_expression,
                                             feature_annotation)
  contrasts <- feature_contrasts(feature_counts, row_index)

  row_rank <- coalesce_col(row, c("global_review_order", "label_gap_rank"),
                           default = row_index)
  file_stem <- paste(
    sprintf("%03d", as.integer(row_index)),
    safe_name(row_rank),
    safe_name(gene),
    safe_name(tx_id),
    sep = "_"
  )
  feature_file <- file.path(feature_dir, paste0(file_stem, "_features.png"))
  if (nrow(feature_counts) && (force || !file.exists(feature_file))) {
    plot_feature_allocation(feature_counts, row, top_feature, feature_file)
  }

  replicate_dt <- summarize_replicate_diagnostics(
    row, row_index, groups, feature_expression, metadata, top_feature
  )
  replicate_summary <- summarize_replicate_groups(replicate_dt, top_feature)
  replicate_file <- file.path(replicate_dir,
                              paste0(file_stem, "_replicates.png"))
  if (nrow(replicate_dt) && (force || !file.exists(replicate_file))) {
    plot_replicate_diagnostics(replicate_dt, row, top_feature, replicate_file)
  }

  flow_file <- file.path(flow_dir, paste0(file_stem, "_rdg_flow.png"))
  if (force || !file.exists(flow_file)) {
    plot_rdg_flow(feature_annotation, rdg_edges, row, tx_id, leader_extension,
                  display_width, top_feature, flow_file)
  }

  seq_dt <- sequence_feature_context(df, tx_region, display_region, row_index,
                                     row, feature_annotation, tx_id)
  flags <- make_review_flags(row, row_index, feature_counts, contrasts,
                             metric_dt, seq_dt, top_feature)

  exact_control_total <- metric_dt[group_id == "exact_control",
                                   sum(display_total_counts, na.rm = TRUE)]
  exact_case_total <- metric_dt[group_id == "exact_case",
                                sum(display_total_counts, na.rm = TRUE)]
  if (!is.finite(exact_control_total)) exact_control_total <- NA_real_
  if (!is.finite(exact_case_total)) exact_case_total <- NA_real_

  manifest <- data.table(
    source_row = row_index,
    global_review_order = coalesce_col(row, c("global_review_order"), default = ""),
    batch_id = coalesce_col(row, c("batch_id"), default = ""),
    batch_role = coalesce_col(row, c("batch_role"), default = ""),
    gene_symbol = gene,
    tx_id = tx_id,
    design_family = coalesce_col(row, c("design_family"), default = ""),
    experiment = exp_name,
    study = coalesce_col(row, c("study"), default = ""),
    tissue = coalesce_col(row, c("TISSUE"), default = ""),
    cell_line = coalesce_col(row, c("CELL_LINE"), default = ""),
    case_condition = coalesce_col(row, c("case_condition"), default = ""),
    control_conditions = coalesce_col(row, c("control_conditions"), default = ""),
    display_width = display_width,
    leader_extension = leader_extension,
    available_run_count = length(available_runs),
    context_group_count = nrow(groups),
    profile_group_count = nrow(plot_groups),
    profile_status = profile_status,
    exact_control_display_total_counts = exact_control_total,
    exact_case_display_total_counts = exact_case_total,
    observed_exact_log2_case_control_display = log2((exact_case_total + 1) /
                                                      (exact_control_total + 1)),
    case_branch_total_counts = first_finite_num(row, c("case_branch_total_counts")),
    control_branch_total_counts = first_finite_num(row, c("control_branch_total_counts")),
    context_count_gate = coalesce_col(row, c("context_count_gate"), default = ""),
    replicate_gate = coalesce_col(row, c("replicate_gate"), default = ""),
    top_branch_class = coalesce_col(row, c("joint_top_branch_class",
                                           "main_branch_class"), default = ""),
    top_branch_delta = first_finite_num(row, c("joint_top_branch_delta",
                                               "main_branch_delta")),
    top_feature_id = top_feature,
    top_feature_type = coalesce_col(row, c("top_shifted_feature_type"), default = ""),
    transcript_qc_class = coalesce_col(row, c("transcript_qc_class"), default = ""),
    cds_exon_qc_status = coalesce_col(row, c("cds_exon_qc_status"), default = ""),
    expected_browser_pattern = coalesce_col(row, c("expected_browser_pattern"),
                                            default = ""),
    review_question = coalesce_col(row, c("review_question"), default = ""),
    profile_png = profile_file,
    feature_allocation_png = feature_file,
    replicate_diagnostics_png = if (nrow(replicate_dt)) replicate_file else "",
    rdg_flow_png = flow_file,
    status = "ok",
    error = ""
  )
  manifest <- cbind(manifest, flags[, !"source_row"])
  list(
    manifest = manifest,
    groups = groups,
    profile_metrics = metric_dt,
    feature_counts = feature_counts,
    feature_contrasts = contrasts,
    replicate_diagnostics = replicate_dt,
    replicate_summary = replicate_summary,
    sequence_features = seq_dt,
    flags = flags
  )
}

if (!file.exists(queue_file)) {
  stop("Missing queue file: ", queue_file, call. = FALSE)
}
if (!file.exists(metadata_file)) {
  stop("Missing metadata file: ", metadata_file, call. = FALSE)
}
if (!file.exists(feature_expression_file)) {
  stop("Missing feature expression file: ", feature_expression_file,
       call. = FALSE)
}
if (!file.exists(feature_annotation_file)) {
  stop("Missing RDG feature annotation file: ", feature_annotation_file,
       call. = FALSE)
}
if (!file.exists(rdg_edges_file)) {
  stop("Missing RDG edge file: ", rdg_edges_file, call. = FALSE)
}

queue <- fread(queue_file, showProgress = FALSE)
if (!"observatory_url" %chin% names(queue)) {
  stop("Queue must contain observatory_url: ", queue_file, call. = FALSE)
}
if (is.finite(row_n)) {
  queue <- queue[row_n]
} else if (is.finite(limit_n)) {
  queue <- head(queue, limit_n)
}
if (!nrow(queue)) stop("No rows selected from queue: ", queue_file,
                       call. = FALSE)

message("Dominant RDG qualified inspection")
message("  queue: ", queue_file)
message("  rows: ", nrow(queue))
message("  output: ", output_dir)
message("  max profile cells: ", max_profile_cells)

metadata <- fread(metadata_file, showProgress = FALSE)
feature_annotation <- fread(feature_annotation_file, showProgress = FALSE)
rdg_edges <- fread(rdg_edges_file, showProgress = FALSE)

feature_expression_cols <- c(
  "gene_symbol", "tx_id", "feature_id", "feature_type", "category",
  "translon_source", "translon_sources_merged", "feature_bases",
  "cds_overlap_bases", "measured_bases", "measurement_scope", "Run",
  "raw_counts", "fpkm_like", "cds_reference_raw_counts",
  "internal_start_window_raw_counts", "internal_feature_core_raw_counts",
  "internal_upstream_window_raw_counts", "internal_start_vs_body_density_ratio",
  "internal_start_vs_upstream_density_ratio"
)
available_cols <- names(fread(feature_expression_file, nrows = 0,
                              showProgress = FALSE))
feature_expression_cols <- intersect(feature_expression_cols, available_cols)
feature_expression <- fread(feature_expression_file,
                            select = feature_expression_cols,
                            showProgress = FALSE)
setkey(feature_expression, gene_symbol, tx_id)

df_cache <- new.env(parent = emptyenv())
results <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  message("[", i, "/", nrow(queue), "] ",
          coalesce_col(queue[i], c("gene_symbol"), default = "unknown"),
          " ", coalesce_col(queue[i], c("tx_id"), default = ""))
  if (stop_on_error) {
    results[[i]] <- render_qualified_row(
      queue[i], i, df_cache, metadata, feature_expression,
      feature_annotation, rdg_edges
    )
  } else {
    results[[i]] <- tryCatch(
      render_qualified_row(
        queue[i], i, df_cache, metadata, feature_expression,
        feature_annotation, rdg_edges
      ),
      error = function(e) {
        gene <- coalesce_col(queue[i], c("gene_symbol"), default = "")
        tx_id <- coalesce_col(queue[i], c("tx_id"), default = "")
        manifest <- data.table(
          source_row = i,
          global_review_order = coalesce_col(queue[i], c("global_review_order"), default = ""),
          batch_id = coalesce_col(queue[i], c("batch_id"), default = ""),
          batch_role = coalesce_col(queue[i], c("batch_role"), default = ""),
          gene_symbol = gene,
          tx_id = tx_id,
          design_family = coalesce_col(queue[i], c("design_family"), default = ""),
          status = "error",
          error = conditionMessage(e)
        )
        list(
          manifest = manifest,
          groups = data.table(),
          profile_metrics = data.table(),
          feature_counts = data.table(),
          feature_contrasts = data.table(),
          replicate_diagnostics = data.table(),
          replicate_summary = data.table(),
          sequence_features = data.table(),
          flags = data.table()
        )
      }
    )
  }
}

run_set_signature <- function(runs) {
  runs <- sort(unique(as.character(runs)))
  runs <- runs[!is.na(runs) & nzchar(runs)]
  paste(runs, collapse = ";")
}

context_run_set_tables <- function(groups) {
  if (!nrow(groups) || !all(c("runs", "requested_runs") %chin% names(groups))) {
    return(list(
      groups = groups,
      run_sets = data.table(),
      run_set_members = data.table()
    ))
  }

  run_sets_raw <- rbindlist(lapply(c("runs", "requested_runs"), function(col) {
    data.table(
      source_column = col,
      source_group_row = seq_len(nrow(groups)),
      signature = vapply(groups[[col]], run_set_signature, character(1)),
      runs = groups[[col]]
    )
  }), fill = TRUE)
  run_sets_raw <- run_sets_raw[nzchar(signature)]
  unique_sets <- run_sets_raw[!duplicated(signature)]
  unique_sets[, run_set_id := sprintf("RS%05d", seq_len(.N))]
  unique_sets[, run_count := lengths(runs)]

  members <- rbindlist(lapply(seq_len(nrow(unique_sets)), function(i) {
    runs <- unique_sets$runs[[i]]
    data.table(
      run_set_id = unique_sets$run_set_id[[i]],
      run_index = seq_along(runs),
      Run = runs
    )
  }), fill = TRUE)

  map <- unique_sets[, .(signature, run_set_id)]
  groups_out <- copy(groups)
  run_sig <- vapply(groups_out$runs, run_set_signature, character(1))
  requested_sig <- vapply(groups_out$requested_runs, run_set_signature,
                          character(1))
  groups_out[, `:=`(
    run_set_id = map$run_set_id[match(run_sig, map$signature)],
    requested_run_set_id = map$run_set_id[match(requested_sig, map$signature)]
  )]
  groups_out[, c("runs", "requested_runs", "run_signature") := NULL]

  list(
    groups = groups_out,
    run_sets = unique_sets[, .(run_set_id, run_count)],
    run_set_members = members
  )
}

make_review_priority <- function(manifest, replicate_summary,
                                 replicate_effects = data.table()) {
  if (!nrow(manifest)) return(data.table())
  dt <- copy(manifest[status == "ok"])
  if (!nrow(dt)) return(data.table())
  readiness_rank <- c(
    strong_visual_review_ready = 1L,
    qualified_visual_review_ready = 2L,
    spike_risk_review_required = 3L,
    not_claimable_from_coverage = 4L,
    incomplete_inputs = 5L,
    weak_count_visual_review_only = 6L
  )
  role_rank <- c(
    candidate_positive = 1L,
    counterfactual_review = 2L,
    context_audit = 3L,
    replication_followup = 4L,
    measurement_review = 5L,
    negative_control = 6L
  )
  dt[, `:=`(
    readiness_rank = readiness_rank[inspection_readiness],
    batch_role_rank = role_rank[batch_role]
  )]
  dt[is.na(readiness_rank), readiness_rank := 99L]
  dt[is.na(batch_role_rank), batch_role_rank := 99L]
  dt[, flag_text := fifelse(is.na(flags), "", as.character(flags))]
  dt[, recommended_action := fcase(
    grepl("top_feature_exact_replicates_all_missing", flag_text),
    "audit_missing_top_feature_measurement",
    inspection_readiness %chin% c("strong_visual_review_ready",
                                  "qualified_visual_review_ready",
                                  "spike_risk_review_required") &
      grepl("top_feature_exact_replicates_", flag_text),
    "audit_missing_replicate_measurements_before_claim",
    inspection_readiness %chin% c("strong_visual_review_ready",
                                  "qualified_visual_review_ready",
                                  "spike_risk_review_required") &
      grepl("aggregate_replicate_top_feature_direction_discordant", flag_text),
    "inspect_replicate_effect_direction_before_claim",
    inspection_readiness == "strong_visual_review_ready",
    "manual_mechanism_review",
    inspection_readiness == "qualified_visual_review_ready",
    "manual_review_with_count_caution",
    inspection_readiness == "spike_risk_review_required",
    "inspect_replicate_plot_before_claim",
    inspection_readiness == "not_claimable_from_coverage",
    "do_not_claim_low_coverage",
    inspection_readiness == "incomplete_inputs",
    "needs_branch_count_fields_or_audit_only",
    default = "review_case_by_case"
  )]

  if (nrow(replicate_summary)) {
    rep_wide <- dcast(
      replicate_summary,
      source_row ~ group_id,
      value.var = c("exact_signal_run_count", "exact_run_count",
                    "exact_measured_run_count", "exact_missing_run_count",
                    "exact_max_run_fraction",
                    "exact_total_top_feature_counts"),
      fill = NA
    )
    dt <- rep_wide[dt, on = "source_row"]
  }
  if (nrow(replicate_effects)) {
    effect_cols <- c(
      "source_row",
      "replicate_mean_prop_delta",
      "replicate_median_prop_delta",
      "replicate_mean_allocation_prop_delta",
      "replicate_median_allocation_prop_delta",
      "replicate_mean_raw_log2_ratio",
      "replicate_mean_prop_all_exact_case",
      "replicate_mean_prop_all_exact_control",
      "replicate_median_prop_all_exact_case",
      "replicate_median_prop_all_exact_control"
    )
    effect_cols <- intersect(effect_cols, names(replicate_effects))
    dt <- replicate_effects[, ..effect_cols][dt, on = "source_row"]
  }

  setorder(dt, readiness_rank, batch_role_rank,
           -display_min_exact_counts, -top_feature_min_exact_counts)
  dt[, review_priority_rank := seq_len(.N)]
  keep_cols <- intersect(c(
    "review_priority_rank", "recommended_action", "source_row",
    "global_review_order", "batch_id", "batch_role", "gene_symbol", "tx_id",
    "design_family", "study", "tissue", "cell_line", "case_condition",
    "control_conditions", "inspection_readiness", "flags", "top_branch_class",
    "top_branch_delta", "top_feature_id", "display_min_exact_counts",
    "top_feature_min_exact_counts", "exact_top_feature_prop_delta",
    "replicate_mean_prop_delta", "replicate_median_prop_delta",
    "replicate_mean_allocation_prop_delta",
    "replicate_median_allocation_prop_delta",
    "replicate_mean_raw_log2_ratio",
    "top_feature_max_run_fraction",
    "replicate_mean_prop_all_exact_case",
    "replicate_mean_prop_all_exact_control",
    "replicate_median_prop_all_exact_case",
    "replicate_median_prop_all_exact_control",
    "exact_signal_run_count_exact_case", "exact_run_count_exact_case",
    "exact_measured_run_count_exact_case",
    "exact_missing_run_count_exact_case",
    "exact_max_run_fraction_exact_case",
    "exact_total_top_feature_counts_exact_case",
    "exact_signal_run_count_exact_control", "exact_run_count_exact_control",
    "exact_measured_run_count_exact_control",
    "exact_missing_run_count_exact_control",
    "exact_max_run_fraction_exact_control",
    "exact_total_top_feature_counts_exact_control",
    "profile_png", "feature_allocation_png", "replicate_diagnostics_png",
    "rdg_flow_png"
  ), names(dt))
  dt[, ..keep_cols]
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

format_review_value <- function(x, digits = 3L) {
  if (!length(x) || is.na(x[[1]])) return("")
  if (is.numeric(x) || is.integer(x)) {
    value <- as.numeric(x[[1]])
    if (!is.finite(value)) return("")
    if (abs(value) >= 100 || identical(value, round(value))) {
      return(format(round(value), trim = TRUE, scientific = FALSE))
    }
    return(format(round(value, digits), trim = TRUE, scientific = FALSE))
  }
  as.character(x[[1]])
}

html_path <- function(path, output_dir) {
  path <- as.character(path)
  if (!length(path) || is.na(path[[1]]) || !nzchar(path[[1]])) return("")
  path <- path[[1]]
  out_norm <- normalizePath(output_dir, mustWork = FALSE)
  path_norm <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(out_norm, .Platform$file.sep)
  if (startsWith(path_norm, prefix)) {
    path <- substring(path_norm, nchar(prefix) + 1L)
  } else {
    path <- path_norm
  }
  path <- gsub("\\\\", "/", path)
  utils::URLencode(path, reserved = FALSE)
}

metric_table <- function(row, fields) {
  rows <- vapply(fields, function(field) {
    if (!field %chin% names(row)) return("")
    value <- format_review_value(row[[field]])
    if (!nzchar(value)) return("")
    paste0(
      "<tr><th>", html_escape(field), "</th><td>",
      html_escape(value), "</td></tr>"
    )
  }, character(1))
  rows <- rows[nzchar(rows)]
  if (!length(rows)) return("")
  paste0("<table class=\"metrics\">", paste(rows, collapse = "\n"), "</table>")
}

image_figure <- function(row, field, label, output_dir) {
  if (!field %chin% names(row)) return("")
  src <- html_path(row[[field]], output_dir)
  if (!nzchar(src)) return("")
  paste0(
    "<figure>",
    "<a href=\"", src, "\"><img loading=\"lazy\" src=\"", src,
    "\" alt=\"", html_escape(label), "\"></a>",
    "<figcaption>", html_escape(label), "</figcaption>",
    "</figure>"
  )
}

write_review_index <- function(review_priority, output_dir, output_file) {
  if (!nrow(review_priority)) return(FALSE)
  dt <- copy(review_priority)
  setorder(dt, review_priority_rank)

  actions <- sort(unique(dt$recommended_action[nzchar(dt$recommended_action)]))
  readiness <- sort(unique(dt$inspection_readiness[nzchar(dt$inspection_readiness)]))
  option_tags <- function(values) {
    paste0(
      "<option value=\"", html_escape(values), "\">",
      html_escape(values), "</option>",
      collapse = "\n"
    )
  }

  overview_rows <- vapply(seq_len(nrow(dt)), function(i) {
    row <- dt[i]
    text <- paste(
      row$review_priority_rank, row$source_row, row$gene_symbol, row$tx_id,
      row$recommended_action, row$inspection_readiness, row$flags,
      collapse = " "
    )
    paste0(
      "<tr data-action=\"", html_escape(row$recommended_action),
      "\" data-readiness=\"", html_escape(row$inspection_readiness),
      "\" data-text=\"", html_escape(tolower(text)), "\">",
      "<td><a href=\"#card-", row$review_priority_rank, "\">",
      row$review_priority_rank, "</a></td>",
      "<td>", html_escape(row$recommended_action), "</td>",
      "<td>", html_escape(row$inspection_readiness), "</td>",
      "<td>", html_escape(row$source_row), "</td>",
      "<td>", html_escape(row$gene_symbol), "</td>",
      "<td>", html_escape(row$tx_id), "</td>",
      "<td>", html_escape(row$top_feature_id), "</td>",
      "<td>", html_escape(format_review_value(row$top_feature_min_exact_counts)), "</td>",
      "<td>", html_escape(format_review_value(row$replicate_mean_prop_delta)), "</td>",
      "</tr>"
    )
  }, character(1))

  metric_fields <- c(
    "source_row", "global_review_order", "batch_id", "batch_role",
    "design_family", "study", "tissue", "cell_line", "case_condition",
    "control_conditions", "top_branch_class", "top_branch_delta",
    "top_feature_id", "display_min_exact_counts",
    "top_feature_min_exact_counts", "exact_top_feature_prop_delta",
    "replicate_mean_prop_delta", "replicate_median_prop_delta",
    "replicate_mean_allocation_prop_delta",
    "replicate_mean_raw_log2_ratio", "top_feature_max_run_fraction",
    "exact_signal_run_count_exact_case", "exact_run_count_exact_case",
    "exact_measured_run_count_exact_case",
    "exact_missing_run_count_exact_case",
    "exact_signal_run_count_exact_control", "exact_run_count_exact_control",
    "exact_measured_run_count_exact_control",
    "exact_missing_run_count_exact_control"
  )

  cards <- vapply(seq_len(nrow(dt)), function(i) {
    row <- dt[i]
    text <- paste(
      row$review_priority_rank, row$source_row, row$gene_symbol, row$tx_id,
      row$recommended_action, row$inspection_readiness, row$flags,
      row$study, row$tissue, row$cell_line, row$case_condition,
      row$control_conditions,
      collapse = " "
    )
    images <- paste(c(
      image_figure(row, "profile_png", "Profile", output_dir),
      image_figure(row, "feature_allocation_png", "Feature allocation", output_dir),
      image_figure(row, "replicate_diagnostics_png", "Exact replicate diagnostics", output_dir),
      image_figure(row, "rdg_flow_png", "RDG flow", output_dir)
    ), collapse = "\n")
    paste0(
      "<section id=\"card-", row$review_priority_rank,
      "\" class=\"card\" data-action=\"", html_escape(row$recommended_action),
      "\" data-readiness=\"", html_escape(row$inspection_readiness),
      "\" data-text=\"", html_escape(tolower(text)), "\">",
      "<header><div><h2>", row$review_priority_rank, ". ",
      html_escape(row$gene_symbol), " | ", html_escape(row$tx_id),
      "</h2><p>", html_escape(row$recommended_action), " | ",
      html_escape(row$inspection_readiness), "</p></div>",
      "<a class=\"toplink\" href=\"#top\">Top</a></header>",
      "<p class=\"flags\">", html_escape(row$flags), "</p>",
      metric_table(row, metric_fields),
      "<div class=\"figures\">", images, "</div>",
      "</section>"
    )
  }, character(1))

  css <- "
body{font-family:system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;background:#f7f7f5;color:#1f2328}
header.page{position:sticky;top:0;z-index:10;background:#ffffff;border-bottom:1px solid #d8d8d4;padding:12px 18px}
h1{margin:0 0 8px;font-size:20px} h2{margin:0;font-size:18px} p{margin:4px 0}
.controls{display:flex;gap:8px;flex-wrap:wrap}.controls input,.controls select{font-size:13px;padding:6px;border:1px solid #bbb;border-radius:4px;background:white}
main{padding:18px}.overview{width:100%;border-collapse:collapse;margin:0 0 18px;background:white;font-size:12px}
.overview th,.overview td{border:1px solid #ddd;padding:5px;text-align:left}.overview th{background:#ececea}
.card{background:white;border:1px solid #d8d8d4;border-radius:6px;margin:0 0 20px;padding:14px}
.card>header{display:flex;justify-content:space-between;gap:12px;align-items:flex-start;border:0;padding:0;margin:0 0 8px}
.toplink{font-size:12px}.flags{font-size:12px;color:#5b342d;overflow-wrap:anywhere}
.metrics{border-collapse:collapse;font-size:12px;margin:8px 0 12px;max-width:1100px}.metrics th,.metrics td{border:1px solid #ddd;padding:4px 6px;text-align:left}.metrics th{background:#f0f0ee;color:#444}
.figures{display:grid;grid-template-columns:repeat(2,minmax(280px,1fr));gap:12px}.figures figure{margin:0;border:1px solid #ddd;background:#fbfbfa}.figures img{display:block;width:100%;height:auto}.figures figcaption{font-size:12px;padding:5px 7px;background:#ececea}
.hidden{display:none}@media(max-width:900px){.figures{grid-template-columns:1fr}.card{padding:10px}main{padding:10px}}
"
  js <- "
function applyFilters(){
  const q=document.getElementById('q').value.toLowerCase();
  const action=document.getElementById('action').value;
  const readiness=document.getElementById('readiness').value;
  document.querySelectorAll('[data-text]').forEach(el=>{
    const okText=!q || el.dataset.text.includes(q);
    const okAction=!action || el.dataset.action===action;
    const okReadiness=!readiness || el.dataset.readiness===readiness;
    el.classList.toggle('hidden', !(okText && okAction && okReadiness));
  });
}
"
  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Dominant RDG Qualified Review Index</title>",
    "<style>", css, "</style></head>",
    "<body><header id=\"top\" class=\"page\"><h1>Dominant RDG Qualified Review Index</h1>",
    "<div class=\"controls\">",
    "<input id=\"q\" type=\"search\" placeholder=\"Search gene, tx, flag, study\" oninput=\"applyFilters()\">",
    "<select id=\"action\" onchange=\"applyFilters()\"><option value=\"\">All actions</option>",
    option_tags(actions), "</select>",
    "<select id=\"readiness\" onchange=\"applyFilters()\"><option value=\"\">All readiness</option>",
    option_tags(readiness), "</select>",
    "</div></header><main>",
    "<table class=\"overview\"><thead><tr><th>Rank</th><th>Action</th><th>Readiness</th><th>Row</th><th>Gene</th><th>Transcript</th><th>Top feature</th><th>Top min</th><th>Rep mean delta</th></tr></thead><tbody>",
    overview_rows,
    "</tbody></table>",
    cards,
    "</main><script>", js, "</script></body></html>"
  )
  writeLines(html, output_file)
  TRUE
}

manifest <- rbindlist(lapply(results, `[[`, "manifest"), fill = TRUE)
groups <- rbindlist(lapply(results, `[[`, "groups"), fill = TRUE)
profile_metrics_dt <- rbindlist(lapply(results, `[[`, "profile_metrics"),
                                fill = TRUE)
feature_counts_dt <- rbindlist(lapply(results, `[[`, "feature_counts"),
                               fill = TRUE)
feature_contrasts_dt <- rbindlist(lapply(results, `[[`, "feature_contrasts"),
                                  fill = TRUE)
replicate_diagnostics_dt <- rbindlist(
  lapply(results, `[[`, "replicate_diagnostics"),
  fill = TRUE
)
replicate_summary_dt <- rbindlist(lapply(results, `[[`, "replicate_summary"),
                                  fill = TRUE)
sequence_features_dt <- rbindlist(lapply(results, `[[`, "sequence_features"),
                                  fill = TRUE)
flags_dt <- rbindlist(lapply(results, `[[`, "flags"), fill = TRUE)
replicate_effects_dt <- summarize_replicate_effects(replicate_diagnostics_dt)
manifest <- append_replicate_measurement_flags(manifest, replicate_summary_dt)
manifest <- append_replicate_effect_flags(manifest, replicate_effects_dt)
if (nrow(flags_dt)) {
  flag_cols <- intersect(names(flags_dt), names(manifest))
  flags_dt <- manifest[, ..flag_cols]
}

manifest_file <- file.path(output_dir, "qualified_inspection_manifest.csv")
groups_file <- file.path(output_dir, "qualified_context_groups.csv")
run_sets_file <- file.path(output_dir, "qualified_run_sets.csv")
run_set_members_file <- file.path(output_dir, "qualified_run_set_members.csv")
profile_metrics_file <- file.path(output_dir, "qualified_profile_metrics.csv")
feature_counts_file <- file.path(output_dir, "qualified_feature_counts.csv")
feature_contrasts_file <- file.path(output_dir, "qualified_feature_contrasts.csv")
replicate_diagnostics_file <- file.path(output_dir,
                                        "qualified_replicate_diagnostics.csv")
replicate_summary_file <- file.path(output_dir,
                                    "qualified_replicate_summary.csv")
replicate_effects_file <- file.path(output_dir,
                                    "qualified_replicate_effects.csv")
sequence_file <- file.path(output_dir, "qualified_sequence_features.csv")
flags_file <- file.path(output_dir, "qualified_review_flags.csv")
review_sheet_file <- file.path(output_dir, "qualified_review_sheet.csv")
review_priority_file <- file.path(output_dir, "qualified_review_priority.csv")
review_index_file <- file.path(output_dir, "qualified_review_index.html")

fwrite(manifest, manifest_file)
run_tables <- context_run_set_tables(groups)
fwrite(run_tables$groups, groups_file)
fwrite(run_tables$run_sets, run_sets_file)
fwrite(run_tables$run_set_members, run_set_members_file)
fwrite(profile_metrics_dt, profile_metrics_file)
fwrite(feature_counts_dt, feature_counts_file)
fwrite(feature_contrasts_dt, feature_contrasts_file)
fwrite(replicate_diagnostics_dt, replicate_diagnostics_file)
fwrite(replicate_summary_dt, replicate_summary_file)
fwrite(replicate_effects_dt, replicate_effects_file)
fwrite(sequence_features_dt, sequence_file)
fwrite(flags_dt, flags_file)

review_priority <- make_review_priority(manifest, replicate_summary_dt,
                                        replicate_effects_dt)
fwrite(review_priority, review_priority_file)
write_review_index(review_priority, output_dir, review_index_file)

review_sheet <- copy(manifest)
review_sheet[, `:=`(
  reviewer_signal_label = "",
  reviewer_mechanism_label = "",
  reviewer_context_label = "",
  reviewer_confidence = "",
  reviewer_notes = "",
  next_action = ""
)]
fwrite(review_sheet, review_sheet_file)

ok_count <- sum(manifest$status == "ok", na.rm = TRUE)
readiness_summary <- if ("inspection_readiness" %chin% names(manifest)) {
  manifest[status == "ok", .N, by = inspection_readiness][order(-N)]
} else {
  data.table()
}
flag_summary <- if (nrow(flags_dt) && "flags" %chin% names(flags_dt)) {
  flag_text <- flags_dt$flags[!is.na(flags_dt$flags) & nzchar(flags_dt$flags)]
  flag_values <- if (length(flag_text)) {
    unlist(strsplit(flag_text, ";", fixed = TRUE), use.names = FALSE)
  } else {
    character()
  }
  if (length(flag_values)) {
    data.table(flag = flag_values)[, .N, by = flag][order(-N)]
  } else {
    data.table(flag = character(), N = integer())
  }
} else {
  data.table(flag = character(), N = integer())
}

summary_file <- file.path(output_dir, "qualified_inspection_summary.md")
summary_lines <- c(
  "# Dominant RDG Qualified Inspection",
  "",
  paste0("Generated by `scripts/dominant_rdg_qualified_inspection.R` on ",
         format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "."),
  "",
  "## Inputs",
  "",
  paste0("- Queue: `", queue_file, "`"),
  paste0("- Rows selected: `", nrow(queue), "`"),
  paste0("- Metadata: `", metadata_file, "`"),
  paste0("- Feature expression: `", feature_expression_file, "`"),
  paste0("- RDG annotations: `", feature_annotation_file, "` and `",
         rdg_edges_file, "`"),
  "",
  "## Outputs",
  "",
  "- `qualified_inspection_manifest.csv`: one row per card with readiness flags and figure paths.",
  "- `qualified_context_groups.csv`: compact exact URL selections plus study/tissue/condition/all-merged comparison groups.",
  "- `qualified_run_sets.csv` and `qualified_run_set_members.csv`: run-set provenance for context groups.",
  "- `qualified_profile_metrics.csv`: displayed-coverage totals and spike metrics for plotted groups.",
  "- `qualified_feature_counts.csv`: per-feature raw counts, proportions, and single-run dominance metrics for every group.",
  "- `qualified_feature_contrasts.csv`: exact and broad-context per-feature deltas.",
  "- `qualified_replicate_diagnostics.csv`: exact case/control per-run top-feature and clean-CDS evidence, preserving missing raw-count measurements.",
  "- `qualified_replicate_summary.csv`: exact-group measured/missing runs, replicate support, and single-run dominance summary.",
  "- `qualified_replicate_effects.csv`: replicate-aware exact case/control top-feature effect estimates.",
  "- `qualified_review_priority.csv`: sorted card-level review queue with recommended actions and replicate support columns.",
  "- `qualified_review_index.html`: static review index with filters, row metrics, and all four PNG layers.",
  "- `qualified_sequence_features.csv`: RDG feature sequence, frame, Kozak, coding-potential, and sequence-availability fields.",
  "- `profiles/*.png`: normalized coverage profiles with all modeled RDG feature classes.",
  "- `feature_allocation/*.png`: per-feature allocation comparison across exact and broad groups.",
  "- `replicate_diagnostics/*.png`: per-run exact case/control top-feature and clean-CDS diagnostics; x marks missing values.",
  "- `rdg_flow/*.png`: RDG topology with scanning, translation, post-termination, competitive-exit, and CDS edges.",
  "",
  "## Render Summary",
  "",
  paste0("- Successful rows: `", ok_count, " / ", nrow(manifest), "`."),
  paste0("- Unique transcripts: `", uniqueN(manifest[status == "ok"]$tx_id), "`."),
  paste0("- Profile PNGs: `", length(list.files(profile_dir, pattern = "\\.png$")), "`."),
  paste0("- Feature-allocation PNGs: `", length(list.files(feature_dir, pattern = "\\.png$")), "`."),
  paste0("- Replicate-diagnostics PNGs: `", length(list.files(replicate_dir, pattern = "\\.png$")), "`."),
  paste0("- RDG-flow PNGs: `", length(list.files(flow_dir, pattern = "\\.png$")), "`.")
)
if (nrow(readiness_summary)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## Readiness Classes",
    "",
    paste0("- `", readiness_summary$inspection_readiness, "`: ",
           readiness_summary$N, collapse = "\n")
  )
}
if (nrow(flag_summary)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## Most Common Flags",
    "",
    paste0("- `", flag_summary$flag, "`: ", flag_summary$N,
           collapse = "\n")
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "## Inspection Notes",
  "",
  "- The profile layer is normalized by displayed counts, so it is for shape and localization, not total-expression claims.",
  "- The feature-count and contrast tables are the primary layer for branch/allocation interpretation.",
  "- The all-merged and tissue-general groups are included in feature-count tables even when they are too large for profile plotting.",
  "- `single_position_spike_risk`, `top3_position_spike_risk`, and `top_feature_single_run_dominated` mark cases that should not be accepted from screenshots alone.",
  "- Missing per-run top-feature counts are marked in replicate PNGs and counted in `qualified_replicate_summary.csv`; they should not be read as true zero support.",
  "- `aggregate_replicate_top_feature_direction_discordant` marks rows where aggregate exact proportions and unweighted replicate-mean proportions disagree in direction.",
  "- Rows with internal ORF, NTE, or NTT architecture need RDG-flow inspection because the allocation-only strip can hide the mechanism."
)
writeLines(summary_lines, summary_file)

message("Done: ", ok_count, " ok / ", nrow(manifest), " rows")
message("Manifest: ", manifest_file)
message("Summary: ", summary_file)
