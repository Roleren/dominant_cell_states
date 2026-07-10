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
    if (identical(arg, "--html")) {
      out$html <- TRUE
    } else if (identical(arg, "--force")) {
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
  "dominant_rdg_validation_label_bridge",
  "dominant_rdg_validation_label_gap_queue.csv"
)
queue_file <- opts$queue %||% Sys.getenv("RDG_COVERAGE_QUEUE",
                                         unset = default_queue)
output_dir <- opts$output_dir %||% Sys.getenv(
  "RDG_COVERAGE_OUTPUT_DIR",
  unset = file.path(analysis_dir, "dominant_rdg_card_coverage_review")
)
limit_raw <- opts$limit %||% Sys.getenv("RDG_COVERAGE_LIMIT", unset = "")
row_raw <- opts$row %||% Sys.getenv("RDG_COVERAGE_ROW", unset = "")
write_html <- isTRUE(opts$html) ||
  as_bool(Sys.getenv("RDG_COVERAGE_HTML", unset = ""), FALSE)
force <- isTRUE(opts$force) ||
  as_bool(Sys.getenv("RDG_COVERAGE_FORCE", unset = ""), FALSE)
stop_on_error <- as_bool(Sys.getenv("RDG_COVERAGE_STOP_ON_ERROR", unset = ""),
                         FALSE)

limit_n <- suppressWarnings(as.integer(limit_raw))
if (length(limit_n) == 0 || is.na(limit_n) || limit_n <= 0L) limit_n <- NA_integer_
row_n <- suppressWarnings(as.integer(row_raw))
if (length(row_n) == 0 || is.na(row_n) || row_n <= 0L) row_n <- NA_integer_

figure_dir <- file.path(output_dir, "figures")
html_dir <- file.path(output_dir, "html")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (write_html) dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)

feature_annotation_file <- file.path(
  analysis_dir,
  "dominant_rdg_outputs",
  "rdg_feature_annotation.csv"
)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

scalar_chr <- function(x, default = "") {
  if (!length(x) || is.na(x[[1]])) return(default)
  as.character(x[[1]])
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

profile_groups <- function(reads, groups) {
  mat <- as.matrix(reads)
  out <- rbindlist(lapply(seq_along(groups), function(i) {
    group <- groups[[i]]
    cols <- intersect(group$runs, colnames(mat))
    counts <- if (length(cols)) {
      matrixStats::rowSums2(mat, cols = match(cols, colnames(mat)), na.rm = TRUE)
    } else {
      rep(0, nrow(mat))
    }
    data.table(
      group_index = i,
      group_id = group$id,
      group_label = group$label,
      group_role = group$role,
      run_count = length(cols),
      position = seq_len(nrow(mat)),
      frame = factor((seq_len(nrow(mat)) - 1L) %% 3L, levels = 0:2),
      count = as.numeric(counts)
    )
  }))
  out[, group_label := factor(group_label, levels = unique(group_label))]
  out
}

feature_class <- function(feature_type, category) {
  out <- rep("other", length(feature_type))
  out[feature_type == "clean_CDS" | category == "clean_CDS"] <- "clean_CDS"
  out[feature_type == "leader_uorf" | category == "uORF"] <- "leader_uORF"
  out[feature_type == "overlapping_uorf" |
        category %chin% c("uoORF", "ouORF")] <- "overlapping_uORF"
  out
}

feature_strip <- function(feature_annotation, row, tx_id, leader_extension,
                          profile_dt) {
  empty <- data.table()
  attr(empty, "feature_qc") <- list(
    feature_rows = 0L,
    feature_rows_outside_display = 0L,
    feature_max_display_x = NA_real_
  )
  if (!nrow(feature_annotation) || !"tx_id" %in% names(feature_annotation)) {
    return(empty)
  }
  keep_tx <- feature_annotation[["tx_id"]] == tx_id
  keep_tx[is.na(keep_tx)] <- FALSE
  fx <- copy(feature_annotation[which(keep_tx)])
  if (!nrow(fx)) return(empty)
  fx[, feature_class := feature_class(feature_type, category)]
  fx <- fx[feature_class %chin% c("clean_CDS", "leader_uORF", "overlapping_uORF")]
  if (!nrow(fx)) return(empty)
  top_feature <- coalesce_col(
    row,
    c("top_shifted_feature_id", "branch_usage_top_feature_id",
      "grouped_branch_usage_top_feature_id", "count_aware_top_feature_id"),
    default = ""
  )
  fx[, is_top_feature := nzchar(top_feature) & feature_id == top_feature]
  fx[, `:=`(
    xmin = as.numeric(tx_start) + leader_extension,
    xmax = as.numeric(tx_end) + leader_extension
  )]
  display_width <- max(profile_dt$position, na.rm = TRUE)
  feature_rows <- nrow(fx)
  feature_rows_outside_display <- fx[
    !is.finite(xmin) | !is.finite(xmax) | xmax < 1 | xmin > display_width,
    .N
  ]
  fx <- fx[is.finite(xmin) & is.finite(xmax) & xmax >= 1 & xmin <= display_width]
  if (!nrow(fx)) {
    attr(empty, "feature_qc") <- list(
      feature_rows = feature_rows,
      feature_rows_outside_display = feature_rows_outside_display,
      feature_max_display_x = NA_real_
    )
    return(empty)
  }
  fx[, `:=`(
    xmin = pmax(1, xmin),
    xmax = pmin(display_width, xmax)
  )]
  feature_max_display_x <- max(fx$xmax, na.rm = TRUE)
  groups <- unique(profile_dt[, .(group_label, group_index)])
  strip <- fx[rep(seq_len(nrow(fx)), nrow(groups))]
  strip[, group_label := rep(groups$group_label, each = nrow(fx))]
  strip[, group_index := rep(groups$group_index, each = nrow(fx))]
  y <- profile_dt[, .(ymax_group = max(count, 1, na.rm = TRUE)), by = group_label]
  strip <- y[strip, on = "group_label"]
  strip[, `:=`(
    ymin = -0.18 * ymax_group,
    ymax = -0.04 * ymax_group,
    label_y = -0.11 * ymax_group
  )]
  attr(strip, "feature_qc") <- list(
    feature_rows = feature_rows,
    feature_rows_outside_display = feature_rows_outside_display,
    feature_max_display_x = feature_max_display_x
  )
  strip
}

coverage_plot <- function(profile_dt, strip_dt, row, state, output_file) {
  gene <- coalesce_col(row, c("gene_symbol"), default = state$browser$gene %||% "")
  tx_id <- coalesce_col(row, c("tx_id"), default = state$browser$tx %||% "")
  design <- coalesce_col(row, c("design_family"), default = "")
  context <- coalesce_col(
    row,
    c("observatory_context_label", "branch_context_label",
      "observatory_link_note"),
    default = ""
  )
  expected <- coalesce_col(row, c("expected_browser_pattern"), default = "")
  question <- coalesce_col(row, c("review_question"), default = "")
  clean_pct <- first_finite_num(
    row,
    c("atlas_clean_cds_pct_change", "half_normalization_cds_pct_change")
  )
  clean_text <- if (is.finite(clean_pct)) {
    paste0("clean CDS ", sprintf("%+.1f%%", clean_pct))
  } else {
    "clean CDS effect not available"
  }
  title <- paste(gene, tx_id, design, sep = " | ")
  subtitle <- paste(c(clean_text, expected, question), collapse = "\n")
  caption <- context

  frame_cols <- c("0" = "#D95F02", "1" = "#1B9E77", "2" = "#7570B3")
  feature_cols <- c(
    clean_CDS = "#4C78A8",
    leader_uORF = "#F58518",
    overlapping_uORF = "#E45756"
  )

  p <- ggplot(profile_dt, aes(position, count)) +
    geom_col(aes(fill = frame), width = 1, linewidth = 0) +
    facet_grid(group_label ~ ., scales = "free_y", switch = "y") +
    scale_y_continuous(expand = expansion(mult = c(0.22, 0.08))) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Transcript coordinate in displayed range",
      y = "P-site coverage",
      fill = "Frame",
      caption = caption
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8.5, lineheight = 0.95),
      plot.caption = element_text(size = 7.2, hjust = 0),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  if (nrow(strip_dt)) {
    p <- p +
      geom_rect(
        data = strip_dt,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
            fill = feature_class),
        inherit.aes = FALSE,
        alpha = 0.70,
        color = "grey25",
        linewidth = 0.2
      ) +
      geom_rect(
        data = strip_dt[is_top_feature == TRUE],
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = "black",
        linewidth = 0.55
      )
  }
  p <- p +
    scale_fill_manual(
      values = c(frame_cols, feature_cols),
      breaks = c(names(frame_cols), names(feature_cols)),
      drop = TRUE
    )

  ggsave(output_file, p, width = 12, height = max(5.2, 2.35 + 1.55 *
                                                    uniqueN(profile_dt$group_label)),
         dpi = 140, limitsize = FALSE)
  p
}

write_plot_html <- function(plot, output_file) {
  if (!requireNamespace("plotly", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    return(FALSE)
  }
  helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
  if (file.exists(helper)) source(helper)
  widget <- suppressWarnings(plotly::ggplotly(plot))
  if (exists("dominant_save_widget", mode = "function")) {
    dominant_save_widget(widget, output_file = output_file,
                         title = tools::file_path_sans_ext(basename(output_file)))
  } else {
    htmlwidgets::saveWidget(widget, output_file, selfcontained = FALSE)
  }
  TRUE
}

render_card <- function(row, row_index, df_cache, feature_annotation) {
  url <- coalesce_col(row, c("observatory_url"), default = "")
  state <- parse_observatory_state(url)
  if (is.null(state)) {
    stop("Could not parse observatory_url")
  }
  tx_id <- coalesce_col(row, c("tx_id"), default = state$browser$tx %||% "")
  gene <- coalesce_col(row, c("gene_symbol"), default = state$browser$gene %||% "")
  if (!nzchar(tx_id)) stop("Missing transcript id")

  groups <- state_selection_groups(state)
  groups <- groups[vapply(groups, function(x) x$role != "all_merged", logical(1))]
  if (!length(groups)) stop("No non-empty control/case selections in URL state")

  exp_name <- scalar_chr(state$exp, "all_samples-Homo_sapiens")
  if (is.null(df_cache[[exp_name]])) {
    df_cache[[exp_name]] <- ORFik::read.experiment(exp_name, validate = FALSE)
  }
  df <- df_cache[[exp_name]]

  tx_region <- ORFik::loadRegion(df, "mrna", names.keep = tx_id)
  if (!length(tx_region)) stop("Transcript not found in experiment: ", tx_id)
  display_region <- display_region_from_state(tx_region, state$browser)
  collection_path <- RiboCrypt:::collection_path_from_exp(
    df, tx_id, grl_all = tx_region
  )
  requested_runs <- unique(unlist(lapply(groups, `[[`, "runs"), use.names = FALSE))
  requested_runs <- requested_runs[!is.na(requested_runs) & nzchar(requested_runs)]
  available_runs <- collection_available_runs(collection_path)
  selected_runs <- intersect(requested_runs, available_runs)
  missing_runs <- setdiff(requested_runs, available_runs)
  if (!length(selected_runs)) {
    stop("No selected runs are present in collection for ", tx_id)
  }
  groups <- lapply(groups, function(group) {
    group$runs <- intersect(group$runs, selected_runs)
    group
  })
  groups <- groups[lengths(lapply(groups, `[[`, "runs")) > 0]
  if (!length(groups)) stop("No selected groups remain after run matching")

  reads <- RiboCrypt:::load_collection(
    collection_path,
    grl = display_region,
    columns = selected_runs
  )
  profile_dt <- profile_groups(reads, groups)
  leader_extension <- suppressWarnings(as.numeric(state$browser$extendLeaders %||% 0))
  if (!is.finite(leader_extension)) leader_extension <- 0
  strip_dt <- feature_strip(feature_annotation, row, tx_id,
                            leader_extension, profile_dt)
  feature_qc <- attr(strip_dt, "feature_qc", exact = TRUE) %||% list(
    feature_rows = 0L,
    feature_rows_outside_display = 0L,
    feature_max_display_x = NA_real_
  )

  row_rank <- coalesce_col(row, c("label_gap_rank", "global_review_order"),
                           default = row_index)
  file_stem <- paste(
    sprintf("%03d", as.integer(row_index)),
    safe_name(row_rank),
    safe_name(gene),
    safe_name(tx_id),
    sep = "_"
  )
  png_file <- file.path(figure_dir, paste0(file_stem, ".png"))
  html_file <- file.path(html_dir, paste0(file_stem, ".html"))

  plot <- NULL
  if (force || !file.exists(png_file)) {
    plot <- coverage_plot(profile_dt, strip_dt, row, state, png_file)
  }
  html_written <- FALSE
  if (write_html && (force || !file.exists(html_file))) {
    if (is.null(plot)) plot <- coverage_plot(profile_dt, strip_dt, row, state, png_file)
    html_written <- write_plot_html(plot, html_file)
  } else if (write_html && file.exists(html_file)) {
    html_written <- TRUE
  }

  group_summary <- profile_dt[, .(
    run_count = unique(run_count),
    display_total_counts = sum(count, na.rm = TRUE),
    max_position_count = max(count, na.rm = TRUE)
  ), by = .(group_label, group_role)]

  control_total <- group_summary[group_role == "control",
                                 sum(display_total_counts, na.rm = TRUE)]
  case_total <- group_summary[group_role == "case",
                              sum(display_total_counts, na.rm = TRUE)]
  observed_log2_case_control <- if (is.finite(control_total) &&
                                   is.finite(case_total)) {
    log2((case_total + 1) / (control_total + 1))
  } else {
    NA_real_
  }

  data.table(
    source_row = row_index,
    label_gap_rank = coalesce_col(row, c("label_gap_rank"), default = ""),
    batch_id = coalesce_col(row, c("batch_id"), default = ""),
    gene_symbol = gene,
    tx_id = tx_id,
    design_family = coalesce_col(row, c("design_family"), default = ""),
    experiment = exp_name,
    selected_group_count = length(groups),
    selected_run_count = length(selected_runs),
    missing_run_count = length(missing_runs),
    missing_runs = paste(missing_runs, collapse = ";"),
    display_width = nrow(profile_dt) / length(groups),
    feature_rows = feature_qc$feature_rows %||% 0L,
    feature_rows_outside_display = feature_qc$feature_rows_outside_display %||% 0L,
    feature_max_display_x = feature_qc$feature_max_display_x %||% NA_real_,
    control_display_total_counts = control_total,
    case_display_total_counts = case_total,
    observed_log2_case_control = observed_log2_case_control,
    expected_clean_cds_pct_change = first_finite_num(
      row,
      c("atlas_clean_cds_pct_change", "half_normalization_cds_pct_change")
    ),
    top_feature_id = coalesce_col(
      row,
      c("top_shifted_feature_id", "branch_usage_top_feature_id",
        "grouped_branch_usage_top_feature_id", "count_aware_top_feature_id"),
      default = ""
    ),
    observatory_context_label = coalesce_col(
      row,
      c("observatory_context_label", "branch_context_label"),
      default = ""
    ),
    png_file = png_file,
    html_file = if (write_html && html_written) html_file else "",
    status = "ok",
    error = ""
  )
}

if (!file.exists(queue_file)) {
  stop("Missing queue file: ", queue_file, call. = FALSE)
}

queue <- fread(queue_file, showProgress = FALSE)
if (!"observatory_url" %in% names(queue)) {
  stop("Queue must contain observatory_url: ", queue_file, call. = FALSE)
}
if (is.finite(row_n)) {
  queue <- queue[row_n]
} else if (is.finite(limit_n)) {
  queue <- head(queue, limit_n)
}
if (!nrow(queue)) stop("No rows selected from queue: ", queue_file, call. = FALSE)

feature_annotation <- if (file.exists(feature_annotation_file)) {
  fread(feature_annotation_file, showProgress = FALSE)
} else {
  data.table()
}

message("Dominant RDG card coverage review")
message("  queue: ", queue_file)
message("  rows: ", nrow(queue))
message("  output: ", output_dir)
message("  html: ", write_html)

df_cache <- new.env(parent = emptyenv())
manifest <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  message("[", i, "/", nrow(queue), "] ",
          coalesce_col(queue[i], c("gene_symbol"), default = "unknown"),
          " ", coalesce_col(queue[i], c("tx_id"), default = ""))
  if (stop_on_error) {
    manifest[[i]] <- render_card(queue[i], i, df_cache, feature_annotation)
  } else {
    manifest[[i]] <- tryCatch(
      render_card(queue[i], i, df_cache, feature_annotation),
      error = function(e) {
        data.table(
          source_row = i,
          label_gap_rank = coalesce_col(queue[i], c("label_gap_rank"), default = ""),
          batch_id = coalesce_col(queue[i], c("batch_id"), default = ""),
          gene_symbol = coalesce_col(queue[i], c("gene_symbol"), default = ""),
          tx_id = coalesce_col(queue[i], c("tx_id"), default = ""),
          design_family = coalesce_col(queue[i], c("design_family"), default = ""),
          experiment = "",
          selected_group_count = NA_integer_,
          selected_run_count = NA_integer_,
          missing_run_count = NA_integer_,
          missing_runs = "",
        display_width = NA_real_,
        feature_rows = NA_integer_,
        feature_rows_outside_display = NA_integer_,
        feature_max_display_x = NA_real_,
        control_display_total_counts = NA_real_,
          case_display_total_counts = NA_real_,
          observed_log2_case_control = NA_real_,
          expected_clean_cds_pct_change = first_finite_num(
            queue[i],
            c("atlas_clean_cds_pct_change", "half_normalization_cds_pct_change")
          ),
          top_feature_id = coalesce_col(
            queue[i],
            c("top_shifted_feature_id", "branch_usage_top_feature_id",
              "grouped_branch_usage_top_feature_id", "count_aware_top_feature_id"),
            default = ""
          ),
          observatory_context_label = coalesce_col(
            queue[i],
            c("observatory_context_label", "branch_context_label"),
            default = ""
          ),
          png_file = "",
          html_file = "",
          status = "error",
          error = conditionMessage(e)
        )
      }
    )
  }
}

manifest <- rbindlist(manifest, fill = TRUE)
manifest_file <- file.path(output_dir, "dominant_rdg_card_coverage_review_manifest.csv")
fwrite(manifest, manifest_file)

review_file <- file.path(output_dir, "dominant_rdg_card_coverage_review_sheet.csv")
review_sheet <- copy(manifest)
review_sheet[, `:=`(
  reviewer_label = "",
  reviewer_signal_grade = "",
  reviewer_notes = "",
  next_action = ""
)]
fwrite(review_sheet, review_file)

readme <- c(
  "Dominant RDG Card Coverage Review",
  "=================================",
  "",
  "Generated by `scripts/dominant_rdg_card_coverage_review.R`.",
  "",
  "Default usage from the RiboCrypt repo root:",
  "",
  "```sh",
  "env -u LC_ALL \\",
  "  XDG_CACHE_HOME=\"$PWD/dominant_cell_states/.runtime/xdg-cache\" \\",
  "  XDG_CONFIG_HOME=\"$PWD/dominant_cell_states/.runtime/xdg-config\" \\",
  "  BFC_CACHE=\"$PWD/dominant_cell_states/.runtime/biocfilecache\" \\",
  "  R --vanilla -q -e 'devtools::load_all(\".\"); source(\"dominant_cell_states/scripts/dominant_rdg_card_coverage_review.R\")'",
  "```",
  "",
  "Useful environment overrides:",
  "- `RDG_COVERAGE_LIMIT=20` renders the first 20 queue rows.",
  "- `RDG_COVERAGE_ROW=7` renders a single queue row.",
  "- `RDG_COVERAGE_HTML=1` also writes interactive HTML widgets.",
  "- `RDG_COVERAGE_FORCE=1` overwrites existing plot files.",
  "",
  "Outputs:",
  "- `figures/*.png`: static control/case coverage panels.",
  "- `dominant_rdg_card_coverage_review_manifest.csv`: render and QC manifest.",
  "- `dominant_rdg_card_coverage_review_sheet.csv`: empty review columns plus manifest fields."
)
writeLines(readme, file.path(output_dir, "README.md"))

ok_count <- sum(manifest$status == "ok", na.rm = TRUE)
message("Done: ", ok_count, " ok / ", nrow(manifest), " rows")
message("Manifest: ", manifest_file)
