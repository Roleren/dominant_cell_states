#!/usr/bin/env Rscript

# Add direct RiboCrypt Observatory links and compact review cards to the
# active-review batches. If the metadata drive is not mounted, links fall back
# to gene/transcript-only Observatory state and are flagged accordingly.

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

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for Observatory state URLs.",
       call. = FALSE)
}

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
output_dir <- file.path(analysis_dir, "dominant_rdg_review_links")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

batch_file <- file.path(
  analysis_dir,
  "dominant_rdg_review_batches",
  "dominant_rdg_review_batches.csv"
)
template_file <- file.path(
  analysis_dir,
  "dominant_rdg_review_batches",
  "dominant_rdg_review_batch_template.csv"
)
metadata_file <- local({
  candidates <- c(
    Sys.getenv("DOMINANT_METADATA_FILE", unset = NA_character_),
    "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv",
    path.expand("~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended_qc.csv")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
})

links_out <- file.path(output_dir, "dominant_rdg_review_batch_links.csv")
template_out <- file.path(output_dir, "dominant_rdg_review_batch_template_with_links.csv")
summary_out <- file.path(output_dir, "dominant_rdg_review_link_summary.csv")
html_out <- file.path(output_dir, "dominant_rdg_review_cards.html")
md_out <- file.path(output_dir, "dominant_rdg_review_cards.md")
report_out <- file.path(output_dir, "dominant_rdg_review_link_report.md")

safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(data.table())
  data.table::fread(path, showProgress = FALSE, nThread = 1, ...)
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

is_informative <- function(x) {
  x <- clean_text(x)
  nzchar(x) & !tolower(x) %chin% c("na", "nan", "null", "none", "missing")
}

scalar_text <- function(x, default = "") {
  x <- clean_text(x)
  if (!length(x) || !is_informative(x[1])) default else x[1]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

html_escape <- function(x) {
  x <- clean_text(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

shorten <- function(x, n = 260) {
  x <- gsub("\\s+", " ", clean_text(x))
  too_long <- nchar(x) > n
  x[too_long] <- paste0(substr(x[too_long], 1, n - 3), "...")
  x
}

split_condition_values <- function(x) {
  x <- clean_text(x)
  if (!length(x)) return(character())
  out <- unlist(strsplit(x, "\\s*(;|\\||,)\\s*", perl = TRUE), use.names = FALSE)
  out <- clean_text(out)
  unique(out[is_informative(out)])
}

compact_selections <- function(selections) {
  ids <- as.character(selections$index %||% names(selections$plot_selections))
  if (!length(ids)) ids <- "3"
  plot_sel <- selections$plot_selections %||% list()
  data_sel <- selections$data_table_selections %||% list()
  labels <- selections$labels %||% list()
  active <- as.character(selections$active_selection_id %||% ids[1])
  if (!active %in% ids) active <- ids[1]

  all_runs <- unique(unlist(c(plot_sel[ids], data_sel[ids]), use.names = FALSE))
  all_runs <- clean_text(all_runs)
  all_runs <- all_runs[is_informative(all_runs)]
  run_to_idx <- stats::setNames(seq_along(all_runs), all_runs)

  encode_idx <- function(x) {
    x <- clean_text(x)
    x <- x[is_informative(x)]
    if (!length(x)) return(integer())
    as.integer(run_to_idx[x])
  }

  list(
    active = active,
    order = ids,
    labels = stats::setNames(lapply(ids, function(id) {
      as.character(labels[[id]] %||% "")
    }), ids),
    runs = all_runs,
    p = stats::setNames(lapply(ids, function(id) encode_idx(plot_sel[[id]])), ids),
    d = stats::setNames(lapply(ids, function(id) encode_idx(data_sel[[id]])), ids)
  )
}

observatory_state <- function(gene_label, tx_id, selections) {
  list(
    v = 1L,
    exp = "all_samples-Homo_sapiens",
    color_by = c("tissue", "cell_line"),
    view = "browser",
    browser = list(
      gene = scalar_text(gene_label),
      tx = scalar_text(tx_id),
      frames_type = "columns",
      kmer = 1L,
      extendLeaders = 30L,
      extendTrailers = 30L,
      viewMode = FALSE,
      other_tx = FALSE,
      collapsed_introns = FALSE,
      collapsed_introns_width = 30L,
      genomic_region = "",
      zoom_range = "",
      customSequence = "",
      go = TRUE
    ),
    selections = compact_selections(selections)
  )
}

state_param <- function(state) {
  json <- jsonlite::toJSON(state, auto_unbox = TRUE, null = "null")
  compressed <- memCompress(charToRaw(json), type = "gzip")
  encoded <- jsonlite::base64_enc(compressed)
  encoded <- gsub("[\r\n]", "", encoded)
  encoded <- chartr("+/", "-_", encoded)
  sub("=+$", "", encoded)
}

load_metadata <- function(path) {
  if (!file.exists(path)) {
    return(list(data = data.table(), path = path, available = FALSE))
  }
  header <- names(data.table::fread(path, nrows = 0, showProgress = FALSE))
  wanted <- c(
    "Run", "study", "CONDITION", "CELL_LINE", "TISSUE", "GENE",
    "INHIBITOR", "FRACTION", "Cancer_type", "Sex", "TIMEPOINT", "AUTHOR"
  )
  cols <- intersect(wanted, header)
  if (!all(c("Run", "CONDITION") %chin% cols)) {
    return(list(data = data.table(), path = path, available = FALSE))
  }
  dt <- data.table::fread(path, select = cols, showProgress = FALSE,
                          nThread = 1)
  for (col in names(dt)) data.table::set(dt, j = col, value = clean_text(dt[[col]]))
  list(data = dt, path = path, available = TRUE)
}

make_gene_label_map <- function(dt) {
  out <- unique(dt[, .(
    tx_id = clean_text(tx_id),
    gene_symbol = clean_text(gene_symbol)
  )])
  out[, observatory_gene_label := gene_symbol]

  if (!requireNamespace("ORFik", quietly = TRUE)) return(out)
  labels <- tryCatch({
    df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
    sy <- ORFik::symbols(df)
    data.table::setDT(sy)
    if (!all(c("ensembl_tx_name", "ensembl_gene_id") %chin% names(sy))) {
      stop("Missing symbol columns", call. = FALSE)
    }
    sy[, gene_label_symbol := clean_text(external_gene_name)]
    sy[!is_informative(gene_label_symbol),
       gene_label_symbol := clean_text(ensembl_gene_id)]
    sy[, observatory_gene_label :=
         paste(gene_label_symbol, clean_text(ensembl_gene_id))]
    sy[, observatory_gene_label :=
         sub(" ", "-", observatory_gene_label, fixed = TRUE)]
    sy[, observatory_gene_label :=
         sub("(^-)|(^NA-)", "", observatory_gene_label, perl = TRUE)]
    unique(sy[
      is_informative(ensembl_tx_name) & is_informative(observatory_gene_label),
      .(tx_id = clean_text(ensembl_tx_name), observatory_gene_label)
    ])
  }, error = function(e) data.table())

  if (!nrow(labels)) return(out)
  out <- merge(out[, !"observatory_gene_label"], labels,
               by = "tx_id", all.x = TRUE, sort = FALSE)
  out[!is_informative(observatory_gene_label),
      observatory_gene_label := gene_symbol]
  out
}

filter_metadata_by_row <- function(meta, row) {
  if (!nrow(meta)) return(meta)
  for (col in c("AUTHOR", "CELL_LINE", "TISSUE", "GENE", "INHIBITOR",
                "FRACTION", "Cancer_type", "Sex", "TIMEPOINT")) {
    if (col %in% names(meta) && col %in% names(row) &&
        is_informative(row[[col]][1])) {
      val <- scalar_text(row[[col]][1])
      meta <- meta[get(col) == val]
    }
  }
  meta
}

condition_runs <- function(meta, values) {
  values <- clean_text(values)
  values <- values[is_informative(values)]
  if (!nrow(meta) || !"Run" %in% names(meta) || !"CONDITION" %in% names(meta) ||
      !length(values)) {
    return(character())
  }
  unique(meta[CONDITION %chin% values, Run])
}

row_subset <- function(row, metadata) {
  study <- scalar_text(row$study, scalar_text(row$best_study))
  case_condition <- scalar_text(row$case_condition,
                                scalar_text(row$best_case_condition))
  control_conditions <- split_condition_values(
    scalar_text(row$control_conditions,
                scalar_text(row$best_control_conditions))
  )
  context_label <- scalar_text(
    row$branch_context_label,
    paste(study, case_condition, sep = " / ")
  )

  if (!nrow(metadata) || !"study" %in% names(metadata)) {
    return(list(
      case_runs = character(),
      control_runs = character(),
      note = "gene_transcript_only_metadata_unavailable",
      mode = "gene_only",
      context_label = context_label,
      case_condition = case_condition,
      control_conditions = paste(control_conditions, collapse = "; "),
      study = study
    ))
  }

  base <- metadata
  if (is_informative(study)) {
    study_value <- study
    base <- base[get("study") == study_value]
  }
  exact <- filter_metadata_by_row(copy(base), row)
  exact_case_runs <- condition_runs(exact, case_condition)
  exact_control_runs <- condition_runs(exact, control_conditions)

  if (length(exact_case_runs) && length(exact_control_runs)) {
    case_runs <- exact_case_runs
    control_runs <- exact_control_runs
    mode <- "exact_context"
    note <- "exact_context_case_control"
  } else {
    study_case_runs <- condition_runs(base, case_condition)
    study_control_runs <- condition_runs(base, control_conditions)
    if (length(study_case_runs) && length(study_control_runs)) {
      case_runs <- study_case_runs
      control_runs <- study_control_runs
      mode <- "study_condition"
      note <- "study_condition_case_control_context_fields_not_matched"
    } else if (length(exact_case_runs) || length(exact_control_runs)) {
      case_runs <- exact_case_runs
      control_runs <- exact_control_runs
      mode <- "exact_context_partial"
      note <- "exact_context_partial_case_or_control_missing"
    } else {
      case_runs <- study_case_runs
      control_runs <- study_control_runs
      mode <- if (length(case_runs) || length(control_runs)) {
        "study_condition_partial"
      } else {
        "gene_only"
      }
      note <- if (mode == "study_condition_partial") {
        "study_condition_partial_case_or_control_missing"
      } else {
        "gene_transcript_only_no_matching_runs"
      }
    }
  }

  list(
    case_runs = case_runs,
    control_runs = control_runs,
    note = note,
    mode = mode,
    context_label = context_label,
    case_condition = case_condition,
    control_conditions = paste(control_conditions, collapse = "; "),
    study = study
  )
}

make_link_row <- function(row, metadata, gene_labels) {
  gene_label <- gene_labels[tx_id == row$tx_id[1], observatory_gene_label][1]
  if (!is_informative(gene_label)) gene_label <- row$gene_symbol[1]
  subset <- row_subset(row, metadata)

  ids <- character()
  plot_selections <- list()
  data_table_selections <- list()
  labels <- list()

  if (length(subset$control_runs)) {
    ids <- c(ids, "1")
    plot_selections[["1"]] <- subset$control_runs
    data_table_selections[["1"]] <- subset$control_runs
    labels[["1"]] <- paste0("Control: ", subset$control_conditions)
  }
  if (length(subset$case_runs)) {
    case_id <- as.character(length(ids) + 1L)
    ids <- c(ids, case_id)
    plot_selections[[case_id]] <- subset$case_runs
    data_table_selections[[case_id]] <- subset$case_runs
    labels[[case_id]] <- paste0("Case: ", subset$case_condition)
  }

  # Keep an empty third selection as a compatibility guard for the current
  # Observatory URL restore behavior, where the final populated group can be
  # lost. It also gives a visible all-merged comparison slot.
  ids <- unique(c(ids, "3"))
  plot_selections[["3"]] <- character()
  data_table_selections[["3"]] <- character()
  labels[["3"]] <- "All merged"

  state <- observatory_state(
    gene_label = gene_label,
    tx_id = row$tx_id[1],
    selections = list(
      index = ids,
      plot_selections = plot_selections,
      data_table_selections = data_table_selections,
      labels = labels,
      active_selection_id = "3"
    )
  )

  url <- paste0("https://ribocrypt.org/#Observatory?obs_state=",
                state_param(state))
  data.table(
    global_review_order = row$global_review_order[1],
    observatory_url = url,
    observatory_link_mode = subset$mode,
    observatory_link_note = paste(
      subset$note,
      paste0("study=", subset$study),
      paste0("context=", subset$context_label),
      paste0("case_n=", length(subset$case_runs)),
      paste0("control_n=", length(subset$control_runs)),
      sep = "; "
    ),
    observatory_gene_label = gene_label,
    observatory_case_n = length(subset$case_runs),
    observatory_control_n = length(subset$control_runs),
    observatory_context_label = subset$context_label,
    observatory_case_condition = subset$case_condition,
    observatory_control_conditions = subset$control_conditions
  )
}

batches <- safe_fread(batch_file)
if (nrow(batches) == 0) {
  stop("Missing review batch table: ", batch_file, call. = FALSE)
}
template <- safe_fread(template_file)
metadata_info <- load_metadata(metadata_file)
gene_labels <- make_gene_label_map(batches)

link_rows <- rbindlist(lapply(seq_len(nrow(batches)), function(i) {
  make_link_row(batches[i], metadata_info$data, gene_labels)
}), fill = TRUE)

linked <- merge(batches, link_rows, by = "global_review_order",
                all.x = TRUE, sort = FALSE)
setorder(linked, global_review_order)
link_cols <- c(
  "global_review_order", "batch_id", "batch_name", "batch_review_order",
  "batch_role", "gene_symbol", "tx_id", "design_family",
  "observatory_url", "observatory_link_mode", "observatory_link_note",
  "observatory_gene_label", "observatory_case_n", "observatory_control_n",
  "observatory_context_label", "review_instruction",
  setdiff(names(linked), c(
    "global_review_order", "batch_id", "batch_name", "batch_review_order",
    "batch_role", "gene_symbol", "tx_id", "design_family",
    "observatory_url", "observatory_link_mode", "observatory_link_note",
    "observatory_gene_label", "observatory_case_n", "observatory_control_n",
    "observatory_context_label", "review_instruction"
  ))
)
setcolorder(linked, intersect(link_cols, names(linked)))
fwrite(linked, links_out)

if (nrow(template)) {
  template_linked <- merge(template, link_rows, by = "global_review_order",
                           all.x = TRUE, sort = FALSE)
  setorder(template_linked, global_review_order)
  first_cols <- c(
    "global_review_order", "batch_id", "batch_name", "batch_review_order",
    "batch_role", "gene_symbol", "tx_id", "design_family",
    "observatory_url", "observatory_link_mode", "observatory_link_note",
    "review_instruction", "expected_browser_pattern",
    "expected_context_pattern", "review_question"
  )
  setcolorder(template_linked, intersect(c(first_cols,
                                           setdiff(names(template_linked),
                                                   first_cols)),
                                         names(template_linked)))
  fwrite(template_linked, template_out)
}

summary_dt <- linked[, .(
  n_rows = .N,
  n_genes = uniqueN(gene_symbol),
  n_exact_context_links = sum(observatory_link_mode == "exact_context"),
  n_exact_context_partial_links = sum(
    observatory_link_mode == "exact_context_partial"
  ),
  n_study_condition_links = sum(observatory_link_mode == "study_condition"),
  n_study_condition_partial_links = sum(
    observatory_link_mode == "study_condition_partial"
  ),
  n_gene_only_links = sum(observatory_link_mode == "gene_only"),
  median_case_n = as.numeric(median(observatory_case_n, na.rm = TRUE)),
  median_control_n = as.numeric(median(observatory_control_n, na.rm = TRUE))
), by = .(batch_id, batch_name, batch_role)]
summary_dt[, metadata_file_available := metadata_info$available]
summary_dt[, metadata_file := metadata_info$path]
fwrite(summary_dt, summary_out)

card_rows <- rbindlist(lapply(split(linked, linked$batch_id), function(x) {
  x <- copy(x)
  setorder(x, batch_review_order)
  x[seq_len(min(nrow(x), 35L))]
}), fill = TRUE)
setorder(card_rows, batch_id, batch_review_order)
cards <- vapply(seq_len(nrow(card_rows)), function(i) {
  row <- card_rows[i]
  paste0(
    "<article class=\"card role-", html_escape(row$batch_role), "\">",
    "<div class=\"meta\">", html_escape(row$batch_id), " / ",
    html_escape(row$batch_role), " / #", row$batch_review_order, "</div>",
    "<h3><a href=\"", html_escape(row$observatory_url),
    "\" target=\"_blank\" rel=\"noopener noreferrer\">",
    html_escape(row$gene_symbol), "</a></h3>",
    "<div class=\"sub\">", html_escape(row$design_family), " | ",
    html_escape(row$tx_id), "</div>",
    "<p><strong>Instruction:</strong> ",
    html_escape(shorten(row$review_instruction, 220)), "</p>",
    "<p><strong>Expected:</strong> ",
    html_escape(shorten(row$expected_browser_pattern, 260)), "</p>",
    "<p><strong>Link:</strong> ", html_escape(row$observatory_link_mode),
    " | case n=", row$observatory_case_n,
    " | control n=", row$observatory_control_n, "</p>",
    "</article>"
  )
}, character(1))

html <- c(
  "<!doctype html>",
  "<html><head><meta charset=\"utf-8\">",
  "<title>RDG Review Cards</title>",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:24px;background:#f6f7f9;color:#17202a}",
  "h1{font-size:26px;margin:0 0 6px 0} .note{color:#52606d;margin-bottom:18px}",
  ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:12px}",
  ".card{background:white;border:1px solid #d9e2ec;border-radius:8px;padding:12px;box-shadow:0 1px 2px rgba(16,24,40,.05)}",
  ".card h3{margin:4px 0 4px 0;font-size:18px}.card a{color:#005ea8;text-decoration:none}.card a:hover{text-decoration:underline}",
  ".meta{font-size:12px;color:#62748a;text-transform:uppercase;letter-spacing:.04em}.sub{font-size:13px;color:#486581;margin-bottom:8px}",
  "p{font-size:13px;line-height:1.35;margin:8px 0}",
  "</style></head><body>",
  "<h1>RDG Active-Review Cards</h1>",
  paste0("<div class=\"note\">Generated ", Sys.time(),
         ". Metadata file available: ", metadata_info$available,
         ". Links include an empty third all-merged group as an Observatory restore guard.</div>"),
  "<div class=\"grid\">",
  cards,
  "</div></body></html>"
)
writeLines(html, html_out, useBytes = TRUE)

md_lines <- c(
  "# RDG Active-Review Cards",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Metadata file available: ", metadata_info$available),
  "",
  "The Observatory links include an empty third `All merged` group as a restore guard.",
  "",
  vapply(seq_len(nrow(card_rows)), function(i) {
    row <- card_rows[i]
    paste0(
      "- ", row$batch_id, " #", row$batch_review_order, " ",
      "[", row$gene_symbol, " | ", row$design_family, "](",
      row$observatory_url, ")",
      " - ", row$observatory_link_mode,
      " - case n=", row$observatory_case_n,
      ", control n=", row$observatory_control_n
    )
  }, character(1))
)
writeLines(md_lines, md_out, useBytes = TRUE)

report_lines <- c(
  "# RDG Review Link Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("- Metadata file: ", metadata_info$path),
  paste0("- Metadata available: ", metadata_info$available),
  paste0("- Review rows linked: ", nrow(linked)),
  paste0("- Exact context links: ",
         sum(linked$observatory_link_mode == "exact_context")),
  paste0("- Exact context partial links: ",
         sum(linked$observatory_link_mode == "exact_context_partial")),
  paste0("- Study/condition links: ",
         sum(linked$observatory_link_mode == "study_condition")),
  paste0("- Study/condition partial links: ",
         sum(linked$observatory_link_mode == "study_condition_partial")),
  paste0("- Gene-only links: ",
         sum(linked$observatory_link_mode == "gene_only")),
  "",
  "Use `dominant_rdg_review_batch_template_with_links.csv` for manual review when direct links are needed."
)
writeLines(report_lines, report_out, useBytes = TRUE)

message("Wrote review links: ", links_out)
message("Wrote review cards: ", html_out)
