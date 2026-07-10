#!/usr/bin/env Rscript

# Render static documentation pages for the RDG Atlas Browser.
# The Shiny app serves these prebuilt HTML files through an iframe, so page
# navigation stays fast and Markdown conversion is kept out of request handling.

dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

if (!requireNamespace("commonmark", quietly = TRUE)) {
  stop("Package 'commonmark' is required to render browser docs.",
       call. = FALSE)
}
if (!requireNamespace("htmltools", quietly = TRUE)) {
  stop("Package 'htmltools' is required to render browser docs.",
       call. = FALSE)
}

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        file.exists(file.path(here, "md",
                              "dominant_cell_state_scientific_note.md"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (file.exists(file.path(candidate,
                              "md",
                              "dominant_cell_state_scientific_note.md"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states from: ", start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
docs_dir <- file.path(
  analysis_dir,
  "inst",
  "shiny",
  "rdg_atlas_browser",
  "docs"
)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

doc_css <- paste(
  "body { font-family: -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif; color: #1f2933; line-height: 1.58; margin: 0; background: #f7fafc; }",
  ".doc-shell { max-width: 980px; margin: 0 auto; padding: 34px 34px 60px; background: #fff; min-height: 100vh; border-left: 1px solid #e5eaf0; border-right: 1px solid #e5eaf0; }",
  "h1, h2, h3 { color: #102a43; line-height: 1.2; }",
  "h1 { margin-top: 0; font-size: 30px; }",
  "h2 { border-top: 1px solid #edf1f5; padding-top: 22px; margin-top: 32px; }",
  "a { color: #0b5cab; }",
  "code { background: #f1f5f9; padding: 2px 4px; border-radius: 4px; }",
  "pre { background: #111827; color: #f9fafb; padding: 14px 16px; border-radius: 7px; overflow-x: auto; }",
  "pre code { background: transparent; padding: 0; color: inherit; }",
  "blockquote { border-left: 4px solid #cbd5e1; margin-left: 0; padding-left: 14px; color: #52606d; }",
  "table { border-collapse: collapse; width: 100%; margin: 12px 0; }",
  "th, td { border: 1px solid #d9dde3; padding: 7px 9px; vertical-align: top; }",
  "th { background: #f1f5f9; }",
  "ul, ol { padding-left: 24px; }",
  "li { margin: 4px 0; }",
  sep = "\n"
)

render_doc <- function(source_rel_path, output_file, title) {
  source_file <- file.path(analysis_dir, source_rel_path)
  if (!file.exists(source_file)) {
    stop("Missing Markdown source: ", source_file, call. = FALSE)
  }
  md <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  body <- commonmark::markdown_html(md, extensions = TRUE)
  page <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", htmltools::htmlEscape(title), "</title>\n<style>\n",
    doc_css,
    "\n</style>\n</head>\n<body>\n<main class=\"doc-shell\">\n",
    body,
    "\n</main>\n</body>\n</html>\n"
  )
  writeLines(page, output_file, useBytes = TRUE)
  data.frame(
    source = source_rel_path,
    output = sub(paste0("^", normalizePath(analysis_dir, mustWork = TRUE),
                        .Platform$file.sep), "",
                 normalizePath(output_file, mustWork = TRUE)),
    source_mtime = as.character(file.info(source_file)$mtime),
    output_mtime = as.character(file.info(output_file)$mtime),
    output_bytes = file.info(output_file)$size,
    stringsAsFactors = FALSE
  )
}

rendered <- do.call(rbind, list(
  render_doc(
    file.path("md", "dominant_rdg_atlas_browser_tutorial.md"),
    file.path(docs_dir, "tutorial.html"),
    "RDG Atlas Browser Tutorial"
  ),
  render_doc(
    file.path("md", "dominant_cell_state_scientific_note.md"),
    file.path(docs_dir, "scientific_note.html"),
    "Dominant Cell States and RDGs Scientific Note"
  )
))

manifest_file <- file.path(docs_dir, "docs_manifest.csv")
utils::write.csv(rendered, manifest_file, row.names = FALSE, quote = TRUE)

message("Rendered RDG Atlas Browser docs:")
for (i in seq_len(nrow(rendered))) {
  message("  ", rendered$output[i], " (", rendered$output_bytes[i], " bytes)")
}
message("  ", sub(paste0("^", normalizePath(analysis_dir, mustWork = TRUE),
                         .Platform$file.sep), "",
                  normalizePath(manifest_file, mustWork = TRUE)))
