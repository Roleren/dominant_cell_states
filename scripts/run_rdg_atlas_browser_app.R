#!/usr/bin/env Rscript

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        file.exists(file.path(here, "DESCRIPTION")) &&
        dir.exists(file.path(here, "inst", "shiny", "rdg_atlas_browser"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (file.exists(file.path(candidate, "DESCRIPTION")) &&
        dir.exists(file.path(candidate, "inst", "shiny", "rdg_atlas_browser"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) {
      stop("Could not find dominant_cell_states analysis root from: ",
           start, call. = FALSE)
    }
    here <- parent
  }
}

analysis_dir <- find_analysis_dir()
for (helper in list.files(file.path(analysis_dir, "R"), pattern = "\\.R$",
                          full.names = TRUE)) {
  sys.source(helper, envir = globalenv())
}

args <- commandArgs(trailingOnly = TRUE)
port <- as.integer(Sys.getenv("RDG_ATLAS_BROWSER_PORT", "3876"))
host <- Sys.getenv("RDG_ATLAS_BROWSER_HOST", "127.0.0.1")
if (length(args) >= 1L && nzchar(args[[1]])) port <- as.integer(args[[1]])
if (length(args) >= 2L && nzchar(args[[2]])) host <- args[[2]]

message("Starting RDG atlas browser at http://", host, ":", port)
run_dcs_app(
  root = analysis_dir,
  host = host,
  port = port,
  launch.browser = FALSE
)
