#!/usr/bin/env Rscript

source(file.path(Sys.getenv("DOMINANT_CELL_STATES_DIR", "."), "scripts",
                 "dominant_script_common.R"), local = TRUE)

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
