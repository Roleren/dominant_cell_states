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

suppressPackageStartupMessages({
  library(data.table)
})

analysis_dir <- find_analysis_dir_local()
repo_root <- analysis_dir
script_dir <- file.path(analysis_dir, "scripts")
source(file.path(script_dir, "dominant_pipeline_cache.R"))
source(file.path(script_dir, "run_dominant_cell_state_pipeline.R"))

steps <- pipeline_steps(analysis_dir)
outputs <- pipeline_step_outputs(analysis_dir)
cache_version <- Sys.getenv(
  "DOMINANT_PIPELINE_CACHE_VERSION",
  "dominant-pipeline-cache-v1"
)
source_helpers <- file.path(
  script_dir,
  c("dominant_state_module_definitions.R", "dominant_coverage_cache.R")
)

timing <- system.time({
  status_dt <- rbindlist(lapply(seq_len(nrow(steps)), function(i) {
    step <- steps$step[[i]]
    deps <- pipeline_step_dependencies(step, analysis_dir)
    status <- pipeline_cache_status(
      step_name = step,
      script_file = steps$script[[i]],
      analysis_dir = analysis_dir,
      outputs = outputs[[step]],
      dependencies = deps,
      cache_version = cache_version,
      extra_source_files = source_helpers
    )
    data.table(
      step = step,
      cache_valid = isTRUE(status$valid),
      reason = status$reason,
      n_outputs = length(outputs[[step]]),
      n_dependencies = length(deps)
    )
  }), fill = TRUE)
})

status_dt[, `:=`(
  checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  cache_check_seconds_total = unname(timing[["elapsed"]])
)]

out_file <- file.path(
  analysis_dir,
  "results",
  "test_reports",
  "pipeline_cache_status.csv"
)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
fwrite(status_dt, out_file)
print(status_dt)
message("Cache status check elapsed seconds: ", round(timing[["elapsed"]], 3))
message("Saved: ", out_file)
