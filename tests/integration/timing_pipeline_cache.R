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
setup_runtime_env(analysis_dir)
load_ribocrypt(find_repo_root(analysis_dir))
source(file.path(script_dir, "dominant_pipeline_runtime_inputs.R"))
runtime_inputs <- pipeline_runtime_input_files()

steps <- pipeline_steps(analysis_dir)
outputs <- pipeline_step_outputs(analysis_dir)
cache_version <- Sys.getenv(
  "DOMINANT_PIPELINE_CACHE_VERSION",
  "dominant-pipeline-cache-v1"
)


timing <- system.time({
  status_dt <- rbindlist(lapply(seq_len(nrow(steps)), function(i) {
    step <- steps$step[[i]]
    deps <- pipeline_step_dependencies(step, analysis_dir, runtime_inputs = runtime_inputs)
    status <- pipeline_cache_status(
      step_name = step,
      script_file = steps$script[[i]],
      analysis_dir = analysis_dir,
      outputs = pipeline_step_cache_outputs(outputs[[step]],
        Sys.getenv("DOMINANT_RDG_PLOT_MODE", "all")),
      dependencies = deps,
      cache_version = cache_version,
      extra_source_files = pipeline_step_source_files(steps$script[[i]], analysis_dir)
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
