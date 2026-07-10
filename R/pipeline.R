dcs_pipeline_runner <- function(root = dcs_project_root()) {
  file.path(root, "scripts", "run_dominant_cell_state_pipeline.R")
}

dcs_load_pipeline <- function(root = dcs_project_root(),
                              envir = parent.frame()) {
  runner <- dcs_pipeline_runner(root)
  if (!file.exists(runner)) {
    stop("Missing pipeline runner: ", runner, call. = FALSE)
  }
  sys.source(runner, envir = envir)
  invisible(envir)
}

dcs_run_pipeline <- function(..., root = dcs_project_root()) {
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = root, envir = env)
  env$run_dominant_cell_state_pipeline(
    ...,
    repo_root = dcs_ribocrypt_repo(required = TRUE),
    analysis_dir = root
  )
}
