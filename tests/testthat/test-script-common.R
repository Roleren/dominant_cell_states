test_that("step scripts find a project root with any directory name", {
  env <- new.env()
  sys.source(file.path(dcs_project_root(), "scripts", "dominant_script_common.R"), envir = env)
  root <- withr::local_tempdir(pattern = "checkout_with_other_name")
  dir.create(file.path(root, "scripts"))
  dir.create(file.path(root, "md"))
  file.create(file.path(root, "md", "dominant_cell_state_scientific_note.md"))
  dir.create(file.path(root, "a", "b"), recursive = TRUE)
  withr::local_envvar(DOMINANT_CELL_STATES_DIR = "")
  expect_identical(env$find_analysis_dir(file.path(root, "a", "b")),
                   normalizePath(root))
  withr::local_envvar(DOMINANT_CELL_STATES_DIR = root)
  expect_identical(env$find_analysis_dir(tempdir()), normalizePath(root))
  withr::local_envvar(DOMINANT_CELL_STATES_DIR = file.path(root, "a"))
  expect_error(env$find_analysis_dir(), "not a dominant_cell_states root")
})

test_that("no tracked step script keeps a private copy of the root finder", {
  root <- dcs_project_root()
  scripts <- list.files(file.path(root, "scripts"), "[.]R$", full.names = TRUE)
  legacy <- vapply(scripts, function(f) {
    any(grepl('basename(here) == "dominant_cell_states"', readLines(f, warn = FALSE),
              fixed = TRUE))
  }, logical(1))
  # The runner keeps its own fallback after dcs_project_root().
  expect_identical(basename(names(legacy)[legacy]), "run_dominant_cell_state_pipeline.R")
})
