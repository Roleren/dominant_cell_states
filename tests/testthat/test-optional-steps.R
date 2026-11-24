# Optional internal steps are registered through an untracked local hook file.

make_hook_root <- function(hook_lines = NULL) {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(file.path(root, "scripts"))
  if (!is.null(hook_lines)) {
    writeLines(hook_lines, file.path(root, "scripts", "local_optional_steps.R"))
  }
  root
}

hook_lines <- c(
  "local_optional_steps <- function(analysis_dir) list(",
  "  list(step = 'extra_step', script = 'extra_step.R', after = 'rdg_graphs',",
  "       coverage_consumer = TRUE,",
  "       outputs = function(p) p('extra_dir', 'extra.csv'),",
  "       dependencies = function(rp, metadata_file) rp('input.csv'),",
  "       cache_check = function(results_dir, cache_state) {",
  "         cache_state$valid <- FALSE; cache_state$reason <- 'hook'; cache_state })",
  ")"
)

test_that("pipeline runs without optional steps when no hook file exists", {
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = dcs_project_root(), envir = env)
  root <- make_hook_root()
  expect_identical(env$pipeline_optional_steps(root), list())
  expect_identical(env$pipeline_scopes(root), env$pipeline_core_scopes())
  expect_false("extra_step" %in% env$pipeline_steps(root)$step)
})

test_that("hook file inserts optional steps after their anchor", {
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = dcs_project_root(), envir = env)
  root <- make_hook_root(hook_lines)
  withr::local_envvar(DOMINANT_CELL_STATES_RESULTS_DIR = file.path(root, "results"))

  steps <- env$pipeline_steps(root)
  idx <- match("extra_step", steps$step)
  expect_identical(steps$step[[idx - 1L]], "rdg_graphs")
  expect_identical(steps$script[[idx]], file.path(root, "scripts", "extra_step.R"))

  scopes <- env$pipeline_scopes(root)
  core <- env$pipeline_core_scopes()
  for (name in names(core)) {
    expect_identical("extra_step" %in% scopes[[name]], "rdg_graphs" %in% core[[name]])
  }
  expect_identical(scopes$rdg[match("extra_step", scopes$rdg) - 1L], "rdg_graphs")

  outputs <- env$pipeline_step_outputs(root)
  expect_identical(basename(outputs$extra_step), "extra.csv")
  expect_identical(env$pipeline_optional_coverage_consumers(root), "extra_step")
  spec <- env$pipeline_optional_step_spec("extra_step", root)
  expect_identical(spec$cache_check(root, list(valid = TRUE))$reason, "hook")
  expect_null(env$pipeline_optional_step_spec("rdg_graphs", root))
})

test_that("optional insertion ignores scopes without the anchor and duplicates", {
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = dcs_project_root(), envir = env)
  spec <- list(list(step = "x", after = "b"))
  expect_identical(env$pipeline_insert_optional_steps(c("a", "b", "c"), spec),
                   c("a", "b", "x", "c"))
  expect_identical(env$pipeline_insert_optional_steps(c("a", "c"), spec), c("a", "c"))
  expect_identical(env$pipeline_insert_optional_steps(c("b", "x"), spec), c("b", "x"))
})
