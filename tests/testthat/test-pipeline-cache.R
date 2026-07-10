test_that("pipeline cache invalidates on output, dependency, and version changes", {
  root <- dcs_project_root()
  script_dir <- file.path(root, "scripts")
  source(file.path(script_dir, "dominant_pipeline_cache.R"))

  tmp <- tempfile("dominant-cache-test-")
  dir.create(tmp, recursive = TRUE)
  dir.create(file.path(tmp, "scripts"), recursive = TRUE)
  file.copy(
    file.path(script_dir, "dominant_pipeline_cache.R"),
    file.path(tmp, "scripts", "dominant_pipeline_cache.R")
  )

  script <- file.path(tmp, "dummy_step.R")
  dep <- file.path(tmp, "input.txt")
  out <- file.path(tmp, "output.csv")
  writeLines("message('dummy')", script)
  writeLines("input-v1", dep)

  status <- pipeline_cache_status(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v1"
  )
  expect_false(status$valid)
  expect_match(status$reason, "missing output")

  writeLines("x,y\n1,2", out)
  pipeline_cache_write(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v1",
    seconds = 0
  )
  status <- pipeline_cache_status(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v1"
  )
  expect_true(status$valid)

  writeLines(c("input-v1", "input-v2"), dep)
  status <- pipeline_cache_status(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v1"
  )
  expect_false(status$valid)

  pipeline_cache_write(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v1",
    seconds = 0
  )
  status <- pipeline_cache_status(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v2"
  )
  expect_false(status$valid)

  unlink(out)
  status <- pipeline_cache_status(
    step_name = "dummy",
    script_file = script,
    analysis_dir = tmp,
    outputs = out,
    dependencies = dep,
    cache_version = "v1"
  )
  expect_false(status$valid)
  expect_match(status$reason, "missing output")
})
