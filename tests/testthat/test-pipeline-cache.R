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

test_that("identical input rewrites reuse cache but altered outputs do not", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  script <- file.path(tmp, "step.R"); dep <- file.path(tmp, "input.csv")
  out <- file.path(tmp, "output.csv")
  writeLines("1 + 1", script); writeLines("a,b", dep); writeLines("x,y", out)
  write <- function() pipeline_cache_write("step", script, tmp, out, dep)
  status <- function() pipeline_cache_status("step", script, tmp, out, dep)
  write()
  Sys.setFileTime(dep, Sys.time() + 10)
  expect_true(status()$valid)
  writeLines("a,c", dep)
  expect_false(status()$valid)
  expect_match(status()$reason, "input.csv")
  write()
  writeLines("x,z", out)
  expect_false(status()$valid)
  expect_match(status()$reason, "output files changed")
})

test_that("helper dependencies are scoped and followed transitively", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  tmp <- tempfile(); dir.create(file.path(tmp, "scripts"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE))
  script <- file.path(tmp, "scripts", "step.R")
  helper <- file.path(tmp, "scripts", "helper.R")
  nested <- file.path(tmp, "scripts", "nested.R")
  writeLines('source("helper.R")', script)
  writeLines('source("nested.R")', helper)
  writeLines('source("helper.R")', nested)
  writeLines('1', file.path(tmp, "scripts", "unused.R"))
  expect_setequal(pipeline_step_source_files(script, tmp), c(helper, nested))
})

test_that("package code sourced from R/ is tracked as a step dependency", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  tmp <- tempfile()
  dir.create(file.path(tmp, "scripts"), recursive = TRUE)
  dir.create(file.path(tmp, "R"))
  on.exit(unlink(tmp, recursive = TRUE))
  script <- file.path(tmp, "scripts", "step.R")
  common <- file.path(tmp, "scripts", "dominant_script_common.R")
  engine <- file.path(tmp, "R", "ridge.R")
  writeLines(c(
    'source(file.path(Sys.getenv("DOMINANT_CELL_STATES_DIR", "."), "scripts", "dominant_script_common.R"), local = TRUE)',
    'source(file.path(analysis_dir, "R", "ridge.R"), local = TRUE)'), script)
  writeLines("1", common)
  writeLines("1", engine)
  expect_setequal(pipeline_step_source_files(script, tmp),
                  normalizePath(c(common, engine)))
})

test_that("megafst index timestamps invalidate step and measurement caches", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_runtime_inputs.R"))
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  index <- file.path(tmp, "coverage_index.fst")
  count <- file.path(tmp, "countTable_cds.qs")
  script <- file.path(tmp, "step.R"); out <- file.path(tmp, "out.csv")
  writeLines("index", index); writeLines("counts", count)
  writeLines("1", script); writeLines("1", out)
  deps <- pipeline_step_runtime_inputs("clean_cds", list(index = index, counts = count))
  pipeline_cache_write("clean_cds", script, tmp, out, deps)
  expect_true(pipeline_cache_status("clean_cds", script, tmp, out, deps)$valid)
  Sys.setFileTime(index, Sys.time() + 100)
  status <- pipeline_cache_status("clean_cds", script, tmp, out, deps)
  expect_false(status$valid)
  expect_match(status$reason, "coverage_index.fst")
  pipeline_cache_write("clean_cds", script, tmp, out, deps)
  Sys.setFileTime(count, Sys.time() + 100)
  expect_false(pipeline_cache_status("clean_cds", script, tmp, out, deps)$valid)
})

test_that("table mode does not require suppressed figures", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  outputs <- c("data.csv", "figures/chart.png", "figures/chart.html", "review.html")
  expect_identical(pipeline_step_cache_outputs(outputs, "tables"), c("data.csv", "review.html"))
  expect_identical(pipeline_step_cache_outputs(outputs, "all"), outputs)
})

test_that("unchanged project binary table rewrites retain downstream caches", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  script <- file.path(tmp, "step.R"); dep <- file.path(tmp, "activity.fst")
  out <- file.path(tmp, "out.csv")
  writeLines("1", script); writeLines("binary-table-fixture", dep); writeLines("1", out)
  pipeline_cache_write("test", script, tmp, out, dep)
  Sys.setFileTime(dep, Sys.time() + 100)
  expect_true(pipeline_cache_status("test", script, tmp, out, dep)$valid)
})

test_that("standalone coverage QC cannot skip a rebuilt megafst collection", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_runtime_inputs.R"))
  inputs <- list(index = "coverage_index.fst", sample_info = "libraries.rds")
  for (step in c("rdg_graphs", "cds_frame_bumpiness_qc", "infected_off_frame_hotspots")) {
    expect_true("coverage_index.fst" %in% pipeline_step_runtime_inputs(step, inputs))
  }
})

test_that("source discovery ignores existence probes and follows helper variables", {
  source(file.path(dcs_project_root(), "scripts", "dominant_pipeline_cache.R"))
  tmp <- tempfile(); dir.create(file.path(tmp,"scripts"),recursive=TRUE)
  on.exit(unlink(tmp,recursive=TRUE))
  script <- file.path(tmp,"scripts","step.R")
  helper <- file.path(tmp,"scripts","helper.R")
  writeLines(c('helper_path <- file.path("scripts", "helper.R")',
    'helper_path <- helper_path[file.exists(helper_path)][1]',
    'file.exists("unused.R")', '# source("also_unused.R")', 'source(helper_path)'), script)
  for (name in c("helper.R","unused.R","also_unused.R")) writeLines("1", file.path(tmp,"scripts",name))
  expect_identical(pipeline_step_source_files(script,tmp), helper)
})

test_that("signatures ignore input ordering and full outputs satisfy a subset", {
  root <- dcs_project_root()
  source(file.path(root,'scripts/dominant_pipeline_cache.R'))
  tmp <- tempfile();dir.create(tmp);dir.create(file.path(tmp,'scripts'))
  on.exit(unlink(tmp,recursive=TRUE))
  file.copy(file.path(root,'scripts/dominant_pipeline_cache.R'),file.path(tmp,'scripts/dominant_pipeline_cache.R'))
  script <- file.path(tmp,'step.R');writeLines('1 + 1',script)
  a <- file.path(tmp,'a.csv');b <- file.path(tmp,'b.csv')
  writeLines('a',a);writeLines('b',b)
  expect_identical(pipeline_file_signature(c(a,b)),pipeline_file_signature(c(b,a)))
  pipeline_cache_write('subset',script,tmp,outputs=c(b,a))
  expect_true(pipeline_cache_status('subset',script,tmp,outputs=a)$valid)
  writeLines('changed',a)
  expect_false(pipeline_cache_status('subset',script,tmp,outputs=a)$valid)
})
