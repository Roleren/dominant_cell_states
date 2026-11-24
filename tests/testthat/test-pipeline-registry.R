test_that("pipeline registry loads and declares package-local outputs", {
  root <- dcs_project_root()
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = root, envir = env)

  steps <- env$pipeline_steps(root)
  expect_true(nrow(steps) > 0)
  expect_true(all(file.exists(steps$script)))
  expect_identical(anyDuplicated(steps$step), 0L)

  outputs <- env$pipeline_step_outputs(root)
  flat_outputs <- unlist(outputs, use.names = FALSE)
  expect_true(length(flat_outputs) > 0)
  normalized <- normalizePath(flat_outputs, mustWork = FALSE)
  result_prefix <- paste0(dcs_results_dir(root = root), .Platform$file.sep)
  expect_true(all(startsWith(normalized, result_prefix)))
  expect_true(any(grepl("merged_replicate_translon_shift_candidate_queue.csv",
                        flat_outputs, fixed = TRUE)))
})

test_that("pipeline step slicing works without generated outputs", {
  root <- dcs_project_root()
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = root, envir = env)

  sliced <- env$slice_steps(
    c("manual_translons", "fst_pages", "clean_cds"),
    from = "fst_pages",
    to = "clean_cds"
  )
  expect_identical(sliced, c("fst_pages", "clean_cds"))
})

test_that("measurement dependencies exclude unused generated summaries", {
  root <- dcs_project_root()
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = root, envir = env)
  clean <- basename(env$pipeline_step_dependencies("clean_cds", root, runtime_inputs = list()))
  uorf <- basename(env$pipeline_step_dependencies("uorf_regulation", root, runtime_inputs = list()))
  expect_false("required_fst_pages_summary.txt" %in% clean)
  expect_false("manual_translon_candidates.csv" %in% clean)
  expect_true("manual_translons_ranges.qs" %in% clean)
  expect_true("human_dominant_cell_states_clean_cds_gene_diagnostics.csv" %in% uorf)
  expect_false("human_dominant_cell_states_clean_cds_gene_expression.csv" %in% uorf)
})

test_that("pipeline thread defaults are bounded and explicit overrides are validated", {
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = dcs_project_root(), envir = env)
  expect_equal(env$pipeline_data_table_threads("", current = 24L), 4L)
  expect_equal(env$pipeline_data_table_threads("", current = 1L), 1L)
  expect_equal(env$pipeline_data_table_threads("24"), 24L)
  for (setting in c("0", "-1", "1.5", "abc")) {
    expect_error(env$pipeline_data_table_threads(setting), "positive integer")
  }
})

test_that("sourcing the runner from Rscript does not launch its CLI", {
  root <- dcs_project_root();tmp <- tempfile(fileext='.R')
  on.exit(unlink(tmp))
  writeLines(c(sprintf('source(%s)',encodeString(file.path(root,'scripts/run_dominant_cell_state_pipeline.R'),quote='"')),
    'cat("SOURCE_ONLY_OK\\n")'),tmp)
  result <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(tmp)),stdout=TRUE,stderr=TRUE)
  expect_true(any(grepl('SOURCE_ONLY_OK',result,fixed=TRUE)))
  expect_false(any(grepl('Dominant-cell-state pipeline',result,fixed=TRUE)))
})

test_that("no pipeline step runs or sources the one-off benchmark scripts", {
  root <- dcs_project_root()
  env <- new.env(parent = globalenv())
  dcs_load_pipeline(root = root, envir = env)
  sys.source(file.path(root, "scripts", "dominant_pipeline_cache.R"), envir = env)
  steps <- env$pipeline_steps(root)
  bench <- basename(list.files(file.path(root, "scripts", "benchmarks")))
  expect_true(length(bench) > 0)
  expect_false(any(basename(steps$script) %in% bench))
  present <- steps$script[file.exists(steps$script)]
  sourced <- unlist(lapply(present, env$pipeline_step_source_files, analysis_dir = root))
  expect_false(any(grepl("/benchmarks/", sourced, fixed = TRUE)))
  expect_false(any(basename(sourced) %in% bench))
})
