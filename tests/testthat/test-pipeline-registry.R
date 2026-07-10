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
