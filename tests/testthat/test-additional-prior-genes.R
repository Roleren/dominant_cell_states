test_that("tracked and legacy additional priorities are combined without scoring", {
  skip_if_not_installed("data.table")
  env <- new.env(parent = globalenv())
  source(dcs_source_path("scripts", "dominant_state_module_definitions.R"),
         local = env)
  root <- tempfile("dcs-prior-inputs-")
  dir.create(file.path(root, "curated_inputs"), recursive = TRUE)
  dir.create(file.path(root, "results", "curated_inputs"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  withr::local_envvar(DOMINANT_CELL_STATES_RESULTS_DIR = NA_character_)
  data.table::fwrite(data.frame(gene_symbol = c("BOD1L1", "FMR1")),
    file.path(root, "curated_inputs", "curated_additional_gene_priorities.csv"))
  data.table::fwrite(data.frame(gene_symbol = c("FMR1", "PELO")),
    file.path(root, "results", "curated_inputs", "curated_additional_gene_priorities.csv"))
  expect_identical(env$load_additional_prior_gene_set(root),
                   c("BOD1L1", "FMR1", "PELO"))
  modules <- env$base_dominant_state_modules()
  expect_true(all(c("BOD1L1", "FMR1", "PELO") %in%
    env$dominant_state_target_genes(modules, root)))
  expect_false(any(c("BOD1L1", "FMR1") %in%
    unlist(c(modules$up_genes, modules$down_genes))))
})
