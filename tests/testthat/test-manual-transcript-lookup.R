test_that("manual transcript lookup preserves measured choices and fills missing genes", {
  skip_if_not_installed("data.table")
  env <- new.env(parent = globalenv())
  env$`%chin%` <- data.table::`%chin%`
  for (name in c("fread", "setnames", "data.table", "as.data.table", "fifelse")) {
    env[[name]] <- getExportedValue("data.table", name)
  }
  code <- parse(dcs_source_path("scripts", "dominant_manual_translons.R"))
  for (expr in code) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% c("normalize_manual_type", "read_manual_additions")) {
      eval(expr, env)
    }
  }
  root <- tempfile("manual-lookup-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  queue <- file.path(root, "queue.txt")
  writeLines(c("gene,start,stop,type", "FMR1,18,98,uORF",
               "PELO,57,113,uORF"), queue)
  data.table::fwrite(data.frame(gene_symbol = "FMR1", tx_id = "measured_fmr1"),
    file.path(root, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"))
  fallback <- data.table::data.table(gene_symbol = c("FMR1", "PELO"),
                                    tx_id = c("other_fmr1", "canonical_pelo"))
  rows <- env$read_manual_additions(queue, root, fallback)
  expect_identical(rows$tx_id[rows$gene_symbol == "FMR1"], "measured_fmr1")
  expect_identical(rows$tx_id[rows$gene_symbol == "PELO"], "canonical_pelo")
  expect_identical(rows$lookup_source[rows$gene_symbol == "PELO"], "canonical_annotation")
  ambiguous <- rbind(fallback, data.table::data.table(gene_symbol = "PELO", tx_id = "other_pelo"))
  expect_error(env$read_manual_additions(queue, root, ambiguous), "Ambiguous canonical")
  unlink(file.path(root, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"))
  expect_equal(nrow(env$read_manual_additions(queue, root, fallback)), 2L)
})
