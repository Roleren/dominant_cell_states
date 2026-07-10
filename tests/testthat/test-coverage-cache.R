test_that("coverage cache reuses transcript windows until index signature changes", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
  })

  root <- dcs_project_root()
  source(file.path(root, "scripts", "dominant_coverage_cache.R"))

  tmp <- tempfile("dominant-coverage-cache-test-")
  dir.create(tmp, recursive = TRUE)
  fst_index <- file.path(tmp, "coverage_index.fst")
  writeLines("index-v1", fst_index)

  display <- GenomicRanges::GRangesList(
    TX1 = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(start = 101L, end = 110L),
      strand = "+"
    )
  )

  calls <- 0L
  fake_reader <- function(display_region, fst_index, columns = NULL) {
    calls <<- calls + 1L
    list(matrix(
      seq_len(20),
      nrow = 10,
      dimnames = list(NULL, c("RUN1", "RUN2"))
    ))
  }

  first <- dominant_cached_coverage_by_transcript(
    display,
    fst_index,
    gene_symbol = "GENE1",
    tx_id = "TX1",
    analysis_dir = tmp,
    reader = fake_reader
  )
  second <- dominant_cached_coverage_by_transcript(
    display,
    fst_index,
    gene_symbol = "GENE1",
    tx_id = "TX1",
    analysis_dir = tmp,
    reader = fake_reader
  )

  expect_identical(calls, 1L)
  expect_identical(first, second)

  writeLines(c("index-v1", "index-v2"), fst_index)
  third <- dominant_cached_coverage_by_transcript(
    display,
    fst_index,
    gene_symbol = "GENE1",
    tx_id = "TX1",
    analysis_dir = tmp,
    reader = fake_reader
  )

  expect_identical(calls, 2L)
  expect_identical(first, third)
})
