test_that("checked transcript mapping rejects unmappable and unpaired ranges", {
  skip_if_not_installed("ORFik")
  suppressWarnings(suppressPackageStartupMessages(library(ORFik)))
  source(file.path(dcs_project_root(), "scripts", "dominant_transcript_mapping.R"))
  plus <- GRanges("chr1", IRanges(c(100, 200), c(109, 219)), strand = "+")
  minus <- GRanges("chr1", IRanges(c(200, 100), c(219, 109)), strand = "-")
  query <- GRangesList(
    GRanges("chr1", IRanges(101, 103), strand = "+"),
    GRanges("chr1", IRanges(205, 208), strand = "-"),
    GRanges("chr1", IRanges(40, 45), strand = "+"),
    GRanges("chr1", IRanges(c(105, 200), c(109, 203)), strand = "+"))
  subject <- GRangesList(plus, minus, plus, plus)
  mapped <- dominant_checked_pmap(query, subject)
  expect_equal(start(mapped[[1]]), 2L)
  expect_equal(end(mapped[[1]]), 4L)
  expect_equal(start(mapped[[2]]), 12L)
  expect_equal(end(mapped[[2]]), 15L)
  expect_length(mapped[[3]], 0L)
  expect_equal(start(mapped[[4]]), 6L)
  expect_equal(end(mapped[[4]]), 14L)
  scalar <- GRangesList(lapply(seq_along(query), function(i) {
    dominant_checked_pmap(query[i], subject[i])[[1]]
  }))
  expect_identical(start(unlist(mapped)), start(unlist(scalar)))
  expect_identical(end(unlist(mapped)), end(unlist(scalar)))
  swapped <- dominant_checked_pmap(
    GRangesList(GRanges("chr1", IRanges(401, 403), strand = "+"),
                GRanges("chr1", IRanges(101, 103), strand = "+")),
    GRangesList(plus, GRanges("chr1", IRanges(400, 410), strand = "+")))
  expect_identical(lengths(swapped), c(0L, 0L))
  zero <- dominant_checked_pmap(
    GRangesList(GRanges("chr1", IRanges(220, width = 0), strand = "+")),
    GRangesList(plus))
  expect_length(zero[[1]], 0L)
  for (script in c("dominant_cell_state_required_fst_pages.R",
      "dominant_cell_states_fst_clean_cds.R", "dominant_uorf_regulation_inference.R",
      "dominant_ribosome_decision_graphs.R")) {
    definition <- Filter(function(expr) is.call(expr) &&
      identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("collapse_to_tx_iranges")),
      as.list(parse(file.path(dcs_project_root(), "scripts", script))))
    scope <- new.env(parent = environment())
    eval(definition[[1]], scope)
    ranges <- scope$collapse_to_tx_iranges(mapped)
    expect_equal(start(ranges)[3], 0L)
    expect_equal(width(ranges)[3], 0L)
  }
})


test_that("FST selection loads checked mapping before categorizing ranges", {
  code <- parse(file.path(dcs_project_root(), "scripts", "dominant_cell_state_required_fst_pages.R"))
  scope <- new.env(parent = globalenv())
  scope$analysis_dir <- dcs_project_root()
  helper <- which(vapply(code, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("source")) &&
      grepl("dominant_transcript_mapping.R", paste(deparse(expr), collapse = " "), fixed = TRUE)
  }, logical(1)))
  first_use <- which(vapply(code, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("annotate_translon_categories"))
  }, logical(1)))
  expect_length(helper, 1L)
  expect_true(helper < first_use)
  eval(code[[helper]], envir = scope)
  expect_true(exists("dominant_checked_pmap", envir = scope, inherits = FALSE))
})
