test_that("adding a gene does not recompute existing measurements", {
  source(file.path(dcs_project_root(), "scripts", "dominant_object_cache.R"))
  tmp <- tempfile(); on.exit(unlink(tmp, recursive = TRUE))
  calls <- 0L
  measure <- function(gene, coords = 1:3, coverage = 1L) {
    dominant_cached_object(tmp, gene, list(coords = coords, coverage = coverage),
      function() { calls <<- calls + 1L; sum(coords) })
  }
  expect_equal(measure("A"), 6)
  expect_equal(measure("B"), 6)
  expect_equal(measure("A"), 6)
  expect_equal(measure("B"), 6)
  expect_equal(measure("C"), 6)
  expect_identical(calls, 3L)
  expect_equal(measure("A", 2:4), 9)
  expect_identical(calls, 4L)
  measure("B", coverage = 2L)
  expect_identical(calls, 5L)
})

test_that("changes in calculation helpers invalidate derived measurements", {
  source(file.path(dcs_project_root(), "scripts", "dominant_object_cache.R"))
  scope <- new.env(parent = baseenv())
  scope$helper <- function(x) x + 1
  scope$calculate <- eval(quote(function(x) helper(x)), scope)
  before <- dominant_computation_signature(scope$calculate)
  scope$unrelated_plot <- function(x) x * 2
  expect_identical(dominant_computation_signature(scope$calculate), before)
  scope$helper <- function(x) x + 2
  expect_false(identical(dominant_computation_signature(scope$calculate), before))
})

test_that("derived measurement caches can be bypassed explicitly", {
  source(file.path(dcs_project_root(), "scripts", "dominant_object_cache.R"))
  tmp <- tempfile(); on.exit(unlink(tmp, recursive = TRUE))
  calls <- 0L
  get <- function(use_cache) dominant_cached_object(tmp, "test", 1,
    function() { calls <<- calls + 1L; calls }, use_cache = use_cache)
  expect_equal(get(TRUE), 1L)
  expect_equal(get(TRUE), 1L)
  expect_equal(get(FALSE), 2L)
  expect_equal(get(TRUE), 1L)
})

test_that("serialized data.table signatures reuse values after a QS round trip", {
  source(file.path(dcs_project_root(), "scripts", "dominant_object_cache.R"))
  tmp <- tempfile(); on.exit(unlink(tmp, recursive = TRUE))
  calls <- 0L
  measure <- function() dominant_cached_object(tmp, "qc",
    list(exon_qc = data.table::data.table(tx_id = character(), score = numeric())),
    function() { calls <<- calls + 1L; calls })
  expect_equal(measure(), 1L)
  expect_equal(measure(), 1L)
  expect_identical(calls, 1L)
})

test_that("biologically empty results can be cached without rerunning coverage", {
  source(file.path(dcs_project_root(), "scripts", "dominant_object_cache.R"))
  tmp <- tempfile(); on.exit(unlink(tmp, recursive = TRUE))
  calls <- 0L
  measure <- function() dominant_cached_object(tmp, "empty", list(translons = integer()),
    function() { calls <<- calls + 1L; NULL })
  expect_null(measure())
  expect_null(measure())
  expect_identical(calls, 1L)
})

test_that("confirmed ORFik fixes share a cache version without masking future updates", {
  expect_identical(dominant_orfik_cache_version("1.29.12"), "1.29.14")
  expect_identical(dominant_orfik_cache_version("1.29.14"), "1.29.14")
  expect_identical(dominant_orfik_cache_version("1.29.15"), "1.29.14")
  expect_identical(dominant_orfik_cache_version("1.29.16"), "1.29.16")
})
