Sys.setenv(ORFIK_BENCHMARK_DEFINITIONS_ONLY = "1")
source("scripts/benchmarks/benchmark_orfik_hotpaths.R")
root <- tempfile("fst-regression-")
dir.create(root)
idx <- rbindlist(lapply(c("1", "2"), function(chr) {
  rbindlist(lapply(c(1L, 11L), function(begin) {
    paths <- file.path(root, paste0(chr, "_", begin, "_", c("plus", "minus"), ".fst"))
    offset <- if (chr == "1") 0L else 100L
    fst::write_fst(data.table(A = offset + begin:(begin + 9L), B = -begin:-(begin + 9L)), paths[1])
    fst::write_fst(data.table(A = offset + 1000L + begin:(begin + 9L), B = 2L * begin:(begin + 9L)), paths[2])
    data.table(chr = chr, start = begin, end = begin + 9L,
               file_forward = paths[1], file_reverse = paths[2])
  }))
}))
index <- file.path(root, "coverage_index.fst")
fst::write_fst(idx, index)
grl <- GRangesList(plus = GRanges("1", IRanges(c(12, 3), c(14, 5)), "+"),
  minus = GRanges("2", IRanges(c(15, 2), c(17, 4)), "-"),
  crossing = GRanges("1", IRanges(8, 13), "+"),
  duplicate = GRanges("1", IRanges(c(3, 3), c(5, 5)), "+"))
for (cols in list(NULL, "B", c("B", "A"))) {
  a <- old$coverageByTranscriptFST(grl, index, cols)
  b <- candidate_fst(grl, index, cols)
  stopifnot(isTRUE(all.equal(a, b, check.attributes = TRUE, tolerance = 0)))
}
set.seed(817L)
for (i in seq_len(100L)) {
  cov <- RleList(lapply(seq_len(sample.int(20L, 1)), function(i) {
    z <- sample(c(0L, 0L, 0L, 2L, 4L), sample.int(100L, 1), replace = TRUE)
    if (i %% 2L) z <- z / 3
    Rle(z)
  }))
  names(cov) <- sample(paste0("tx", seq_along(cov)))
  for (keep in c(TRUE, FALSE)) for (frames in c(TRUE, FALSE)) for (drop in c(TRUE, FALSE)) {
    args <- list(coverage = cov, drop.zero.dt = drop, keep.names = keep,
                 withFrames = frames, fraction = 28L)
    stopifnot(isTRUE(all.equal(suppressWarnings(do.call(old$coverage_to_dt, args)),
      suppressWarnings(do.call(candidate_dt, args)), check.attributes = TRUE, tolerance = 0)))
  }
}
unlink(root, recursive = TRUE)
cat("Candidate equivalence passed: multi-page/strand/column FST inputs and 800 randomized sparse/dense profile combinations.\n")
