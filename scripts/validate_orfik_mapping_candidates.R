Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1")
source("scripts/benchmark_orfik_accessors_mapping.R")
set.seed(510L)
checks <- 0L
for (direction in c("+", "-")) for (iteration in seq_len(100L)) {
  n <- sample.int(25L, 1L)
  lengths <- sample.int(50L, n, replace = TRUE)
  starts <- cumsum(c(100L, head(lengths, -1L) + sample.int(40L, n - 1L)))
  ends <- starts + lengths - 1L
  if (direction == "-") { starts <- rev(starts); ends <- rev(ends); lengths <- rev(lengths) }
  picks <- sample.int(n, 50L, replace = TRUE)
  qs <- starts[picks] + sample.int(10L, 50L, replace = TRUE) %% lengths[picks]
  qe <- pmin(qs + sample.int(10L, 50L, replace = TRUE) - 1L, ends[picks])
  args <- list(xStart = qs, xEnd = qe, transcriptStart = starts,
    transcriptEnd = ends, indices = rep.int(1L, n), direction = direction,
    nExonCumSum = rep.int(0L, 50L))
  a <- candidate_pmap_to(qs, qe, starts, ends, rep.int(1L, n), direction, rep.int(0L, 50L))
  prefix <- c(0L, head(cumsum(lengths), -1L))[picks]
  expected <- if (direction == "+") list(prefix + qs - starts[picks] + 1L,
    prefix + qe - starts[picks] + 1L) else list(prefix + ends[picks] - qe + 1L,
    prefix + ends[picks] - qs + 1L)
  stopifnot(identical(a$ranges, expected))
  stopifnot(identical(a, do.call(old$pmapToTranscriptsCPP, args)))
  for (drop in c(FALSE, TRUE)) {
    # Independently enumerate transcript bases, then recover exon-specific slices.
    first <- sample.int(sum(lengths), 1L)
    last <- first + sample.int(sum(lengths) - first + 1L, 1L) - 1L
    positions <- lapply(seq_len(n), function(j) {
      z <- seq.int(starts[j], ends[j]); if (direction == "-") rev(z) else z
    })
    bounds <- c(0L, cumsum(lengths))
    s <- en <- integer(n)
    for (j in seq_len(n)) {
      use <- intersect(first:last, seq.int(bounds[j] + 1L, bounds[j + 1L]))
      if (length(use)) {
        z <- positions[[j]][use - bounds[j]]
        s[j] <- min(z); en[j] <- max(z)
      }
    }
    keep <- if (drop) s > 0L else rep(TRUE, n)
    a <- candidate_pmap_from(first, last, starts, ends, rep.int(1L, n), direction, drop)
    stopifnot(identical(a$ranges, list(s[keep], en[keep])))
    stopifnot(identical(a, old$pmapFromTranscriptsCPP(first, last, starts, ends,
      rep.int(1L, n), direction, drop)))
    checks <- checks + 1L
  }
  checks <- checks + 1L
}
# Never call the old unsafe implementation on unmappable/malformed ranges.
bad <- candidate_pmap_to(c(205L, 400L, 150L, 0L, NA_integer_),
  c(210L, 405L, 155L, 1L, NA_integer_), c(100L, 200L), c(120L, 220L),
  c(1L, 2L), "+", rep.int(0L, 5L))
stopifnot(identical(bad$ranges, list(integer(5L), rep.int(-1L, 5L))))
empty <- candidate_pmap_to(integer(), integer(), integer(), integer(), integer(), "+", integer())
stopifnot(identical(empty$ranges, list(integer(), integer())))
stopifnot(inherits(try(candidate_pmap_from(1L, 100L, 100L, 120L, 1L, "+", TRUE), silent = TRUE), "try-error"))
cat(checks, "independent randomized native mapping checks passed, plus bounded/unmapped/empty checks.\n")
