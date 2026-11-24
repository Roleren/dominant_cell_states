# Coverage-table stress cases, independent of the external megafst timings.
Sys.setenv(ORFIK_BENCHMARK_DEFINITIONS_ONLY = "1")
source("scripts/benchmark_orfik_hotpaths.R")
out <- file.path(out, "coverage_tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(421L)
for (density in c(.01, .5)) {
  cov <- RleList(lapply(seq_len(100L), function(i) {
    z <- integer(5000L); k <- sample.int(length(z), length(z) * density)
    z[k] <- sample.int(20L, length(k), replace = TRUE); Rle(z)
  }))
  names(cov) <- paste0("tx", seq_along(cov))
  bench(paste0("100x5000_density_", density), "coverage_to_dt",
    old$coverage_to_dt, candidate_dt, list(coverage = cov, drop.zero.dt = TRUE,
                                        withFrames = TRUE))
  variable <- RleList(lapply(seq_along(cov), function(i) cov[[i]][seq_len(5000L - i)]))
  bench(paste0("100_variable_windows_density_", density), "coverage_to_dt",
    old$coverage_to_dt, candidate_dt, list(coverage = variable, drop.zero.dt = FALSE,
                                        withFrames = TRUE))
}
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
