Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1")
source("scripts/benchmark_orfik_accessors_mapping.R")
out <- file.path(out, "accessor_crossover")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (n in c(2048L, 4096L, 8192L)) {
  tx <- GRangesList(GRanges("1", IRanges(seq.int(1L, by = 30L, length.out = 33L), width = 20L), "+"))
  grl <- rep(tx, n)
  direct <- function(grl, keep.names = FALSE) as.character(grl@unlistData@strand[grl@partitioning@end])
  measure(paste0("strand_", n, "x33_exons"), old$strandPerGroup, direct,
          list(grl = grl, keep.names = FALSE), batch = 100L, n = 5L)
}
