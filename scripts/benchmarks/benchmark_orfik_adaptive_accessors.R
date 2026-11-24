Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1")
source("scripts/benchmarks/benchmark_orfik_accessors_mapping.R")
out <- file.path(out, "adaptive_accessors")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
make_adaptive <- function(field) {
  f <- function(grl, keep.names = TRUE) {
    validGRL(class(grl))
    positions <- grl@partitioning@end
    values <- grl@unlistData@strand
    if (!keep.names) {
      ng <- length(positions)
      ne <- if (ng) positions[ng] else 0L
      if (ne > 65536L && ne > 32 * ng) return(as.character(values[positions]))
      return(as.character(values)[positions])
    }
    values[positions]
  }
  if (field == "seqnames") {
    txt <- paste(deparse(body(f)), collapse = "\n")
    body(f) <- parse(text = gsub("@strand", "@seqnames", txt, fixed = TRUE))[[1]]
  }
  environment(f) <- ns
  f
}
for (exons in c(17L, 200L)) {
  tx <- GRangesList(GRanges("1", IRanges(seq.int(1L, by = 30L, length.out = exons), width = 20L), "+"))
  grl <- rep(tx, if (exons == 17L) 1000L else 5000L)
  for (name in c("strandPerGroup", "seqnamesPerGroup")) {
    f <- make_adaptive(if (name == "strandPerGroup") "strand" else "seqnames")
    measure(paste0(name, "_", length(grl), "x", exons, "exons"), old[[name]],
      f, list(grl = grl, keep.names = FALSE), batch = 100L, n = 5L)
  }
}
