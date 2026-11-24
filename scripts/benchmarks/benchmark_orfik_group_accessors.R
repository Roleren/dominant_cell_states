Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1")
source("scripts/benchmarks/benchmark_orfik_accessors_mapping.R")
out <- file.path(out, "group_accessors")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fast_strand <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (!keep.names) return(as.character(grl@unlistData@strand[grl@partitioning@end]))
  grl@unlistData@strand[grl@partitioning@end]
}
fast_seqnames <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (!keep.names) return(as.character(grl@unlistData@seqnames[grl@partitioning@end]))
  grl@unlistData@seqnames[grl@partitioning@end]
}
environment(fast_strand) <- environment(fast_seqnames) <- ns
stopifnot(!identical(body(fast_strand), body(old$strandPerGroup)))
fwrite(data.table(function_name = c("strandPerGroup", "seqnamesPerGroup"),
  earlier_prototype_changed = c(!identical(body(e$strandPerGroup), body(old$strandPerGroup)),
    !identical(body(e$seqnamesPerGroup), body(old$seqnamesPerGroup)))),
  file.path(out, "prototype_transform_check.csv"))
for (exons in c(17L, 200L)) {
  tx <- GRangesList(GRanges("1", IRanges(seq.int(1L, by = 30L, length.out = exons), width = 20L), "+"))
  grl <- rep(tx, if (exons == 17L) 1000L else 5000L)
  for (name in c("strandPerGroup", "seqnamesPerGroup")) {
    improved <- if (name == "strandPerGroup") fast_strand else fast_seqnames
    measure(paste0(name, "_", length(grl), "x", exons, "exons"), old[[name]],
      improved, list(grl = grl, keep.names = FALSE), batch = 100L)
  }
}
