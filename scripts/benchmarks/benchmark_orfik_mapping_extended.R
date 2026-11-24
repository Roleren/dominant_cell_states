Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1")
source("scripts/benchmarks/benchmark_orfik_accessors_mapping.R")
out <- file.path(out, "extended")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Large annotation accessor inputs, with both short and long Rle runs.
for (genes in c(5000L)) {
  tx <- GRangesList(GRanges("1", IRanges(seq.int(1L, by = 30L, length.out = 200L), width = 20L), "+"))
  grl <- rep(tx, genes)
  for (name in c("strandPerGroup", "seqnamesPerGroup")) {
    improved <- if (name == "strandPerGroup") fast_strand else fast_seqnames
    measure(paste0(name, "_5000x200exons"), old[[name]], improved, list(grl = grl), batch = 10L)
  }
}
for (exons in c(3L, 20L, 200L)) for (direction in c("+", "-")) {
  ng <- 10000L
  s <- seq.int(100L, by = 30L, length.out = exons)
  if (direction == "-") s <- rev(s)
  for (drop in c(FALSE, TRUE)) {
    args <- list(xStart = rep.int((exons - 1L) * 20L + 2L, ng),
      xEnd = rep.int(exons * 20L - 1L, ng),
      transcriptStart = rep.int(s, ng), transcriptEnd = rep.int(s + 19L, ng),
      indices = rep(seq_len(ng), each = exons), direction = direction, removeEmpty = drop)
    measure(paste0("from_native_10k_", exons, "exons_", direction, "_drop", drop),
      old$pmapFromTranscriptsCPP, candidate_pmap_from, args, batch = if (exons == 3L) 3L else 1L, n = 5L)
  }
}
e$pmapFromTranscriptsCPP <- candidate_pmap_from
from <- old$pmapFromTranscriptF; environment(from) <- e
for (n in c(1L, 100L, 2000L)) {
  tx <- rep(suppressWarnings(do.call(c, real)), n)
  widths <- widthPerGroup(tx, FALSE)
  x <- IRanges(pmax(1L, widths - 10L), widths)
  measure(paste0("from_real_", length(tx), "_dropTRUE"), old$pmapFromTranscriptF,
    from, list(x = x, transcripts = tx, removeEmpty = TRUE), n = 7L)
}
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
