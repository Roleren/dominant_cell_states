Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1", ORFIK_MAPPING_BASELINE_CACHE = "rcpp-final")
source("scripts/benchmarks/benchmark_orfik_accessors_mapping.R")
out <- file.path(out, "source_verification")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (n in c(1L, 100L, 2000L)) {
  tx <- rep(suppressWarnings(do.call(c, real)), n)
  measure(paste0("to_real_exons_", length(tx)), old$pmapToTranscriptF,
          ORFik::pmapToTranscriptF, list(x = tx, transcripts = tx),
          batch = if (n == 1L) 5L else 1L, n = 5L)
  widths <- widthPerGroup(tx, FALSE)
  x <- IRanges(pmax(1L, widths - 10L), widths)
  measure(paste0("from_real_", length(tx), "_dropTRUE"), old$pmapFromTranscriptF,
          ORFik::pmapFromTranscriptF, list(x = x, transcripts = tx, removeEmpty = TRUE), n = 5L)
}
for (exons in c(3L, 20L, 200L)) {
  starts <- seq.int(100L, by = 30L, length.out = exons)
  for (direction in c("+", "-")) {
    s <- if (direction == "+") starts else rev(starts)
    args <- list(xStart = rep.int(tail(s, 1L) + 1L, 100000L),
      xEnd = rep.int(tail(s, 1L) + 5L, 100000L), transcriptStart = s,
      transcriptEnd = s + 19L, indices = rep.int(1L, exons), direction = direction,
      nExonCumSum = rep.int(0L, 100000L))
    measure(paste0("to_native_100k_", exons, "exons_", direction),
      original_native$to, ORFik:::pmapToTranscriptsCPP, args, batch = 30L, n = 5L)
  }
  args <- list(xStart = rep.int((exons - 1L) * 20L + 2L, 10000L),
    xEnd = rep.int(exons * 20L - 1L, 10000L),
    transcriptStart = rep.int(starts, 10000L), transcriptEnd = rep.int(starts + 19L, 10000L),
    indices = rep(seq_len(10000L), each = exons), direction = "+", removeEmpty = TRUE)
  measure(paste0("from_native_10k_", exons, "exons_dropTRUE"),
    original_native$from, ORFik:::pmapFromTranscriptsCPP, args,
    batch = if (exons == 3L) 30L else 3L, n = 5L)
}
for (exons in c(17L, 200L)) {
  tx <- GRangesList(GRanges("1", IRanges(seq.int(1L, by = 30L, length.out = exons), width = 20L), "+"))
  grl <- rep(tx, if (exons == 17L) 1000L else 5000L)
  for (name in c("strandPerGroup", "seqnamesPerGroup")) {
    measure(paste0(name, "_", length(grl), "x", exons, "exons_FALSE"), old[[name]],
      get(name, ns), list(grl = grl, keep.names = FALSE), batch = 100L, n = 5L)
  }
}
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
cat("Edited-source comparisons passed against independently compiled original native code.\n")
