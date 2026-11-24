# Screen the other coverage/mapping calls used by RiboCrypt and this pipeline.
Sys.setenv(ORFIK_BENCHMARK_DEFINITIONS_ONLY = "1")
source("scripts/benchmark_orfik_hotpaths.R")
records <- lapply(c("ATF4", "FMR1", "BOD1L1"), function(gene) {
  x <- qs2::qs_read(list.files(".coverage_cache", pattern = paste0("^", gene, "__.*\\.qs$"), full.names = TRUE)[1])
  s <- x$signature
  gr <- GRanges(s$seqnames, IRanges(s$start, s$end), s$strand)
  p <- unlist(lapply(seq_along(gr), function(i) {
    z <- seq.int(start(gr)[i], end(gr)[i]); if (as.character(strand(gr)[i]) == "-") rev(z) else z
  }))
  reads <- GRanges(rep(as.character(seqnames(gr)), width(gr)), IRanges(p, width = 1L),
                   rep(as.character(strand(gr)), width(gr)))
  reads$score <- rowSums(as.matrix(x$coverage))
  names(gr) <- NULL
  tx <- GRangesList(gr); names(tx) <- s$transcript_names
  list(tx = tx, reads = reads)
})
tx <- do.call(c, lapply(records, `[[`, "tx"))
reads <- do.call(c, lapply(records, `[[`, "reads"))
seqlengths(reads) <- vapply(seqlevels(reads), function(chr) max(end(reads)[as.character(seqnames(reads)) == chr]), 0L)
seqinfo(tx) <- seqinfo(reads)
cov <- covRleFromGR(reads, weight = "score")
tasks <- list(
  coverageByTranscriptC = list(x = cov, transcripts = tx),
  coverageByTranscriptSum = list(x = cov, transcripts = tx),
  pmapToTranscriptF = list(x = tx, transcripts = tx),
  coveragePerTiling_covRle = list(grl = tx, reads = cov, is.sorted = TRUE,
                                 as.data.table = TRUE, withFrames = TRUE),
  coveragePerTiling_weighted_reads = list(grl = tx, reads = reads, is.sorted = TRUE,
                                        as.data.table = TRUE, withFrames = TRUE))
results <- list()
if (Sys.getenv("ORFIK_SCREEN_ONLY_CANDIDATE") != "1") {
for (size in c(1L, 100L)) for (label in names(tasks)) {
  args <- tasks[[label]]
  if (size > 1L) {
    if (!is.null(args$transcripts)) args$transcripts <- rep(args$transcripts, size)
    if (label == "pmapToTranscriptF") args$x <- rep(args$x, size)
    if (!is.null(args$grl)) args$grl <- rep(args$grl, size)
  }
  f <- if (grepl("coveragePerTiling", label)) ORFik::coveragePerTiling else get(label, ns)
  invisible(do.call(f, args))
  profile <- file.path(out, paste0("screen_", label, "_", size, ".out"))
  Rprof(profile, interval = .01)
  for (i in seq_len(10L)) invisible(do.call(f, args))
  Rprof(NULL)
  fwrite(as.data.table(summaryRprof(profile)$by.self, keep.rownames = "call"), paste0(profile, ".csv"))
  for (iteration in seq_len(7L)) {
    gc()
    elapsed <- system.time(for (i in seq_len(5L)) invisible(do.call(f, args)))[["elapsed"]] / 5L
    results[[length(results) + 1L]] <- data.table(function_name = label,
      transcripts = length(tx) * size, iteration = iteration, seconds = elapsed)
  }
  fwrite(rbindlist(results), file.path(out, "screen_timings.csv"))
}
}
out <- file.path(out, "coverage_extraction")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (size in c(1L, 100L, 2000L)) {
  bench(paste0("real_coordinates_", 3L * size, "_transcripts"), "coverageByTranscriptC",
    old$coverageByTranscriptC, candidate_cov,
    list(x = cov, transcripts = rep(tx, size)), n = 7L, batch = 3L)
}
if (length(results)) print(rbindlist(results)[, .(median_seconds = median(seconds)), by = .(function_name, transcripts)])
