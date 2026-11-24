# Accessor and native mapping candidates, before source edits.
devtools::load_all(path.expand("~/livemount/forks/ORFik"), compile = FALSE)
devtools::load_all(path.expand("~/livemount/forks/RiboCrypt"))
library(data.table)
setDTthreads(4L)
out <- Sys.getenv("ORFIK_MAPPING_BENCH_DIR", "results/orfik_mapping_2026-10-10")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("ORFik")
file <- file.path(out, "baseline.qs")
if (!file.exists(file)) {
  if (Sys.getenv("ORFIK_MAPPING_CAPTURE_BASELINE") != "1") stop("Missing pre-edit mapping snapshot")
  qs2::qs_save(mget(c("pmapToTranscriptF", "pmapFromTranscriptF",
  "pmapToTranscriptsCPP", "pmapFromTranscriptsCPP", "groupings", "widthPerGroup",
  "strandPerGroup", "seqnamesPerGroup"), ns), file)
}
old <- qs2::qs_read(file)
for (i in seq_along(old)) environment(old[[i]]) <- ns
source("scripts/benchmarks/orfik_mapping_benchmark_baseline.R")
original_native <- load_original_mapping_native(out)
old$pmapToTranscriptsCPP <- original_native$to
old$pmapFromTranscriptsCPP <- original_native$from
b <- new.env(parent = ns)
for (name in names(old)) b[[name]] <- old[[name]]
b$strandBool <- ORFik::strandBool
environment(b$strandBool) <- b
environment(old$pmapToTranscriptF) <- environment(old$pmapFromTranscriptF) <- b
if (Sys.getenv("ORFIK_MAPPING_SKIP_COMPILE") != "1")
  Rcpp::sourceCpp("scripts/benchmarks/benchmark_orfik_mapping_candidate.cpp", cacheDir = file.path(out, "rcpp"))

fast_groupings <- function(grl) {
  ends <- grl@partitioning@end
  rep.int(seq_along(ends), diff(c(0L, ends)))
}
fast_strand <- old$strandPerGroup
body(fast_strand) <- parse(text = gsub(
  "as.character(grl@unlistData@strand)[grl@partitioning@end]",
  "as.character(grl@unlistData@strand[grl@partitioning@end])",
  paste(deparse(body(fast_strand)), collapse = "\n"), fixed = TRUE))[[1]]
fast_seqnames <- old$seqnamesPerGroup
body(fast_seqnames) <- parse(text = gsub(
  "as.character(grl@unlistData@seqnames)[grl@partitioning@end]",
  "as.character(grl@unlistData@seqnames[grl@partitioning@end])",
  paste(deparse(body(fast_seqnames)), collapse = "\n"), fixed = TRUE))[[1]]
fast_widths <- old$widthPerGroup
text <- paste(deparse(body(fast_widths)), collapse = "\n")
text <- gsub("end(grl@partitioning)", "grl@partitioning@end", text, fixed = TRUE)
text <- gsub("names(res) <- names(grl)", "names(res) <- grl@partitioning@NAMES", text, fixed = TRUE)
body(fast_widths) <- parse(text = text)[[1]]
environment(fast_groupings) <- ns

e <- new.env(parent = ns)
e$pmapToTranscriptsCPP <- candidate_pmap_to
e$groupings <- fast_groupings; e$strandPerGroup <- fast_strand
e$seqnamesPerGroup <- fast_seqnames; e$widthPerGroup <- fast_widths
candidate <- old$pmapToTranscriptF
environment(candidate) <- e

timings <- list()
measure <- function(label, baseline, improved, args, batch = 1L, n = 7L) {
  a <- do.call(baseline, args); b <- do.call(improved, args)
  stopifnot(isTRUE(all.equal(a, b, check.attributes = TRUE, tolerance = 0)))
  for (iteration in seq_len(n)) for (method in if (iteration %% 2) c("baseline", "candidate") else c("candidate", "baseline")) {
    f <- if (method == "baseline") baseline else improved
    gc()
    seconds <- system.time(for (j in seq_len(batch)) invisible(do.call(f, args)),
                           gcFirst = FALSE)[["elapsed"]] / batch
    timings[[length(timings) + 1L]] <<- data.table(input = label, method = method,
      iteration = iteration, batch = batch, seconds = seconds)
  }
  d <- rbindlist(timings)
  fwrite(d, file.path(out, "timings.csv"))
  s <- d[, .(seconds = median(seconds), min = min(seconds), max = max(seconds)), by = .(input, method)]
  fwrite(s, file.path(out, "summary.csv")); print(s[input == label])
}

real <- lapply(c("ATF4", "FMR1", "BOD1L1"), function(gene) {
  x <- qs2::qs_read(list.files(".coverage_cache", pattern = paste0("^", gene, "__.*\\.qs$"), full.names = TRUE)[1])
  s <- x$signature
  grl <- GRangesList(GRanges(s$seqnames, IRanges(s$start, s$end), s$strand))
  names(grl) <- s$transcript_names
  grl
})

if (Sys.getenv("ORFIK_MAPPING_DEFINITIONS_ONLY") != "1") {
for (n in c(1L, 100L, 2000L)) {
  tx <- rep(suppressWarnings(do.call(c, real)), n)
  measure(paste0("to_real_exons_", length(tx)), old$pmapToTranscriptF, candidate,
          list(x = tx, transcripts = tx), batch = if (n == 1L) 5L else 1L)
}
grl <- rep(real[[2]], 1000L)
for (name in c("groupings", "strandPerGroup", "seqnamesPerGroup", "widthPerGroup")) {
  improved <- switch(name, groupings = fast_groupings, strandPerGroup = fast_strand,
    seqnamesPerGroup = fast_seqnames, widthPerGroup = fast_widths)
  measure(paste0(name, "_1000x17exons"), old[[name]], improved, list(grl = grl), batch = 100L)
}
# Equal-result accessor expressions, with sufficient calls for microsecond timing.
measure("names_1000_transcripts", function(x) names(x), function(x) x@partitioning@NAMES,
        list(x = grl), batch = 100000L, n = 5L)
measure("lengths_1000_transcripts", function(x) lengths(x, use.names = FALSE),
        function(x) diff(c(0L, x@partitioning@end)), list(x = grl), batch = 5000L, n = 5L)
for (exons in c(3L, 20L, 200L)) {
  starts <- seq.int(100L, by = 30L, length.out = exons)
  for (direction in c("+", "-")) {
    s <- if (direction == "+") starts else rev(starts)
    args <- list(xStart = rep.int(tail(starts, 1L) + 1L, 100000L),
      xEnd = rep.int(tail(starts, 1L) + 5L, 100000L),
      transcriptStart = s, transcriptEnd = s + 19L,
      indices = rep.int(1L, exons), direction = direction,
      nExonCumSum = rep.int(0L, 100000L))
    # Last exon in transcript order, for both strands.
    args$xStart[] <- tail(s, 1L) + 1L; args$xEnd[] <- tail(s, 1L) + 5L
    measure(paste0("native_100k_queries_", exons, "exons_", direction),
      old$pmapToTranscriptsCPP, candidate_pmap_to, args, batch = 3L)
  }
}
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
}
