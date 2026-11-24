Sys.setenv(ORFIK_BENCHMARK_DEFINITIONS_ONLY = "1")
source("scripts/benchmarks/benchmark_orfik_hotpaths.R")
answer <- list()
for (gene in c("ATF4", "FMR1", "BOD1L1")) {
  x <- qs2::qs_read(list.files(".coverage_cache", pattern = paste0("^", gene, "__.*\\.qs$"), full.names = TRUE)[1])
  s <- x$signature
  grl <- GRangesList(GRanges(s$seqnames, IRanges(s$start, s$end), s$strand))
  names(grl) <- s$transcript_names
  for (method in c("baseline", "candidate")) {
    calls <- list()
    e <- new.env(parent = ns)
    e$read_fst <- function(path, columns = NULL, from = 1, to = NULL, ...) {
      if (!identical(path, s$fst_index)) calls[[length(calls) + 1L]] <<-
        data.table(path = path, from = from, to = to, rows = to - from + 1L)
      fst::read_fst(path, columns = columns, from = from, to = to, ...)
    }
    f <- if (method == "baseline") old$coverageByTranscriptFST else candidate_fst
    if (method == "candidate") {
      e$read_fst_tiles <- read_fst_tiles
      environment(e$read_fst_tiles) <- e
    }
    environment(f) <- e
    value <- f(grl, s$fst_index, head(names(x$coverage), 4L))
    stopifnot(isTRUE(all.equal(value[[1]], x$coverage[, 1:4, with = FALSE], check.attributes = TRUE, tolerance = 0)))
    d <- rbindlist(calls)
    answer[[length(answer) + 1L]] <- data.table(gene = gene, method = method,
      page_reads = nrow(d), rows_read = sum(d$rows), largest_read_rows = max(d$rows))
  }
}
fwrite(rbindlist(answer), file.path(out, "fst_read_audit.csv"))
print(rbindlist(answer))
