# Recheck the edited package against the saved unmodified functions in a fresh
# session. Run only after the candidate benchmarks have demonstrated speedups.
base_dir <- Sys.getenv("ORFIK_BENCHMARK_DIR", "results/orfik_performance_2026-10-09")
final_dir <- file.path(base_dir, "source_verification")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(file.path(base_dir, "baseline_functions.qs")))
file.copy(file.path(base_dir, "baseline_functions.qs"),
          file.path(final_dir, "baseline_functions.qs"), overwrite = TRUE)
Sys.setenv(ORFIK_BENCHMARK_DEFINITIONS_ONLY = "1", ORFIK_BENCHMARK_DIR = final_dir)
source("scripts/benchmark_orfik_hotpaths.R")
for (gene in c("ATF4", "FMR1", "BOD1L1")) {
  x <- qs2::qs_read(list.files(".coverage_cache", pattern = paste0("^", gene, "__.*\\.qs$"), full.names = TRUE)[1])
  s <- x$signature
  grl <- GRangesList(GRanges(s$seqnames, IRanges(s$start, s$end), s$strand))
  names(grl) <- s$transcript_names
  a <- ORFik::coverageByTranscriptFST(grl, s$fst_index)[[1]]
  stopifnot(isTRUE(all.equal(a, x$coverage, check.attributes = TRUE, tolerance = 0)))
  bench(paste0(gene, "_4_libraries"), "coverageByTranscriptFST",
    old$coverageByTranscriptFST, ORFik::coverageByTranscriptFST,
    list(grl = grl, fst_index = s$fst_index, columns = head(names(x$coverage), 4L)), n = 7L, batch = 3L)
  bench(paste0(gene, "_all_libraries"), "coverageByTranscriptFST",
    old$coverageByTranscriptFST, ORFik::coverageByTranscriptFST,
    list(grl = grl, fst_index = s$fst_index), n = 5L)
  cov <- RleList(rowSums(as.matrix(x$coverage))); names(cov) <- names(grl)
  bench(paste0(gene, "_summed_profile"), "coverage_to_dt", old$coverage_to_dt,
    ORFik:::coverage_to_dt, list(coverage = cov, drop.zero.dt = TRUE,
                               withFrames = TRUE), n = 7L, batch = 20L)
}
# A diverse panel rather than duplicate coordinates: first 30 cached genes in
# filename order, crossing chromosomes, strands and exon structures.
files <- head(sort(list.files(".coverage_cache", pattern = "\\.qs$", full.names = TRUE)), 30L)
objects <- lapply(files, qs2::qs_read)
panel <- do.call(c, lapply(objects, function(x) {
  s <- x$signature
  z <- GRangesList(GRanges(s$seqnames, IRanges(s$start, s$end), s$strand))
  names(z) <- s$transcript_names; z
}))
columns <- head(names(objects[[1]]$coverage), 4L)
bench("30_gene_panel_4_libraries", "coverageByTranscriptFST",
  old$coverageByTranscriptFST, ORFik::coverageByTranscriptFST,
  list(grl = panel, fst_index = objects[[1]]$signature$fst_index, columns = columns), n = 7L)
fwrite(data.table(gene = vapply(objects, `[[`, "", "gene_symbol"),
  transcript = names(panel), exons = lengths(panel),
  nucleotides = vapply(panel, function(z) sum(width(z)), 0L)), file.path(out, "panel_inputs.csv"))
set.seed(421L)
for (density in c(.01, .5)) {
  cov <- RleList(lapply(seq_len(100L), function(i) {
    z <- integer(5000L); k <- sample.int(length(z), length(z) * density)
    z[k] <- sample.int(20L, length(k), replace = TRUE); Rle(z)
  }))
  names(cov) <- paste0("tx", seq_along(cov))
  bench(paste0("100x5000_density_", density), "coverage_to_dt",
    old$coverage_to_dt, ORFik:::coverage_to_dt,
    list(coverage = cov, drop.zero.dt = TRUE, withFrames = TRUE), n = 7L)
}
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
cat("Edited source reproduced saved production coverage and passed every timed output comparison.\n")
