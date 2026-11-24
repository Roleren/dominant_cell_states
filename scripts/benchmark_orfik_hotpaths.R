# Compare unmodified ORFik with in-memory candidates before editing its source.
# Run from dominant_cell_states; external coverage inputs are read-only.
devtools::load_all(path.expand("~/livemount/forks/ORFik"))
devtools::load_all(path.expand("~/livemount/forks/RiboCrypt"))
library(data.table)
setDTthreads(4L)
out <- Sys.getenv("ORFIK_BENCHMARK_DIR", "results/orfik_performance_2026-10-09")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("ORFik")
baseline_file <- file.path(out, "baseline_functions.qs")
if (!file.exists(baseline_file)) {
  if (Sys.getenv("ORFIK_BENCHMARK_CAPTURE_BASELINE") != "1")
    stop("Missing original-function snapshot. Before edits, run with ORFIK_BENCHMARK_CAPTURE_BASELINE=1.")
  qs2::qs_save(mget(c("coverageByTranscriptFST", "coverage_to_dt",
                     "coverageByTranscriptC", "coverageByTranscriptSum",
                     "pmapToTranscriptF"), ns), baseline_file)
}
old <- qs2::qs_read(baseline_file)
for (i in seq_along(old)) environment(old[[i]]) <- ns

# Batch the page/exon join, retaining input exon order and within-exon page order.
candidate_fst <- function(grl, fst_index, columns = NULL) {
  if (!file.exists(fst_index)) stop("No valid fst index file at location: ", fst_index)
  index <- read_fst(fst_index, as.data.table = TRUE)
  if (dirname(index$file_forward[1]) != dirname(fst_index)) {
    index[, file_forward := file.path(dirname(fst_index), basename(file_forward))]
    index[, file_reverse := file.path(dirname(fst_index), basename(file_reverse))]
  }
  if (!file.exists(index$file_forward[1])) stop("First FST coverage page missing")
  gr <- unlist(grl, use.names = FALSE)
  stopifnot(all(as.character(seqnames(gr)) %in% index$chr))
  query <- data.table(chr = as.character(seqnames(gr)), start = start(gr),
    end = end(gr), query_id = seq_along(gr),
    transcript_id = rep(seq_along(grl), lengths(grl)),
    strand = as.character(strand(gr)))
  if (!nrow(query)) return(lapply(seq_along(grl), function(i) data.table()))
  setkey(query, chr, start, end)
  tiles <- foverlaps(index, query, by.x = c("chr", "start", "end"), nomatch = 0L)
  tiles[, start_segment := pmax(start - i.start + 1, 1)]
  tiles[, end_segment := pmin(end - i.start + 1, i.end - i.start + 1)]
  tiles[, file := ifelse(strand == "+", file_forward, file_reverse)]
  setorder(tiles, transcript_id, query_id, i.start, i.end)
  groups <- split(seq_len(nrow(tiles)), factor(tiles$transcript_id,
                                            levels = seq_along(grl)))
  fst::threads_fst(5, reset_after_fork = FALSE)
  answer <- lapply(groups, function(rows) rbindlist(lapply(rows, function(i) {
    d <- read_fst(tiles$file[i], columns, from = tiles$start_segment[i],
                  to = tiles$end_segment[i], as.data.table = TRUE)
    if (tiles$strand[i] == "-") d <- d[rev(seq_len(nrow(d)))]
    d
  })))
  names(answer) <- names(grl)
  answer
}
environment(candidate_fst) <- ns

# Bound transient reads while sharing overlapping/nearby slices of each page.
# No retained cache: every call reads the current index and requested pages.
read_fst_tiles <- function(tiles, columns) {
  pieces <- vector("list", nrow(tiles))
  read_block <- function(rows) {
    first <- min(tiles$start_segment[rows])
    last <- max(tiles$end_segment[rows])
    d <- read_fst(tiles$file[rows[1L]], columns, from = first, to = last,
                  as.data.table = TRUE)
    for (i in rows) {
      offsets <- seq.int(tiles$start_segment[i] - first + 1L,
                         tiles$end_segment[i] - first + 1L)
      if (tiles$strand[i] == "-") offsets <- rev(offsets)
      pieces[[i]] <<- d[offsets]
    }
  }
  for (file in unique(tiles$file)) {
    rows <- which(tiles$file == file)
    rows <- rows[order(tiles$start_segment[rows], tiles$end_segment[rows])]
    block <- integer(); first <- last <- total <- 0
    for (i in rows) {
      s <- tiles$start_segment[i]; e <- tiles$end_segment[i]
      span <- max(last, e) - first + 1L
      if (length(block) && (span > 8192L || span > 8 * (total + e - s + 1L))) {
        read_block(block); block <- integer()
      }
      if (!length(block)) { first <- s; last <- e; total <- 0 }
      block <- c(block, i); last <- max(last, e); total <- total + e - s + 1L
    }
    if (length(block)) read_block(block)
  }
  pieces
}
candidate_env <- new.env(parent = ns)
candidate_env$read_fst_tiles <- read_fst_tiles
environment(candidate_env$read_fst_tiles) <- ns
txt <- paste(deparse(body(candidate_fst)), collapse = "\n")
begin <- regexpr("answer <- lapply", txt, fixed = TRUE)[1]
finish <- regexpr("names(answer)", txt, fixed = TRUE)[1]
stopifnot(begin > 0L, finish > begin)
txt <- paste0(substr(txt, 1L, begin - 1L),
  "pieces <- read_fst_tiles(tiles, columns)\n    answer <- lapply(groups, function(rows) rbindlist(pieces[rows]))\n    ",
  substr(txt, finish, nchar(txt)))
body(candidate_fst) <- parse(text = txt)[[1]]
environment(candidate_fst) <- candidate_env

# Replace tile()/unlist() with expansion of run widths, starts, and group names.
candidate_dt <- old$coverage_to_dt
txt <- paste(deparse(body(candidate_dt)), collapse = "\n")
txt <- sub("ir <- unlist(tile(ir, width = 1))",
  "run_widths <- width(ir)\n    positions <- rep.int(start(ir) - 1L, run_widths) + sequence(run_widths)\n    gene_names <- rep.int(names(ir), run_widths)", txt, fixed = TRUE)
txt <- gsub("factor(names(ir))", "factor(gene_names)", txt, fixed = TRUE)
txt <- gsub("as.integer(names(ir))", "as.integer(gene_names)", txt, fixed = TRUE)
txt <- gsub("position = start(ir)", "position = positions", txt, fixed = TRUE)
body(candidate_dt) <- parse(text = txt)[[1]]
environment(candidate_dt) <- ns

# ORFik already uses chromosome Views for >5,000 transcripts. Test the same
# extraction path on smaller batches instead of repeated GRanges dispatch.
candidate_cov <- old$coverageByTranscriptC
txt <- paste(deparse(body(candidate_cov)), collapse = "\n")
txt <- sub("if \\(length\\(transcripts\\) > [^)]+\\)", "if (TRUE)", txt)
body(candidate_cov) <- parse(text = txt)[[1]]
environment(candidate_cov) <- ns

rows <- list()
bench <- function(label, fun, reference, candidate, args, n = 9L, batch = 1L) {
  a <- do.call(reference, args); b <- do.call(candidate, args)
  stopifnot(isTRUE(all.equal(a, b, check.attributes = TRUE, tolerance = 0)))
  profile <- file.path(out, paste0("profile_", label, ".out"))
  Rprof(profile, interval = .01)
  for (j in seq_len(batch)) invisible(do.call(reference, args))
  Rprof(NULL)
  p <- summaryRprof(profile)
  fwrite(as.data.table(p$by.self, keep.rownames = "call"),
         paste0(profile, ".csv"))
  for (iteration in seq_len(n)) {
    order <- if (iteration %% 2L) c("baseline", "candidate") else c("candidate", "baseline")
    for (method in order) {
      f <- if (method == "baseline") reference else candidate
      gc()
      elapsed <- system.time(for (j in seq_len(batch)) invisible(do.call(f, args)))[["elapsed"]] / batch
      rows[[length(rows) + 1L]] <<- data.table(input = label, function_name = fun,
        method = method, iteration = iteration, batch = batch, seconds = elapsed)
    }
  }
  tab <- rbindlist(rows)
  fwrite(tab, file.path(out, "timings.csv"))
  s <- tab[, .(median_seconds = median(seconds), min_seconds = min(seconds),
               max_seconds = max(seconds)), by = .(input, function_name, method)]
  fwrite(s, file.path(out, "summary.csv"))
  print(s[input == label])
}

# Actual production display coordinates and all-library cached coverage.
if (Sys.getenv("ORFIK_BENCHMARK_DEFINITIONS_ONLY") != "1") {
real <- lapply(c("ATF4", "FMR1", "BOD1L1"), function(gene) {
  f <- list.files(".coverage_cache", pattern = paste0("^", gene, "__.*\\.qs$"), full.names = TRUE)
  stopifnot(length(f) >= 1L)
  x <- qs2::qs_read(f[1])
  s <- x$signature
  grl <- GRangesList(GRanges(s$seqnames, IRanges(s$start, s$end), s$strand))
  names(grl) <- s$transcript_names
  list(gene = gene, grl = grl, index = s$fst_index, coverage = x$coverage)
})
manifest <- rbindlist(lapply(real, function(x) data.table(gene = x$gene,
  transcript = names(x$grl), exons = lengths(x$grl),
  nucleotides = sum(width(unlist(x$grl))), libraries = ncol(x$coverage),
  index = x$index, index_mtime = as.character(file.info(x$index)$mtime))))
fwrite(manifest, file.path(out, "real_inputs.csv"))
for (x in real) {
  bench(paste0(x$gene, "_4_libraries"), "coverageByTranscriptFST",
    old$coverageByTranscriptFST, candidate_fst,
    list(grl = x$grl, fst_index = x$index, columns = head(names(x$coverage), 4L)), batch = 3L)
  bench(paste0(x$gene, "_all_libraries"), "coverageByTranscriptFST",
    old$coverageByTranscriptFST, candidate_fst,
    list(grl = x$grl, fst_index = x$index), n = 5L)
  # Real summed-library profile: retain fractional coverage and actual zero runs.
  cov <- RleList(rowSums(as.matrix(x$coverage)))
  names(cov) <- names(x$grl)
  bench(paste0(x$gene, "_summed_profile"), "coverage_to_dt",
    old$coverage_to_dt, candidate_dt, list(coverage = cov, drop.zero.dt = TRUE,
                                        withFrames = TRUE), batch = 20L)
}
# Synthetic stress inputs have deterministic seeds and varied run structure.
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
  bench(paste0("100_variable_windows_density_", density), "coverage_to_dt",
    old$coverage_to_dt, candidate_dt,
    list(coverage = RleList(lapply(seq_along(cov), function(i) cov[[i]][seq_len(5000L - i)])), drop.zero.dt = FALSE,
         withFrames = TRUE))
}
# Multi-transcript page-join workload from repeated real coordinates (no fake I/O).
x <- real[[1L]]
bench("ATF4_100_transcripts_4_libraries", "coverageByTranscriptFST",
  old$coverageByTranscriptFST, candidate_fst,
  list(grl = rep(x$grl, 100L), fst_index = x$index,
       columns = head(names(x$coverage), 4L)), n = 5L)
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
qs2::qs_save(list(fst = candidate_fst, dt = candidate_dt), file.path(out, "verified_candidates.qs"))
}
