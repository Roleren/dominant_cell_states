dominant_coverage_cache_dir <- function(analysis_dir) {
  cache_dir <- file.path(analysis_dir, ".coverage_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

dominant_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

dominant_hash_object <- function(x) {
  tmp <- tempfile("dominant-hash-")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  unname(tools::md5sum(tmp))
}

dominant_display_signature <- function(display_region, fst_index) {
  gr <- unlist(display_region, use.names = FALSE)
  info <- file.info(fst_index, extra_cols = FALSE)
  list(
    transcript_names = names(display_region),
    seqnames = as.character(seqnames(gr)),
    start = as.integer(start(gr)),
    end = as.integer(end(gr)),
    strand = as.character(strand(gr)),
    fst_index = normalizePath(fst_index, mustWork = FALSE),
    fst_index_exists = file.exists(fst_index),
    fst_index_size = if (!is.na(info$size)) as.numeric(info$size) else NA_real_,
    fst_index_mtime = if (!is.na(info$mtime)) as.numeric(info$mtime) else NA_real_
  )
}

dominant_coverage_cache_file <- function(analysis_dir, gene_symbol, tx_id,
                                         signature) {
  prefix <- paste(
    dominant_safe_name(gene_symbol),
    dominant_safe_name(tx_id),
    sep = "__"
  )
  file.path(
    dominant_coverage_cache_dir(analysis_dir),
    paste0(prefix, "__", dominant_hash_object(signature), ".rds")
  )
}

dominant_cached_coverage_by_transcript <- function(display_region, fst_index,
                                                   gene_symbol, tx_id,
                                                   analysis_dir,
                                                   reader = coverageByTranscriptFST,
                                                   force = FALSE) {
  signature <- dominant_display_signature(display_region, fst_index)
  cache_file <- dominant_coverage_cache_file(
    analysis_dir = analysis_dir,
    gene_symbol = gene_symbol,
    tx_id = tx_id,
    signature = signature
  )

  if (!force && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$signature, signature)) {
      return(cached$coverage)
    }
  }

  coverage <- reader(display_region, fst_index, columns = NULL)[[1]]
  saveRDS(
    list(
      generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      gene_symbol = gene_symbol,
      tx_id = tx_id,
      signature = signature,
      coverage = coverage
    ),
    cache_file
  )
  coverage
}

