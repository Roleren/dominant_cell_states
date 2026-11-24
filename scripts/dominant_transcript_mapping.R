# Guard the fast native mapper against ranges outside the paired transcript.
# ORFik's native exon loop can read beyond its array for unmappable ranges.
dominant_checked_pmap <- function(query, subject) {
  if (!length(query)) return(query)
  if (length(query) != length(subject)) {
    if (length(query) == 1L) query <- rep(query, length(subject))
    else if (length(subject) == 1L) subject <- rep(subject, length(query))
    else stop("Query and transcript lengths must match or recycle one element.")
  }
  mapped <- GenomicRanges::GRangesList(rep(list(GenomicRanges::GRanges()), length(query)))
  names(mapped) <- names(query)
  if (!length(subject)) return(mapped)
  q <- unlist(query, use.names = FALSE)
  s <- unlist(subject, use.names = FALSE)
  query_groups <- rep(seq_along(query), lengths(query))
  subject_groups <- rep(seq_along(subject), lengths(subject))
  covered <- rep(FALSE, length(q))
  hits <- GenomicRanges::findOverlaps(q, s, type = "within", ignore.strand = FALSE)
  paired <- query_groups[S4Vectors::queryHits(hits)] ==
    subject_groups[S4Vectors::subjectHits(hits)]
  covered[S4Vectors::queryHits(hits)[paired]] <- TRUE
  covered <- covered & width(q) > 0L
  valid <- lengths(query) > 0L
  valid[unique(query_groups[!covered])] <- FALSE
  if (any(valid)) {
    mapped[valid] <- suppressWarnings(ORFik::pmapToTranscriptF(query[valid], subject[valid]))
    widths <- vapply(subject, function(gr) sum(width(gr)), numeric(1))
    bounded <- vapply(seq_along(mapped), function(i) {
      gr <- mapped[[i]]
      !length(gr) || all(start(gr) >= 1L & end(gr) <= widths[i] & width(gr) > 0L)
    }, logical(1))
    mapped[!bounded] <- GenomicRanges::GRangesList(
      rep(list(GenomicRanges::GRanges()), sum(!bounded)))
  }
  mapped
}
