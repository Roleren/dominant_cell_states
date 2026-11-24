Sys.setenv(ORFIK_MAPPING_DEFINITIONS_ONLY = "1")
source("scripts/benchmark_orfik_accessors_mapping.R")
out <- file.path(out, "wrapper")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Replace expressions structurally; deparse line wrapping must not determine
# whether a candidate transformation is actually applied.
replace_expr <- function(expr, from, to) {
  if (identical(expr, from)) return(to)
  if (is.call(expr)) for (i in seq_along(expr)) expr[i] <- list(replace_expr(expr[[i]], from, to))
  expr
}
slot_candidate <- old$pmapToTranscriptF
body(slot_candidate) <- replace_expr(body(slot_candidate),
  quote(as.character(unlist(strand(xOriginal), use.names = FALSE))),
  quote(as.character(xOriginal@unlistData@strand)))
body(slot_candidate) <- replace_expr(body(slot_candidate),
  quote(width(xOriginal@partitioning)),
  quote(xOriginal@partitioning@end - c(0L, head(xOriginal@partitioning@end, -1L))))
body(slot_candidate) <- replace_expr(body(slot_candidate),
  quote(width(transcripts@partitioning)),
  quote(transcripts@partitioning@end - c(0L, head(transcripts@partitioning@end, -1L))))
stopifnot(!identical(body(slot_candidate), body(old$pmapToTranscriptF)))
environment(slot_candidate) <- e
for (n in c(1L, 100L, 2000L)) {
  tx <- rep(suppressWarnings(do.call(c, real)), n)
  measure(paste0("to_slots_real_exons_", length(tx)), old$pmapToTranscriptF,
          slot_candidate, list(x = tx, transcripts = tx), batch = if (n == 1L) 5L else 1L)
}
# Isolate the most expensive accessor replacement, including its object creation.
tx <- rep(real[[2]], 1000L)
measure("strand_unlist_1000x17exons", function(x) as.character(unlist(strand(x), use.names = FALSE)),
        function(x) as.character(x@unlistData@strand), list(x = tx), batch = 100L)
qs2::qs_save(slot_candidate, file.path(out, "verified_wrapper.qs"))
