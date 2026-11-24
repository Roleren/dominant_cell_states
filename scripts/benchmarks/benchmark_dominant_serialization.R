# Small real-transcript serialization benchmark; production caches remain untouched.
stopifnot(requireNamespace("qs2", quietly = TRUE))
output_dir <- file.path("results", "serialization_benchmark_2026-10-09")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
genes <- c("FMR1", "BOD1L1", "FN1")
threads <- c(1L, 2L, 4L, 5L, 8L)
repeats <- 3L
rows <- list()
for (gene in genes) {
  paths <- list.files(".coverage_cache", paste0("^", gene, "__.*[.](rds|qs)$"), full.names = TRUE)
  stopifnot(length(paths) > 0L)
  input <- paths[which.max(as.numeric(file.info(paths)$mtime))]
  object <- if (grepl("[.]rds$", input)) readRDS(input) else qs2::qs_read(input, nthreads=4L)
  configurations <- c("rds_gzip", paste0("qs2_", threads))
  for (repeat_id in seq_len(repeats)) {
    # Rotate order to reduce systematic first-format bias.
    order <- ((seq_along(configurations) + repeat_id - 2L) %% length(configurations)) + 1L
    for (index in order) {
      format <- configurations[index]
      nthreads <- if (format == "rds_gzip") NA_integer_ else threads[index - 1L]
      path <- file.path(output_dir, paste0(gene, "_", format, if (format == "rds_gzip") ".rds" else ".qs"))
      gc(verbose = FALSE)
      write_seconds <- system.time({
        if (format == "rds_gzip") saveRDS(object, path)
        else qs2::qs_save(object, path, nthreads = nthreads)
      })[["elapsed"]]
      gc(verbose = FALSE)
      read_seconds <- system.time({
        restored <- if (format == "rds_gzip") readRDS(path)
        else qs2::qs_read(path, nthreads = nthreads)
      })[["elapsed"]]
      # Content and attributes must survive each round trip.
      equal <- isTRUE(all.equal(object, restored, check.attributes = TRUE))
      stopifnot(equal)
      rows[[length(rows) + 1L]] <- data.frame(gene = gene, input = input,
        format = format, threads = nthreads, repeat_id = repeat_id,
        write_seconds = write_seconds, read_seconds = read_seconds,
        file_bytes = file.info(path)$size, object_bytes = as.numeric(object.size(object)),
        roundtrip_equal = equal)
    }
  }
  message("Finished ", gene)
}
measurements <- do.call(rbind, rows)
write.csv(measurements, file.path(output_dir, "measurements.csv"), row.names = FALSE)
medians <- aggregate(cbind(write_seconds, read_seconds, file_bytes, object_bytes) ~ gene + format,
                     measurements, median)
write.csv(medians, file.path(output_dir, "per_transcript_medians.csv"), row.names = FALSE)
summary <- aggregate(cbind(write_seconds, read_seconds, file_bytes) ~ format, medians, sum)
summary <- summary[match(c("rds_gzip", paste0("qs2_", threads)), summary$format), ]
summary$read_speedup <- summary$read_seconds[1] / summary$read_seconds
summary$write_speedup <- summary$write_seconds[1] / summary$write_seconds
summary$file_mib <- summary$file_bytes / 1024^2
write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
sink(file.path(output_dir, "session_info.txt")); print(sessionInfo()); print(setNames(lapply(c("compress_level", "shuffle", "nthreads"), qs2::qopt), c("compress_level", "shuffle", "nthreads"))); sink()
