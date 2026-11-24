# Project-owned cached objects use qs2 (the .qs extension does not imply qs v1).
dominant_serialization_threads <- function() {
  value <- suppressWarnings(as.numeric(Sys.getenv("DOMINANT_QS_THREADS", "4")))
  if (length(value) != 1L || !is.finite(value) || value < 1 ||
      value != floor(value) || value > .Machine$integer.max) {
    stop("DOMINANT_QS_THREADS must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

dominant_read_object <- function(path) {
  qs2::qs_read(path, nthreads = dominant_serialization_threads())
}

dominant_save_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("qs2-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  qs2::qs_save(object, tmp, nthreads = dominant_serialization_threads())
  if (!file.rename(tmp, path)) stop("Could not save cached object: ", path)
  invisible(path)
}

# External ORFik inputs may still be RDS. Prefer a QS sibling only when current.
dominant_preferred_object_file <- function(path) {
  if (!grepl("[.]rds$", path)) return(path)
  qs <- sub("[.]rds$", ".qs", path)
  if (file.exists(qs) && (!file.exists(path) ||
      file.info(qs)$mtime >= file.info(path)$mtime)) qs else path
}

dominant_read_external_object <- function(path) {
  ORFik::read_RDSQS(dominant_preferred_object_file(path),
                  nthread = dominant_serialization_threads())
}
