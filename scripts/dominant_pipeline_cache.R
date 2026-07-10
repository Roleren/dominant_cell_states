pipeline_cache_dir <- function(analysis_dir) {
  cache_dir <- file.path(analysis_dir, ".pipeline_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

pipeline_cache_manifest_file <- function(step_name, analysis_dir) {
  file.path(pipeline_cache_dir(analysis_dir), paste0(step_name, ".rds"))
}

pipeline_normalize_paths <- function(paths, base_dir = NULL) {
  if (is.null(paths) || length(paths) == 0) return(character())
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (!length(paths)) return(character())
  if (!is.null(base_dir)) {
    is_abs <- grepl("^(/|~)", paths)
    paths[!is_abs] <- file.path(base_dir, paths[!is_abs])
  }
  unique(normalizePath(paths, mustWork = FALSE))
}

pipeline_md5_file <- function(path) {
  path <- pipeline_normalize_paths(path)
  if (length(path) != 1L || !file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

pipeline_file_signature <- function(paths, base_dir = NULL, include_md5 = FALSE) {
  paths <- pipeline_normalize_paths(paths, base_dir = base_dir)
  if (!length(paths)) {
    return(data.frame(
      path = character(),
      exists = logical(),
      is_dir = logical(),
      size = numeric(),
      mtime = numeric(),
      md5 = character(),
      stringsAsFactors = FALSE
    ))
  }

  info <- file.info(paths, extra_cols = FALSE)
  exists <- !is.na(info$isdir)
  out <- data.frame(
    path = paths,
    exists = exists,
    is_dir = ifelse(exists, info$isdir, NA),
    size = ifelse(exists, as.numeric(info$size), NA_real_),
    mtime = ifelse(exists, as.numeric(info$mtime), NA_real_),
    md5 = NA_character_,
    stringsAsFactors = FALSE
  )

  if (include_md5) {
    file_rows <- out$exists & !out$is_dir
    out$md5[file_rows] <- vapply(out$path[file_rows], pipeline_md5_file,
                                 character(1))
  }
  out[order(out$path), , drop = FALSE]
}

pipeline_step_fingerprint <- function(step_name, script_file, analysis_dir,
                                      dependencies = character(),
                                      cache_version = "dominant-pipeline-cache-v1",
                                      extra_source_files = character()) {
  helper_file <- file.path(analysis_dir, "scripts", "dominant_pipeline_cache.R")
  source_files <- unique(c(script_file, helper_file, extra_source_files))
  list(
    step = step_name,
    cache_version = cache_version,
    sources = pipeline_file_signature(source_files, include_md5 = TRUE),
    dependencies = pipeline_file_signature(dependencies, include_md5 = FALSE)
  )
}

pipeline_missing_outputs <- function(outputs) {
  sig <- pipeline_file_signature(outputs)
  sig$path[!sig$exists]
}

pipeline_cache_status <- function(step_name, script_file, analysis_dir,
                                  outputs,
                                  dependencies = character(),
                                  cache_version = "dominant-pipeline-cache-v1",
                                  extra_source_files = character()) {
  outputs <- pipeline_normalize_paths(outputs)
  missing_outputs <- pipeline_missing_outputs(outputs)
  if (!length(outputs)) {
    return(list(
      valid = FALSE,
      reason = "no cache outputs declared",
      fingerprint = NULL,
      manifest = NULL
    ))
  }
  if (length(missing_outputs)) {
    return(list(
      valid = FALSE,
      reason = paste0("missing output(s): ",
                      paste(basename(missing_outputs), collapse = ", ")),
      fingerprint = NULL,
      manifest = NULL
    ))
  }

  manifest_file <- pipeline_cache_manifest_file(step_name, analysis_dir)
  if (!file.exists(manifest_file)) {
    return(list(
      valid = FALSE,
      reason = "cache manifest missing",
      fingerprint = NULL,
      manifest = NULL
    ))
  }

  manifest <- tryCatch(readRDS(manifest_file), error = function(e) NULL)
  if (is.null(manifest) || is.null(manifest$fingerprint)) {
    return(list(
      valid = FALSE,
      reason = "cache manifest unreadable",
      fingerprint = NULL,
      manifest = manifest
    ))
  }

  fingerprint <- pipeline_step_fingerprint(
    step_name = step_name,
    script_file = script_file,
    analysis_dir = analysis_dir,
    dependencies = dependencies,
    cache_version = cache_version,
    extra_source_files = extra_source_files
  )

  if (!identical(fingerprint, manifest$fingerprint)) {
    return(list(
      valid = FALSE,
      reason = "inputs, source files, or cache version changed",
      fingerprint = fingerprint,
      manifest = manifest
    ))
  }

  list(
    valid = TRUE,
    reason = "cache valid",
    fingerprint = fingerprint,
    manifest = manifest
  )
}

pipeline_cache_write <- function(step_name, script_file, analysis_dir,
                                 outputs,
                                 dependencies = character(),
                                 cache_version = "dominant-pipeline-cache-v1",
                                 extra_source_files = character(),
                                 seconds = NA_real_) {
  manifest_file <- pipeline_cache_manifest_file(step_name, analysis_dir)
  manifest <- list(
    step = step_name,
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    seconds = seconds,
    outputs = pipeline_file_signature(outputs, include_md5 = FALSE),
    fingerprint = pipeline_step_fingerprint(
      step_name = step_name,
      script_file = script_file,
      analysis_dir = analysis_dir,
      dependencies = dependencies,
      cache_version = cache_version,
      extra_source_files = extra_source_files
    )
  )
  saveRDS(manifest, manifest_file)
  manifest_file
}
