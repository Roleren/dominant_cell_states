serialization_helper <- c(
  file.path("scripts", "dominant_serialization.R"),
  file.path("dominant_cell_states", "scripts", "dominant_serialization.R"),
  file.path("..", "..", "scripts", "dominant_serialization.R")
)
serialization_helper <- serialization_helper[file.exists(serialization_helper)][1]
if (is.na(serialization_helper)) stop("Missing dominant_serialization.R")
source(serialization_helper, local = TRUE)

pipeline_cache_dir <- function(analysis_dir) {
  cache_dir <- file.path(analysis_dir, ".pipeline_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

pipeline_cache_manifest_file <- function(step_name, analysis_dir) {
  file.path(pipeline_cache_dir(analysis_dir), paste0(step_name, ".qs"))
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

.pipeline_digest_memo <- new.env(parent = emptyenv())

pipeline_md5_file <- function(path) {
  path <- pipeline_normalize_paths(path)
  if (length(path) != 1L || !file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  info <- file.info(path)
  key <- paste(info$size, as.numeric(info$mtime), as.numeric(info$ctime))
  previous <- .pipeline_digest_memo[[path]]
  if (!is.null(previous) && identical(previous$key, key)) return(previous$md5)
  value <- unname(tools::md5sum(path))
  .pipeline_digest_memo[[path]] <- list(key = key, md5 = value)
  value
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

  if (isTRUE(include_md5) || identical(include_md5, "auto")) {
    file_rows <- out$exists & !out$is_dir
    if (identical(include_md5, "auto")) {
      file_rows <- file_rows & grepl("\\.(R|r|csv|txt|md|rds|fst|qs)$", out$path) &
        out$size <= 256 * 1024^2 &
        !grepl("/Bio_data/|countTable_|totalCounts_|coverage_index[.]fst$", out$path)
    }
    out$md5[file_rows] <- vapply(out$path[file_rows], pipeline_md5_file,
                                 character(1))
  }
  # Rewriting identical content must not invalidate downstream computations.
  out$mtime[!is.na(out$md5)] <- NA_real_
  out <- out[order(out$path), , drop = FALSE]
  rownames(out) <- NULL
  out
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
    dependencies = pipeline_file_signature(dependencies, include_md5 = "auto")
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

  manifest <- tryCatch(dominant_read_object(manifest_file), error = function(e) NULL)
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
      reason = pipeline_fingerprint_changes(manifest$fingerprint, fingerprint),
      fingerprint = fingerprint,
      manifest = manifest
    ))
  }

  current_outputs <- pipeline_file_signature(outputs, include_md5 = "auto")
  recorded_outputs <- manifest$outputs[match(current_outputs$path, manifest$outputs$path), , drop = FALSE]
  rownames(recorded_outputs) <- NULL
  if (!identical(current_outputs, recorded_outputs)) {
    return(list(valid = FALSE, reason = "output files changed since completion",
                fingerprint = fingerprint, manifest = manifest))
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
    outputs = pipeline_file_signature(outputs, include_md5 = "auto"),
    fingerprint = pipeline_step_fingerprint(
      step_name = step_name,
      script_file = script_file,
      analysis_dir = analysis_dir,
      dependencies = dependencies,
      cache_version = cache_version,
      extra_source_files = extra_source_files
    )
  )
  dominant_save_object(manifest, manifest_file)
  manifest_file
}

# Identify the inputs responsible for a miss, including added/removed inputs.
pipeline_fingerprint_changes <- function(previous, current) {
  changed <- character()
  if (!identical(previous$cache_version, current$cache_version)) {
    changed <- "cache version"
  }
  for (kind in c("sources", "dependencies")) {
    old <- previous[[kind]]
    new <- current[[kind]]
    for (path in union(old$path, new$path)) {
      a <- old[old$path == path, , drop = FALSE]
      b <- new[new$path == path, , drop = FALSE]
      rownames(a) <- rownames(b) <- NULL
      if (!identical(a, b)) changed <- c(changed, paste0(kind, ": ", basename(path)))
    }
  }
  paste0("changed: ", paste(unique(changed), collapse = ", "))
}

# Follow actual source() calls, including filenames assigned to helper variables.
# File-existence probes and filenames in comments are not source dependencies.
.pipeline_source_memo <- new.env(parent = emptyenv())
pipeline_direct_source_references <- function(path) {
  info <- file.info(path)
  key <- paste(info$size, as.numeric(info$mtime), as.numeric(info$ctime))
  previous <- .pipeline_source_memo[[path]]
  if (!is.null(previous) && identical(previous$key, key)) return(previous$references)
  assignments <- list()
  targets <- list()
  strings <- function(expr) {
    if (missing(expr)) return(character())
    if (is.character(expr)) return(expr)
    if (!is.call(expr) && !is.pairlist(expr) && !is.expression(expr)) return(character())
    unlist(lapply(as.list(expr), strings), use.names = FALSE)
  }
  walk <- function(expr) {
    if (missing(expr)) return(invisible(NULL))
    if (!is.call(expr) && !is.expression(expr)) return(invisible(NULL))
    if (is.call(expr)) {
      head <- as.character(expr[[1]])[1]
      if (head %in% c("<-", "=")) {
        if (is.symbol(expr[[2]]) && !(is.call(expr[[3]]) &&
            identical(expr[[3]][[1]], as.name("function")))) {
          name <- as.character(expr[[2]])
          assignments[[name]] <<- c(assignments[[name]], list(expr[[3]]))
        }
      }
      if (head %in% c("source", "sys.source") && length(expr) >= 2L) {
        targets[[length(targets) + 1L]] <<- expr[[2]]
      }
    }
    for (child in as.list(expr)) walk(child)
  }
  walk(parse(path))
  resolve <- function(expr, seen = character()) {
    if (missing(expr)) return(character())
    literals <- strings(expr)
    paths <- literals[grepl("[.]R$", literals)]
    variables <- setdiff(all.names(expr), seen)
    for (variable in intersect(variables, names(assignments))) {
      for (value in assignments[[variable]]) {
        paths <- c(paths, resolve(value, c(seen, variable)))
      }
    }
    paths
  }
  references <- unique(basename(as.character(unlist(lapply(targets, resolve)))))
  .pipeline_source_memo[[path]] <- list(key = key, references = references)
  references
}

pipeline_step_source_files <- function(script_file, analysis_dir) {
  visited <- character()
  visit <- function(path) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!file.exists(path) || path %in% visited) return(invisible(NULL))
    visited <<- c(visited, path)
    for (reference in pipeline_direct_source_references(path)) {
      # Step scripts source helpers from scripts/ and package code from R/.
      candidates <- file.path(analysis_dir, c("scripts", "R"), reference)
      existing <- candidates[file.exists(candidates)]
      visit(if (length(existing)) existing[[1]] else candidates[[1]])
    }
  }
  visit(script_file)
  # The runner injects this helper into every step environment.
  visit(file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R"))
  setdiff(visited, normalizePath(script_file, mustWork = FALSE))
}

pipeline_step_cache_outputs <- function(outputs, rdg_plot_mode = "all") {
  if (identical(rdg_plot_mode, "tables")) {
    outputs <- outputs[!grepl("\\.(png|pdf|svg)$|(^|/)figures/.*\\.html$", outputs)]
  }
  outputs
}
