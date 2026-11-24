serialization_helper <- c(
  file.path("scripts", "dominant_serialization.R"),
  file.path("dominant_cell_states", "scripts", "dominant_serialization.R"),
  file.path("..", "..", "scripts", "dominant_serialization.R")
)
serialization_helper <- serialization_helper[file.exists(serialization_helper)][1]
if (is.na(serialization_helper)) stop("Missing dominant_serialization.R")
source(serialization_helper, local = TRUE)

# Small derived measurements, keyed by biological inputs rather than gene-panel size.
dominant_input_files <- function(paths) {
  paths <- sort(unique(normalizePath(paths, mustWork = FALSE)))
  info <- file.info(paths)
  data.frame(path = paths, size = info$size, mtime = as.numeric(info$mtime))
}

dominant_function_signature <- function(fun) {
  list(formals = formals(fun), body = body(fun))
}

dominant_cached_object <- function(cache_dir, namespace, signature, compute,
                                   cacheable = function(value) TRUE,
                                   use_cache = Sys.getenv("DOMINANT_MEASUREMENT_CACHE", "on") != "off") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key <- digest::digest(signature, algo = "md5")
  path <- file.path(cache_dir, paste0(namespace, "__", key, ".qs"))
  cached <- if (use_cache && file.exists(path)) tryCatch(dominant_read_object(path), error = function(e) NULL)
  # Compare serialized values: data.table self-reference pointers are transient
  # and differ after deserialization, even when the biological inputs are identical.
  if (!is.null(cached) && identical(serialize(cached$signature, NULL),
                                   serialize(signature, NULL))) {
    message("Measurement cache hit: ", namespace)
    return(cached$value)
  }
  message("Measurement cache miss: ", namespace)
  value <- compute()
  if (use_cache && cacheable(value)) {
    dominant_save_object(list(signature = signature, value = value), path)
  }
  value
}

# Track the calculation and its project-local helper functions, rather than
# unrelated plotting or model code in the same producer script.
dominant_computation_signature <- function(fun) {
  signatures <- list()
  visit <- function(name, current) {
    if (name %in% names(signatures)) return(invisible(NULL))
    signatures[[name]] <<- dominant_function_signature(current)
    scope <- environment(current)
    if (is.null(scope) || isNamespace(scope)) return(invisible(NULL))
    for (called in setdiff(codetools::findGlobals(current, merge = FALSE)$functions,
      c("file.path", "dominant_cached_coverage_by_transcript",
        "dominant_read_object", "dominant_save_object", "dominant_serialization_threads"))) {
      helper <- if (exists(called, scope, inherits = FALSE)) {
        get(called, scope, inherits = FALSE)
      } else if (exists(called, globalenv(), inherits = FALSE)) {
        get(called, globalenv(), inherits = FALSE)
      } else NULL
      if (is.function(helper) && !is.primitive(helper)) visit(called, helper)
    }
  }
  visit("calculation", fun)
  signatures[order(names(signatures))]
}

# Explicit user-confirmed compatibility for unrelated development fixes on
# 2026-10-09. Other versions retain normal invalidation behavior.
dominant_orfik_cache_version <- function(version = as.character(utils::packageVersion("ORFik"))) {
  if (length(version) != 1L || is.na(version)) stop("Specify one ORFik version.")
  if (version %in% c("1.29.12", "1.29.14", "1.29.15")) "1.29.14" else version
}
