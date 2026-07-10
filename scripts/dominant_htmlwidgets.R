dominant_find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "scripts", "run_dominant_cell_state_pipeline.R")) &&
        file.exists(file.path(here, "scripts", "dominant_cell_states.R"))) {
      return(here)
    }
    if (basename(here) == "scripts" &&
        file.exists(file.path(here, "run_dominant_cell_state_pipeline.R")) &&
        file.exists(file.path(here, "dominant_cell_states.R"))) {
      return(dirname(here))
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        file.exists(file.path(candidate, "scripts", "run_dominant_cell_state_pipeline.R"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

dominant_htmlwidget_libdir <- function(analysis_dir = dominant_find_analysis_dir()) {
  libdir <- file.path(analysis_dir, "htmlwidget_libs")
  dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(libdir, mustWork = TRUE)
}

dominant_relative_path <- function(from_dir, to_dir) {
  from_dir <- normalizePath(from_dir, mustWork = TRUE)
  to_dir <- normalizePath(to_dir, mustWork = TRUE)
  from_parts <- strsplit(from_dir, "/", fixed = TRUE)[[1]]
  to_parts <- strsplit(to_dir, "/", fixed = TRUE)[[1]]
  common <- 0L
  max_common <- min(length(from_parts), length(to_parts))
  while (common < max_common &&
         identical(from_parts[[common + 1L]], to_parts[[common + 1L]])) {
    common <- common + 1L
  }
  up <- rep("..", length(from_parts) - common)
  down <- if (common < length(to_parts)) {
    to_parts[(common + 1L):length(to_parts)]
  } else {
    character()
  }
  rel <- c(up, down)
  rel <- rel[nzchar(rel)]
  if (length(rel) == 0) "." else do.call(file.path, as.list(rel))
}

dominant_widget_dependency_signature <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, all.files = TRUE,
                        full.names = TRUE, no.. = TRUE)
    files <- files[file.exists(files) & !dir.exists(files)]
    if (!length(files)) return("EMPTY_DIR")
    rel <- substring(files, nchar(path) + 2L)
    ord <- order(rel)
    files <- files[ord]
    rel <- rel[ord]
    hashes <- unname(tools::md5sum(files))
    paste(paste(rel, hashes, sep = "="), collapse = "\n")
  } else {
    paste0(basename(path), "=", unname(tools::md5sum(path)))
  }
}

dominant_copy_widget_libs <- function(private_libdir, shared_libdir) {
  if (!dir.exists(private_libdir)) return(invisible(character()))
  deps <- list.files(private_libdir, all.files = FALSE, full.names = TRUE,
                     no.. = TRUE)
  copied <- character()
  for (dep in deps) {
    target <- file.path(shared_libdir, basename(dep))
    sig <- dominant_widget_dependency_signature(dep)
    sig_file <- paste0(target, ".dominant_widget_signature")
    if ((dir.exists(target) || file.exists(target))) {
      old_sig <- if (file.exists(sig_file)) {
        paste(readLines(sig_file, warn = FALSE), collapse = "\n")
      } else {
        tryCatch(dominant_widget_dependency_signature(target),
                 error = function(e) NA_character_)
      }
      if (identical(old_sig, sig)) {
        if (!file.exists(sig_file)) writeLines(sig, sig_file, useBytes = TRUE)
        copied <- c(copied, target)
        next
      }
    }
    if (dir.exists(target) || file.exists(target)) {
      unlink(target, recursive = TRUE, force = TRUE)
    }
    if (file.exists(sig_file)) unlink(sig_file, force = TRUE)
    ok <- file.copy(dep, shared_libdir, recursive = TRUE, copy.date = TRUE)
    if (!isTRUE(ok)) {
      stop("Could not copy htmlwidget dependency to shared libdir: ",
           dep, call. = FALSE)
    }
    writeLines(sig, sig_file, useBytes = TRUE)
    copied <- c(copied, target)
  }
  invisible(copied)
}

dominant_save_widget <- function(widget, output_file,
                                 title = class(widget)[[1]],
                                 background = "white") {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("Package 'htmlwidgets' is required to save interactive HTML output.",
         call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  local_libdir <- paste0(tools::file_path_sans_ext(basename(output_file)),
                         "_files")
  local_libpath <- file.path(dirname(output_file), local_libdir)
  htmlwidgets::saveWidget(
    widget,
    file = output_file,
    selfcontained = FALSE,
    libdir = local_libdir,
    background = background,
    title = title
  )
  shared_libdir <- dominant_htmlwidget_libdir()
  dominant_copy_widget_libs(local_libpath, shared_libdir)
  rel_libdir <- dominant_relative_path(dirname(output_file), shared_libdir)
  html <- readLines(output_file, warn = FALSE)
  html <- gsub(paste0(local_libdir, "/"), paste0(rel_libdir, "/"),
               html, fixed = TRUE)
  writeLines(html, output_file, useBytes = TRUE)
  unlink(local_libpath, recursive = TRUE, force = TRUE)
  invisible(output_file)
}

dominant_save_ggplotly <- function(plot, output_file,
                                   title = tools::file_path_sans_ext(
                                     basename(output_file)
                                   )) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required to save interactive plotly output.",
         call. = FALSE)
  }
  widget <- suppressWarnings(plotly::ggplotly(plot))
  dominant_save_widget(widget, output_file = output_file, title = title)
}

dominant_remove_private_widget_libs <- function(
    analysis_dir = dominant_find_analysis_dir()) {
  dirs <- list.dirs(analysis_dir, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[grepl("_files$", basename(dirs))]
  dirs <- dirs[basename(dirs) != "htmlwidget_libs"]
  if (length(dirs) > 0) unlink(dirs, recursive = TRUE, force = TRUE)
  invisible(dirs)
}
