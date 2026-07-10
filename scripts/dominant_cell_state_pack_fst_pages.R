suppressPackageStartupMessages({
  library(data.table)
})

# Run this on RStudio Server from the folder containing required_fst_pages.csv.
#
# Typical use:
#   source("scripts/dominant_cell_state_pack_fst_pages.R")
#   pack_required_fst_pages(
#     page_source_dir = "/path/to/all_samples-Homo_sapiens/collection_tables_indexed"
#   )
#
# If the paths in fst_file_from_index are valid on the server, page_source_dir
# can be omitted.

pack_required_fst_pages <- function(
    required_csv = "required_fst_pages.csv",
    page_source_dir = NULL,
    output_dir = "required_fst_pages_export",
    archive_prefix = "required_fst_pages",
    link_mode = c("hardlink", "copy"),
    chunk_gb = Inf,
    compression = c("none", "gzip")) {
  link_mode <- match.arg(link_mode)
  compression <- match.arg(compression)

  if (!file.exists(required_csv)) {
    stop("Missing required page table: ", required_csv)
  }

  pages <- fread(required_csv)
  needed_cols <- c("fst_basename", "fst_file_from_index")
  if (!all(needed_cols %chin% names(pages))) {
    stop("Required page table is missing columns: ",
         paste(setdiff(needed_cols, names(pages)), collapse = ", "))
  }

  pages <- unique(pages[, .(fst_basename, fst_file_from_index)])
  if (!is.null(page_source_dir)) {
    pages[, source_path := file.path(path.expand(page_source_dir), fst_basename)]
  } else {
    pages[, source_path := path.expand(fst_file_from_index)]
  }

  pages[, exists := file.exists(source_path)]
  if (any(!pages$exists)) {
    missing <- pages[exists == FALSE]
    fwrite(missing, "missing_required_fst_pages.csv")
    stop(
      "Missing ", nrow(missing), " required FST pages. ",
      "I wrote missing_required_fst_pages.csv. ",
      "If the files are in another directory, rerun with page_source_dir = '/that/folder'."
    )
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  stage_dir <- file.path(output_dir, "pages")
  if (!dir.exists(stage_dir)) dir.create(stage_dir, recursive = TRUE)

  pages[, dest_path := file.path(stage_dir, fst_basename)]
  existing_dest <- file.exists(pages$dest_path)
  if (any(existing_dest)) {
    file.remove(pages$dest_path[existing_dest])
  }

  if (link_mode == "hardlink") {
    linked <- file.link(pages$source_path, pages$dest_path)
    if (any(!linked)) {
      message("Hardlink failed for ", sum(!linked), " files; copying those files instead.")
      copied <- file.copy(
        pages$source_path[!linked],
        pages$dest_path[!linked],
        overwrite = TRUE,
        copy.mode = FALSE,
        copy.date = FALSE
      )
      if (any(!copied)) {
        stop("Could not stage ", sum(!copied), " files by hardlink or copy.")
      }
    }
  } else {
    copied <- file.copy(
      pages$source_path,
      pages$dest_path,
      overwrite = TRUE,
      copy.mode = FALSE,
      copy.date = FALSE
    )
    if (any(!copied)) {
      stop("Could not copy ", sum(!copied), " required FST pages.")
    }
  }

  pages[, size_bytes := file.info(dest_path)$size]
  fwrite(pages, file.path(output_dir, "manifest.csv"))

  files <- pages[order(fst_basename), .(fst_basename, size_bytes)]
  if (is.finite(chunk_gb)) {
    chunk_bytes <- chunk_gb * 1024^3
    chunk_ids <- integer(nrow(files))
    current_chunk <- 1L
    current_size <- 0
    for (i in seq_len(nrow(files))) {
      if (current_size > 0 && current_size + files$size_bytes[i] > chunk_bytes) {
        current_chunk <- current_chunk + 1L
        current_size <- 0
      }
      chunk_ids[i] <- current_chunk
      current_size <- current_size + files$size_bytes[i]
    }
    files[, chunk := chunk_ids]
  } else {
    files[, chunk := 1L]
  }

  archive_ext <- if (compression == "gzip") ".tar.gz" else ".tar"
  archive_paths <- character()
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(stage_dir)

  for (chunk_id in sort(unique(files$chunk))) {
    chunk_files <- files[chunk == chunk_id, fst_basename]
    suffix <- if (max(files$chunk) == 1L) "" else sprintf("_part%02d", chunk_id)
    archive_path <- file.path(
      old_wd,
      output_dir,
      paste0(archive_prefix, suffix, archive_ext)
    )
    utils::tar(
      tarfile = archive_path,
      files = chunk_files,
      compression = compression
    )
    archive_paths <- c(archive_paths, archive_path)
  }

  setwd(old_wd)
  archives <- data.table(
    archive = archive_paths,
    size_bytes = file.info(archive_paths)$size
  )
  fwrite(archives, file.path(output_dir, "archives.csv"))

  message("Staged ", nrow(pages), " required FST pages in: ", normalizePath(stage_dir))
  message("Wrote manifest: ", normalizePath(file.path(output_dir, "manifest.csv")))
  message("Wrote ", nrow(archives), " archive(s):")
  for (archive in archives$archive) {
    message("  ", normalizePath(archive))
  }

  invisible(list(pages = pages, archives = archives))
}

message("Loaded pack_required_fst_pages(). Run it with page_source_dir if the indexed paths are not valid on this server.")
