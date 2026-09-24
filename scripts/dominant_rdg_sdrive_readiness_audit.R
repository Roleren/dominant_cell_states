#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) source(dominant_loader)

find_repo_root <- function(start = getwd()) {
  if (exists("dominant_find_ribocrypt_repo")) {
    repo <- dominant_find_ribocrypt_repo(start)
    if (nzchar(repo)) return(repo)
  }
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        identical(dominant_description_package(here), "RiboCrypt")) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find RiboCrypt repo root from: ", start, call. = FALSE)
}

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

setup_runtime_env <- function(analysis_dir) {
  # Capture ORFik's already-correct config (resolved via the default,
  # un-redirected BiocFileCache) BEFORE isolating XDG/BFC state below.
  # Once isolated, ORFik::config() can no longer see that cached entry and
  # regenerates one from its own hardcoded, environment-naive default
  # ("~/Bio_data/..."), which is wrong on any server where Bio_data lives
  # elsewhere (e.g. "~/livemount/Bio_data" here). Verified directly: this
  # reproduces with nothing more than pointing XDG_CACHE_HOME/
  # XDG_CONFIG_HOME at a fresh empty directory.
  orfik_config <- if (requireNamespace("ORFik", quietly = TRUE)) {
    tryCatch(ORFik::config(), error = function(e) NULL)
  } else NULL

  runtime_root <- file.path(analysis_dir, ".runtime")
  cache_dir <- file.path(runtime_root, "xdg-cache")
  config_dir <- file.path(runtime_root, "xdg-config")
  bfc_dir <- file.path(runtime_root, "biocfilecache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(bfc_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(
    XDG_CACHE_HOME = cache_dir,
    XDG_CONFIG_HOME = config_dir,
    BFC_CACHE = bfc_dir
  )

  # Re-seed the now-isolated BiocFileCache with the real config captured
  # above, so ORFik::config() inside this isolated environment still
  # resolves to the correct project paths instead of a wrong default.
  if (!is.null(orfik_config) && length(orfik_config) == 4) {
    conf_dt <- data.frame(
      type = names(orfik_config),
      directory = unname(orfik_config),
      stringsAsFactors = FALSE
    )
    tryCatch(ORFik:::config.save(conf = conf_dt), error = function(e) NULL)
  }
}

repo_root <- find_repo_root()
analysis_dir <- find_analysis_dir()
setup_runtime_env(analysis_dir)

if (!"RiboCrypt" %in% loadedNamespaces()) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package 'devtools' is required to load local RiboCrypt.",
         call. = FALSE)
  }
  devtools::load_all(repo_root, quiet = TRUE)
}

safe_fread <- function(path) {
  if (file.exists(path)) {
    return(fread(path, showProgress = FALSE))
  }
  data.table()
}

scalar_text <- function(x, default = "") {
  if (!length(x) || is.na(x[[1]])) return(default)
  trimws(as.character(x[[1]]))
}

as_bool_vec <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  value <- tolower(trimws(as.character(x)))
  value %chin% c("1", "true", "t", "yes", "y")
}

max_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

sum_or_zero <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(0)
  sum(x)
}

file_size <- function(path) {
  if (!length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size[[1]])
}

file_mtime <- function(path) {
  if (!length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return("")
  }
  format(file.info(path)$mtime[[1]], "%Y-%m-%d %H:%M:%S %Z")
}

existing_path_count <- function(paths) {
  paths <- as.character(paths)
  paths[is.na(paths)] <- ""
  sum(file.exists(paths[nzchar(paths)]))
}

normalize_existing_path <- function(path) {
  path <- scalar_text(path)
  if (!nzchar(path)) return("")
  normalizePath(path, mustWork = FALSE)
}

extract_fafile_member <- function(fa, member) {
  if (is.null(fa) || inherits(fa, "error")) return("")
  out <- tryCatch({
    if (exists("path", mode = "function") && identical(member, "path")) {
      return(scalar_text(get("path", mode = "function")(fa)))
    }
    xdata <- slot(fa, ".xData")
    if (exists(member, envir = xdata, inherits = FALSE)) {
      return(scalar_text(get(member, envir = xdata, inherits = FALSE)))
    }
    ""
  }, error = function(e) "")
  scalar_text(out)
}

output_dir <- Sys.getenv(
  "RDG_SDRIVE_READINESS_OUTPUT_DIR",
  unset = file.path(analysis_dir, "dominant_rdg_sdrive_readiness")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sdrive_root <- Sys.getenv("RDG_SDRIVE_ROOT", unset = "/media/roler/S")
experiment_name <- Sys.getenv("RDG_SDRIVE_EXPERIMENT",
                              unset = "all_merged-Homo_sapiens")
# Two layouts are in play: the S-drive convention nests an extra "Ribo-seq/"
# path segment under processed_data/, while the full server (this host)
# keeps "all_samples-Homo_sapiens/collection_tables_indexed/" directly under
# processed_data/ with no such segment. Try the server-relative path (derived
# from ORFik's own processed_data root, not hardcoded) first, then fall back
# to the legacy S-drive path.
all_samples_collection_dir_candidates <- c(
  Sys.getenv("RDG_ALL_SAMPLES_COLLECTION_DIR", unset = NA_character_),
  {
    # "all_samples-Homo_sapiens" is a union experiment spanning many
    # per-study libFolders, so ORFik::libFolder() on it resolves to
    # whichever individual study's row comes first, not a shared
    # "all_samples" directory. The processed_data root from ORFik::config()
    # is the reliable way to locate the real collection_tables_indexed dir.
    server_processed_root <- tryCatch(
      path.expand(unname(ORFik::config()["bam"])),
      error = function(e) NA_character_
    )
    if (!is.na(server_processed_root)) {
      file.path(server_processed_root, "all_samples-Homo_sapiens",
               "collection_tables_indexed")
    } else NA_character_
  },
  file.path(sdrive_root,
           "data/Bio_data/processed_data/Ribo-seq/all_samples-Homo_sapiens",
           "collection_tables_indexed")
)
all_samples_collection_dir_candidates <-
  all_samples_collection_dir_candidates[!is.na(all_samples_collection_dir_candidates)]
existing_collection_dirs <-
  all_samples_collection_dir_candidates[dir.exists(all_samples_collection_dir_candidates)]
all_samples_collection_dir <- if (length(existing_collection_dirs)) {
  existing_collection_dirs[[1]]
} else {
  all_samples_collection_dir_candidates[[1]]
}
fst_index_file <- file.path(all_samples_collection_dir, "coverage_index.fst")

status_rows <- list()
add_status <- function(check, ok, path = "", detail = "", value = "",
                       status = NULL) {
  path <- scalar_text(path)
  status_rows[[length(status_rows) + 1L]] <<- data.table(
    check = check,
    status = status %||% if (isTRUE(ok)) "ok" else "missing",
    ok = isTRUE(ok),
    value = scalar_text(value),
    path = path,
    detail = scalar_text(detail),
    size_bytes = file_size(path),
    mtime = file_mtime(path)
  )
}

mount_ok <- dir.exists(sdrive_root)
add_status("sdrive_mount", mount_ok, sdrive_root,
           "Configured S-drive root is reachable.")

experiment <- NULL
experiment_error <- ""
# Do not gate this on mount_ok: ORFik::read.experiment() resolves its own
# paths (ORFik::config()-based on this server) and does not require the
# configured S-drive root to be mounted. Gating here produced a false
# "error" status (rather than attempting the load) whenever sdrive_root
# genuinely does not apply to the current machine.
experiment <- tryCatch(
  ORFik::read.experiment(experiment_name, validate = FALSE),
  error = function(e) {
    experiment_error <<- conditionMessage(e)
    NULL
  }
)
experiment_ok <- !is.null(experiment)
add_status(
  "orfik_experiment_load",
  experiment_ok,
  "",
  experiment_error,
  experiment_name,
  status = if (experiment_ok) "ok" else "error"
)

lib_folder <- ""
ref_folder <- ""
fasta_path <- ""
fasta_index_path <- ""
fasta_gzindex_path <- ""
txdb_path <- ""
txdb_error <- ""
ofst_path <- ""
if (experiment_ok) {
  lib_folder <- tryCatch(scalar_text(ORFik::libFolder(experiment)),
                         error = function(e) "")
  ref_folder <- tryCatch(scalar_text(ORFik::refFolder(experiment)),
                         error = function(e) "")
  fa <- tryCatch(ORFik::findFa(experiment), error = function(e) e)
  fasta_path <- extract_fafile_member(fa, "path")
  fasta_index_path <- extract_fafile_member(fa, "index")
  fasta_gzindex_path <- extract_fafile_member(fa, "gzindex")
  # ORFik::txdbFile() does not exist in the installed ORFik version on this
  # server (neither exported nor internal); the experiment object itself
  # carries the resolved txdb path directly on its @txdb slot.
  txdb_path <- tryCatch(scalar_text(experiment@txdb),
                        error = function(e) {
                          txdb_error <<- conditionMessage(e)
                          ""
                        })
  if ("filepath" %in% names(experiment)) {
    ofst_candidates <- as.character(experiment$filepath)
    ofst_candidates <- ofst_candidates[grepl("\\.ofst$", ofst_candidates)]
    if (length(ofst_candidates)) ofst_path <- ofst_candidates[1]
  }
}
if (!nzchar(ofst_path) && nzchar(lib_folder)) {
  found_ofst <- list.files(lib_folder, pattern = "\\.ofst$",
                           full.names = TRUE, recursive = FALSE)
  if (length(found_ofst)) ofst_path <- found_ofst[1]
}

add_status("orfik_lib_folder", dir.exists(lib_folder), lib_folder,
           "All-merged ORFik library folder.")
add_status("orfik_ref_folder", dir.exists(ref_folder), ref_folder,
           "Reference folder returned by ORFik.")
add_status("reference_fasta", file.exists(fasta_path), fasta_path,
           "Primary genome FASTA returned by ORFik::findFa().")
add_status("reference_fasta_index", file.exists(fasta_index_path),
           fasta_index_path, "FASTA .fai index.")
# .gzi is only meaningful for a bgzip-compressed FASTA (Rsamtools random
# access); a plain, uncompressed .fa has no .gzi and needs none, so only
# flag this as missing when the FASTA itself is actually gzipped.
fasta_is_gzipped <- grepl("\\.(gz|bgz)$", fasta_path)
add_status(
  "reference_fasta_gzindex",
  !fasta_is_gzipped || file.exists(fasta_gzindex_path),
  fasta_gzindex_path,
  if (fasta_is_gzipped) "FASTA .gzi index if present." else
    "Not applicable: reference FASTA is not bgzip-compressed."
)
add_status(
  "reference_txdb",
  file.exists(txdb_path),
  txdb_path,
  txdb_error,
  status = if (nzchar(txdb_error)) "not_available" else
    if (file.exists(txdb_path)) "ok" else "missing"
)
add_status("all_merged_ofst", file.exists(ofst_path), ofst_path,
           "All-merged ORFik offset library.")

bigwig_candidates <- character()
if (nzchar(lib_folder)) {
  bigwig_candidates <- c(
    file.path(lib_folder, "all_in_one_l50_forward.bigWig"),
    file.path(lib_folder, "all_in_one_l50_reverse.bigWig"),
    file.path(lib_folder, "bigwig", "all_in_onel50_forward.bigWig"),
    file.path(lib_folder, "bigwig", "all_in_onel50_reverse.bigWig"),
    # Naming convention seen on the full server (all_merged-Homo_sapiens):
    # "all_<strand>.bigWig" directly under bigwig/, not "all_in_one_l50_*".
    file.path(lib_folder, "bigwig", "all_forward.bigWig"),
    file.path(lib_folder, "bigwig", "all_reverse.bigWig")
  )
}
forward_bigwigs <- bigwig_candidates[grepl("forward\\.bigWig$", bigwig_candidates)]
reverse_bigwigs <- bigwig_candidates[grepl("reverse\\.bigWig$", bigwig_candidates)]
add_status(
  "all_merged_forward_bigwig",
  any(file.exists(forward_bigwigs)),
  scalar_text(forward_bigwigs[file.exists(forward_bigwigs)]),
  paste(bigwig_candidates, collapse = "; ")
)
add_status(
  "all_merged_reverse_bigwig",
  any(file.exists(reverse_bigwigs)),
  scalar_text(reverse_bigwigs[file.exists(reverse_bigwigs)]),
  paste(bigwig_candidates, collapse = "; ")
)

count_table_kinds <- c("cds", "leaders", "mrna", "trailers", "uorfs")
# The server's ORFik count tables are saved as .qs (qs::qsave), not .rds;
# check both, preferring whichever actually exists on disk.
count_table_path_for <- function(kind) {
  if (!nzchar(lib_folder)) return("")
  candidates <- file.path(lib_folder, "QC_STATS",
                          paste0("countTable_", kind, c(".qs", ".rds")))
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
}
count_table_paths <- vapply(count_table_kinds, count_table_path_for, character(1))
for (i in seq_along(count_table_kinds)) {
  add_status(
    paste0("all_merged_count_table_", count_table_kinds[[i]]),
    file.exists(count_table_paths[[i]]),
    count_table_paths[[i]],
    "All-merged ORFik count table."
  )
}

add_status("all_samples_fst_index", file.exists(fst_index_file),
           fst_index_file, "Indexed all-sample coverage page manifest.")

required_pages <- safe_fread(file.path(analysis_dir, "required_fst_pages.csv"))
required_pages_by_gene <- safe_fread(file.path(
  analysis_dir, "required_fst_pages_by_gene.csv"
))
manifest <- safe_fread(file.path(
  analysis_dir, "dominant_rdg_qualified_inspection",
  "qualified_inspection_manifest.csv"
))
sequence_features <- safe_fread(file.path(
  analysis_dir, "dominant_rdg_qualified_inspection",
  "qualified_sequence_features.csv"
))

page_readiness <- copy(required_pages)
if (!nrow(page_readiness)) {
  page_readiness <- data.table(
    fst_basename = character(),
    expected_local_path_after_download = character()
  )
}
if (!"expected_local_path_after_download" %in% names(page_readiness)) {
  page_readiness[, expected_local_path_after_download := ""]
}
page_readiness[, local_path := normalizePath(
  expected_local_path_after_download, mustWork = FALSE
)]
page_readiness[, local_file_exists_previous := as_bool_vec(
  if ("local_file_exists" %in% names(page_readiness)) {
    local_file_exists
  } else {
    FALSE
  }
)]
page_readiness[, local_file_exists_now := file.exists(local_path)]
page_readiness[, file_size_bytes := vapply(local_path, file_size, numeric(1))]
page_readiness[, file_size_mb := round(file_size_bytes / 1024^2, 2)]
page_readiness[, file_mtime := vapply(local_path, file_mtime, character(1))]
page_readiness[, readiness_status := fifelse(
  local_file_exists_now,
  "present",
  "missing"
)]

required_page_total <- nrow(page_readiness)
required_page_present <- sum(page_readiness$local_file_exists_now)
required_page_missing <- required_page_total - required_page_present
required_page_size <- sum_or_zero(page_readiness$file_size_bytes)
add_status(
  "required_fst_pages",
  required_page_total > 0 && required_page_missing == 0,
  "",
  paste0(required_page_present, "/", required_page_total,
         " pages present; total size ",
         round(required_page_size / 1024^3, 2), " GB."),
  required_page_total,
  status = if (!required_page_total) "missing_manifest" else
    if (!required_page_missing) "ok" else "incomplete"
)

pages_by_tx <- data.table(
  gene_symbol = character(),
  tx_id = character(),
  required_fst_rows = integer(),
  required_fst_page_count = integer(),
  required_fst_pages_present = integer(),
  required_fst_pages_missing = integer(),
  required_fst_page_basenames = character(),
  missing_required_fst_page_basenames = character()
)
if (nrow(required_pages_by_gene)) {
  by_gene <- copy(required_pages_by_gene)
  if (!"expected_local_path_after_download" %in% names(by_gene)) {
    by_gene[, expected_local_path_after_download := ""]
  }
  by_gene[, local_path := normalizePath(
    expected_local_path_after_download, mustWork = FALSE
  )]
  by_gene[, local_file_exists_now := file.exists(local_path)]
  pages_by_tx <- by_gene[, .(
    required_fst_rows = .N,
    required_fst_page_count = uniqueN(fst_basename),
    required_fst_pages_present = uniqueN(
      fst_basename[local_file_exists_now]
    ),
    required_fst_pages_missing = uniqueN(
      fst_basename[!local_file_exists_now]
    ),
    required_fst_page_basenames = paste(sort(unique(fst_basename)),
                                        collapse = ";"),
    missing_required_fst_page_basenames = paste(
      sort(unique(fst_basename[!local_file_exists_now])),
      collapse = ";"
    )
  ), by = .(gene_symbol, tx_id)]
}

sequence_summary <- data.table(source_row = integer())
if (nrow(sequence_features) && "source_row" %in% names(sequence_features)) {
  sf <- copy(sequence_features)
  sequence_summary <- sf[, .(
    rdg_feature_count = .N,
    leader_uorf_features = sum(feature_type %chin% "leader_uorf",
                               na.rm = TRUE),
    overlapping_uorf_features = sum(feature_type %chin% "overlapping_uorf",
                                    na.rm = TRUE),
    internal_orf_features = sum(feature_type %chin% "internal_orf",
                                na.rm = TRUE),
    nte_extension_features = sum(feature_type %chin% "nte_extension",
                                 na.rm = TRUE),
    strong_kozak_features = sum(kozak_strength %chin% "strong",
                                na.rm = TRUE),
    sequence_mrna_available = any(as_bool_vec(mrna_sequence_available)),
    sequence_display_available = any(as_bool_vec(display_sequence_available)),
    sequence_mrna_length = max_or_na(mrna_sequence_length),
    sequence_display_length = max_or_na(display_sequence_length),
    has_nonpositive_feature_coordinates = any(
      suppressWarnings(as.numeric(tx_start)) <= 0 |
        suppressWarnings(as.numeric(tx_end)) <= 0,
      na.rm = TRUE
    )
  ), by = source_row]
}

transcript_readiness <- copy(manifest)
if (!nrow(transcript_readiness)) {
  transcript_readiness <- data.table(
    source_row = integer(),
    gene_symbol = character(),
    tx_id = character(),
    design_family = character()
  )
}
if (!"source_row" %in% names(transcript_readiness)) {
  transcript_readiness[, source_row := seq_len(.N)]
}
transcript_readiness <- merge(
  transcript_readiness,
  sequence_summary,
  by = "source_row",
  all.x = TRUE,
  sort = FALSE
)
if (all(c("gene_symbol", "tx_id") %chin% names(transcript_readiness))) {
  transcript_readiness <- merge(
    transcript_readiness,
    pages_by_tx,
    by = c("gene_symbol", "tx_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
for (column in c(
  "required_fst_rows", "required_fst_page_count",
  "required_fst_pages_present", "required_fst_pages_missing",
  "rdg_feature_count", "leader_uorf_features", "overlapping_uorf_features",
  "internal_orf_features", "nte_extension_features",
  "strong_kozak_features"
)) {
  if (!column %in% names(transcript_readiness)) {
    transcript_readiness[, (column) := 0]
  }
  transcript_readiness[is.na(get(column)), (column) := 0]
}
for (column in c(
  "required_fst_page_basenames", "missing_required_fst_page_basenames",
  # inspection_readiness/status come from qualified_inspection_manifest.csv;
  # on a fresh checkout (that step has not run yet), transcript_readiness
  # falls back to an empty table without these columns, and the
  # sdrive_readiness_class/review_support_class fcase() below reads them
  # unconditionally -- default them so this audit can run before that step.
  "inspection_readiness", "status"
)) {
  if (!column %in% names(transcript_readiness)) {
    transcript_readiness[, (column) := ""]
  }
  transcript_readiness[is.na(get(column)), (column) := ""]
}
for (column in c(
  "sequence_mrna_available", "sequence_display_available",
  "has_nonpositive_feature_coordinates"
)) {
  if (!column %in% names(transcript_readiness)) {
    transcript_readiness[, (column) := FALSE]
  }
  transcript_readiness[is.na(get(column)), (column) := FALSE]
}
manifest_mrna_available <- if ("mrna_sequence_available" %in%
                               names(transcript_readiness)) {
  as_bool_vec(transcript_readiness$mrna_sequence_available)
} else {
  rep(FALSE, nrow(transcript_readiness))
}
manifest_display_available <- if ("display_sequence_available" %in%
                                  names(transcript_readiness)) {
  as_bool_vec(transcript_readiness$display_sequence_available)
} else {
  rep(FALSE, nrow(transcript_readiness))
}
transcript_readiness[, mrna_sequence_available :=
                       manifest_mrna_available | sequence_mrna_available]
transcript_readiness[, display_sequence_available :=
                       manifest_display_available | sequence_display_available]
asset_cols <- intersect(
  c("profile_png", "feature_allocation_png", "replicate_diagnostics_png",
    "rdg_flow_png"),
  names(transcript_readiness)
)
if (length(asset_cols)) {
  transcript_readiness[, static_asset_count := rowSums(as.data.frame(
    lapply(.SD, function(x) file.exists(as.character(x)))
  )), .SDcols = asset_cols]
  transcript_readiness[, static_asset_total := length(asset_cols)]
  transcript_readiness[, static_assets_present :=
                         static_asset_count == static_asset_total]
} else {
  transcript_readiness[, `:=`(
    static_asset_count = 0L,
    static_asset_total = 0L,
    static_assets_present = FALSE
  )]
}

all_merged_core_ready <- file.exists(ofst_path) &&
  any(file.exists(forward_bigwigs)) &&
  any(file.exists(reverse_bigwigs)) &&
  all(file.exists(count_table_paths))
reference_ready <- dir.exists(ref_folder) && file.exists(fasta_path) &&
  file.exists(fasta_index_path)
transcript_readiness[, all_merged_core_ready := all_merged_core_ready]
transcript_readiness[, reference_ready := reference_ready]
transcript_readiness[, sdrive_mount_ready := mount_ok]
transcript_readiness[, fst_pages_ready :=
                       required_fst_page_count > 0 &
                         required_fst_pages_missing == 0]
transcript_readiness[, sequence_ready :=
                       mrna_sequence_available &
                         display_sequence_available]
if (!"status" %in% names(transcript_readiness)) {
  transcript_readiness[, status := ""]
}
transcript_readiness[, static_inspection_status := fifelse(
  status %chin% "ok" & static_assets_present,
  "static_assets_ready",
  fifelse(status %chin% "ok", "static_assets_missing",
          "inspection_packet_failed")
)]
transcript_readiness[, sdrive_readiness_class := fcase(
  !sdrive_mount_ready, "sdrive_not_mounted",
  !(status %chin% "ok"), "inspection_packet_failed",
  !static_assets_present, "missing_static_assets",
  required_fst_page_count == 0, "no_required_fst_mapping",
  required_fst_pages_missing > 0, "missing_required_fst_pages",
  !reference_ready, "reference_unavailable",
  !all_merged_core_ready, "all_merged_unavailable",
  !sequence_ready, "sequence_unavailable",
  default = "sdrive_ready"
)]
transcript_readiness[, review_support_class := fcase(
  inspection_readiness %chin% c("strong_visual_review_ready",
                                "qualified_visual_review_ready"),
  "coverage_reviewable",
  inspection_readiness %chin% "spike_risk_review_required",
  "spike_risk_review_required",
  inspection_readiness %chin% "not_claimable_from_coverage",
  "not_claimable_static_coverage",
  inspection_readiness %chin% "incomplete_inputs",
  "incomplete_static_inputs",
  default = "needs_manual_triage"
)]

priority_cols <- c(
  "source_row", "global_review_order", "batch_id", "batch_role",
  "gene_symbol", "tx_id", "design_family", "study", "tissue", "cell_line",
  "case_condition", "control_conditions", "inspection_readiness",
  "review_support_class", "sdrive_readiness_class",
  "static_inspection_status", "fst_pages_ready", "required_fst_page_count",
  "required_fst_pages_present", "required_fst_pages_missing",
  "mrna_sequence_available", "display_sequence_available",
  "reference_ready", "all_merged_core_ready", "rdg_feature_count",
  "leader_uorf_features", "overlapping_uorf_features",
  "internal_orf_features", "nte_extension_features",
  "strong_kozak_features", "has_nonpositive_feature_coordinates",
  "branch_min_case_control_counts", "display_min_exact_counts",
  "top_feature_min_exact_counts", "max_exact_position_fraction",
  "top_feature_max_run_fraction", "top_feature_in_rdg_annotation",
  "has_internal_or_nte_architecture", "flags",
  "required_fst_page_basenames", "missing_required_fst_page_basenames",
  "profile_png", "feature_allocation_png", "replicate_diagnostics_png",
  "rdg_flow_png"
)
priority_cols <- intersect(priority_cols, names(transcript_readiness))
transcript_readiness <- transcript_readiness[, ..priority_cols]
if ("source_row" %in% names(transcript_readiness)) {
  setorder(transcript_readiness, source_row)
}

status <- rbindlist(status_rows, fill = TRUE)
status[, generated_at := format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")]

class_counts <- if (nrow(transcript_readiness)) {
  transcript_readiness[, .N, by = .(sdrive_readiness_class)][
    order(-N, sdrive_readiness_class)
  ]
} else {
  data.table(sdrive_readiness_class = character(), N = integer())
}

metric_rows <- list(
  data.table(metric = "sdrive_mount_ready", value = as.character(mount_ok),
             detail = sdrive_root),
  data.table(metric = "experiment_loaded",
             value = as.character(experiment_ok), detail = experiment_name),
  data.table(metric = "required_fst_pages",
             value = as.character(required_page_total), detail = ""),
  data.table(metric = "required_fst_pages_present",
             value = as.character(required_page_present), detail = ""),
  data.table(metric = "required_fst_pages_missing",
             value = as.character(required_page_missing), detail = ""),
  data.table(metric = "required_fst_pages_size_gb",
             value = as.character(round(required_page_size / 1024^3, 3)),
             detail = ""),
  data.table(metric = "qualified_source_rows",
             value = as.character(nrow(transcript_readiness)), detail = ""),
  data.table(metric = "qualified_transcripts",
             value = as.character(uniqueN(
               paste(transcript_readiness$gene_symbol,
                     transcript_readiness$tx_id, sep = "\r")
             )), detail = ""),
  data.table(metric = "sdrive_ready_source_rows",
             value = as.character(sum(
               transcript_readiness$sdrive_readiness_class %chin%
                 "sdrive_ready"
             )), detail = ""),
  data.table(metric = "coverage_reviewable_source_rows",
             value = as.character(sum(
               transcript_readiness$review_support_class %chin%
                 "coverage_reviewable"
             )), detail = ""),
  data.table(metric = "spike_risk_source_rows",
             value = as.character(sum(
               transcript_readiness$review_support_class %chin%
                 "spike_risk_review_required"
             )), detail = ""),
  data.table(metric = "not_claimable_static_source_rows",
             value = as.character(sum(
               transcript_readiness$review_support_class %chin%
                 "not_claimable_static_coverage"
             )), detail = "")
)
if (nrow(class_counts)) {
  metric_rows <- c(metric_rows, lapply(seq_len(nrow(class_counts)), function(i) {
    data.table(
      metric = paste0("sdrive_class_", class_counts$sdrive_readiness_class[[i]]),
      value = as.character(class_counts$N[[i]]),
      detail = ""
    )
  }))
}
summary_metrics <- rbindlist(metric_rows, fill = TRUE)
summary_metrics[, generated_at := format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")]

fwrite(status, file.path(output_dir, "sdrive_experiment_status.csv"))
fwrite(page_readiness, file.path(output_dir, "sdrive_fst_page_readiness.csv"))
fwrite(transcript_readiness,
       file.path(output_dir, "sdrive_transcript_readiness.csv"))
fwrite(summary_metrics,
       file.path(output_dir, "sdrive_readiness_summary_metrics.csv"))
fwrite(class_counts,
       file.path(output_dir, "sdrive_readiness_class_counts.csv"))

summary_lines <- c(
  "# S-Drive Readiness Audit",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("S-drive root:", sdrive_root),
  paste("ORFik experiment:", experiment_name),
  paste("Experiment loaded:", experiment_ok),
  paste("Reference ready:", reference_ready),
  paste("All-merged core ready:", all_merged_core_ready),
  "",
  "## Indexed Coverage Pages",
  "",
  paste("- Required FST pages:", required_page_total),
  paste("- Present FST pages:", required_page_present),
  paste("- Missing FST pages:", required_page_missing),
  paste("- Present page size:", round(required_page_size / 1024^3, 2), "GB"),
  "",
  "## Qualified Inspection Sources",
  "",
  paste("- Source rows:", nrow(transcript_readiness)),
  paste("- Distinct transcript keys:", uniqueN(
    paste(transcript_readiness$gene_symbol,
          transcript_readiness$tx_id, sep = "\r")
  )),
  paste("- S-drive ready rows:", sum(
    transcript_readiness$sdrive_readiness_class %chin% "sdrive_ready"
  )),
  paste("- Coverage-reviewable rows:", sum(
    transcript_readiness$review_support_class %chin% "coverage_reviewable"
  )),
  paste("- Spike-risk rows:", sum(
    transcript_readiness$review_support_class %chin%
      "spike_risk_review_required"
  )),
  paste("- Not-claimable static coverage rows:", sum(
    transcript_readiness$review_support_class %chin%
      "not_claimable_static_coverage"
  )),
  "",
  "## Readiness Classes",
  "",
  if (nrow(class_counts)) {
    paste0("- ", class_counts$sdrive_readiness_class, ": ", class_counts$N)
  } else {
    "- No qualified inspection rows were available."
  },
  "",
  "## Interpretation",
  "",
  "This audit separates backing-data availability from biological claim quality.",
  "`sdrive_ready` means the mounted reference, all-merged ORFik files, required indexed FST pages, mRNA/display sequence flags, and static inspection assets are all present.",
  "It does not mean the coverage pattern is biologically convincing; use `review_support_class`, the Inspection Assets page, and Observatory review for that decision."
)
writeLines(summary_lines,
           file.path(output_dir, "sdrive_readiness_summary.md"))

message("Wrote S-drive readiness audit to: ", output_dir)
message("Required FST pages present: ", required_page_present, "/",
        required_page_total)
message("S-drive ready qualified rows: ",
        sum(transcript_readiness$sdrive_readiness_class %chin% "sdrive_ready"),
        "/", nrow(transcript_readiness))
