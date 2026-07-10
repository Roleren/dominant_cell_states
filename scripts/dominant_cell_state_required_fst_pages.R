dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

# Avoid probing unavailable server-backed paths, but fill the readiness column
# automatically once the local page directory is mounted/downloaded. Override
# with DOMINANT_CHECK_LOCAL_FSTS=true/false if needed.
check_local_downloaded_files <- tolower(Sys.getenv(
  "DOMINANT_CHECK_LOCAL_FSTS",
  "auto"
))
if (!check_local_downloaded_files %chin% c("auto", "true", "false", "1", "0",
                                           "yes", "no")) {
  stop("DOMINANT_CHECK_LOCAL_FSTS must be auto, true, or false.",
       call. = FALSE)
}

message("Dominant-state atlas FST page discovery:")
message("  1. Select all marker/prior genes and representative transcripts used by the atlas")
message("  2. Build each transcript leader+CDS display region")
message("  3. Match those regions to coverage_index.fst page intervals")
message("  4. Save one complete page-list product and a gene-to-page audit table")

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."
results_dir <- Sys.getenv(
  "DOMINANT_CELL_STATES_RESULTS_DIR",
  unset = file.path(analysis_dir, "results")
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- function(...) file.path(results_dir, ...)

output_unique_file <- result_file("required_fst_pages.csv")
output_detail_file <- result_file("required_fst_pages_by_gene.csv")
output_selected_file <- result_file("required_fst_selected_transcripts.csv")
output_download_list_file <- result_file("required_fst_pages.txt")
output_basename_list_file <- result_file("required_fst_page_basenames.txt")
output_summary_file <- result_file("required_fst_pages_summary.txt")

module_definition_file <- file.path(
  analysis_dir,
  "scripts",
  "dominant_state_module_definitions.R"
)
if (!file.exists(module_definition_file)) {
  stop("Missing dominant-state module definition file: ", module_definition_file)
}
source(module_definition_file)
dominant_state_modules <- load_dominant_state_modules(analysis_dir)

module_membership <- rbindlist(lapply(seq_len(nrow(dominant_state_modules)), function(i) {
  module <- dominant_state_modules[i]
  rbindlist(list(
    data.table(
      gene_symbol = toupper(module$up_genes[[1]]),
      dominant_state = module$dominant_state,
      model_direction = if (module$dominant_state == "Baseline (control)") {
        "baseline"
      } else {
        "up"
      }
    ),
    data.table(
      gene_symbol = toupper(module$down_genes[[1]]),
      dominant_state = module$dominant_state,
      model_direction = "down"
    )
  ))
}))

gene_membership <- module_membership[
  ,
  .(
    dominant_states = paste(unique(dominant_state), collapse = "; "),
    model_directions = paste(unique(model_direction), collapse = "; ")
  ),
  by = gene_symbol
]

range_positions <- function(grl) {
  gr <- unlist(grl, use.names = FALSE)
  if (length(gr) == 0) return(integer())
  unique(unlist(Map(seq.int, start(gr), end(gr)), use.names = FALSE))
}

safe_pmap_positions <- function(query, subject) {
  if (length(query) == 0 || length(subject) == 0) return(integer())
  mapped <- tryCatch(
    suppressWarnings(pmapToTranscriptF(query, subject)),
    error = function(e) GRangesList()
  )
  range_positions(mapped)
}

make_display_region <- function(leader_tx, cds_tx, tx_id) {
  parts <- if (length(leader_tx) > 0 && widthPerGroup(leader_tx) > 0) {
    c(leader_tx, cds_tx)
  } else {
    cds_tx
  }
  display <- GRangesList(reduce(unlistGrl(parts)))
  display <- sortPerGroup(display, quick.rev = TRUE)
  names(display) <- tx_id
  display
}

load_translon_source <- function(paths, source) {
  path <- unique(paths)
  path <- path[file.exists(path)][1]
  if (is.na(path)) return(GRangesList())
  translons <- read_RDSQS(path)
  mcols(translons)$translon_source <- source
  translons
}

load_predicted_translons <- function(df) {
  t_paths <- c(
    file.path(refFolder(df), "predicted_translons",
              "predicted_translons_with_sequence_ranges.rds"),
    file.path(dirname(df@fafile), "predicted_translons",
              "predicted_translons_with_sequence_ranges.rds")
  )
  tc_paths <- c(
    file.path(refFolder(df), "predicted_translons", "TransCode",
              "predicted_translons_ranges.qs"),
    file.path(dirname(df@fafile), "predicted_translons", "TransCode",
              "predicted_translons_ranges.qs")
  )
  m_paths <- c(
    file.path("dominant_cell_states", "manual_translons", "manual_translons_ranges.rds"),
    file.path("manual_translons", "manual_translons_ranges.rds")
  )
  t <- load_translon_source(t_paths, "T")
  tc <- load_translon_source(tc_paths, "TC")
  m <- load_translon_source(m_paths, "M")
  if (length(t) == 0 && length(tc) == 0 && length(m) == 0) {
    expected <- unique(c(t_paths, tc_paths, m_paths))
    warning("Predicted translon range file not found. Checked: ",
            paste(expected, collapse = "; "))
    return(GRangesList())
  }
  message("Loaded translon sources: T=", length(t), ", TC=", length(tc),
          ", M=", length(m))
  c(t, tc, m)
}

load_symbol_table <- function(df) {
  sym <- tryCatch(
    suppressMessages(symbols(df)),
    error = function(e) {
      warning("symbols(df) failed: ", conditionMessage(e))
      data.table()
    }
  )
  if (nrow(sym) > 0 &&
      all(c("external_gene_name", "ensembl_tx_name") %chin% names(sym))) {
    return(as.data.table(sym))
  }

  table_path <- file.path(refFolder(df), "predicted_translons",
                          "predicted_translons_with_sequence.fst")
  if (file.exists(table_path)) {
    tt <- fst::read_fst(table_path, as.data.table = TRUE)
    needed <- c("external_gene_name", "ensembl_tx_name")
    if (all(needed %chin% names(tt))) {
      return(unique(tt[!is.na(external_gene_name) & !is.na(ensembl_tx_name),
                       ..needed]))
    }
  }

  stop("Could not load transcript-to-gene symbols. Tried symbols(df) and ",
       table_path,
       ". Check that the reference folder is mounted and contains gene symbols ",
       "or predicted_translons_with_sequence.fst.")
}

collapse_to_tx_iranges <- function(mapped) {
  starts <- vapply(mapped, function(x) {
    if (length(x) == 0) return(NA_integer_)
    min(start(x))
  }, integer(1))
  ends <- vapply(mapped, function(x) {
    if (length(x) == 0) return(NA_integer_)
    max(end(x))
  }, integer(1))
  IRanges(start = starts, end = ends)
}

safe_pmap_translons_to_transcript <- function(translons, mrna) {
  mapped_list <- lapply(seq_along(translons), function(i) {
    tx_id <- names(translons)[i]
    tryCatch(
      suppressWarnings(pmapToTranscriptF(translons[i], mrna[tx_id])[[1]]),
      error = function(e) GRanges()
    )
  })
  mapped <- GRangesList(mapped_list)
  names(mapped) <- names(translons)
  mapped
}

annotate_translon_categories <- function(translons, cds, mrna) {
  if (length(translons) == 0) return(translons)

  keep <- names(translons) %chin% names(mrna)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " predicted translons without mRNA ranges.")
    translons <- translons[keep]
  }
  if (length(translons) == 0) return(translons)

  mapped <- safe_pmap_translons_to_transcript(translons, mrna)
  tx_iranges <- collapse_to_tx_iranges(mapped)
  valid <- !is.na(start(tx_iranges)) & !is.na(end(tx_iranges)) &
    end(tx_iranges) >= start(tx_iranges)
  if (!all(valid)) {
    warning("Dropping ", sum(!valid),
            " predicted translons with invalid transcript-coordinate ranges.")
    translons <- translons[valid]
    tx_iranges <- tx_iranges[valid]
  }
  if (length(translons) == 0) return(translons)
  categories <- rep(NA_character_, length(translons))

  valid <- !is.na(start(tx_iranges)) & !is.na(end(tx_iranges))
  if (any(valid)) {
    groupings <- match(names(translons)[valid], names(mrna))
    categories[valid] <- ORFik:::categorize_ORFs(
      tx_iranges[valid],
      groupings = as.integer(groupings),
      cds = cds,
      mrna = mrna,
      verbose = FALSE
    )
  }

  mcols(translons)$category <- categories
  mcols(translons)$tx_start <- start(tx_iranges)
  mcols(translons)$tx_end <- end(tx_iranges)
  translons
}

deduplicate_translons_by_tx_coordinates <- function(translons) {
  if (length(translons) == 0) return(translons)
  if (!all(c("tx_start", "tx_end") %chin% names(mcols(translons)))) {
    stop("Translons must be transcript-coordinate annotated before deduplication.")
  }
  key <- paste(names(translons), mcols(translons)$tx_start,
               mcols(translons)$tx_end, sep = "|")
  source <- as.character(mcols(translons)$translon_source)
  source[is.na(source)] <- "unknown"
  merged_source <- vapply(split(source, key), function(x) {
    paste(unique(x), collapse = "+")
  }, character(1))
  mcols(translons)$translon_sources_merged <- unname(merged_source[key])
  keep <- !duplicated(key)
  dropped <- sum(!keep)
  if (dropped > 0) {
    message("Dropped ", dropped,
            " duplicate translons with identical transcript start/stop.")
  }
  translons[keep]
}

is_cds_masking_uorf <- function(category) {
  as.character(category) %chin% c("uORF", "uoORF")
}

overlapping_uorf_positions <- function(translons, tx_id, leader_tx, cds_tx,
                                       display_region, cds_positions) {
  tx_translons <- translons[names(translons) == tx_id]
  if (length(tx_translons) == 0 || length(leader_tx) == 0) {
    return(list(
      positions = integer(),
      n_overlapping_uorfs = 0L,
      n_full_cds_spanning_translons = 0L,
      n_ignored_overlap_translons = 0L,
      ignored_overlap_categories = NA_character_,
      masking_uorf_categories = NA_character_
    ))
  }
  if (!"category" %in% names(mcols(tx_translons))) {
    stop("Predicted translons must be categorized before clean-CDS masking.")
  }

  leader_gr <- unlist(leader_tx, use.names = FALSE)
  cds_gr <- unlist(cds_tx, use.names = FALSE)
  excluded <- integer()
  n_overlapping <- 0L
  n_full_cds_spanning <- 0L
  n_ignored <- 0L
  ignored_categories <- character()
  masking_categories <- character()

  for (i in seq_along(tx_translons)) {
    translon_gr <- tx_translons[[i]]
    overlaps_leader <- length(findOverlaps(translon_gr, leader_gr,
                                           ignore.strand = FALSE)) > 0
    overlaps_cds <- length(findOverlaps(translon_gr, cds_gr,
                                        ignore.strand = FALSE)) > 0
    if (!overlaps_leader || !overlaps_cds) next

    category <- as.character(mcols(tx_translons)$category[i])
    if (is.na(category)) {
      stop("Could not categorize a CDS-overlapping predicted translon for ", tx_id)
    }
    if (!is_cds_masking_uorf(category)) {
      n_ignored <- n_ignored + 1L
      ignored_categories <- c(ignored_categories, category)
      next
    }

    translon_pos <- safe_pmap_positions(GRangesList(translon_gr), display_region)
    excluded_this <- intersect(cds_positions, translon_pos)
    if (length(excluded_this) >= length(cds_positions)) {
      n_full_cds_spanning <- n_full_cds_spanning + 1L
    }
    n_overlapping <- n_overlapping + 1L
    masking_categories <- c(masking_categories, category)
    excluded <- c(excluded, excluded_this)
  }

  list(
    positions = unique(excluded),
    n_overlapping_uorfs = n_overlapping,
    n_full_cds_spanning_translons = n_full_cds_spanning,
    n_ignored_overlap_translons = n_ignored,
    ignored_overlap_categories = if (length(ignored_categories) == 0) {
      NA_character_
    } else {
      paste(sort(unique(ignored_categories)), collapse = ";")
    },
    masking_uorf_categories = if (length(masking_categories) == 0) {
      NA_character_
    } else {
      paste(sort(unique(masking_categories)), collapse = ";")
    }
  )
}

choose_gene_transcript <- function(gene_symbol, tx_ids, cds_all, leader_all,
                                   translons) {
  candidates <- rbindlist(lapply(tx_ids, function(tx_id) {
    if (!tx_id %in% names(cds_all)) return(NULL)
    cds_tx <- cds_all[tx_id]
    leader_tx <- if (tx_id %in% names(leader_all)) leader_all[tx_id] else GRangesList()
    display_region <- make_display_region(leader_tx, cds_tx, tx_id)
    cds_positions <- safe_pmap_positions(cds_tx, display_region)
    excluded <- overlapping_uorf_positions(
      translons, tx_id, leader_tx, cds_tx, display_region, cds_positions
    )
    clean_positions <- setdiff(cds_positions, excluded$positions)

    data.table(
      gene_symbol = gene_symbol,
      tx_id = tx_id,
      cds_bases = length(cds_positions),
      leader_bases = if (length(leader_tx) > 0) widthPerGroup(leader_tx) else 0L,
      clean_cds_bases = length(clean_positions),
      excluded_cds_bases = length(excluded$positions),
      n_overlapping_uorfs = excluded$n_overlapping_uorfs,
      n_full_cds_spanning_translons = excluded$n_full_cds_spanning_translons,
      n_ignored_overlap_translons = excluded$n_ignored_overlap_translons,
      ignored_overlap_categories = excluded$ignored_overlap_categories,
      masking_uorf_categories = excluded$masking_uorf_categories,
      is_canonical_isoform = tx_id %chin% canonical_tx
    )
  }), fill = TRUE)

  if (nrow(candidates) == 0) return(NULL)
  candidates[, has_overlapping_uorf := n_overlapping_uorfs > 0]
  candidates[, has_clean_cds := clean_cds_bases > 0]
  setorder(candidates, -has_clean_cds, -has_overlapping_uorf,
           -is_canonical_isoform, -clean_cds_bases, -leader_bases, tx_id)
  candidates[, has_overlapping_uorf := NULL]
  candidates[, has_clean_cds := NULL]
  candidates[1]
}

query_fst_pages_for_transcript <- function(selected, cds_all, leader_all,
                                           index, fst_index) {
  tx_id <- selected$tx_id
  if (is.na(tx_id) || !tx_id %in% names(cds_all)) return(NULL)

  cds_tx <- cds_all[tx_id]
  leader_tx <- if (tx_id %in% names(leader_all)) leader_all[tx_id] else GRangesList()
  display_region <- make_display_region(leader_tx, cds_tx, tx_id)
  gr <- unlist(display_region, use.names = FALSE)
  if (length(gr) == 0) return(NULL)

  strand_value <- unique(as.character(strand(gr)))
  strand_value <- strand_value[strand_value %chin% c("+", "-")]
  if (length(strand_value) != 1) {
    stop("Expected one concrete strand for ", selected$gene_symbol,
         " (", tx_id, "), found: ", paste(unique(as.character(strand(gr))), collapse = ", "))
  }
  strand_direction <- if (strand_value == "+") "forward" else "reverse"

  out <- rbindlist(lapply(seq_along(gr), function(i) {
    query_chr <- as.character(seqnames(gr)[i])
    query_start <- start(gr)[i]
    query_end <- end(gr)[i]
    hits <- index[
      chr == query_chr &
        page_start <= query_end &
        page_end >= query_start
    ]
    if (nrow(hits) == 0) {
      stop("No FST index page matched ", selected$gene_symbol, " (", tx_id,
           ") at ", query_chr, ":", query_start, "-", query_end)
    }
    hits[, `:=`(
      query_id = i,
      query_chr = query_chr,
      query_start = query_start,
      query_end = query_end,
      query_strand = strand_value,
      direction = strand_direction,
      start_segment = pmax(page_start - query_start + 1L, 1L),
      end_segment = pmin(page_end - query_start + 1L, query_end - query_start + 1L)
    )]
    hits
  }), fill = TRUE)

  selected_file_column <- if (strand_direction == "forward") "file_forward" else "file_reverse"
  out[, fst_file_from_index := get(selected_file_column)]
  out[, expected_local_path_after_download :=
        file.path(dirname(fst_index), basename(fst_file_from_index))]
  out[, fst_basename := basename(fst_file_from_index)]
  out[, `:=`(
    gene_symbol = selected$gene_symbol,
    tx_id = tx_id,
    cds_bases = selected$cds_bases,
    leader_bases = selected$leader_bases,
    clean_cds_bases = selected$clean_cds_bases,
    excluded_cds_bases = selected$excluded_cds_bases,
    n_overlapping_uorfs = selected$n_overlapping_uorfs,
    n_full_cds_spanning_translons = selected$n_full_cds_spanning_translons,
    n_ignored_overlap_translons = selected$n_ignored_overlap_translons,
    ignored_overlap_categories = selected$ignored_overlap_categories,
    masking_uorf_categories = selected$masking_uorf_categories,
    is_canonical_isoform = selected$is_canonical_isoform
  )]
  out[]
}

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
symbols <- load_symbol_table(df)
translons <- load_predicted_translons(df)
canonical_tx <- canonical_isoforms(df)
fst_index <- file.path(collection_dir_from_exp(df), "coverage_index.fst")
if (!file.exists(fst_index)) stop("Missing FST coverage index: ", fst_index)

index <- fst::read_fst(fst_index, as.data.table = TRUE)
required_index_columns <- c("chr", "part", "start", "end", "file_forward", "file_reverse")
if (!all(required_index_columns %chin% names(index))) {
  stop("coverage_index.fst is missing columns: ",
       paste(setdiff(required_index_columns, names(index)), collapse = ", "))
}
setnames(index, c("start", "end"), c("page_start", "page_end"))
index[, part := as.integer(part)]

target_genes <- dominant_state_target_genes(
  dominant_state_modules,
  analysis_dir,
  include_postviral_auxiliary = TRUE
)
target_map <- unique(symbols[
  toupper(external_gene_name) %chin% target_genes,
  .(gene_symbol = toupper(external_gene_name), tx_id = ensembl_tx_name)
])

target_tx <- unique(target_map$tx_id)
cds_all <- loadRegion(df, part = "cds", names.keep = target_tx)
leader_all <- loadRegion(df, part = "leaders", names.keep = target_tx)
mrna_all <- loadRegion(df, part = "mrna", names.keep = target_tx)
translons <- translons[names(translons) %chin% target_tx]
translons <- annotate_translon_categories(translons, cds_all, mrna_all)
translons <- deduplicate_translons_by_tx_coordinates(translons)
message("Categorized predicted translons for target transcripts: ",
        paste(names(table(mcols(translons)$category)),
              as.integer(table(mcols(translons)$category)),
              sep = "=", collapse = ", "))
message("Translon source composition after deduplication: ",
        paste(names(table(mcols(translons)$translon_sources_merged)),
              as.integer(table(mcols(translons)$translon_sources_merged)),
              sep = "=", collapse = ", "))

selected_transcripts <- rbindlist(lapply(target_genes, function(gene) {
  tx_ids <- target_map[gene_symbol == gene, tx_id]
  if (length(tx_ids) == 0) {
    return(data.table(
      gene_symbol = gene,
      tx_id = NA_character_,
      cds_bases = NA_integer_,
      leader_bases = NA_integer_,
      clean_cds_bases = NA_integer_,
      excluded_cds_bases = NA_integer_,
      n_overlapping_uorfs = NA_integer_,
      n_full_cds_spanning_translons = NA_integer_,
      n_ignored_overlap_translons = NA_integer_,
      ignored_overlap_categories = NA_character_,
      masking_uorf_categories = NA_character_,
      is_canonical_isoform = NA
    ))
  }
  choose_gene_transcript(gene, tx_ids, cds_all, leader_all, translons)
}), fill = TRUE)

selected_transcripts <- merge(
  selected_transcripts,
  gene_membership,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
selected_transcripts[, selection_rank := seq_len(.N)]

message("Selected ", sum(!is.na(selected_transcripts$tx_id)),
        " transcripts for ", nrow(selected_transcripts), " marker genes.")

page_detail <- rbindlist(lapply(seq_len(nrow(selected_transcripts)), function(i) {
  selected <- selected_transcripts[i]
  query_fst_pages_for_transcript(selected, cds_all, leader_all, copy(index), fst_index)
}), fill = TRUE)

page_detail <- merge(
  page_detail,
  selected_transcripts[, .(gene_symbol, selection_rank, dominant_states, model_directions)],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

setcolorder(page_detail, c(
  "gene_symbol", "selection_rank", "dominant_states", "model_directions", "tx_id",
  "query_id", "query_chr", "query_strand", "query_start", "query_end",
  "chr", "part", "page_start", "page_end", "direction", "start_segment",
  "end_segment", "fst_basename", "fst_file_from_index",
  "expected_local_path_after_download", "cds_bases", "leader_bases",
  "clean_cds_bases", "excluded_cds_bases", "n_overlapping_uorfs"
  , "n_full_cds_spanning_translons", "n_ignored_overlap_translons",
  "ignored_overlap_categories", "masking_uorf_categories",
  "is_canonical_isoform"
))

unique_pages <- unique(page_detail[, .(
  direction,
  chr,
  part,
  page_start,
  page_end,
  fst_basename,
  fst_file_from_index,
  expected_local_path_after_download
)])
page_usage <- page_detail[
  ,
  .(
    n_genes = uniqueN(gene_symbol),
    genes = paste(sort(unique(gene_symbol)), collapse = ", "),
    n_transcripts = uniqueN(tx_id)
  ),
  by = fst_basename
]
unique_pages <- merge(unique_pages, page_usage, by = "fst_basename", all.x = TRUE, sort = FALSE)
local_check_enabled <- check_local_downloaded_files %chin% c("true", "1", "yes")
if (check_local_downloaded_files == "auto") {
  local_dirs <- unique(dirname(unique_pages$expected_local_path_after_download))
  local_check_enabled <- any(dir.exists(local_dirs))
}
if (local_check_enabled) {
  unique_pages[, local_file_exists := file.exists(expected_local_path_after_download)]
} else {
  unique_pages[, local_file_exists := NA]
}
setorder(unique_pages, chr, part, direction, fst_basename)

setorder(page_detail, selection_rank, query_id, part, direction, fst_basename)
setorder(selected_transcripts, selection_rank)

fwrite(unique_pages, output_unique_file)
fwrite(page_detail, output_detail_file)
fwrite(selected_transcripts, output_selected_file)
writeLines(unique_pages$fst_file_from_index, output_download_list_file)
writeLines(unique_pages$fst_basename, output_basename_list_file)
writeLines(c(
  "Dominant-state atlas required FST page summary",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("FST index: ", fst_index),
  "",
  paste0("Target genes: ", length(target_genes)),
  paste0("Genes resolved to at least one transcript: ", uniqueN(target_map$gene_symbol)),
  paste0("Selected transcripts: ", sum(!is.na(selected_transcripts$tx_id))),
  paste0("Unique required FST pages: ", nrow(unique_pages)),
  paste0("Forward pages: ", unique_pages[direction == "forward", .N]),
  paste0("Reverse pages: ", unique_pages[direction == "reverse", .N]),
  paste0("Chromosomes covered: ", paste(sort(unique(unique_pages$chr)), collapse = ", ")),
  "",
  "Canonical files:",
  paste0("- ", output_unique_file),
  paste0("- ", output_download_list_file),
  paste0("- ", output_basename_list_file),
  paste0("- ", output_detail_file),
  paste0("- ", output_selected_file)
), output_summary_file)

message("Saved unique page table: ", output_unique_file)
message("Saved gene-to-page table: ", output_detail_file)
message("Saved selected transcript table: ", output_selected_file)
message("Saved download path list: ", output_download_list_file)
message("Saved basename list: ", output_basename_list_file)
message("Saved summary: ", output_summary_file)
message("Unique page FST files: ", nrow(unique_pages))
message("  Forward: ", unique_pages[direction == "forward", .N])
message("  Reverse: ", unique_pages[direction == "reverse", .N])
message("Chromosomes covered: ", paste(sort(unique(unique_pages$chr)), collapse = ", "))
