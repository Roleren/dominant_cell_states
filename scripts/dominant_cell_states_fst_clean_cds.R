dominant_loader <- c(
  file.path("scripts", "dominant_ribocrypt_loader.R"),
  file.path("dominant_cell_states", "scripts", "dominant_ribocrypt_loader.R")
)
dominant_loader <- dominant_loader[file.exists(dominant_loader)][1]
if (!is.na(dominant_loader)) {
  source(dominant_loader)
  dominant_load_ribocrypt_if_available()
}

library(data.table)
library(GenomicRanges)
library(IRanges)

analysis_dir <- if (dir.exists("dominant_cell_states")) "dominant_cell_states" else "."

output_state_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds.csv"
)
output_expr_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds_gene_expression.csv"
)
output_diagnostics_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"
)
output_cds_exon_qc_file <- file.path(
  analysis_dir, "human_dominant_cell_states_cds_exon_coverage_qc.csv"
)
allow_count_table_fallback <- FALSE

message("Dominant-state QC annotation from FST coverage:")
message("  1. Load per-base coverage through coverageByTranscriptFST()")
message("  2. Use leader+CDS transcript coordinates to identify CDS bases")
message("  3. Remove CDS bases overlapped by predicted uORFs that span leader and CDS")
message("  4. Sum clean CDS coverage and convert to FPKM-like coverage")
message("  5. Score modules as log2(clean CDS FPKM + 1) minus baseline genes")

module_definition_file <- file.path(
  analysis_dir,
  "scripts",
  "dominant_state_module_definitions.R"
)
if (!file.exists(module_definition_file)) {
  stop("Missing dominant-state module definition file: ", module_definition_file)
}
source(module_definition_file)
coverage_cache_file <- file.path(analysis_dir, "scripts", "dominant_coverage_cache.R")
if (!file.exists(coverage_cache_file)) {
  stop("Missing dominant coverage cache helper: ", coverage_cache_file)
}
source(coverage_cache_file)
dominant_state_modules <- load_dominant_state_modules(analysis_dir)

gene_list_text <- function(x) {
  vapply(x, paste, collapse = ", ", FUN.VALUE = character(1))
}

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

get_sample_names <- function(df) {
  sample <- bamVarName(df, skip.replicate = FALSE, skip.condition = FALSE,
                       skip.stage = TRUE, skip.fraction = FALSE,
                       skip.experiment = TRUE)
  sample <- sub("_PRJ[A-Z0-9]+", "", sample)
  sample <- sub("_GSE[0-9]+", "", sample)
  sample
}

library_sizes_by_run <- function(df) {
  lib_sizes <- readRDS(get_lib_sizes_file(df))
  if (length(lib_sizes) == nrow(df)) {
    names(lib_sizes) <- runIDs(df)
  }
  lib_sizes
}

score_module <- function(gene_z, up_genes, down_genes = character()) {
  up_present <- intersect(up_genes, rownames(gene_z))
  down_present <- intersect(down_genes, rownames(gene_z))

  up_score <- if (length(up_present) == 0) {
    rep(NA_real_, ncol(gene_z))
  } else {
    colMeans(gene_z[up_present, , drop = FALSE], na.rm = TRUE)
  }

  down_score <- if (length(down_genes) == 0) {
    rep(0, ncol(gene_z))
  } else if (length(down_present) == 0) {
    rep(NA_real_, ncol(gene_z))
  } else {
    colMeans(gene_z[down_present, , drop = FALSE], na.rm = TRUE)
  }

  up_score - down_score
}

strip_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

get_canonical_isoform_info <- function(df) {
  canonical_raw <- tryCatch(
    canonical_isoforms(df),
    error = function(e) {
      warning("canonical_isoforms(df) unavailable: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(canonical_raw) || length(canonical_raw) == 0) {
    return(list(
      tx = character(),
      source = "fallback_no_canonical_isoforms",
      note = "canonical_isoforms(df) unavailable or empty; fallback transcript ranking used"
    ))
  }
  canonical_tx <- unique(strip_version(canonical_raw))
  canonical_tx <- canonical_tx[!is.na(canonical_tx) & nzchar(canonical_tx)]
  if (length(canonical_tx) == 0) {
    return(list(
      tx = character(),
      source = "fallback_empty_canonical_isoforms",
      note = "canonical_isoforms(df) returned no usable transcript ids; fallback transcript ranking used"
    ))
  }
  list(
    tx = canonical_tx,
    source = "canonical_isoforms(df)",
    note = NA_character_
  )
}

empty_cds_exon_coverage_qc <- function() {
  data.table(
    transcript_id = character(),
    transcript_index = integer(),
    exon_id = character(),
    exon_rank = integer(),
    exon_width = integer(),
    exon_counts = numeric(),
    exon_mean_coverage = numeric(),
    transcript_exon_mean_coverage = numeric(),
    exon_vs_transcript_mean = numeric(),
    transcript_counts = numeric(),
    n_exons = integer(),
    exon_qc_measurable = logical(),
    low_coverage_exon = logical(),
    transcript_has_low_coverage_exon = logical()
  )
}

measure_aggregate_cds_exon_qc <- function(cds_all) {
  aggregate_df <- tryCatch(
    read.experiment("human_all_merged_l50", validate = FALSE),
    error = function(e) {
      warning("Could not load human_all_merged_l50 for CDS exon QC: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(aggregate_df)) return(empty_cds_exon_coverage_qc())

  aggregate_reads <- tryCatch(filepath(aggregate_df, "bigwig"),
                              error = function(e) NULL)
  aggregate_paths <- unlist(aggregate_reads, use.names = FALSE)
  if (length(aggregate_paths) == 0 || !all(file.exists(aggregate_paths))) {
    warning("Skipping CDS exon QC because human_all_merged_l50 bigWig paths ",
            "are unavailable.")
    return(empty_cds_exon_coverage_qc())
  }

  tryCatch(
    cdsExonCoverageQC(
      cds_all,
      aggregate_paths,
      low_exon_fraction = 0.1,
      min_transcript_counts = 10
    ),
    error = function(e) {
      warning("Could not compute aggregate CDS exon QC: ", conditionMessage(e))
      empty_cds_exon_coverage_qc()
    }
  )
}

summarize_transcript_exon_qc <- function(exon_qc) {
  if (nrow(exon_qc) == 0) {
    return(data.table(
      tx_id = character(),
      cds_exon_qc_status = character(),
      cds_exon_qc_measurable = logical(),
      cds_exon_qc_transcript_counts = numeric(),
      cds_exon_qc_n_exons = integer(),
      cds_exon_qc_min_ratio = numeric(),
      cds_exon_qc_low_exon_ranks = character()
    ))
  }
  exon_qc[
    ,
    .(
      cds_exon_qc_status = fifelse(
        any(transcript_has_low_coverage_exon, na.rm = TRUE),
        "low_aggregate_cds_exon",
        fifelse(any(exon_qc_measurable, na.rm = TRUE),
                "pass", "not_measurable")
      ),
      cds_exon_qc_measurable = any(exon_qc_measurable, na.rm = TRUE),
      cds_exon_qc_transcript_counts = max(transcript_counts, na.rm = TRUE),
      cds_exon_qc_n_exons = max(n_exons, na.rm = TRUE),
      cds_exon_qc_min_ratio = if (all(!is.finite(exon_vs_transcript_mean))) {
        NA_real_
      } else {
        min(exon_vs_transcript_mean, na.rm = TRUE)
      },
      cds_exon_qc_low_exon_ranks = {
        ranks <- exon_rank[low_coverage_exon == TRUE]
        if (length(ranks) == 0) NA_character_ else paste(ranks, collapse = ";")
      }
    ),
    by = .(tx_id = transcript_id)
  ]
}

count_table_gene_matrix <- function(df, symbols_dt, target_genes) {
  counts <- countTable(df, "cds", "fpkm")
  count_dt <- as.data.table(counts)
  count_matrix <- as.matrix(count_dt)
  storage.mode(count_matrix) <- "numeric"
  if (ncol(count_matrix) != nrow(df)) {
    stop("Fallback count table columns (", ncol(count_matrix),
         ") do not match experiment rows (", nrow(df), ").")
  }
  colnames(count_matrix) <- runIDs(df)

  feature_id <- rownames(counts)
  if (is.null(feature_id) || length(feature_id) != nrow(count_matrix)) {
    stop("Fallback count table has no usable transcript row names.")
  }

  symbol_map <- unique(as.data.table(symbols_dt)[
    toupper(external_gene_name) %chin% target_genes,
    .(
      feature_id = as.character(ensembl_tx_name),
      feature_id_no_version = strip_version(ensembl_tx_name),
      gene_symbol = toupper(external_gene_name)
    )
  ])

  feature_map <- data.table(
    row_index = seq_along(feature_id),
    feature_id = as.character(feature_id),
    feature_id_no_version = strip_version(feature_id)
  )
  feature_map[, gene_symbol := symbol_map$gene_symbol[
    match(feature_id, symbol_map$feature_id)
  ]]
  missing_symbol <- is.na(feature_map$gene_symbol)
  feature_map[missing_symbol, gene_symbol := symbol_map$gene_symbol[
    match(feature_id_no_version, symbol_map$feature_id_no_version)
  ]]
  feature_map <- feature_map[!is.na(gene_symbol)]

  out <- matrix(NA_real_, nrow = length(target_genes), ncol = ncol(count_matrix),
                dimnames = list(target_genes, runIDs(df)))
  for (gene in target_genes) {
    rows <- feature_map[gene_symbol == gene, row_index]
    if (length(rows) == 0) next
    out[gene, ] <- colMeans(count_matrix[rows, , drop = FALSE], na.rm = TRUE)
  }
  out
}

choose_gene_transcript <- function(gene_symbol, tx_ids, cds_all, leader_all,
                                   translons, transcript_exon_qc) {
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

    target_tx_id <- tx_id
    exon_qc <- transcript_exon_qc[tx_id == target_tx_id]
    if (nrow(exon_qc) == 0) {
      exon_qc <- data.table(
        cds_exon_qc_status = "not_checked",
        cds_exon_qc_measurable = FALSE,
        cds_exon_qc_transcript_counts = NA_real_,
        cds_exon_qc_n_exons = length(cds_tx[[1]]),
        cds_exon_qc_min_ratio = NA_real_,
        cds_exon_qc_low_exon_ranks = NA_character_
      )
    }
    raw_canonical <- strip_version(tx_id) %chin% canonical_tx
    canonical_qc_priority <- raw_canonical &&
      !identical(exon_qc$cds_exon_qc_status[[1]], "low_aggregate_cds_exon")

    cbind(data.table(
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
      is_canonical_isoform = raw_canonical,
      canonical_isoform_qc_priority = canonical_qc_priority
    ), exon_qc[, setdiff(names(exon_qc), "tx_id"), with = FALSE])
  }), fill = TRUE)

  if (nrow(candidates) == 0) return(NULL)
  candidates[, has_overlapping_uorf := n_overlapping_uorfs > 0]
  candidates[, has_clean_cds := clean_cds_bases > 0]
  candidates[, usable_cds_exon_qc_candidate :=
               cds_exon_qc_status != "low_aggregate_cds_exon"]
  setorder(candidates, -has_clean_cds, -usable_cds_exon_qc_candidate,
           -has_overlapping_uorf,
           -canonical_isoform_qc_priority,
           -clean_cds_bases, -leader_bases, tx_id)
  candidates[, has_overlapping_uorf := NULL]
  candidates[, has_clean_cds := NULL]
  candidates[, usable_cds_exon_qc_candidate := NULL]
  candidates[1]
}

measure_gene_clean_cds <- function(selected, cds_all, leader_all, translons,
                                   fst_index, lib_sizes) {
  tx_id <- selected$tx_id
  cds_tx <- cds_all[tx_id]
  leader_tx <- if (tx_id %in% names(leader_all)) leader_all[tx_id] else GRangesList()
  display_region <- make_display_region(leader_tx, cds_tx, tx_id)
  cds_positions <- safe_pmap_positions(cds_tx, display_region)
  excluded <- overlapping_uorf_positions(
    translons, tx_id, leader_tx, cds_tx, display_region, cds_positions
  )
  clean_positions <- setdiff(cds_positions, excluded$positions)

  if (length(clean_positions) == 0) {
    return(rep(NA_real_, length(lib_sizes)))
  }

  coverage <- tryCatch(
    dominant_cached_coverage_by_transcript(
      display_region = display_region,
      fst_index = fst_index,
      gene_symbol = selected$gene_symbol,
      tx_id = tx_id,
      analysis_dir = analysis_dir
    ),
    error = function(e) {
      warning("Could not read FST coverage for ", selected$gene_symbol,
              " (", tx_id, "): ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(coverage)) {
    out <- rep(NA_real_, length(lib_sizes))
    names(out) <- names(lib_sizes)
    return(out)
  }
  coverage_matrix <- as.matrix(coverage)
  storage.mode(coverage_matrix) <- "numeric"

  run_order <- colnames(coverage_matrix)
  if (!all(run_order %in% names(lib_sizes))) {
    stop("Coverage FST columns are not all present in library size vector.")
  }

  clean_sum <- colSums(coverage_matrix[clean_positions, , drop = FALSE],
                       na.rm = TRUE)
  clean_kb <- length(clean_positions) / 1000
  lib_millions <- as.numeric(lib_sizes[run_order]) / 1e6
  fpkm_like <- (clean_sum / clean_kb) / lib_millions
  names(fpkm_like) <- run_order
  fpkm_like
}

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
symbols <- load_symbol_table(df)
translons <- load_predicted_translons(df)
canonical_info <- get_canonical_isoform_info(df)
canonical_tx <- canonical_info$tx
fst_index <- file.path(collection_dir_from_exp(df), "coverage_index.fst")
if (!file.exists(fst_index)) stop("Missing FST coverage index: ", fst_index)
lib_sizes <- library_sizes_by_run(df)

if (length(canonical_tx) > 0) {
  message("Canonical isoform prioritization enabled for ",
          length(canonical_tx), " transcript IDs.")
} else {
  message("Canonical isoform prioritization unavailable. ",
          "Fallback transcript ranking is active.")
}

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
cds_exon_qc <- measure_aggregate_cds_exon_qc(cds_all)
fwrite(cds_exon_qc, output_cds_exon_qc_file)
transcript_exon_qc <- summarize_transcript_exon_qc(cds_exon_qc)
message("Aggregate CDS exon QC measured ", nrow(cds_exon_qc),
        " exons across ", uniqueN(cds_exon_qc$transcript_id),
        " candidate transcripts; flagged ",
        transcript_exon_qc[cds_exon_qc_status == "low_aggregate_cds_exon", .N],
        " transcript(s) with <10% exon support.")
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
  choose_gene_transcript(gene, tx_ids, cds_all, leader_all, translons,
                         transcript_exon_qc)
}), fill = TRUE)

module_membership <- rbindlist(lapply(seq_len(nrow(dominant_state_modules)), function(i) {
  module <- dominant_state_modules[i]
  rbindlist(list(
    data.table(
      gene_symbol = toupper(module$up_genes[[1]]),
      dominant_state_membership = module$dominant_state,
      model_direction = if (module$dominant_state == "Baseline (control)") {
        "baseline"
      } else {
        "up"
      }
    ),
    data.table(
      gene_symbol = toupper(module$down_genes[[1]]),
      dominant_state_membership = module$dominant_state,
      model_direction = "down"
    )
  ))
}), fill = TRUE)
module_membership <- module_membership[
  ,
  .(
    dominant_states = paste(unique(dominant_state_membership), collapse = "; "),
    model_directions = paste(unique(model_direction), collapse = "; ")
  ),
  by = gene_symbol
]
fatigue_membership <- load_postviral_fatigue_gene_set(analysis_dir)[
  ,
  .(
    postviral_modules = paste(unique(Module), collapse = "; "),
    postviral_tiers = paste(unique(Tier), collapse = "; "),
    postviral_directions = paste(unique(Direction), collapse = "; ")
  ),
  by = .(gene_symbol = Gene)
]
selected_transcripts <- merge(
  selected_transcripts,
  module_membership,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
selected_transcripts <- merge(
  selected_transcripts,
  fatigue_membership,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
selected_transcripts[, selection_rank := seq_len(.N)]
selected_transcripts[, canonical_isoform_source := canonical_info$source]
selected_transcripts[, canonical_isoform_note := canonical_info$note]

message("Selected ", sum(!is.na(selected_transcripts$tx_id)),
        " transcripts for ", nrow(selected_transcripts), " marker genes.")

expr_list <- lapply(seq_len(nrow(selected_transcripts)), function(i) {
  selected <- selected_transcripts[i]
  if (is.na(selected$tx_id)) {
    out <- rep(NA_real_, length(lib_sizes))
    names(out) <- names(lib_sizes)
    return(out)
  }
  message("Measuring ", selected$gene_symbol, " (", selected$tx_id, ")")
  measure_gene_clean_cds(selected, cds_all, leader_all, translons,
                         fst_index, lib_sizes)
})

gene_expr <- do.call(rbind, expr_list)
rownames(gene_expr) <- selected_transcripts$gene_symbol

run_order <- colnames(gene_expr)
if (is.null(run_order)) run_order <- names(expr_list[[1]])
if (!identical(run_order, runIDs(df))) {
  gene_expr <- gene_expr[, runIDs(df), drop = FALSE]
}
stopifnot(identical(colnames(gene_expr), runIDs(df)))

fst_measured <- rowSums(is.finite(gene_expr)) > 0
zero_clean_cds <- selected_transcripts$clean_cds_bases == 0
selected_transcripts[, expression_source := fifelse(
  fst_measured,
  "fst_clean_cds",
  fifelse(zero_clean_cds, "zero_clean_cds_after_uorf_exclusion", "missing_fst")
)]

zero_clean_genes <- selected_transcripts[
  expression_source == "zero_clean_cds_after_uorf_exclusion",
  gene_symbol
]
if (length(zero_clean_genes) > 0) {
  message("No non-uORF-overlapping CDS bases remain for ",
          length(zero_clean_genes), " marker genes: ",
          paste(zero_clean_genes, collapse = ", "),
          ". These genes are kept as NA and omitted from module means.")
}

missing_fst_genes <- selected_transcripts[expression_source == "missing_fst", gene_symbol]
if (length(missing_fst_genes) > 0) {
  message("FST coverage unavailable for ", length(missing_fst_genes),
          " marker genes: ", paste(missing_fst_genes, collapse = ", "))
  if (!allow_count_table_fallback) {
    stop("Pure FST scoring cannot continue without complete marker coverage.")
  }

  message("Filling unavailable marker genes from countTable(df, 'cds', 'fpkm') fallback.")
  fallback_expr <- count_table_gene_matrix(df, symbols, target_genes)
  gene_expr[missing_fst_genes, ] <- fallback_expr[missing_fst_genes, , drop = FALSE]
  selected_transcripts[
    gene_symbol %chin% missing_fst_genes,
    expression_source := "count_table_fpkm_fallback"
  ]
}

expr_norm <- log2(gene_expr + 1)
baseline_genes <- dominant_state_modules[
  dominant_state == "Baseline (control)", up_genes[[1]]
]
baseline_present <- intersect(baseline_genes, rownames(expr_norm))
if (length(baseline_present) < 5) {
  stop("Too few baseline genes were found in the FST-derived matrix: ",
       length(baseline_present), " / ", length(baseline_genes))
}

baseline <- colMeans(expr_norm[baseline_present, , drop = FALSE], na.rm = TRUE)
gene_z <- sweep(expr_norm, 2, baseline, "-")

sample_runs <- data.table(
  Run = runIDs(df),
  sample = get_sample_names(df)
)

dominant_state_scores <- rbindlist(lapply(seq_len(nrow(dominant_state_modules)), function(i) {
  state <- dominant_state_modules[i]
  score <- score_module(gene_z, state$up_genes[[1]], state$down_genes[[1]])
  data.table(
    Run = runIDs(df),
    sample = sample_runs$sample,
    dominant_state = state$dominant_state,
    score = as.numeric(score),
    up_genes_present = length(intersect(state$up_genes[[1]], rownames(gene_z))),
    up_genes_total = length(state$up_genes[[1]]),
    down_genes_present = length(intersect(state$down_genes[[1]], rownames(gene_z))),
    down_genes_total = length(state$down_genes[[1]])
  )
}))

dominant_state_qc_annotations <- dominant_state_scores[
  dominant_state != "Baseline (control)",
  .SD[which.max(score)],
  by = Run
][
  ,
  .(Run, sample, dominant_state, dominant_state_score = score)
]

dominant_state_scores_wide <- dcast(
  dominant_state_scores[, .(Run, dominant_state, score)],
  Run ~ dominant_state,
  value.var = "score"
)

dominant_state_qc_annotations <- merge(
  dominant_state_qc_annotations,
  dominant_state_scores_wide,
  by = "Run",
  sort = FALSE
)
dominant_state_qc_annotations <- merge(
  sample_runs[, .(Run, run_order = .I)],
  dominant_state_qc_annotations,
  by = "Run",
  sort = FALSE
)
setorder(dominant_state_qc_annotations, run_order)
dominant_state_qc_annotations[, run_order := NULL]
setcolorder(
  dominant_state_qc_annotations,
  c("Run", setdiff(names(dominant_state_qc_annotations), "Run"))
)

marker_coverage <- dominant_state_modules[
  ,
  .(
    dominant_state,
    up_genes = gene_list_text(up_genes),
    down_genes = gene_list_text(down_genes),
    up_genes_present = vapply(up_genes, function(x) {
      length(intersect(x, rownames(gene_z)))
    }, integer(1)),
    up_genes_total = lengths(up_genes),
    down_genes_present = vapply(down_genes, function(x) {
      length(intersect(x, rownames(gene_z)))
    }, integer(1)),
    down_genes_total = lengths(down_genes)
  )
]

gene_expression_dt <- data.table(gene_symbol = rownames(gene_expr), gene_expr)
fwrite(dominant_state_qc_annotations, output_state_file)
fwrite(gene_expression_dt, output_expr_file)
fwrite(selected_transcripts, output_diagnostics_file)

stopifnot(nrow(dominant_state_qc_annotations) == nrow(df))
stopifnot(identical(dominant_state_qc_annotations$Run, runIDs(df)))
stopifnot(all(marker_coverage[
  dominant_state != "Baseline (control)", up_genes_present
] > 0))
stopifnot(all(marker_coverage[
  down_genes_total > 0, down_genes_present
] > 0))

message("\nMarker definition / coverage table:")
print(marker_coverage)

message("\nSelected transcript diagnostics:")
print(selected_transcripts)

message("\nDominant-state QC annotation result table:")
print(dominant_state_qc_annotations)

message("\nSaved: ", output_state_file)
message("Saved: ", output_expr_file)
message("Saved: ", output_diagnostics_file)

invisible(list(
  modules = dominant_state_modules,
  marker_coverage = marker_coverage,
  selected_transcripts = selected_transcripts,
  gene_expression = gene_expr,
  expr_norm = expr_norm,
  gene_z = gene_z,
  scores = dominant_state_scores,
  qc_annotations = dominant_state_qc_annotations
))
