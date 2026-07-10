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
  library(ggplot2)
  library(grid)
})

analysis_dir <- if (file.exists("dominant_cell_states/dominant_uorf_feature_expression.csv")) {
  "dominant_cell_states"
} else {
  "."
}

feature_expression_file <- file.path(analysis_dir, "dominant_uorf_feature_expression.csv")
state_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds.csv")
diagnostics_file <- file.path(analysis_dir, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv")
enrichment_file <- file.path(analysis_dir, "dominant_state_clean_cds_enrichment_all_main.csv")
broad_rdg_terms_file <- file.path(analysis_dir, "dominant_state_broad_metadata_rdg_terms.csv")
glmnet_rdg_terms_file <- file.path(analysis_dir, "dominant_state_glmnet_metadata_rdg_terms.csv")
selected_fields_file <- file.path(analysis_dir, "dominant_state_broad_metadata_selected_fields.csv")
metadata_file <- "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv"
output_dir <- file.path(analysis_dir, "dominant_rdg_outputs")

metadata_fields <- c("CELL_LINE", "TISSUE", "GENE", "CONDITION")
min_branch_samples <- 2L
min_enrichment_branch_samples <- 3L
max_enrichment_branch_terms <- Inf
max_enrichment_branch_terms_per_state <- 4L
max_genes_to_plot <- Inf
max_branch_panels_per_gene <- 5L
min_clean_cds_raw_counts <- 30
max_matrix_plot_rows <- Inf
rdg_plot_mode <- if (exists("dominant_rdg_plot_mode", inherits = TRUE)) {
  get("dominant_rdg_plot_mode", inherits = TRUE)
} else {
  Sys.getenv("DOMINANT_RDG_PLOT_MODE", "all")
}
rdg_plot_mode <- tolower(trimws(as.character(rdg_plot_mode)[1]))
if (!nzchar(rdg_plot_mode) || is.na(rdg_plot_mode)) rdg_plot_mode <- "all"
if (rdg_plot_mode %in% c("full", "plots", "plot")) rdg_plot_mode <- "all"
if (rdg_plot_mode %in% c("table", "none", "skip", "skip_plots",
                         "no_plots", "no-plot", "no-plots")) {
  rdg_plot_mode <- "tables"
}
if (!rdg_plot_mode %in% c("all", "tables")) {
  stop("Unknown RDG plot mode: ", rdg_plot_mode,
       ". Use 'all' or 'tables'.")
}
generate_rdg_plots <- identical(rdg_plot_mode, "all")

message("Dominant-state ribosome decision graph prototype:")
message("  1. Build transcript-level RDG topology from predicted translons + clean CDS")
message("  2. Use dominant states and significant metadata enrichments as branch overlays")
message("  3. Estimate relative translon use and approximate branch probabilities")
message("  4. Save RDG node/edge tables, branch score tables, and ggplot figures")
message("  RDG plot mode: ", rdg_plot_mode)

if (file.exists(selected_fields_file)) {
  selected_profile <- fread(selected_fields_file)
  if (all(c("field", "used_for_broad_model") %chin% names(selected_profile))) {
    metadata_fields <- unique(c(
      metadata_fields,
      selected_profile[used_for_broad_model == TRUE, field]
    ))
  }
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures_by_state"), showWarnings = FALSE,
           recursive = TRUE)

clean_level <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none")] <- NA_character_
  x
}

file_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  substr(x, 1L, 160L)
}

single_field_level <- function(field, level) {
  level <- as.character(level)
  prefix <- paste0(field, "=")
  has_prefix <- startsWith(level, prefix)
  ifelse(has_prefix, substring(level, nchar(prefix) + 1L), level)
}

parse_conjunction_level <- function(field, level) {
  fields <- strsplit(as.character(field), ":", fixed = TRUE)[[1]]
  parts <- strsplit(as.character(level), " | ", fixed = TRUE)[[1]]
  values <- rep(NA_character_, length(fields))
  names(values) <- fields
  for (part in parts) {
    eq <- regexpr("=", part, fixed = TRUE)[1]
    if (!is.finite(eq) || eq < 1L) next
    key <- substr(part, 1L, eq - 1L)
    value <- substr(part, eq + 1L, nchar(part))
    if (key %chin% fields) values[[key]] <- value
  }
  values
}

branch_label_text <- function(branch_source, field, level, dominant_state) {
  if (branch_source == "sample_dominant_call") {
    return(paste0("dominant: ", dominant_state))
  }
  branch_text <- if (grepl(":", field, fixed = TRUE)) {
    gsub(" | ", "\n", level, fixed = TRUE)
  } else {
    paste0(field, "=", single_field_level(field, level))
  }
  paste0(branch_text, "\n", dominant_state)
}

short_state_label <- function(x) {
  x <- as.character(x)
  x <- sub("Post-viral fatigue / ribosome stress", "Post-viral fatigue", x, fixed = TRUE)
  x <- sub(" / Cell cycle$", "", x)
  x <- sub(" / HIF1A program$", "", x)
  x <- sub(" / antiviral$", "", x)
  x <- sub(" / Mesenchymal shift$", "", x)
  x <- sub(" \\(ATF4 axis\\)$", "", x)
  x <- sub(" \\(mitochondrial\\)$", "", x)
  x <- sub("mTOR / Translation capacity \\(TOP program\\)", "mTOR/TOP", x)
  x
}

truncate_label <- function(x, max_chars = 46L) {
  x <- as.character(x)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 3L), "...")
  x
}

compact_branch_split_label <- function(field, level) {
  field <- as.character(field)
  level <- as.character(level)
  if (grepl(":", field, fixed = TRUE)) {
    values <- parse_conjunction_level(field, level)
    keep <- names(values)[!is.na(values) & nzchar(values)]
    keep <- keep[seq_len(min(length(keep), 3L))]
    parts <- paste0(keep, "=", truncate_label(values[keep], 14L))
    extra <- sum(!is.na(values) & nzchar(values)) - length(keep)
    if (extra > 0L) parts <- c(parts, paste0("+", extra, " more"))
    return(paste(parts, collapse = " | "))
  }
  paste0(field, "=", truncate_label(single_field_level(field, level), 30L))
}

sum_or_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
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
  mcols(translons)$translon_source_path <- path
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
    file.path(analysis_dir, "manual_translons", "manual_translons_ranges.rds"),
    file.path("dominant_cell_states", "manual_translons", "manual_translons_ranges.rds"),
    file.path("manual_translons", "manual_translons_ranges.rds")
  )
  t <- load_translon_source(t_paths, "T")
  tc <- load_translon_source(tc_paths, "TC")
  m <- load_translon_source(m_paths, "M")
  if (length(t) == 0 && length(tc) == 0 && length(m) == 0) {
    stop("Predicted translon range files not found. Checked: ",
         paste(unique(c(
           t_paths,
           tc_paths,
           m_paths
         )), collapse = "; "))
  }
  message("Loaded translon sources: T=", length(t), ", TC=", length(tc),
          ", M=", length(m))
  c(t, tc, m)
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

annotate_translon_categories <- function(translons, cds, mrna) {
  if (length(translons) == 0) return(translons)

  keep <- names(translons) %chin% names(mrna)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " predicted translons without mRNA ranges.")
    translons <- translons[keep]
  }
  if (length(translons) == 0) return(translons)

  mapped <- suppressWarnings(pmapToTranscriptF(translons, mrna[names(translons)]))
  tx_iranges <- collapse_to_tx_iranges(mapped)
  valid <- !is.na(start(tx_iranges)) & !is.na(end(tx_iranges))
  categories <- rep(NA_character_, length(translons))

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

is_model_translon_category <- function(category) {
  as.character(category) %chin% c("uORF", "uoORF", "NTE", "NTT", "internal")
}

feature_prefix_from_category <- function(category) {
  category <- as.character(category)
  if (category == "uORF") return("uORF")
  if (category == "uoORF") return("uORF")
  if (category == "NTE") return("NTE")
  if (category == "NTT") return("NTT")
  if (category == "internal") return("iORF")
  "translon"
}

is_cds_allocation_feature_type <- function(feature_type) {
  as.character(feature_type) %chin% c(
    "leader_uorf",
    "overlapping_uorf",
    "clean_CDS"
  )
}

classify_translon <- function(translon_gr, leader_gr, cds_gr) {
  overlaps_leader <- length(findOverlaps(translon_gr, leader_gr,
                                         ignore.strand = FALSE)) > 0
  overlaps_cds <- length(findOverlaps(translon_gr, cds_gr,
                                      ignore.strand = FALSE)) > 0
  if (overlaps_leader && overlaps_cds) return("overlapping_uorf")
  if (overlaps_leader) return("leader_uorf")
  if (overlaps_cds) return("cds_translon")
  "other_translon"
}

classify_feature <- function(translon_gr, leader_gr, cds_gr, positions,
                             cds_positions, category) {
  feature_type <- classify_translon(translon_gr, leader_gr, cds_gr)
  if (category == "uORF") feature_type <- "leader_uorf"
  if (category == "uoORF") feature_type <- "overlapping_uorf"
  if (category == "NTE") feature_type <- "nte_extension"
  if (category == "NTT") feature_type <- "ntt_truncation"
  if (category == "internal") feature_type <- "internal_orf"
  cds_overlap_bases <- length(intersect(cds_positions, positions))
  list(feature_type = feature_type, cds_overlap_bases = cds_overlap_bases)
}

annotate_feature_sequence <- function(info, display, cds_positions, df) {
  if (nrow(info) == 0 || !"positions" %in% names(info)) return(info)
  display_seq <- tryCatch(
    as.character(txSeqsFromFa(display, df, TRUE)[[1]]),
    error = function(e) NA_character_
  )
  if (is.na(display_seq) || !nzchar(display_seq)) {
    info[, `:=`(
      start_codon = NA_character_,
      stop_codon = NA_character_,
      start_codon_class = NA_character_,
      frame_relative_to_cds = NA_integer_,
      peptide_length_aa = NA_integer_,
      kozak_minus3 = NA_character_,
      kozak_plus4 = NA_character_,
      kozak_strength = NA_character_,
      coding_potential_proxy = NA
    )]
    return(info)
  }

  seq_chars <- strsplit(display_seq, "", fixed = TRUE)[[1]]
  cds_start <- if (length(cds_positions) > 0) min(cds_positions) else NA_integer_
  stop_codons <- c("TAA", "TAG", "TGA")
  near_cognate <- c("CTG", "GTG", "TTG", "ACG", "ATT", "ATC", "ATA")

  seq_info <- rbindlist(lapply(seq_len(nrow(info)), function(i) {
    positions <- info$positions[[i]]
    positions <- positions[is.finite(positions) & positions >= 1L &
                             positions <= length(seq_chars)]
    if (length(positions) == 0) {
      return(data.table(
        start_codon = NA_character_,
        stop_codon = NA_character_,
        start_codon_class = NA_character_,
        frame_relative_to_cds = NA_integer_,
        peptide_length_aa = NA_integer_,
        kozak_minus3 = NA_character_,
        kozak_plus4 = NA_character_,
        kozak_strength = NA_character_,
        coding_potential_proxy = NA
      ))
    }
    feature_seq <- paste(seq_chars[positions], collapse = "")
    feature_len <- nchar(feature_seq)
    tx_start <- min(positions)
    start_codon <- if (feature_len >= 3L) substr(feature_seq, 1L, 3L) else NA_character_
    stop_codon <- if (feature_len >= 3L) {
      substr(feature_seq, feature_len - 2L, feature_len)
    } else {
      NA_character_
    }
    start_class <- if (is.na(start_codon)) {
      NA_character_
    } else if (start_codon == "ATG") {
      "ATG"
    } else if (start_codon %chin% near_cognate) {
      "near_cognate"
    } else {
      "other"
    }
    minus3 <- if (tx_start > 3L) seq_chars[[tx_start - 3L]] else NA_character_
    plus4 <- if (tx_start + 3L <= length(seq_chars)) {
      seq_chars[[tx_start + 3L]]
    } else {
      NA_character_
    }
    kozak_strength <- if (is.na(minus3) || is.na(plus4)) {
      NA_character_
    } else if (minus3 %chin% c("A", "G") && plus4 == "G") {
      "strong"
    } else if (minus3 %chin% c("A", "G") || plus4 == "G") {
      "moderate"
    } else {
      "weak"
    }
    has_stop <- !is.na(stop_codon) && stop_codon %chin% stop_codons
    peptide_length <- if (feature_len >= 3L) {
      floor(feature_len / 3L) - as.integer(has_stop)
    } else {
      NA_integer_
    }
    data.table(
      start_codon = start_codon,
      stop_codon = stop_codon,
      start_codon_class = start_class,
      frame_relative_to_cds = if (is.finite(cds_start)) {
        as.integer((tx_start - cds_start) %% 3L)
      } else {
        NA_integer_
      },
      peptide_length_aa = as.integer(peptide_length),
      kozak_minus3 = minus3,
      kozak_plus4 = plus4,
      kozak_strength = kozak_strength,
      coding_potential_proxy = start_class %chin% c("ATG", "near_cognate") &&
        has_stop && is.finite(peptide_length) && peptide_length >= 2L
    )
  }), fill = TRUE)
  cbind(info, seq_info)
}

assign_track_y <- function(feature_type, feature_order) {
  ifelse(
    feature_type == "clean_CDS", -0.55,
    ifelse(
      feature_type == "overlapping_uorf", 0.35,
      ifelse(
        feature_type == "nte_extension", -0.95,
        ifelse(
          feature_type == "ntt_truncation", -1.18,
          ifelse(
            feature_type == "internal_orf", -0.82,
            ifelse(
              feature_type == "cds_translon", -0.95,
              0.85 + ((feature_order - 1L) %% 3L) * 0.28
            )
          )
        )
      )
    )
  )
}

build_feature_annotation_one <- function(df, translons, gene, tx_id) {
  cds_tx <- loadRegion(df, "cds", names.keep = tx_id)
  if (length(cds_tx) == 0 || !tx_id %in% names(cds_tx)) return(NULL)
  leader_tx <- loadRegion(df, "leaders", names.keep = tx_id)
  if (!tx_id %in% names(leader_tx)) leader_tx <- GRangesList()

  display <- make_display_region(leader_tx, cds_tx, tx_id)
  cds_positions <- safe_pmap_positions(cds_tx, display)
  leader_gr <- unlist(leader_tx, use.names = FALSE)
  cds_gr <- unlist(cds_tx, use.names = FALSE)
  tx_translons <- translons[names(translons) == tx_id]
  if (length(tx_translons) > 0 && !"category" %in% names(mcols(tx_translons))) {
    stop("Predicted translons must be categorized before RDG annotation.")
  }

  uorf_info <- rbindlist(lapply(seq_along(tx_translons), function(i) {
    category <- as.character(mcols(tx_translons)$category[i])
    if (is.na(category)) {
      stop("Could not categorize predicted translon ", i, " for ", tx_id)
    }
    if (!is_model_translon_category(category)) return(NULL)

    translon_gr <- tx_translons[[i]]
    positions <- safe_pmap_positions(GRangesList(translon_gr), display)
    if (length(positions) == 0) return(NULL)
    feature_class <- classify_feature(
      translon_gr, leader_gr, cds_gr, positions, cds_positions, category
    )
    feature_prefix <- feature_prefix_from_category(category)
    data.table(
      feature_id = paste0(feature_prefix, i),
      feature_prefix = feature_prefix,
      feature_type = feature_class$feature_type,
      category = category,
      translon_source = as.character(mcols(tx_translons)$translon_source[i]),
      translon_sources_merged = as.character(mcols(tx_translons)$translon_sources_merged[i]),
      manual_start_codon = if ("manual_start_codon" %chin% names(mcols(tx_translons))) {
        as.character(mcols(tx_translons)$manual_start_codon[i])
      } else {
        NA_character_
      },
      manual_stop_codon = if ("manual_stop_codon" %chin% names(mcols(tx_translons))) {
        as.character(mcols(tx_translons)$manual_stop_codon[i])
      } else {
        NA_character_
      },
      tx_start = min(positions),
      tx_end = max(positions),
      feature_bases = length(positions),
      cds_overlap_bases = feature_class$cds_overlap_bases,
      positions = list(positions)
    )
  }), fill = TRUE)

  if (nrow(uorf_info) > 0) {
    setorder(uorf_info, tx_start, tx_end)
    uorf_info[, feature_seq := seq_len(.N), by = feature_prefix]
    uorf_info[, feature_id := paste0(feature_prefix, feature_seq)]
    uorf_info[, c("feature_prefix", "feature_seq") := NULL]
  }

  overlapping_positions <- if (nrow(uorf_info) > 0) {
    unique(unlist(uorf_info[feature_type == "overlapping_uorf", positions],
                  use.names = FALSE))
  } else {
    integer()
  }
  clean_cds_positions <- setdiff(cds_positions, intersect(cds_positions, overlapping_positions))
  clean_cds <- data.table(
    feature_id = "clean_CDS",
    feature_type = "clean_CDS",
    category = "clean_CDS",
    tx_start = if (length(clean_cds_positions) > 0) min(clean_cds_positions) else min(cds_positions),
    tx_end = if (length(clean_cds_positions) > 0) max(clean_cds_positions) else max(cds_positions),
    feature_bases = length(clean_cds_positions),
    positions = list(clean_cds_positions)
  )

  info <- rbind(
    if (nrow(uorf_info) > 0) copy(uorf_info) else NULL,
    clean_cds,
    fill = TRUE
  )
  info <- annotate_feature_sequence(info, display, cds_positions, df)
  manual_start <- "manual_start_codon" %chin% names(info) &
    !is.na(info$manual_start_codon) & nzchar(info$manual_start_codon)
  manual_stop <- "manual_stop_codon" %chin% names(info) &
    !is.na(info$manual_stop_codon) & nzchar(info$manual_stop_codon)
  if (any(manual_start)) {
    info[manual_start, start_codon := manual_start_codon]
  }
  if (any(manual_stop)) {
    info[manual_stop, stop_codon := manual_stop_codon]
  }
  near_cognate <- c("CTG", "GTG", "TTG", "ACG", "ATT", "ATC", "ATA")
  stop_codons <- c("TAA", "TAG", "TGA")
  info[, start_codon_class := fifelse(
    is.na(start_codon),
    NA_character_,
    fifelse(
      start_codon == "ATG",
      "ATG",
      fifelse(start_codon %chin% near_cognate, "near_cognate", "other")
    )
  )]
  info[, peptide_length_aa := fifelse(
    is.finite(feature_bases) & feature_bases >= 3L,
    as.integer(floor(feature_bases / 3L) -
                 as.integer(!is.na(stop_codon) & stop_codon %chin% stop_codons)),
    NA_integer_
  )]
  info[, coding_potential_proxy := start_codon_class %chin% c("ATG", "near_cognate") &
         !is.na(stop_codon) & stop_codon %chin% stop_codons &
         is.finite(peptide_length_aa) & peptide_length_aa >= 2L]
  info[, positions := NULL]
  info[, clean_order := fifelse(feature_id == "clean_CDS", 1L, 0L)]
  setorder(info, tx_start, clean_order, tx_end)
  info[, clean_order := NULL]
  info[, `:=`(
    gene_symbol = gene,
    tx_id = tx_id,
    transcript_length = max(c(tx_end, cds_positions), na.rm = TRUE),
    feature_order = seq_len(.N)
  )]
  info[, track_y := assign_track_y(feature_type, feature_order)]
  info[]
}

build_feature_annotation <- function(df, translons, feature_dt) {
  tx_by_gene <- unique(feature_dt[, .(gene_symbol, tx_id)])
  rbindlist(lapply(seq_len(nrow(tx_by_gene)), function(i) {
    row <- tx_by_gene[i]
    message("Building RDG feature annotation for ", row$gene_symbol,
            " (", row$tx_id, ")")
    build_feature_annotation_one(df, translons, row$gene_symbol, row$tx_id)
  }), fill = TRUE)
}

build_rdg_topology_one <- function(ann) {
  gene <- ann$gene_symbol[[1]]
  tx_id <- ann$tx_id[[1]]
  transcript_length <- max(ann$transcript_length, na.rm = TRUE)
  start_nodes <- ann[, .(
    node_id = paste(gene, feature_id, "start", sep = "|"),
    gene_symbol,
    tx_id,
    feature_id,
    feature_type,
    node_type = fifelse(feature_id == "clean_CDS", "cds_start", "branch_start"),
    x = tx_start,
    y = 0
  )]
  stop_nodes <- ann[, .(
    node_id = paste(gene, feature_id, "stop", sep = "|"),
    gene_symbol,
    tx_id,
    feature_id,
    feature_type,
    node_type = "termination",
    x = tx_end,
    y = track_y
  )]
  nodes <- rbind(
    data.table(
      node_id = paste(gene, "cap", sep = "|"),
      gene_symbol = gene,
      tx_id = tx_id,
      feature_id = NA_character_,
      feature_type = "cap",
      node_type = "cap",
      x = 1,
      y = 0
    ),
    start_nodes,
    stop_nodes,
    data.table(
      node_id = paste(gene, "end", sep = "|"),
      gene_symbol = gene,
      tx_id = tx_id,
      feature_id = NA_character_,
      feature_type = "end",
      node_type = "transcript_end",
      x = transcript_length,
      y = 0
    ),
    fill = TRUE
  )

  branch_x <- sort(unique(c(1L, ann$tx_start, transcript_length)))
  scan_edges <- rbindlist(lapply(seq_len(length(branch_x) - 1L), function(i) {
    data.table(
      edge_id = paste(gene, "scan", i, sep = "|"),
      gene_symbol = gene,
      tx_id = tx_id,
      source = if (i == 1L) paste(gene, "cap", sep = "|") else paste(gene, branch_x[i], "scan", sep = "|"),
      target = if (i == length(branch_x) - 1L) paste(gene, "end", sep = "|") else paste(gene, branch_x[i + 1L], "scan", sep = "|"),
      edge_type = "scanning",
      feature_id = NA_character_,
      x = branch_x[i],
      xend = branch_x[i + 1L],
      y = 0,
      yend = 0
    )
  }))

  translation_edges <- ann[, .(
    edge_id = paste(gene_symbol, feature_id, "translate", sep = "|"),
    gene_symbol,
    tx_id,
    source = paste(gene_symbol, feature_id, "start", sep = "|"),
    target = paste(gene_symbol, feature_id, "stop", sep = "|"),
    edge_type = fifelse(feature_id == "clean_CDS", "cds_translation", "translon_translation"),
    feature_id,
    x = tx_start,
    xend = tx_end,
    y = track_y,
    yend = track_y
  )]

  rejoin_edges <- ann[feature_id != "clean_CDS", .(
    edge_id = paste(gene_symbol, feature_id, "postterm", sep = "|"),
    gene_symbol,
    tx_id,
    source = paste(gene_symbol, feature_id, "stop", sep = "|"),
    target = fifelse(
      feature_type == "overlapping_uorf",
      paste(gene_symbol, "exit_after_overlapping_uorf", sep = "|"),
      paste(gene_symbol, tx_end, "scan", sep = "|")
    ),
    edge_type = fifelse(feature_type == "overlapping_uorf",
                       "competitive_exit",
                       "posttermination_scanning"),
    feature_id,
    x = tx_end,
    xend = tx_end,
    y = track_y,
    yend = 0
  )]

  list(nodes = nodes, edges = rbind(scan_edges, translation_edges, rejoin_edges, fill = TRUE))
}

build_rdg_topology <- function(feature_annotation) {
  split_ann <- split(feature_annotation, feature_annotation$gene_symbol)
  pieces <- lapply(split_ann, build_rdg_topology_one)
  list(
    nodes = rbindlist(lapply(pieces, `[[`, "nodes"), fill = TRUE),
    edges = rbindlist(lapply(pieces, `[[`, "edges"), fill = TRUE)
  )
}

build_sample_context <- function(state_dt, metadata_file) {
  state_cols <- setdiff(names(state_dt), c("Run", "sample", "dominant_state",
                                           "dominant_state_score"))
  context <- state_dt[, c("Run", "sample", "dominant_state", state_cols), with = FALSE]
  if (file.exists(metadata_file)) {
    m <- fread(metadata_file)
    keep <- intersect(c("Run", metadata_fields), names(m))
    context <- merge(context, m[, ..keep], by = "Run", all.x = TRUE, sort = FALSE)
  } else {
    warning("Metadata file not found: ", metadata_file)
  }
  context[]
}

read_model_branch_terms <- function(model_file, branch_source_default) {
  out <- data.table()
  if (!is.null(model_file) && file.exists(model_file)) {
    model_terms <- fread(model_file)
    required <- c("field", "level", "dominant_state", "adjusted_effect_z",
                  "adjusted_wilcox_greater_p_adj", "n_in")
    if (all(required %chin% names(model_terms))) {
      out <- model_terms[
        is.finite(adjusted_effect_z) &
          adjusted_effect_z > 0 &
          is.finite(adjusted_wilcox_greater_p_adj) &
          adjusted_wilcox_greater_p_adj < 0.1 &
          n_in >= min_enrichment_branch_samples
      ][order(adjusted_wilcox_greater_p_adj, -adjusted_effect_z)]
      out <- out[
        ,
        head(.SD, max_enrichment_branch_terms_per_state),
        by = dominant_state
      ]
      if (!"branch_source" %chin% names(out)) {
        out[, branch_source := branch_source_default]
      }
      out <- out[, .(
        branch_source, field, level, dominant_state, adjusted_effect_z,
        adjusted_wilcox_greater_p_adj, n_in
      )]
    }
  }
  out
}

build_branch_terms <- function(state_dt, enrichment_file, broad_file = NULL,
                               glmnet_file = NULL) {
  state_cols <- setdiff(names(state_dt), c("Run", "sample", "dominant_state",
                                           "dominant_state_score"))
  state_cols <- state_cols[state_cols != "Baseline (control)"]
  dominant_terms <- data.table(
    branch_source = "sample_dominant_call",
    field = "dominant_state_call",
    level = state_cols,
    dominant_state = state_cols,
    adjusted_effect_z = NA_real_,
    adjusted_wilcox_greater_p_adj = NA_real_,
    n_in = NA_integer_
  )

  enrichment_terms <- data.table()
  if (file.exists(enrichment_file)) {
    enrichment <- fread(enrichment_file)
    required <- c("field", "level", "dominant_state", "adjusted_effect_z",
                  "adjusted_wilcox_greater_p_adj", "n_in")
    if (all(required %chin% names(enrichment))) {
      enrichment_terms <- enrichment[
        is.finite(adjusted_effect_z) &
          adjusted_effect_z > 0 &
          is.finite(adjusted_wilcox_greater_p_adj) &
          adjusted_wilcox_greater_p_adj < 0.05 &
          n_in >= min_enrichment_branch_samples
      ][order(adjusted_wilcox_greater_p_adj, -adjusted_effect_z)]
      enrichment_terms <- enrichment_terms[
        ,
        head(.SD, max_enrichment_branch_terms_per_state),
        by = dominant_state
      ]
      setorder(enrichment_terms, adjusted_wilcox_greater_p_adj, -adjusted_effect_z)
      if (is.finite(max_enrichment_branch_terms)) {
        enrichment_terms <- enrichment_terms[seq_len(min(.N, max_enrichment_branch_terms))]
      }
      enrichment_terms[, branch_source := "metadata_enrichment"]
      enrichment_terms <- enrichment_terms[, .(
        branch_source, field, level, dominant_state, adjusted_effect_z,
        adjusted_wilcox_greater_p_adj, n_in
      )]
    }
  }

  broad_terms <- read_model_branch_terms(broad_file, "broad_metadata_nwise")
  glmnet_terms <- read_model_branch_terms(glmnet_file, "glmnet_regularized_metadata")

  terms <- rbind(dominant_terms, enrichment_terms, broad_terms, glmnet_terms,
                 fill = TRUE)
  terms[, branch_id := paste(branch_source, file_safe(field), file_safe(level),
                             file_safe(dominant_state), sep = "|")]
  terms[, branch_label := mapply(
    branch_label_text,
    branch_source,
    field,
    level,
    dominant_state,
    USE.NAMES = FALSE
  )]
  unique(terms, by = "branch_id")
}

branch_runs <- function(context, branch) {
  if (branch$field == "dominant_state_call") {
    in_branch <- context$dominant_state == branch$dominant_state
  } else if (grepl(":", branch$field, fixed = TRUE)) {
    values <- parse_conjunction_level(branch$field, branch$level)
    if (!all(names(values) %chin% names(context)) || anyNA(values)) {
      in_branch <- rep(FALSE, nrow(context))
    } else {
      in_branch <- rep(TRUE, nrow(context))
      for (field in names(values)) {
        in_branch <- in_branch &
          clean_level(context[[field]]) == single_field_level(field, values[[field]])
      }
    }
  } else if (branch$field %chin% names(context)) {
    in_branch <- clean_level(context[[branch$field]]) ==
      single_field_level(branch$field, branch$level)
  } else {
    in_branch <- rep(FALSE, nrow(context))
  }
  in_branch[is.na(in_branch)] <- FALSE
  setNames(
    list(context[in_branch, Run], context[!in_branch, Run]),
    c("in", "out")
  )
}

score_features_for_branches <- function(feature_dt, feature_annotation,
                                        context, branch_terms) {
  if (!"raw_counts" %chin% names(feature_dt)) {
    stop(
      "Missing raw_counts in ", feature_expression_file,
      ". Rerun dominant_uorf_regulation_inference.R before RDG scoring."
    )
  }
  feature_dt[, raw_counts := as.numeric(raw_counts)]
  feature_dt <- merge(
    feature_dt,
    feature_annotation[, .(
      gene_symbol, tx_id, feature_id, tx_start, tx_end, transcript_length,
      track_y, feature_order, annotation_feature_type = feature_type,
      annotation_feature_bases = feature_bases
    )],
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
  feature_dt[is.na(feature_type), feature_type := annotation_feature_type]
  feature_dt[is.na(feature_bases), feature_bases := annotation_feature_bases]

  rbindlist(lapply(seq_len(nrow(branch_terms)), function(i) {
    branch <- branch_terms[i]
    runs <- branch_runs(context, branch)
    rbindlist(lapply(names(runs), function(group_name) {
      group_runs <- runs[[group_name]]
      if (length(group_runs) < min_branch_samples) return(NULL)
      scores <- feature_dt[
        Run %chin% group_runs,
	        .(
	          n_samples = uniqueN(Run),
	          sum_raw_counts = sum_or_na(raw_counts),
	          mean_fpkm_like = mean(fpkm_like, na.rm = TRUE),
	          median_fpkm_like = median(fpkm_like, na.rm = TRUE)
	        ),
        by = .(
          gene_symbol, tx_id, feature_id, feature_type, feature_bases,
          tx_start, tx_end, transcript_length, track_y, feature_order
        )
      ]
      scores[, branch_group := group_name]
      scores
    }), fill = TRUE)[
      ,
      `:=`(
        branch_id = branch$branch_id,
        branch_source = branch$branch_source,
        field = branch$field,
        level = branch$level,
        dominant_state = branch$dominant_state,
        branch_label = branch$branch_label,
        adjusted_effect_z = branch$adjusted_effect_z,
        adjusted_wilcox_greater_p_adj = branch$adjusted_wilcox_greater_p_adj
      )
    ]
	  }), fill = TRUE)[
	    ,
	    allocation_feature := is_cds_allocation_feature_type(feature_type)
	  ][
	    ,
	    total_mean_fpkm_like := sum(mean_fpkm_like, na.rm = TRUE),
	    by = .(branch_id, branch_group, gene_symbol)
	  ][
	    ,
	    allocation_total_mean_fpkm_like :=
	      sum(mean_fpkm_like[allocation_feature], na.rm = TRUE),
	    by = .(branch_id, branch_group, gene_symbol)
	  ][
	    ,
	    clean_cds_raw_counts := sum_raw_counts[feature_id == "clean_CDS"][1],
	    by = .(branch_id, branch_group, gene_symbol)
	  ][
	    ,
	    passes_clean_cds_count_cutoff :=
	      is.finite(clean_cds_raw_counts) &
	      clean_cds_raw_counts >= min_clean_cds_raw_counts
	  ][
	    ,
	    relative_translation := fifelse(
	      passes_clean_cds_count_cutoff &
	        is.finite(allocation_total_mean_fpkm_like) &
	        allocation_total_mean_fpkm_like > 0,
	      mean_fpkm_like / allocation_total_mean_fpkm_like,
	      NA_real_
	    )
  ][
    order(branch_id, branch_group, gene_symbol, tx_start, tx_end)
  ][
	    ,
	    remaining_mean_fpkm_like := rev(cumsum(rev(mean_fpkm_like))),
	    by = .(branch_id, branch_group, gene_symbol)
	  ][
	    passes_clean_cds_count_cutoff == FALSE,
	    remaining_mean_fpkm_like := NA_real_
	  ][
	    ,
	    initiation_probability_proxy := fifelse(
	      passes_clean_cds_count_cutoff &
	        is.finite(remaining_mean_fpkm_like) & remaining_mean_fpkm_like > 0,
	      pmin(pmax(mean_fpkm_like / remaining_mean_fpkm_like, 0), 1),
	      NA_real_
	    )
	  ][
	    ,
	    clean_cds_count_filter := fifelse(
	      passes_clean_cds_count_cutoff,
	      NA_character_,
	      paste0("clean_CDS raw counts < ", min_clean_cds_raw_counts)
	    )
	  ][]
}

contrast_branch_scores <- function(scores) {
  in_scores <- scores[branch_group == "in"]
  out_scores <- scores[branch_group == "out"]
  key_cols <- c("branch_id", "branch_source", "field", "level", "dominant_state",
                "branch_label", "gene_symbol", "tx_id", "feature_id", "feature_type")
	  merged <- merge(
	    in_scores[, c(key_cols, "n_samples", "mean_fpkm_like",
	                  "sum_raw_counts", "clean_cds_raw_counts",
	                  "passes_clean_cds_count_cutoff", "relative_translation",
	                  "initiation_probability_proxy"), with = FALSE],
	    out_scores[, c(key_cols, "n_samples", "mean_fpkm_like",
	                   "sum_raw_counts", "clean_cds_raw_counts",
	                   "passes_clean_cds_count_cutoff", "relative_translation",
	                   "initiation_probability_proxy"), with = FALSE],
	    by = key_cols,
	    suffixes = c("_in", "_out"),
	    all = FALSE,
    sort = FALSE
	  )
	  merged[, `:=`(
	    comparison_passes_clean_cds_count_cutoff =
	      passes_clean_cds_count_cutoff_in &
	      passes_clean_cds_count_cutoff_out,
	    log2_mean_fpkm_like_in_vs_out = log2((mean_fpkm_like_in + 1) /
	                                           (mean_fpkm_like_out + 1)),
	    delta_relative_translation = relative_translation_in - relative_translation_out,
	    delta_initiation_probability_proxy =
	      initiation_probability_proxy_in - initiation_probability_proxy_out
	  )]
	  merged[
	    comparison_passes_clean_cds_count_cutoff == FALSE,
	    `:=`(
	      log2_mean_fpkm_like_in_vs_out = NA_real_,
	      delta_relative_translation = NA_real_,
	      delta_initiation_probability_proxy = NA_real_
	    )
	  ]
	  merged[]
	}

select_genes_to_plot <- function(feature_annotation, contrasts) {
  multi_feature_genes <- feature_annotation[
    ,
    .(n_features = uniqueN(feature_id),
      n_uorfs = uniqueN(feature_id[feature_id != "clean_CDS"])),
    by = gene_symbol
  ][n_uorfs > 0][order(-n_features), gene_symbol]

  ranked <- contrasts[
    feature_id == "clean_CDS" & is.finite(log2_mean_fpkm_like_in_vs_out)
  ][
    ,
    .(max_abs_cds_shift = max(abs(log2_mean_fpkm_like_in_vs_out), na.rm = TRUE)),
    by = gene_symbol
  ][order(-max_abs_cds_shift), gene_symbol]

  genes <- unique(c("ATF4", ranked, multi_feature_genes))
  genes <- genes[genes %chin% feature_annotation$gene_symbol]
  genes[seq_len(min(length(genes), max_genes_to_plot))]
}

select_branch_panels_for_gene <- function(gene, branch_terms, contrasts) {
  gene_ranked <- contrasts[
    gene_symbol == gene &
      branch_id %chin% branch_terms$branch_id &
      feature_id == "clean_CDS" &
      is.finite(log2_mean_fpkm_like_in_vs_out)
  ][
    order(-abs(log2_mean_fpkm_like_in_vs_out)),
    .SD[1],
    by = dominant_state
  ][
    order(-abs(log2_mean_fpkm_like_in_vs_out)),
    branch_id
  ]
  selected <- unique(c(
    gene_ranked,
    branch_terms[branch_source == "sample_dominant_call", branch_id],
    branch_terms[branch_source != "sample_dominant_call", branch_id]
  ))
  selected[seq_len(min(length(selected), max_branch_panels_per_gene))]
}

plot_rdg_gene <- function(gene, feature_annotation, scores, branch_terms,
                          state_filter = NULL) {
  gene_ann <- feature_annotation[gene_symbol == gene]
  if (nrow(gene_ann) == 0) return(NULL)
  plot_branch_terms <- if (!is.null(state_filter)) {
    branch_terms[dominant_state %chin% state_filter]
  } else {
    branch_terms
  }
  if (nrow(plot_branch_terms) == 0) return(NULL)

  plot_branch_ids <- select_branch_panels_for_gene(
    gene,
    plot_branch_terms,
    contrast_branch_scores(scores[gene_symbol == gene])
  )
  plot_dt <- scores[gene_symbol == gene & branch_id %chin% plot_branch_ids]
  if (nrow(plot_dt) == 0) return(NULL)

  branch_meta <- unique(plot_dt[, .(
    branch_id, branch_source, field, level, dominant_state, branch_label,
    adjusted_effect_z, adjusted_wilcox_greater_p_adj
  )])
  branch_meta <- branch_meta[match(plot_branch_ids, branch_id)]
  branch_meta <- branch_meta[!is.na(branch_id)]
  branch_meta[, branch_rank := seq_len(.N)]
  branch_meta[, row_y := rev(seq_len(.N)) * 1.52]
  branch_meta[, section := fifelse(
	    branch_source == "sample_dominant_call",
	    "Dominant-state winner split",
    paste0("Metadata split enriched for ", short_state_label(dominant_state))
  )]
  branch_meta[, split_label := mapply(
    compact_branch_split_label,
    field,
    level,
    USE.NAMES = FALSE
  )]
  branch_meta[, display_label := fifelse(
    branch_source == "sample_dominant_call",
    paste0("IN: top score is\n", short_state_label(dominant_state)),
    paste0("IN: ", split_label,
           "\nscore: ", short_state_label(dominant_state))
  )]
  branch_meta[, effect_label := fifelse(
    is.finite(adjusted_effect_z),
    paste0("z=", sprintf("%.2f", adjusted_effect_z),
           "  q=", format(adjusted_wilcox_greater_p_adj, digits = 2, scientific = TRUE)),
    ""
  )]

  plot_dt <- merge(
    plot_dt,
    branch_meta[, .(branch_id, branch_rank, row_y, section, display_label,
                    effect_label)],
    by = "branch_id",
    all.x = TRUE,
    sort = FALSE
  )
  plot_dt[, group_offset := fifelse(branch_group == "out", 0.26, -0.26)]
  plot_dt[, track_y_plot := row_y + group_offset]
  plot_dt[, feature_lane_offset := fifelse(
    feature_id == "clean_CDS", -0.14,
    fifelse(
      feature_type == "overlapping_uorf", 0.18,
      fifelse(feature_type == "internal_orf", -0.20, 0.13)
    )
  )]
  plot_dt[, feature_y_plot := track_y_plot + feature_lane_offset]
  plot_dt[, group_label := fifelse(branch_group == "in", "IN", "OUT")]

  max_tx <- max(plot_dt$transcript_length, na.rm = TRUE)
  label_space <- max(620, max_tx * 0.42)
  x_left <- -label_space
  x_group <- -max(48, label_space * 0.07)
  x_axis_start <- 1
  x_axis_end <- max_tx

  row_bg <- unique(branch_meta[, .(
    branch_id, branch_rank, row_y, section, display_label, effect_label
  )])
  row_bg[, ymin := row_y - 0.66]
  row_bg[, ymax := row_y + 0.66]
  row_bg[, bg_fill := fifelse(branch_rank %% 2L == 0L, "#f8fafc", "#ffffff")]

	  group_labels <- unique(plot_dt[, .(
	    branch_id, branch_source, field, level, dominant_state, branch_group,
	    group_label, n_samples, clean_cds_raw_counts,
	    passes_clean_cds_count_cutoff, track_y_plot
	  )])
	  group_labels[, clean_cds_label := fifelse(
	    is.finite(clean_cds_raw_counts),
	    format(round(clean_cds_raw_counts), big.mark = ",", scientific = FALSE),
	    "NA"
	  )]
	  group_labels[, label := fifelse(
	    branch_source == "sample_dominant_call",
	    fifelse(
	      branch_group == "in",
	      paste0("IN\nn=", n_samples, "\nCDS=", clean_cds_label),
	      paste0("OUT\nn=", n_samples, "\nCDS=", clean_cds_label)
	    ),
	    fifelse(
	      branch_group == "in",
	      paste0("IN\nn=", n_samples, "\nCDS=", clean_cds_label),
	      paste0("OUT\nn=", n_samples, "\nCDS=", clean_cds_label)
	    )
	  )]

  branch_points <- plot_dt[
    ,
    .SD[which.max(mean_fpkm_like)],
    by = .(branch_id, branch_group, feature_id)
  ]
  branch_points <- branch_points[is.finite(tx_start)]

  feature_labels <- unique(plot_dt[
    ,
    .SD[which.max(mean_fpkm_like)],
    by = .(feature_id)
  ][
    ,
    .(feature_id, feature_order, tx_start, tx_end)
  ])
  feature_labels[, label_x := (tx_start + tx_end) / 2]
  feature_labels[, feature_y := max(row_bg$ymax) + 0.34 +
                   ((feature_order - 1L) %% 3L) * 0.16]
  section_labels <- row_bg[
    ,
    .(section_y = max(ymax) + 0.18),
    by = section
  ]

  p <- ggplot() +
    geom_rect(
      data = row_bg,
      aes(xmin = x_left, xmax = x_axis_end * 1.02, ymin = ymin, ymax = ymax),
      fill = row_bg$bg_fill,
      color = NA,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = row_bg,
      aes(x = x_left, xend = x_axis_end * 1.02, y = ymin, yend = ymin),
      linewidth = 0.22,
      color = "#e5e7eb",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = section_labels,
      aes(x = x_left, y = section_y, label = section),
      hjust = 0,
      vjust = 0,
      size = 3.3,
      fontface = "bold",
      color = "#374151",
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = unique(plot_dt[, .(branch_id, branch_group, track_y_plot, transcript_length)]),
      aes(x = x_axis_start, xend = transcript_length, y = track_y_plot, yend = track_y_plot),
      linewidth = 0.42,
      color = "#4b5563",
      arrow = arrow(length = unit(0.065, "inches"), type = "closed"),
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = plot_dt[feature_id != "clean_CDS" & is.finite(tx_end) & is.finite(track_y)],
      aes(x = tx_end, xend = tx_end, y = feature_y_plot, yend = track_y_plot),
      linewidth = 0.23,
      color = "#94a3b8",
      inherit.aes = FALSE
    ) +
    geom_rect(
      data = plot_dt[is.finite(tx_start) & is.finite(tx_end) & is.finite(track_y)],
      aes(
        xmin = tx_start,
        xmax = pmax(tx_end, tx_start + 1),
        ymin = feature_y_plot - 0.06,
        ymax = feature_y_plot + 0.06,
        fill = relative_translation
      ),
      color = NA,
      linewidth = 0,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = branch_points,
      aes(
        x = tx_start,
        y = track_y_plot,
        size = initiation_probability_proxy,
        fill = initiation_probability_proxy
      ),
      shape = 21,
      color = "#111827",
      stroke = 0.28,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = row_bg,
      aes(x = x_left, y = row_y + 0.16, label = display_label),
      hjust = 0,
      vjust = 0.5,
      lineheight = 0.92,
      size = 3.05,
      fontface = "bold",
      color = "#111827",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = row_bg[nchar(effect_label) > 0],
      aes(x = x_left, y = row_y - 0.42, label = effect_label),
      hjust = 0,
      vjust = 0.5,
      size = 2.65,
      color = "#6b7280",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = group_labels,
      aes(x = x_group, y = track_y_plot, label = label),
      hjust = 1,
      vjust = 0.5,
      lineheight = 0.9,
      size = 2.65,
      color = "#374151",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = feature_labels,
      aes(x = label_x, y = feature_y, label = feature_id),
      angle = 0,
      hjust = 0.5,
      vjust = 0,
      size = 2.7,
      color = "#111827",
      check_overlap = TRUE,
      inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffff",
      high = "#8b0000",
      midpoint = 0.5,
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "#d1d5db",
      name = "relative use\n/ branch proxy"
    ) +
    scale_size_continuous(
      range = c(1.4, 5.4),
      limits = c(0, 1),
      name = "branch proxy"
    ) +
    scale_x_continuous(
      limits = c(x_left, x_axis_end * 1.03),
      breaks = pretty(c(0, x_axis_end), n = 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(
      ylim = c(min(row_bg$ymin) - 0.20, max(row_bg$ymax) + 0.78),
      clip = "off"
    ) +
	    labs(
	      title = paste0(
	        gene, " ribosome decision graph",
	        if (!is.null(state_filter)) {
	          paste0(" - ", paste(state_filter, collapse = ", "))
	        } else {
	          ""
	        }
	      ),
	      subtitle = paste0(
	        "Winner split: IN = tested state is the sample's top module score.\n",
	        "Metadata split: IN = metadata level enriched for tested score, not necessarily top-state calls. ",
	        "z = adjusted effect; q = BH-adjusted Wilcoxon P.\n",
	        "CDS/uORF color uses the CDS-allocation denominator; NTE/NTT/iORF color is branch evidence. ",
	        "Circle size is the branch proxy. ",
	        "Relative-use/proxy values are NA when branch clean-CDS raw counts < ",
	        min_clean_cds_raw_counts, "."
	      ),
      x = "transcript coordinate in leader+CDS display",
      y = NULL
    ) +
    theme_minimal(base_size = 11.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(12, 20, 12, 18),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.22, "inches"),
      legend.key.height = unit(0.36, "inches"),
      legend.spacing.y = unit(0.08, "inches"),
      legend.margin = margin(0, 0, 0, 2),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(color = "#374151", size = 9.4, lineheight = 1.12),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.text.x = element_text(color = "#4b5563")
    )

  p
}

save_plot_pair <- function(plot, path_base, width = 11, height = 8) {
  png_path <- paste0(path_base, ".png")
  pdf_path <- paste0(path_base, ".pdf")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_path, plot, width = width, height = height, dpi = 180,
           device = ragg::agg_png, bg = "white")
  } else {
    ggsave(png_path, plot, width = width, height = height, dpi = 180,
           bg = "white")
  }
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  c(png = png_path, pdf = pdf_path)
}

read_existing_relative_usage_plot_manifest <- function(out_dir, state,
                                                       state_safe) {
  files <- list.files(
    out_dir,
    pattern = paste0("^relative_usage_matrix_", state_safe, "\\.(png|pdf)$"),
    full.names = TRUE
  )
  if (!length(files)) {
    return(data.table(
      dominant_state = state,
      file_type = "plots_skipped",
      path = NA_character_
    ))
  }
  data.table(
    dominant_state = state,
    file_type = tools::file_ext(files),
    path = files
  )
}

read_existing_rdg_figure_manifest <- function(output_dir, genes) {
  files <- list.files(
    file.path(output_dir, "figures"),
    pattern = "^rdg_.*\\.(png|pdf)$",
    full.names = TRUE
  )
  if (!length(files)) {
    return(data.table(
      gene_symbol = character(),
      file_type = character(),
      path = character()
    ))
  }
  gene_map <- setNames(genes, file_safe(genes))
  safe_gene <- sub("^rdg_", "", tools::file_path_sans_ext(basename(files)))
  data.table(
    gene_symbol = unname(gene_map[safe_gene]),
    file_type = tools::file_ext(files),
    path = files
  )[!is.na(gene_symbol)]
}

read_existing_state_rdg_figure_manifest <- function(output_dir, genes, states) {
  files <- list.files(
    file.path(output_dir, "figures_by_state"),
    pattern = "^rdg_.*\\.(png|pdf)$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (!length(files)) {
    return(data.table(
      dominant_state = character(),
      gene_symbol = character(),
      file_type = character(),
      path = character()
    ))
  }
  gene_map <- setNames(genes, file_safe(genes))
  state_map <- setNames(states, file_safe(states))
  state_safe <- basename(dirname(files))
  prefix <- paste0("rdg_")
  suffix <- paste0("_", state_safe)
  stem <- tools::file_path_sans_ext(basename(files))
  safe_gene <- mapply(function(x, suff) {
    x <- sub(paste0("^", prefix), "", x)
    sub(paste0(suff, "$"), "", x)
  }, stem, suffix, USE.NAMES = FALSE)
  data.table(
    dominant_state = unname(state_map[state_safe]),
    gene_symbol = unname(gene_map[safe_gene]),
    file_type = tools::file_ext(files),
    path = files
  )[!is.na(dominant_state) & !is.na(gene_symbol)]
}

ordered_feature_columns <- function(feature_ids) {
  parse_group <- function(prefix) {
    ids <- feature_ids[grepl(paste0("^", prefix, "[0-9]+$"), feature_ids)]
    ids[order(as.integer(sub(paste0("^", prefix), "", ids)))]
  }
  uorfs <- parse_group("uORF")
  ntes <- parse_group("NTE")
  ntts <- parse_group("NTT")
  iorfs <- parse_group("iORF")
  clean <- intersect("clean_CDS", feature_ids)
  used <- c(uorfs, ntes, ntts, iorfs, clean)
  unique(c(uorfs, ntes, ntts, iorfs, clean, setdiff(sort(feature_ids), used)))
}

save_relative_usage_matrices <- function(scores) {
  matrix_dir <- file.path(output_dir, "relative_usage_matrices")
  state_dir <- file.path(matrix_dir, "by_state")
  dir.create(matrix_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(state_dir, showWarnings = FALSE, recursive = TRUE)

  in_scores <- copy(scores[branch_group == "in"])
  if (nrow(in_scores) == 0) return(data.table())
  in_scores[, low_clean_cds_count := !passes_clean_cds_count_cutoff]
  long <- in_scores[, .(
    branch_id, branch_source, field, level, dominant_state, branch_label,
    gene_symbol, tx_id, feature_id, feature_type, feature_order,
    n_samples, sum_raw_counts, clean_cds_raw_counts,
    passes_clean_cds_count_cutoff, low_clean_cds_count,
	    mean_fpkm_like, allocation_feature, allocation_total_mean_fpkm_like,
	    relative_translation, initiation_probability_proxy
  )]
  setorder(long, dominant_state, field, level, gene_symbol, feature_order)
  fwrite(long, file.path(matrix_dir, "relative_usage_in_long.csv"))

  id_cols <- c("branch_id", "branch_source", "field", "level",
               "dominant_state", "branch_label", "gene_symbol", "tx_id",
               "n_samples", "clean_cds_raw_counts",
               "passes_clean_cds_count_cutoff")
  wide <- dcast(
    long,
    branch_id + branch_source + field + level + dominant_state + branch_label +
      gene_symbol + tx_id + n_samples + clean_cds_raw_counts +
      passes_clean_cds_count_cutoff ~ feature_id,
    value.var = "relative_translation"
  )
  feature_cols <- ordered_feature_columns(setdiff(names(wide), id_cols))
  setcolorder(wide, c(id_cols, feature_cols[feature_cols %chin% names(wide)]))
  setorder(wide, dominant_state, field, level, gene_symbol)
  fwrite(wide, file.path(matrix_dir, "relative_usage_in_wide.csv"))

  manifest <- rbindlist(lapply(sort(unique(long$dominant_state)), function(state) {
    state_safe <- file_safe(state)
    out_dir <- file.path(state_dir, state_safe)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    state_long <- long[dominant_state == state]
    state_wide <- wide[dominant_state == state]
    fwrite(state_long, file.path(out_dir, "relative_usage_in_long.csv"))
    fwrite(state_wide, file.path(out_dir, "relative_usage_in_wide.csv"))
    if (!generate_rdg_plots) {
      return(read_existing_relative_usage_plot_manifest(out_dir, state, state_safe))
    }

    plot_dt <- copy(state_long)
    plot_dt[, feature_id := factor(
      feature_id,
      levels = ordered_feature_columns(unique(feature_id))
    )]
    plot_dt[, branch_label_compact := gsub("\\n", " / ", branch_label)]
    plot_dt[, row_label := paste0(gene_symbol, " | ", branch_label_compact)]

    row_order <- plot_dt[
      ,
      .(
        pass = any(passes_clean_cds_count_cutoff),
        clean_cds = max(clean_cds_raw_counts, na.rm = TRUE),
        cds_relative = suppressWarnings(max(
          relative_translation[as.character(feature_id) == "clean_CDS"],
          na.rm = TRUE
        ))
      ),
      by = row_label
    ]
    row_order[!is.finite(clean_cds), clean_cds := NA_real_]
    row_order[!is.finite(cds_relative), cds_relative := NA_real_]
    setorder(row_order, -pass, -cds_relative, -clean_cds, row_label)
    if (is.finite(max_matrix_plot_rows) && nrow(row_order) > max_matrix_plot_rows) {
      row_order <- row_order[seq_len(max_matrix_plot_rows)]
    }
    plot_dt <- plot_dt[row_label %chin% row_order$row_label]
    plot_dt[, row_label := factor(row_label, levels = rev(row_order$row_label))]
    if (nrow(plot_dt) == 0) return(NULL)

    p <- ggplot(
      plot_dt,
      aes(x = feature_id, y = row_label, fill = relative_translation)
    ) +
      geom_tile(color = "#f3f4f6", linewidth = 0.18) +
      scale_fill_gradient2(
        low = "#2166ac",
        mid = "#ffffff",
        high = "#8b0000",
        midpoint = 0.5,
        limits = c(0, 1),
        oob = scales::squish,
        na.value = "#d1d5db",
        name = "relative use"
      ) +
      labs(
        title = paste0("IN-group relative usage matrix: ", state),
        subtitle = paste0(
          "Rows are gene x branch terms. CDS and uORFs use the CDS-allocation denominator; ",
          "NTE/NTT/iORF values are branch-evidence ratios. Grey = clean-CDS raw counts < ",
          min_clean_cds_raw_counts,
          " or unavailable."
        ),
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10.5) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "#111827"),
        axis.text.y = element_text(size = 6.8, color = "#111827"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#374151", size = 9.2),
        legend.position = "right",
        legend.key.width = unit(0.24, "inches"),
        legend.key.height = unit(0.36, "inches")
      )

    height <- min(26, max(5, 0.17 * uniqueN(plot_dt$row_label) + 2.4))
    width <- min(13, max(7.5, 0.38 * uniqueN(plot_dt$feature_id) + 6.5))
    paths <- save_plot_pair(
      p,
      file.path(out_dir, paste0("relative_usage_matrix_", state_safe)),
      width = width,
      height = height
    )
    data.table(
      dominant_state = state,
      file_type = names(paths),
      path = unname(paths)
    )
  }), fill = TRUE)
  fwrite(manifest, file.path(matrix_dir, "relative_usage_matrix_manifest.csv"))
  manifest[]
}

write_readme <- function() {
  text <- c(
    "Dominant-state ribosome decision graph outputs",
    "================================================",
    "",
    "This directory is produced by dominant_ribosome_decision_graphs.R.",
    "",
    "The script follows the RDG abstraction from Tierney et al. Genome Research 2024:",
    "- ribosome paths include scanned and translated segments;",
    "- translated segments are represented as translons;",
    "- starts/reinitiation-like events are treated as branch points;",
    "- quantitative Ribo-seq data are used to color translons and approximate branch probabilities.",
    "",
    "Inputs expected from the current dominant-state workflow:",
    "- dominant_uorf_feature_expression.csv",
    "- human_dominant_cell_states_clean_cds.csv",
    "- human_dominant_cell_states_clean_cds_gene_diagnostics.csv",
    "- dominant_state_clean_cds_enrichment_all_main.csv",
    "- dominant_state_broad_metadata_rdg_terms.csv, when available",
    "- dominant_state_glmnet_metadata_rdg_terms.csv, when available",
    "- metadata_done_samples_extended_qc.csv",
    "",
	    "Important caveat:",
	    "The branch probabilities here are proxies derived from relative regional coverage.",
	    "They are not yet mechanistic initiation probabilities. They become candidate branch",
	    "parameters to test with better path models after all required FST pages are local.",
	    "",
	    paste0("Low-count filter: relative-use and branch-proxy values are set to NA when"),
	    paste0("the branch-level clean_CDS raw-count total is < ", min_clean_cds_raw_counts, "."),
	    "",
	    "Additional compact views:",
	    "- relative_usage_matrices/relative_usage_in_long.csv",
	    "- relative_usage_matrices/relative_usage_in_wide.csv",
	    "- relative_usage_matrices/by_state/<dominant_state>/"
	  )
  writeLines(text, file.path(output_dir, "README.txt"))
}

if (!file.exists(feature_expression_file)) {
  stop("Missing ", feature_expression_file,
       ". Run dominant_uorf_regulation_inference.R after the required FST pages are available.")
}
if (!file.exists(state_file)) stop("Missing ", state_file)

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
translons <- load_predicted_translons(df)
feature_dt <- fread(feature_expression_file)
state_dt <- fread(state_file)

feature_dt[, gene_symbol := toupper(gene_symbol)]
feature_dt <- feature_dt[is.finite(fpkm_like)]
if (!"raw_counts" %chin% names(feature_dt)) {
  stop(
    "Missing raw_counts in ", feature_expression_file,
    ". Rerun dominant_uorf_regulation_inference.R before RDG scoring."
  )
}
if (nrow(feature_dt) == 0) stop("No finite feature-expression rows in ", feature_expression_file)

target_tx <- unique(feature_dt$tx_id)
cds_all <- loadRegion(df, "cds", names.keep = target_tx)
mrna_all <- loadRegion(df, "mrna", names.keep = target_tx)
translons <- translons[names(translons) %chin% target_tx]
translons <- annotate_translon_categories(translons, cds_all, mrna_all)
translons <- deduplicate_translons_by_tx_coordinates(translons)
message("Categorized predicted translons for RDG transcripts: ",
        paste(names(table(mcols(translons)$category)),
              as.integer(table(mcols(translons)$category)),
              sep = "=", collapse = ", "))
message("Translon source composition after deduplication: ",
        paste(names(table(mcols(translons)$translon_sources_merged)),
              as.integer(table(mcols(translons)$translon_sources_merged)),
              sep = "=", collapse = ", "))

feature_annotation <- build_feature_annotation(df, translons, feature_dt)
feature_annotation <- feature_annotation[feature_id %chin% unique(feature_dt$feature_id)]
if (nrow(feature_annotation) == 0) {
  stop("Could not build feature annotation for rows in ", feature_expression_file)
}

topology <- build_rdg_topology(feature_annotation)
sample_context <- build_sample_context(state_dt, metadata_file)
branch_terms <- build_branch_terms(
  state_dt,
  enrichment_file,
  broad_rdg_terms_file,
  glmnet_rdg_terms_file
)
scores <- score_features_for_branches(feature_dt, feature_annotation,
                                      sample_context, branch_terms)
contrasts <- contrast_branch_scores(scores)
genes_to_plot <- select_genes_to_plot(feature_annotation, contrasts)

fwrite(feature_annotation, file.path(output_dir, "rdg_feature_annotation.csv"))
fwrite(topology$nodes, file.path(output_dir, "rdg_nodes.csv"))
fwrite(topology$edges, file.path(output_dir, "rdg_edges.csv"))
fwrite(branch_terms, file.path(output_dir, "rdg_branch_terms.csv"))
fwrite(scores, file.path(output_dir, "rdg_branch_feature_scores.csv"))
fwrite(contrasts, file.path(output_dir, "rdg_branch_feature_contrasts.csv"))
writeLines(genes_to_plot, file.path(output_dir, "rdg_genes_plotted.txt"))
matrix_manifest <- save_relative_usage_matrices(scores)
write_readme()

if (generate_rdg_plots) {
  figure_manifest <- rbindlist(lapply(genes_to_plot, function(gene) {
    message("Plotting RDG overlays for ", gene)
    p <- plot_rdg_gene(gene, feature_annotation, scores, branch_terms)
    if (is.null(p)) return(NULL)
    base <- file.path(output_dir, "figures", paste0("rdg_", file_safe(gene)))
    paths <- save_plot_pair(
      p,
      base,
      width = 12.8,
      height = 1.45 * max_branch_panels_per_gene + 2
    )
    data.table(gene_symbol = gene, file_type = names(paths), path = unname(paths))
  }), fill = TRUE)
  if (is.null(figure_manifest)) figure_manifest <- data.table()
} else {
  message("Skipping RDG overlay figures; table outputs were regenerated.")
  figure_manifest <- read_existing_rdg_figure_manifest(output_dir, genes_to_plot)
}
fwrite(figure_manifest, file.path(output_dir, "rdg_figure_manifest.csv"))

if (generate_rdg_plots) {
  state_figure_manifest <- rbindlist(lapply(sort(unique(branch_terms$dominant_state)), function(state) {
    state_genes <- contrasts[
      dominant_state == state &
        feature_id == "clean_CDS" &
        is.finite(log2_mean_fpkm_like_in_vs_out),
      unique(gene_symbol)
    ]
    if (length(state_genes) == 0) return(NULL)
    state_safe <- file_safe(state)
    state_dir <- file.path(output_dir, "figures_by_state", state_safe)
    dir.create(state_dir, showWarnings = FALSE, recursive = TRUE)
    rbindlist(lapply(state_genes, function(gene) {
      message("Plotting state-specific RDG overlay for ", gene, " / ", state)
      p <- plot_rdg_gene(gene, feature_annotation, scores, branch_terms,
                         state_filter = state)
      if (is.null(p)) return(NULL)
      base <- file.path(state_dir, paste0("rdg_", file_safe(gene), "_", state_safe))
      paths <- save_plot_pair(
        p,
        base,
        width = 12.8,
        height = 1.45 * max_branch_panels_per_gene + 2
      )
      data.table(
        dominant_state = state,
        gene_symbol = gene,
        file_type = names(paths),
        path = unname(paths)
      )
    }), fill = TRUE)
  }), fill = TRUE)
  if (is.null(state_figure_manifest)) state_figure_manifest <- data.table()
} else {
  state_figure_manifest <- read_existing_state_rdg_figure_manifest(
    output_dir,
    genes_to_plot,
    sort(unique(branch_terms$dominant_state))
  )
}
fwrite(state_figure_manifest, file.path(output_dir, "rdg_state_figure_manifest.csv"))

message("Saved RDG outputs in: ", normalizePath(output_dir))
message("Genes plotted: ", paste(genes_to_plot, collapse = ", "))
message(if (generate_rdg_plots) {
  "Figures written:"
} else {
  "Existing figures referenced in manifest:"
})
if (nrow(figure_manifest) > 0) {
  for (path in figure_manifest$path) message("  ", path)
} else {
  message("  none")
}

invisible(list(
  feature_annotation = feature_annotation,
  nodes = topology$nodes,
  edges = topology$edges,
  branch_terms = branch_terms,
	    branch_feature_scores = scores,
	    branch_feature_contrasts = contrasts,
	    figure_manifest = figure_manifest,
	    state_figure_manifest = state_figure_manifest,
	    matrix_manifest = matrix_manifest
	  ))
