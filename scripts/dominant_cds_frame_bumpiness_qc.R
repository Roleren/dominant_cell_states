#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
})

find_analysis_dir <- function(start = getwd()) {
  here <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (basename(here) == "dominant_cell_states" &&
        dir.exists(file.path(here, "scripts"))) {
      return(here)
    }
    candidate <- file.path(here, "dominant_cell_states")
    if (dir.exists(candidate) &&
        dir.exists(file.path(candidate, "scripts"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("Could not find dominant_cell_states analysis directory from: ",
       start, call. = FALSE)
}

analysis_dir <- find_analysis_dir()
output_dir <- file.path(analysis_dir, "dominant_cds_frame_bumpiness_qc")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)
coverage_cache_file <- file.path(analysis_dir, "scripts", "dominant_coverage_cache.R")
if (!file.exists(coverage_cache_file)) {
  stop("Missing dominant coverage cache helper: ", coverage_cache_file)
}
source(coverage_cache_file)

message("CDS frame and bumpiness QC")
message("  1. Load FST transcript coverage for all modeled marker CDSs")
message("  2. Estimate per-run CDS-frame fractions in middle CDS regions")
message("  3. Test whether RRM1/MYC/EGLN3 viral runs are globally mis-shifted")
message("  4. Compare CDS bumpiness against coverage depth and viral controls")

target_genes <- c("RRM1", "MYC", "EGLN3", "DDIT3")
candidate_genes <- c("RRM1", "MYC", "EGLN3")
control_genes <- c("DDIT3")
housekeeping_candidate_genes <- c(
  "ACTB", "GAPDH", "HPRT1", "TBP", "PPIA", "UBC", "YWHAZ",
  "EEF1A1", "PSMB4", "TPT1", "RPL27", "RPL30", "RPS18", "RPS27A"
)
metadata_file <- local({
  candidates <- c(
    Sys.getenv("DOMINANT_METADATA_FILE", unset = NA_character_),
    "/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv",
    path.expand("~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended_qc.csv")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
})
diagnostics_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"
)
feature_contrasts_file <- file.path(
  analysis_dir, "dominant_next_model",
  "dominant_next_model_matched_control_feature_contrasts.csv"
)

if (!file.exists(diagnostics_file)) stop("Missing diagnostics: ", diagnostics_file)
if (!file.exists(feature_contrasts_file)) {
  stop("Missing feature contrasts: ", feature_contrasts_file)
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

add_metadata_frame_qc <- function(dt) {
  usage_cols <- paste0("Frame_usage_", 0:2)
  out_cols <- c(
    paste0("metadata_frame", 0:2, "_fraction"),
    "metadata_canonical_frame_fraction",
    "metadata_dominant_frame",
    "metadata_possible_misshift"
  )
  if (!all(usage_cols %in% names(dt))) {
    for (col in out_cols) {
      if (!col %in% names(dt)) dt[, (col) := NA]
    }
    return(dt)
  }

  raw0 <- suppressWarnings(as.numeric(dt[["Frame_usage_0"]]))
  raw1 <- suppressWarnings(as.numeric(dt[["Frame_usage_1"]]))
  raw2 <- suppressWarnings(as.numeric(dt[["Frame_usage_2"]]))
  raw_total <- raw0 + raw1 + raw2
  usage_is_percent <- is.finite(raw_total) & raw_total > 1.5
  frame0 <- ifelse(usage_is_percent, raw0 / 100, raw0)
  frame1 <- ifelse(usage_is_percent, raw1 / 100, raw1)
  frame2 <- ifelse(usage_is_percent, raw2 / 100, raw2)
  available <- is.finite(frame0) & is.finite(frame1) & is.finite(frame2) &
    is.finite(raw_total) & raw_total > 0

  frame_mat <- cbind(frame0, frame1, frame2)
  dominant_idx <- max.col(frame_mat, ties.method = "first") - 1L
  dominant_frame <- paste0("frame", dominant_idx)
  dominant_frame[!available] <- NA_character_

  dt[, metadata_frame0_fraction := frame0]
  dt[, metadata_frame1_fraction := frame1]
  dt[, metadata_frame2_fraction := frame2]
  dt[!available, c("metadata_frame0_fraction", "metadata_frame1_fraction",
                   "metadata_frame2_fraction") := .(NA_real_, NA_real_,
                                                    NA_real_)]
  dt[, metadata_dominant_frame := dominant_frame]
  dt[, metadata_canonical_frame_fraction := metadata_frame0_fraction]
  dt[, metadata_possible_misshift :=
       is.finite(metadata_canonical_frame_fraction) &
       metadata_canonical_frame_fraction < 0.45 &
       metadata_dominant_frame != "frame0"]
  dt[is.na(metadata_possible_misshift), metadata_possible_misshift := FALSE]
  dt
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

parse_run_signature <- function(x) {
  x <- paste(as.character(x), collapse = ";")
  runs <- regmatches(
    x,
    gregexpr("(SRR|ERR|DRR)[0-9]+", x, perl = TRUE)
  )[[1]]
  unique(runs[nzchar(runs)])
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, y)$p.value)
}

collapse_text <- function(x, max_n = 5L) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) return(NA_character_)
  suffix <- if (length(x) > max_n) ";..." else ""
  paste0(paste(head(x, max_n), collapse = ";"), suffix)
}

make_viral_design_run_map <- function(viral_dt) {
  if (!nrow(viral_dt)) return(data.table())
  rbindlist(lapply(seq_len(nrow(viral_dt)), function(i) {
    row <- viral_dt[i]
    case_runs <- parse_run_signature(row$case_run_signature)
    control_runs <- parse_run_signature(row$control_run_signature)
    design_id <- paste(
      row$gene_symbol, row$study, row$CELL_LINE, row$case_condition,
      i, sep = "|"
    )
    common <- list(
      gene_symbol = row$gene_symbol,
      viral_design_id = design_id,
      viral_study = row$study,
      viral_AUTHOR = row$AUTHOR,
      viral_CELL_LINE = row$CELL_LINE,
      viral_TISSUE = row$TISSUE,
      viral_case_condition = row$case_condition,
      viral_control_conditions = row$control_conditions,
      viral_case_n_runs = length(case_runs),
      viral_control_n_runs = length(control_runs)
    )
    rbind(
      as.data.table(c(common, list(
        Run = case_runs,
        viral_design_role = "case"
      ))),
      as.data.table(c(common, list(
        Run = control_runs,
        viral_design_role = "control"
      )))
    )
  }), fill = TRUE)
}

collapse_viral_run_map <- function(run_map, role) {
  if (!nrow(run_map)) return(data.table())
  run_map[
    viral_design_role == role,
    .(
      viral_design_id = collapse_text(viral_design_id, max_n = 3L),
      viral_study = collapse_text(viral_study, max_n = 3L),
      viral_AUTHOR = collapse_text(viral_AUTHOR, max_n = 3L),
      viral_CELL_LINE = collapse_text(viral_CELL_LINE, max_n = 3L),
      viral_TISSUE = collapse_text(viral_TISSUE, max_n = 3L),
      viral_case_condition = collapse_text(viral_case_condition, max_n = 3L),
      viral_control_conditions =
        collapse_text(viral_control_conditions, max_n = 3L),
      viral_case_n_runs = max(viral_case_n_runs, na.rm = TRUE),
      viral_control_n_runs = max(viral_control_n_runs, na.rm = TRUE)
    ),
    by = .(target_gene = gene_symbol, Run)
  ]
}

middle_cds_positions <- function(cds_positions, min_edge = 45L,
                                 fraction_edge = 0.10) {
  cds_positions <- sort(unique(as.integer(cds_positions)))
  n <- length(cds_positions)
  if (n < 120L) return(cds_positions)
  edge <- max(as.integer(min_edge), floor(n * fraction_edge))
  if (n <= 2L * edge + 30L) {
    edge <- max(15L, floor(n * 0.05))
  }
  cds_positions[(edge + 1L):(n - edge)]
}

safe_col_quantile <- function(mat, prob) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colQuantiles(mat, probs = prob, drop = TRUE,
                                     na.rm = TRUE))
  }
  apply(mat, 2, stats::quantile, probs = prob, na.rm = TRUE, names = FALSE)
}

safe_col_max <- function(mat) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colMaxs(mat, na.rm = TRUE))
  }
  apply(mat, 2, max, na.rm = TRUE)
}

coverage_metrics <- function(mat, cds_positions, use_positions, gene, tx_id,
                             region = "middle_cds") {
  use_positions <- sort(unique(as.integer(use_positions)))
  use_positions <- use_positions[
    is.finite(use_positions) & use_positions >= 1L & use_positions <= nrow(mat)
  ]
  runs <- colnames(mat)
  if (length(use_positions) == 0L) {
    return(data.table(
      gene_symbol = gene,
      tx_id = tx_id,
      Run = runs,
      region = region,
      cds_positions_used = 0L
    ))
  }
  sub <- mat[use_positions, , drop = FALSE]
  storage.mode(sub) <- "numeric"
  total <- as.numeric(colSums(sub, na.rm = TRUE))
  phase <- as.integer((use_positions - min(cds_positions)) %% 3L)
  phase_counts <- lapply(0:2, function(frame) {
    as.numeric(colSums(sub[phase == frame, , drop = FALSE], na.rm = TRUE))
  })
  phase_counts <- do.call(cbind, phase_counts)
  colnames(phase_counts) <- paste0("frame", 0:2, "_counts")
  frame_fraction <- phase_counts / pmax(total, 1)
  frame_fraction[total <= 0, ] <- NA_real_
  colnames(frame_fraction) <- paste0("frame", 0:2, "_fraction")

  n_pos <- length(use_positions)
  mean_cov <- total / n_pos
  sumsq <- as.numeric(colSums(sub^2, na.rm = TRUE))
  variance <- (sumsq - n_pos * mean_cov^2) / pmax(n_pos - 1L, 1L)
  variance[variance < 0] <- 0
  cv <- sqrt(variance) / mean_cov
  cv[!is.finite(cv)] <- NA_real_
  max_count <- safe_col_max(sub)
  max_share <- max_count / total
  max_share[!is.finite(max_share)] <- NA_real_
  zero_fraction <- as.numeric(colMeans(sub == 0, na.rm = TRUE))

  q99 <- safe_col_quantile(sub, 0.99)
  top_mask <- sweep(sub, 2, q99, `>=`)
  top_signal <- as.numeric(colSums(sub * top_mask, na.rm = TRUE))
  top1pct_share <- top_signal / total
  top1pct_share[!is.finite(top1pct_share)] <- NA_real_
  top_phase_counts <- lapply(0:2, function(frame) {
    as.numeric(colSums(sub[phase == frame, , drop = FALSE] *
                         top_mask[phase == frame, , drop = FALSE],
                       na.rm = TRUE))
  })
  top_phase_counts <- do.call(cbind, top_phase_counts)
  top_phase_fraction <- top_phase_counts / pmax(top_signal, 1)
  top_phase_fraction[top_signal <= 0, ] <- NA_real_
  colnames(top_phase_fraction) <- paste0("top1pct_frame", 0:2, "_fraction")

  out <- data.table(
    gene_symbol = gene,
    tx_id = tx_id,
    Run = runs,
    region = region,
    cds_positions_used = n_pos,
    cds_total_counts = total,
    cds_mean_coverage = mean_cov,
    cds_coverage_cv = cv,
    cds_zero_fraction = zero_fraction,
    cds_max_position_share = max_share,
    cds_top1pct_share = top1pct_share,
    cds_top1pct_signal = top_signal
  )
  out <- cbind(out, as.data.table(phase_counts), as.data.table(frame_fraction),
               as.data.table(top_phase_fraction))
  out[, canonical_frame_fraction := frame0_fraction]
  out[, off_frame_fraction := frame1_fraction + frame2_fraction]
  out[, dominant_cds_frame := {
    frac <- as.matrix(.SD)
    paste0("frame", max.col(frac, ties.method = "first") - 1L)
  }, .SDcols = paste0("frame", 0:2, "_fraction")]
  out[, dominant_cds_frame_fraction := pmax(frame0_fraction, frame1_fraction,
                                            frame2_fraction, na.rm = TRUE)]
  out[, top1pct_canonical_fraction := top1pct_frame0_fraction]
  out[, top1pct_off_frame_fraction :=
        top1pct_frame1_fraction + top1pct_frame2_fraction]
  out[, possible_library_misshift_for_gene :=
        cds_total_counts >= 100 &
        is.finite(canonical_frame_fraction) &
        canonical_frame_fraction < 0.45 &
        dominant_cds_frame != "frame0"]
  out[]
}

aggregate_profile <- function(mat, runs, positions, cds_positions, gene, group) {
  runs <- intersect(runs, colnames(mat))
  positions <- sort(unique(as.integer(positions)))
  positions <- positions[
    is.finite(positions) & positions >= 1L & positions <= nrow(mat)
  ]
  if (!length(runs) || !length(positions)) return(data.table())
  counts <- as.numeric(rowSums(mat[positions, runs, drop = FALSE], na.rm = TRUE))
  total <- sum(counts, na.rm = TRUE)
  data.table(
    gene_symbol = gene,
    profile_group = group,
    n_runs = length(runs),
    cds_relative_position = seq_along(positions),
    transcript_position = positions,
    cds_frame = paste0("frame", (positions - min(cds_positions)) %% 3L),
    counts = counts,
    normalized_density = if (is.finite(total) && total > 0) {
      counts / total * 1000
    } else {
      NA_real_
    }
  )
}

cds_hotspot_profile <- function(mat, runs, cds_positions, gene, group) {
  runs <- intersect(runs, colnames(mat))
  positions <- sort(unique(as.integer(cds_positions)))
  positions <- positions[
    is.finite(positions) & positions >= 1L & positions <= nrow(mat)
  ]
  if (!length(runs) || !length(positions)) return(data.table())
  counts <- as.numeric(rowSums(mat[positions, runs, drop = FALSE],
                               na.rm = TRUE))
  total <- sum(counts, na.rm = TRUE)
  rel_pos <- seq_along(positions)
  out <- data.table(
    gene_symbol = gene,
    hotspot_group = group,
    n_runs = length(runs),
    cds_relative_position = rel_pos,
    codon_index = floor((rel_pos - 1L) / 3L) + 1L,
    transcript_position = positions,
    cds_frame = paste0("frame", (positions - min(positions)) %% 3L),
    counts = counts,
    normalized_density = if (is.finite(total) && total > 0) {
      counts / total * 1000
    } else {
      NA_real_
    }
  )
  if (is.finite(total) && total > 0) {
    setorder(out, -normalized_density, cds_relative_position)
    out[, hotspot_rank := seq_len(.N)]
    setorder(out, cds_relative_position)
  } else {
    out[, hotspot_rank := NA_integer_]
  }
  out[]
}

per_run_cds_hotspots <- function(mat, runs, cds_positions, gene,
                                 group_lookup,
                                 top_fraction = 0.01,
                                 max_positions_per_run = 20L) {
  runs <- intersect(runs, colnames(mat))
  positions <- sort(unique(as.integer(cds_positions)))
  positions <- positions[
    is.finite(positions) & positions >= 1L & positions <= nrow(mat)
  ]
  if (!length(runs) || !length(positions)) return(data.table())
  rbindlist(lapply(runs, function(run) {
    counts <- as.numeric(mat[positions, run, drop = TRUE])
    total <- sum(counts, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) return(data.table())
    n_top <- max(1L, ceiling(length(positions) * top_fraction))
    ord <- order(counts, decreasing = TRUE, na.last = TRUE)
    ord <- ord[seq_len(min(length(ord), max(n_top, max_positions_per_run)))]
    ord <- ord[is.finite(counts[ord]) & counts[ord] > 0]
    if (!length(ord)) return(data.table())
    ord <- ord[seq_len(min(length(ord), max_positions_per_run))]
    rel_pos <- ord
    data.table(
      gene_symbol = gene,
      Run = run,
      infected_library_qc_class = group_lookup[run],
      cds_relative_position = rel_pos,
      codon_index = floor((rel_pos - 1L) / 3L) + 1L,
      transcript_position = positions[ord],
      cds_frame = paste0("frame", (positions[ord] - min(positions)) %% 3L),
      position_count = counts[ord],
      position_count_share = counts[ord] / total,
      hotspot_rank_in_run = seq_along(ord),
      run_cds_total_counts = total
    )
  }), fill = TRUE)
}

summarize_recurrent_hotspots <- function(per_run_hotspots, group_name,
                                         group_runs) {
  group_runs <- unique(group_runs)
  dt <- per_run_hotspots[Run %chin% group_runs]
  if (!length(group_runs) || nrow(dt) == 0) return(data.table())
  dt[
    ,
    .(
      group_n_runs = length(group_runs),
      n_runs_with_hotspot = uniqueN(Run),
      fraction_runs_with_hotspot = uniqueN(Run) / length(group_runs),
      total_position_counts = sum(position_count, na.rm = TRUE),
      median_position_count_share =
        median(position_count_share, na.rm = TRUE),
      best_rank_in_run = min(hotspot_rank_in_run, na.rm = TRUE)
    ),
    by = .(gene_symbol, cds_relative_position, codon_index,
           transcript_position, cds_frame)
  ][
    ,
    hotspot_group := group_name
  ][
    order(-fraction_runs_with_hotspot, -total_position_counts,
          best_rank_in_run)
  ]
}

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
run_order <- runIDs(df)
# collection_dir_from_exp() returns RiboCrypt's generic "collection_tables"
# dir; this project's own indexed coverage-page cache (see
# dominant_cell_state_pack_fst_pages.R's documented page_source_dir) is one
# level further, at "<that>_indexed".
fst_index <- file.path(paste0(collection_dir_from_exp(df), "_indexed"), "coverage_index.fst")
if (!file.exists(fst_index)) stop("Missing FST coverage index: ", fst_index)

diagnostics <- fread(diagnostics_file)
diagnostics <- diagnostics[
  expression_source == "fst_clean_cds" & !is.na(tx_id) & nzchar(tx_id)
]
diagnostics[, gene_symbol := toupper(gene_symbol)]
setorder(diagnostics, gene_symbol)

metadata <- if (file.exists(metadata_file)) {
  fread(metadata_file, showProgress = FALSE)
} else {
  warning("Metadata file unavailable: ", metadata_file)
  data.table(Run = run_order)
}
metadata <- unique(metadata, by = "Run")
metadata <- add_metadata_frame_qc(metadata)

feature_contrasts <- fread(feature_contrasts_file, showProgress = FALSE)
feature_contrasts[, gene_symbol := toupper(gene_symbol)]
viral_clean <- feature_contrasts[
  design_family == "viral_infection" &
    feature_id == "clean_CDS" &
    gene_symbol %chin% target_genes
]
viral_case_union <- unique(unlist(lapply(
  viral_clean$case_run_signature,
  parse_run_signature
), use.names = FALSE))
viral_control_union <- unique(unlist(lapply(
  viral_clean$control_run_signature,
  parse_run_signature
), use.names = FALSE))

target_tx <- unique(diagnostics$tx_id)
cds_all <- loadRegion(df, part = "cds", names.keep = target_tx)
leader_all <- loadRegion(df, part = "leaders", names.keep = target_tx)

global_frame_counts <- matrix(0, nrow = length(run_order), ncol = 3,
                              dimnames = list(run_order, paste0("frame", 0:2)))
run_metrics <- vector("list", nrow(diagnostics))
target_profiles <- list()
hotspot_genes <- c("MYC", "RRM1")
hotspot_mats <- list()
hotspot_cds_positions <- list()
hotspot_tx_ids <- list()

for (i in seq_len(nrow(diagnostics))) {
  row <- diagnostics[i]
  gene <- row$gene_symbol
  tx_id <- row$tx_id
  if (!tx_id %in% names(cds_all)) next
  cds_tx <- cds_all[tx_id]
  leader_tx <- if (tx_id %in% names(leader_all)) {
    leader_all[tx_id]
  } else {
    GRangesList()
  }
  display_region <- make_display_region(leader_tx, cds_tx, tx_id)
  cds_positions <- safe_pmap_positions(cds_tx, display_region)
  mid_positions <- middle_cds_positions(cds_positions)
  if (!length(mid_positions)) next
  message("Measuring CDS frame/bumpiness for ", gene, " (", tx_id, ")")
  coverage <- dominant_cached_coverage_by_transcript(
    display_region = display_region,
    fst_index = fst_index,
    gene_symbol = gene,
    tx_id = tx_id,
    analysis_dir = analysis_dir
  )
  mat <- as.matrix(coverage)
  storage.mode(mat) <- "numeric"
  if (!identical(colnames(mat), run_order)) {
    mat <- mat[, run_order, drop = FALSE]
  }
  metric <- coverage_metrics(mat, cds_positions, mid_positions, gene, tx_id)
  run_metrics[[i]] <- metric

  phase <- as.integer((mid_positions - min(cds_positions)) %% 3L)
  sub <- mat[mid_positions, , drop = FALSE]
  for (frame in 0:2) {
    global_frame_counts[, frame + 1L] <-
      global_frame_counts[, frame + 1L] +
      as.numeric(colSums(sub[phase == frame, , drop = FALSE], na.rm = TRUE))
  }

  if (gene %chin% target_genes) {
    if (gene %chin% hotspot_genes) {
      hotspot_mats[[gene]] <- mat
      hotspot_cds_positions[[gene]] <- cds_positions
      hotspot_tx_ids[[gene]] <- tx_id
    }
    gene_viral <- viral_clean[gene_symbol == gene]
    viral_case_runs <- unique(unlist(lapply(
      gene_viral$case_run_signature,
      parse_run_signature
    ), use.names = FALSE))
    viral_control_runs <- unique(unlist(lapply(
      gene_viral$control_run_signature,
      parse_run_signature
    ), use.names = FALSE))
    target_profiles[[gene]] <- rbindlist(list(
      aggregate_profile(mat, viral_case_runs, cds_positions, cds_positions,
                        gene, "viral_case_union"),
      aggregate_profile(mat, viral_control_runs, cds_positions, cds_positions,
                        gene, "viral_control_union"),
      aggregate_profile(mat, run_order, cds_positions, cds_positions,
                        gene, "all_3857_runs")
    ), fill = TRUE)
  }
}

run_metrics <- rbindlist(run_metrics, fill = TRUE)
global_totals <- rowSums(global_frame_counts, na.rm = TRUE)
global_qc <- data.table(
  Run = rownames(global_frame_counts),
  global_cds_middle_counts = global_totals,
  global_frame0_counts = global_frame_counts[, 1],
  global_frame1_counts = global_frame_counts[, 2],
  global_frame2_counts = global_frame_counts[, 3],
  global_frame0_fraction = global_frame_counts[, 1] / pmax(global_totals, 1),
  global_frame1_fraction = global_frame_counts[, 2] / pmax(global_totals, 1),
  global_frame2_fraction = global_frame_counts[, 3] / pmax(global_totals, 1)
)
global_qc[global_cds_middle_counts <= 0,
          c("global_frame0_fraction", "global_frame1_fraction",
            "global_frame2_fraction") := .(NA_real_, NA_real_, NA_real_)]
global_qc[, global_dominant_frame := {
  frac <- as.matrix(.SD)
  paste0("frame", max.col(frac, ties.method = "first") - 1L)
}, .SDcols = paste0("global_frame", 0:2, "_fraction")]
global_qc[, global_canonical_frame_fraction := global_frame0_fraction]
global_qc[, global_possible_misshift :=
            global_cds_middle_counts >= 1000 &
            is.finite(global_canonical_frame_fraction) &
            global_canonical_frame_fraction < 0.45 &
            global_dominant_frame != "frame0"]

run_metrics <- merge(run_metrics, global_qc, by = "Run", all.x = TRUE,
                     sort = FALSE)
run_metrics <- merge(run_metrics, metadata, by = "Run", all.x = TRUE,
                     sort = FALSE)

run_metrics[, log10_cds_total_counts := log10(cds_total_counts + 1)]
eligible_background <- run_metrics[
  cds_total_counts >= 30 &
    is.finite(cds_top1pct_share) &
    is.finite(log10_cds_total_counts)
]
if (nrow(eligible_background) > 100) {
  coverage_breaks <- unique(quantile(
    eligible_background$log10_cds_total_counts,
    probs = seq(0, 1, by = 0.1),
    na.rm = TRUE,
    names = FALSE
  ))
  eligible_background[, coverage_bin := cut(
    log10_cds_total_counts,
    breaks = coverage_breaks,
    include.lowest = TRUE
  )]
  bin_ref <- eligible_background[
    ,
    .(
      bin_n = .N,
      bin_median_top1pct_share = median(cds_top1pct_share, na.rm = TRUE),
      bin_mad_top1pct_share = mad(cds_top1pct_share, na.rm = TRUE),
      bin_median_cv = median(cds_coverage_cv, na.rm = TRUE),
      bin_mad_cv = mad(cds_coverage_cv, na.rm = TRUE)
    ),
    by = coverage_bin
  ]
  run_metrics[, coverage_bin := cut(
    log10_cds_total_counts,
    breaks = coverage_breaks,
    include.lowest = TRUE
  )]
  run_metrics <- merge(run_metrics, bin_ref, by = "coverage_bin",
                       all.x = TRUE, sort = FALSE)
  run_metrics[, top1pct_bumpiness_z_vs_depth :=
                (cds_top1pct_share - bin_median_top1pct_share) /
                fifelse(is.finite(bin_mad_top1pct_share) &
                          bin_mad_top1pct_share > 0,
                        bin_mad_top1pct_share, NA_real_)]
  run_metrics[, cv_bumpiness_z_vs_depth :=
                (cds_coverage_cv - bin_median_cv) /
                fifelse(is.finite(bin_mad_cv) & bin_mad_cv > 0,
                        bin_mad_cv, NA_real_)]
}

housekeeping_gene_summary <- run_metrics[
  gene_symbol %chin% housekeeping_candidate_genes,
  .(
    tx_id = tx_id[1L],
    n_runs_total = .N,
    n_runs_ge_30 = sum(cds_total_counts >= 30, na.rm = TRUE),
    frac_runs_ge_30 = sum(cds_total_counts >= 30, na.rm = TRUE) / .N,
    median_cds_total_counts = median(cds_total_counts, na.rm = TRUE),
    median_canonical_frame_fraction =
      median(canonical_frame_fraction, na.rm = TRUE),
    median_top1pct_share = median(cds_top1pct_share, na.rm = TRUE),
    median_top1pct_off_frame_fraction =
      median(top1pct_off_frame_fraction, na.rm = TRUE),
    median_top1pct_bumpiness_z_vs_depth =
      median(top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
    n_global_misshift_candidates = sum(global_possible_misshift == TRUE,
                                       na.rm = TRUE),
    n_gene_misshift_candidates =
      sum(possible_library_misshift_for_gene == TRUE, na.rm = TRUE)
  ),
  by = gene_symbol
]
housekeeping_gene_summary[
  ,
  selected_for_housekeeping_baseline :=
    frac_runs_ge_30 >= 0.50 & median_cds_total_counts >= 30
]
if (sum(housekeeping_gene_summary$selected_for_housekeeping_baseline,
        na.rm = TRUE) < 6) {
  setorder(housekeeping_gene_summary, -n_runs_ge_30,
           -median_cds_total_counts)
  housekeeping_gene_summary[
    seq_len(min(.N, 8L)),
    selected_for_housekeeping_baseline := TRUE
  ]
}
housekeeping_genes <- housekeeping_gene_summary[
  selected_for_housekeeping_baseline == TRUE,
  gene_symbol
]
housekeeping_gene_summary[, baseline_role := fifelse(
  selected_for_housekeeping_baseline == TRUE,
  "used_housekeeping_baseline",
  "candidate_not_used"
)]

housekeeping_run_metrics <- run_metrics[
  gene_symbol %chin% housekeeping_genes &
    cds_total_counts >= 30
]
housekeeping_qc <- housekeeping_run_metrics[
  ,
  .(
    housekeeping_n_genes = .N,
    housekeeping_genes = paste(sort(unique(gene_symbol)), collapse = ";"),
    housekeeping_total_counts = sum(cds_total_counts, na.rm = TRUE),
    housekeeping_frame0_counts = sum(frame0_counts, na.rm = TRUE),
    housekeeping_frame1_counts = sum(frame1_counts, na.rm = TRUE),
    housekeeping_frame2_counts = sum(frame2_counts, na.rm = TRUE),
    housekeeping_median_canonical_frame_fraction =
      median(canonical_frame_fraction, na.rm = TRUE),
    housekeeping_median_top1pct_share =
      median(cds_top1pct_share, na.rm = TRUE),
    housekeeping_median_top1pct_off_frame_fraction =
      median(top1pct_off_frame_fraction, na.rm = TRUE),
    housekeeping_median_top1pct_bumpiness_z_vs_depth =
      median(top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
    housekeeping_median_coverage_cv = median(cds_coverage_cv, na.rm = TRUE),
    housekeeping_n_gene_misshift_candidates =
      sum(possible_library_misshift_for_gene == TRUE, na.rm = TRUE)
  ),
  by = Run
]
if (nrow(housekeeping_qc) > 0) {
  housekeeping_qc[
    ,
    `:=`(
      housekeeping_frame0_fraction =
        housekeeping_frame0_counts / pmax(housekeeping_total_counts, 1),
      housekeeping_frame1_fraction =
        housekeeping_frame1_counts / pmax(housekeeping_total_counts, 1),
      housekeeping_frame2_fraction =
        housekeeping_frame2_counts / pmax(housekeeping_total_counts, 1)
    )
  ]
  housekeeping_qc[
    housekeeping_total_counts <= 0,
    c("housekeeping_frame0_fraction", "housekeeping_frame1_fraction",
      "housekeeping_frame2_fraction") := .(NA_real_, NA_real_, NA_real_)
  ]
  housekeeping_qc[
    ,
    housekeeping_dominant_frame := {
      frac <- as.matrix(.SD)
      paste0("frame", max.col(frac, ties.method = "first") - 1L)
    },
    .SDcols = paste0("housekeeping_frame", 0:2, "_fraction")
  ]
  housekeeping_qc[
    ,
    housekeeping_canonical_frame_fraction := housekeeping_frame0_fraction
  ]
  housekeeping_qc[
    ,
    housekeeping_possible_misshift :=
      housekeeping_n_genes >= 4 &
      housekeeping_total_counts >= 1000 &
      is.finite(housekeeping_canonical_frame_fraction) &
      housekeeping_canonical_frame_fraction < 0.45 &
      housekeeping_dominant_frame != "frame0"
  ]
  housekeeping_qc <- merge(
    housekeeping_qc,
    global_qc[, .(Run, global_canonical_frame_fraction,
                  global_possible_misshift)],
    by = "Run",
    all.x = TRUE,
    sort = FALSE
  )
  if (nrow(metadata) > 0) {
    metadata_cols <- intersect(
      c("Run", "study", "STUDY", "BioProject", "PRJ", "CELL_LINE", "TISSUE",
        "CONDITION", "GENE", "AUTHOR", "metadata_frame0_fraction",
        "metadata_frame1_fraction", "metadata_frame2_fraction",
        "metadata_canonical_frame_fraction", "metadata_dominant_frame",
        "metadata_possible_misshift"),
      names(metadata)
    )
    housekeeping_qc <- merge(
      housekeeping_qc,
      metadata[, ..metadata_cols],
      by = "Run",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

target_run_metrics <- run_metrics[gene_symbol %chin% target_genes]

viral_design_summaries <- rbindlist(lapply(seq_len(nrow(viral_clean)), function(i) {
  row <- viral_clean[i]
  case_runs <- parse_run_signature(row$case_run_signature)
  control_runs <- parse_run_signature(row$control_run_signature)
  dt <- target_run_metrics[gene_symbol == row$gene_symbol]
  case_dt <- dt[Run %chin% case_runs]
  control_dt <- dt[Run %chin% control_runs]
  summarize_group <- function(x, prefix) {
    setNames(
      list(
        uniqueN(x$Run),
        sum(x$cds_total_counts, na.rm = TRUE),
        median(x$canonical_frame_fraction, na.rm = TRUE),
        median(x$off_frame_fraction, na.rm = TRUE),
        median(x$global_canonical_frame_fraction, na.rm = TRUE),
        median(x$cds_top1pct_share, na.rm = TRUE),
        median(x$top1pct_off_frame_fraction, na.rm = TRUE),
        median(x$cds_coverage_cv, na.rm = TRUE),
        median(x$top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
        sum(x$global_possible_misshift == TRUE, na.rm = TRUE),
        sum(x$possible_library_misshift_for_gene == TRUE, na.rm = TRUE)
      ),
      paste0(prefix, c(
        "_n_runs", "_cds_total_counts", "_median_canonical_frame_fraction",
        "_median_off_frame_fraction", "_median_global_canonical_frame_fraction",
        "_median_top1pct_share", "_median_top1pct_off_frame_fraction",
        "_median_coverage_cv", "_median_top1pct_bumpiness_z_vs_depth",
        "_n_global_misshift_candidates", "_n_gene_misshift_candidates"
      ))
    )
  }
  out <- data.table(
    gene_symbol = row$gene_symbol,
    tx_id = row$tx_id,
    gene_role = fifelse(row$gene_symbol %chin% candidate_genes,
                        "candidate_bumpy_cds",
                        fifelse(row$gene_symbol %chin% control_genes,
                                "browser_clean_control", "other")),
    study = row$study,
    AUTHOR = row$AUTHOR,
    CELL_LINE = row$CELL_LINE,
    TISSUE = row$TISSUE,
    case_condition = row$case_condition,
    control_conditions = row$control_conditions
  )
  out <- cbind(out, as.data.table(summarize_group(case_dt, "case")),
               as.data.table(summarize_group(control_dt, "control")))
  out[, delta_median_canonical_frame_fraction :=
        case_median_canonical_frame_fraction -
        control_median_canonical_frame_fraction]
  out[, delta_median_global_canonical_frame_fraction :=
        case_median_global_canonical_frame_fraction -
        control_median_global_canonical_frame_fraction]
  out[, delta_median_top1pct_share :=
        case_median_top1pct_share - control_median_top1pct_share]
  out[, delta_median_top1pct_off_frame_fraction :=
        case_median_top1pct_off_frame_fraction -
        control_median_top1pct_off_frame_fraction]
  out[, delta_median_top1pct_bumpiness_z_vs_depth :=
        case_median_top1pct_bumpiness_z_vs_depth -
        control_median_top1pct_bumpiness_z_vs_depth]
  out[, case_vs_control_top1pct_wilcox_p := if (
    nrow(case_dt) >= 2 && nrow(control_dt) >= 2
  ) {
    suppressWarnings(wilcox.test(case_dt$cds_top1pct_share,
                                 control_dt$cds_top1pct_share)$p.value)
  } else NA_real_]
  out[, case_vs_control_canonical_frame_wilcox_p := if (
    nrow(case_dt) >= 2 && nrow(control_dt) >= 2
  ) {
    suppressWarnings(wilcox.test(case_dt$canonical_frame_fraction,
                                 control_dt$canonical_frame_fraction)$p.value)
  } else NA_real_]
  out
}), fill = TRUE)
viral_design_summaries[, case_vs_control_top1pct_wilcox_q :=
                         p.adjust(case_vs_control_top1pct_wilcox_p, "BH")]
viral_design_summaries[, case_vs_control_canonical_frame_wilcox_q :=
                         p.adjust(case_vs_control_canonical_frame_wilcox_p, "BH")]

target_summary <- target_run_metrics[
  cds_total_counts >= 30,
  .(
    n_runs = .N,
    median_cds_total_counts = median(cds_total_counts, na.rm = TRUE),
    median_canonical_frame_fraction =
      median(canonical_frame_fraction, na.rm = TRUE),
    median_global_canonical_frame_fraction =
      median(global_canonical_frame_fraction, na.rm = TRUE),
    median_top1pct_share = median(cds_top1pct_share, na.rm = TRUE),
    median_top1pct_off_frame_fraction =
      median(top1pct_off_frame_fraction, na.rm = TRUE),
    median_top1pct_bumpiness_z_vs_depth =
      median(top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
    n_global_misshift_candidates = sum(global_possible_misshift == TRUE,
                                       na.rm = TRUE),
    n_gene_misshift_candidates = sum(possible_library_misshift_for_gene == TRUE,
                                     na.rm = TRUE)
  ),
  by = gene_symbol
]
target_summary[, gene_role := fifelse(
  gene_symbol %chin% candidate_genes,
  "candidate_bumpy_cds",
  fifelse(gene_symbol %chin% control_genes, "browser_clean_control", "other")
)]

shared_viral_libraries <- dcast(
  target_run_metrics[
    gene_symbol %chin% target_genes &
      Run %chin% viral_case_union &
      cds_total_counts >= 30,
    .(
      Run,
      gene_symbol,
      cds_total_counts,
      canonical_frame_fraction,
      cds_top1pct_share,
      top1pct_off_frame_fraction,
      top1pct_bumpiness_z_vs_depth,
      global_canonical_frame_fraction,
      global_possible_misshift
    )
  ],
  Run ~ gene_symbol,
  value.var = c(
    "cds_total_counts", "canonical_frame_fraction",
    "cds_top1pct_share", "top1pct_off_frame_fraction",
    "top1pct_bumpiness_z_vs_depth"
  )
)
if (nrow(shared_viral_libraries) > 0) {
  shared_global <- unique(target_run_metrics[
    Run %chin% shared_viral_libraries$Run,
    .(Run, global_canonical_frame_fraction, global_possible_misshift)
  ], by = "Run")
  shared_viral_libraries <- merge(shared_viral_libraries, shared_global,
                                  by = "Run", all.x = TRUE, sort = FALSE)
}

shared_ddit3_controls <- rbindlist(lapply(candidate_genes, function(gene) {
  candidate_dt <- target_run_metrics[
    gene_symbol == gene &
      Run %chin% viral_case_union &
      cds_total_counts >= 30,
    .(
      Run,
      candidate_cds_total_counts = cds_total_counts,
      candidate_canonical_frame_fraction = canonical_frame_fraction,
      candidate_top1pct_share = cds_top1pct_share,
      candidate_top1pct_off_frame_fraction = top1pct_off_frame_fraction,
      candidate_top1pct_bumpiness_z_vs_depth =
        top1pct_bumpiness_z_vs_depth,
      candidate_possible_library_misshift_for_gene =
        possible_library_misshift_for_gene,
      global_canonical_frame_fraction,
      global_possible_misshift
    )
  ]
  ddit3_dt <- target_run_metrics[
    gene_symbol == "DDIT3" &
      Run %chin% viral_case_union &
      cds_total_counts >= 30,
    .(
      Run,
      ddit3_cds_total_counts = cds_total_counts,
      ddit3_canonical_frame_fraction = canonical_frame_fraction,
      ddit3_top1pct_share = cds_top1pct_share,
      ddit3_top1pct_off_frame_fraction = top1pct_off_frame_fraction,
      ddit3_top1pct_bumpiness_z_vs_depth =
        top1pct_bumpiness_z_vs_depth,
      ddit3_possible_library_misshift_for_gene =
        possible_library_misshift_for_gene
    )
  ]
  out <- merge(candidate_dt, ddit3_dt, by = "Run", all = FALSE,
               sort = FALSE)
  if (nrow(out) == 0) return(data.table())
  out[, candidate_gene := gene]
  out[, candidate_minus_ddit3_top1pct_share :=
        candidate_top1pct_share - ddit3_top1pct_share]
  out[, candidate_minus_ddit3_top1pct_off_frame_fraction :=
        candidate_top1pct_off_frame_fraction -
        ddit3_top1pct_off_frame_fraction]
  out[, candidate_minus_ddit3_top1pct_bumpiness_z_vs_depth :=
        candidate_top1pct_bumpiness_z_vs_depth -
        ddit3_top1pct_bumpiness_z_vs_depth]
  out[, candidate_minus_ddit3_canonical_frame_fraction :=
        candidate_canonical_frame_fraction -
        ddit3_canonical_frame_fraction]
  setcolorder(out, c("candidate_gene", setdiff(names(out), "candidate_gene")))
  out[]
}), fill = TRUE)

if (nrow(shared_ddit3_controls) > 0 && nrow(metadata) > 0) {
  metadata_cols <- intersect(
    c("Run", "study", "STUDY", "BioProject", "PRJ", "CELL_LINE", "TISSUE",
      "CONDITION", "GENE", "AUTHOR"),
    names(metadata)
  )
  shared_ddit3_controls <- merge(
    shared_ddit3_controls,
    metadata[, ..metadata_cols],
    by = "Run",
    all.x = TRUE,
    sort = FALSE
  )
  setcolorder(shared_ddit3_controls,
              c("candidate_gene", "Run",
                setdiff(names(shared_ddit3_controls),
                        c("candidate_gene", "Run"))))
}

shared_ddit3_summary <- if (nrow(shared_ddit3_controls) > 0) {
  shared_ddit3_controls[
    ,
    .(
      shared_n_runs = .N,
      median_candidate_cds_total_counts =
        median(candidate_cds_total_counts, na.rm = TRUE),
      median_ddit3_cds_total_counts =
        median(ddit3_cds_total_counts, na.rm = TRUE),
      median_candidate_canonical_frame_fraction =
        median(candidate_canonical_frame_fraction, na.rm = TRUE),
      median_ddit3_canonical_frame_fraction =
        median(ddit3_canonical_frame_fraction, na.rm = TRUE),
      median_candidate_top1pct_share =
        median(candidate_top1pct_share, na.rm = TRUE),
      median_ddit3_top1pct_share =
        median(ddit3_top1pct_share, na.rm = TRUE),
      median_candidate_top1pct_off_frame_fraction =
        median(candidate_top1pct_off_frame_fraction, na.rm = TRUE),
      median_ddit3_top1pct_off_frame_fraction =
        median(ddit3_top1pct_off_frame_fraction, na.rm = TRUE),
      median_candidate_top1pct_bumpiness_z_vs_depth =
        median(candidate_top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
      median_ddit3_top1pct_bumpiness_z_vs_depth =
        median(ddit3_top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
      median_candidate_minus_ddit3_top1pct_bumpiness_z_vs_depth =
        median(candidate_minus_ddit3_top1pct_bumpiness_z_vs_depth,
               na.rm = TRUE),
      n_candidate_bumpiness_z_ge_1 =
        sum(candidate_top1pct_bumpiness_z_vs_depth >= 1, na.rm = TRUE),
      n_ddit3_bumpiness_z_ge_1 =
        sum(ddit3_top1pct_bumpiness_z_vs_depth >= 1, na.rm = TRUE),
      n_candidate_outpaces_ddit3_by_z_ge_1 =
        sum(candidate_minus_ddit3_top1pct_bumpiness_z_vs_depth >= 1,
            na.rm = TRUE),
      n_global_misshift_candidates =
        sum(global_possible_misshift == TRUE, na.rm = TRUE),
      median_global_canonical_frame_fraction =
        median(global_canonical_frame_fraction, na.rm = TRUE)
    ),
    by = candidate_gene
  ]
} else {
  data.table()
}

make_target_housekeeping_controls <- function(run_filter = NULL) {
  if (nrow(housekeeping_qc) == 0) return(data.table())
  rbindlist(lapply(target_genes, function(gene) {
    target_dt <- target_run_metrics[
      gene_symbol == gene &
        cds_total_counts >= 30,
      .(
        Run,
        target_gene = gene_symbol,
        target_cds_total_counts = cds_total_counts,
        target_frame0_counts = frame0_counts,
        target_frame1_counts = frame1_counts,
        target_frame2_counts = frame2_counts,
        target_frame0_fraction = frame0_fraction,
        target_frame1_fraction = frame1_fraction,
        target_frame2_fraction = frame2_fraction,
        target_canonical_frame_fraction = canonical_frame_fraction,
        target_off_frame_fraction = off_frame_fraction,
        target_dominant_cds_frame = dominant_cds_frame,
        target_dominant_cds_frame_fraction = dominant_cds_frame_fraction,
        target_top1pct_share = cds_top1pct_share,
        target_top1pct_off_frame_fraction = top1pct_off_frame_fraction,
        target_top1pct_bumpiness_z_vs_depth =
          top1pct_bumpiness_z_vs_depth,
        target_possible_library_misshift_for_gene =
          possible_library_misshift_for_gene,
        metadata_canonical_frame_fraction,
        metadata_dominant_frame,
        metadata_possible_misshift,
        global_canonical_frame_fraction,
        global_possible_misshift
      )
    ]
    if (!is.null(run_filter)) {
      target_dt <- target_dt[Run %chin% run_filter]
    }
    hk_dt <- housekeeping_qc[
      housekeeping_n_genes >= 4 &
        housekeeping_total_counts >= 1000,
      .(
        Run,
        housekeeping_n_genes,
        housekeeping_total_counts,
        housekeeping_genes,
        housekeeping_canonical_frame_fraction,
        housekeeping_median_canonical_frame_fraction,
        housekeeping_median_top1pct_share,
        housekeeping_median_top1pct_off_frame_fraction,
        housekeeping_median_top1pct_bumpiness_z_vs_depth,
        housekeeping_median_coverage_cv,
        housekeeping_n_gene_misshift_candidates,
        housekeeping_possible_misshift
      )
    ]
    out <- merge(target_dt, hk_dt, by = "Run", all = FALSE, sort = FALSE)
    if (nrow(out) == 0) return(data.table())
    out[, target_minus_housekeeping_canonical_frame_fraction :=
          target_canonical_frame_fraction -
          housekeeping_canonical_frame_fraction]
    out[, target_minus_housekeeping_median_canonical_frame_fraction :=
          target_canonical_frame_fraction -
          housekeeping_median_canonical_frame_fraction]
    out[, target_minus_housekeeping_top1pct_share :=
          target_top1pct_share - housekeeping_median_top1pct_share]
    out[, target_minus_housekeeping_top1pct_off_frame_fraction :=
          target_top1pct_off_frame_fraction -
          housekeeping_median_top1pct_off_frame_fraction]
    out[, target_minus_housekeeping_top1pct_bumpiness_z_vs_depth :=
          target_top1pct_bumpiness_z_vs_depth -
          housekeeping_median_top1pct_bumpiness_z_vs_depth]
    out[]
  }), fill = TRUE)
}

target_housekeeping_all_controls <- make_target_housekeeping_controls()

if (nrow(target_housekeeping_all_controls) > 0 && nrow(metadata) > 0) {
  metadata_cols <- intersect(
    c("Run", "study", "STUDY", "BioProject", "PRJ", "CELL_LINE", "TISSUE",
      "CONDITION", "GENE", "AUTHOR"),
    names(metadata)
  )
  target_housekeeping_all_controls <- merge(
    target_housekeeping_all_controls,
    metadata[, ..metadata_cols],
    by = "Run",
    all.x = TRUE,
    sort = FALSE
  )
  setcolorder(target_housekeeping_all_controls,
              c("target_gene", "Run",
                setdiff(names(target_housekeeping_all_controls),
                        c("target_gene", "Run"))))
}

target_housekeeping_controls <- target_housekeeping_all_controls[
  Run %chin% viral_case_union
]

target_housekeeping_summary <- if (nrow(target_housekeeping_controls) > 0) {
  target_housekeeping_controls[
    ,
    .(
      shared_n_runs = .N,
      median_target_cds_total_counts =
        median(target_cds_total_counts, na.rm = TRUE),
      median_housekeeping_n_genes =
        median(housekeeping_n_genes, na.rm = TRUE),
      median_housekeeping_total_counts =
        median(housekeeping_total_counts, na.rm = TRUE),
      median_target_canonical_frame_fraction =
        median(target_canonical_frame_fraction, na.rm = TRUE),
      median_housekeeping_canonical_frame_fraction =
        median(housekeeping_canonical_frame_fraction, na.rm = TRUE),
      median_target_top1pct_share =
        median(target_top1pct_share, na.rm = TRUE),
      median_housekeeping_top1pct_share =
        median(housekeeping_median_top1pct_share, na.rm = TRUE),
      median_target_top1pct_off_frame_fraction =
        median(target_top1pct_off_frame_fraction, na.rm = TRUE),
      median_housekeeping_top1pct_off_frame_fraction =
        median(housekeeping_median_top1pct_off_frame_fraction, na.rm = TRUE),
      median_target_top1pct_bumpiness_z_vs_depth =
        median(target_top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
      median_housekeeping_top1pct_bumpiness_z_vs_depth =
        median(housekeeping_median_top1pct_bumpiness_z_vs_depth,
               na.rm = TRUE),
      median_target_minus_housekeeping_top1pct_bumpiness_z_vs_depth =
        median(target_minus_housekeeping_top1pct_bumpiness_z_vs_depth,
               na.rm = TRUE),
      n_target_bumpiness_z_ge_1 =
        sum(target_top1pct_bumpiness_z_vs_depth >= 1, na.rm = TRUE),
      n_housekeeping_median_bumpiness_z_ge_1 =
        sum(housekeeping_median_top1pct_bumpiness_z_vs_depth >= 1,
            na.rm = TRUE),
      n_target_outpaces_housekeeping_by_z_ge_1 =
        sum(target_minus_housekeeping_top1pct_bumpiness_z_vs_depth >= 1,
            na.rm = TRUE),
      n_global_misshift_candidates =
        sum(global_possible_misshift == TRUE, na.rm = TRUE),
      n_housekeeping_misshift_candidates =
        sum(housekeeping_possible_misshift == TRUE, na.rm = TRUE),
      n_target_gene_misshift_candidates =
        sum(target_possible_library_misshift_for_gene == TRUE,
            na.rm = TRUE)
    ),
    by = target_gene
  ]
} else {
  data.table()
}

add_ddit3_and_library_qc <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  ddit3_by_run <- target_run_metrics[
    gene_symbol == "DDIT3",
    .(
      Run,
      ddit3_cds_total_counts = cds_total_counts,
      ddit3_canonical_frame_fraction = canonical_frame_fraction,
      ddit3_top1pct_share = cds_top1pct_share,
      ddit3_top1pct_off_frame_fraction = top1pct_off_frame_fraction,
      ddit3_top1pct_bumpiness_z_vs_depth =
        top1pct_bumpiness_z_vs_depth,
      ddit3_possible_library_misshift_for_gene =
        possible_library_misshift_for_gene
    )
  ]
  dt <- merge(
    dt,
    ddit3_by_run,
    by = "Run",
    all.x = TRUE,
    sort = FALSE
  )
  dt[
    ,
    has_ddit3_fair_control :=
      !is.na(ddit3_cds_total_counts) & ddit3_cds_total_counts >= 30
  ]
  dt[
    ,
    target_minus_ddit3_top1pct_bumpiness_z_vs_depth :=
      target_top1pct_bumpiness_z_vs_depth -
      ddit3_top1pct_bumpiness_z_vs_depth
  ]
  dt[
    ,
    target_minus_ddit3_canonical_frame_fraction :=
      target_canonical_frame_fraction - ddit3_canonical_frame_fraction
  ]
  dt[
    ,
    target_outpaces_ddit3_by_z_ge_1 :=
      has_ddit3_fair_control == TRUE &
      is.finite(target_minus_ddit3_top1pct_bumpiness_z_vs_depth) &
      target_minus_ddit3_top1pct_bumpiness_z_vs_depth >= 1
  ]
  dt[
    ,
    technical_frame_bias_possible :=
      metadata_possible_misshift == TRUE |
      global_possible_misshift == TRUE |
      housekeeping_possible_misshift == TRUE
  ]
  dt[
    is.na(technical_frame_bias_possible),
    technical_frame_bias_possible := FALSE
  ]
  dt[
    ,
    target_outpaces_housekeeping_by_z_ge_1 :=
      is.finite(target_minus_housekeeping_top1pct_bumpiness_z_vs_depth) &
      target_minus_housekeeping_top1pct_bumpiness_z_vs_depth >= 1
  ]
  dt[
    ,
    candidate_specific_bumpy_no_hk_artifact :=
      target_top1pct_bumpiness_z_vs_depth >= 1 &
      target_outpaces_housekeeping_by_z_ge_1 == TRUE &
      technical_frame_bias_possible == FALSE
  ]
  dt[
    ,
    candidate_specific_frame_shift_no_hk_artifact :=
      target_canonical_frame_fraction < 0.45 &
      housekeeping_canonical_frame_fraction >= 0.50 &
      technical_frame_bias_possible == FALSE
  ]
  dt[
    ,
    infected_library_qc_class := fcase(
      technical_frame_bias_possible == TRUE,
      "technical_frame_bias_possible",
      candidate_specific_bumpy_no_hk_artifact == TRUE &
        target_outpaces_ddit3_by_z_ge_1 == TRUE,
      "candidate_bumpy_hk_and_ddit3_supported",
      candidate_specific_bumpy_no_hk_artifact == TRUE,
      "candidate_bumpy_vs_housekeeping",
      candidate_specific_frame_shift_no_hk_artifact == TRUE,
      "candidate_frame_shift_vs_housekeeping",
      default = "no_strong_target_specific_signal"
    )
  ]
  dt[]
}

myc_rrm1_all_library_metrics <- add_ddit3_and_library_qc(
  target_housekeeping_all_controls[target_gene %chin% c("MYC", "RRM1")]
)

myc_rrm1_viral_design_run_map <- make_viral_design_run_map(
  viral_clean[gene_symbol %chin% c("MYC", "RRM1")]
)
myc_rrm1_case_design_run_map <- collapse_viral_run_map(
  myc_rrm1_viral_design_run_map, "case"
)
myc_rrm1_control_design_run_map <- collapse_viral_run_map(
  myc_rrm1_viral_design_run_map, "control"
)

myc_rrm1_infected_library_metrics <- myc_rrm1_all_library_metrics[
  Run %chin% viral_case_union
]
if (nrow(myc_rrm1_infected_library_metrics) > 0) {
  if (nrow(myc_rrm1_case_design_run_map) > 0) {
    myc_rrm1_infected_library_metrics <- merge(
      myc_rrm1_infected_library_metrics,
      myc_rrm1_case_design_run_map,
      by = c("target_gene", "Run"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  metric_cols <- c(
    "target_gene", "Run", "infected_library_qc_class",
    "target_cds_total_counts", "target_frame0_counts",
    "target_frame1_counts", "target_frame2_counts",
    "target_frame0_fraction", "target_frame1_fraction",
    "target_frame2_fraction", "target_canonical_frame_fraction",
    "target_off_frame_fraction", "target_dominant_cds_frame",
    "target_dominant_cds_frame_fraction", "target_top1pct_share",
    "target_top1pct_off_frame_fraction",
    "target_top1pct_bumpiness_z_vs_depth",
    "housekeeping_n_genes", "housekeeping_total_counts",
    "housekeeping_canonical_frame_fraction",
    "housekeeping_median_canonical_frame_fraction",
    "housekeeping_median_top1pct_share",
    "housekeeping_median_top1pct_off_frame_fraction",
    "housekeeping_median_top1pct_bumpiness_z_vs_depth",
    "target_minus_housekeeping_median_canonical_frame_fraction",
    "target_minus_housekeeping_top1pct_share",
    "target_minus_housekeeping_top1pct_off_frame_fraction",
    "target_minus_housekeeping_top1pct_bumpiness_z_vs_depth",
    "metadata_canonical_frame_fraction", "metadata_dominant_frame",
    "metadata_possible_misshift",
    "global_canonical_frame_fraction", "global_possible_misshift",
    "housekeeping_possible_misshift",
    "target_possible_library_misshift_for_gene",
    "technical_frame_bias_possible",
    "candidate_specific_bumpy_no_hk_artifact",
    "candidate_specific_frame_shift_no_hk_artifact",
    "has_ddit3_fair_control", "ddit3_cds_total_counts",
    "ddit3_canonical_frame_fraction",
    "ddit3_top1pct_bumpiness_z_vs_depth",
    "target_minus_ddit3_top1pct_bumpiness_z_vs_depth",
    "target_minus_ddit3_canonical_frame_fraction",
    "target_outpaces_ddit3_by_z_ge_1",
    "viral_study", "viral_AUTHOR", "viral_CELL_LINE", "viral_TISSUE",
    "viral_case_condition", "viral_control_conditions",
    "CELL_LINE", "TISSUE", "CONDITION", "GENE", "AUTHOR"
  )
  metric_cols <- intersect(metric_cols,
                           names(myc_rrm1_infected_library_metrics))
  setcolorder(myc_rrm1_infected_library_metrics,
              c(metric_cols,
                setdiff(names(myc_rrm1_infected_library_metrics),
                        metric_cols)))
  setorder(myc_rrm1_infected_library_metrics, target_gene,
           -candidate_specific_bumpy_no_hk_artifact,
           -target_top1pct_bumpiness_z_vs_depth)
}

myc_rrm1_infected_library_summary <- if (
  nrow(myc_rrm1_infected_library_metrics) > 0
) {
  myc_rrm1_infected_library_metrics[
    ,
    .(
      infected_n_libraries = .N,
      median_target_counts =
        median(target_cds_total_counts, na.rm = TRUE),
      median_target_canonical_frame_fraction =
        median(target_canonical_frame_fraction, na.rm = TRUE),
      median_housekeeping_canonical_frame_fraction =
        median(housekeeping_canonical_frame_fraction, na.rm = TRUE),
      median_target_top1pct_bumpiness_z_vs_depth =
        median(target_top1pct_bumpiness_z_vs_depth, na.rm = TRUE),
      median_housekeeping_top1pct_bumpiness_z_vs_depth =
        median(housekeeping_median_top1pct_bumpiness_z_vs_depth,
               na.rm = TRUE),
      median_target_minus_housekeeping_top1pct_bumpiness_z_vs_depth =
        median(target_minus_housekeeping_top1pct_bumpiness_z_vs_depth,
               na.rm = TRUE),
      n_technical_frame_bias_possible =
        sum(technical_frame_bias_possible == TRUE, na.rm = TRUE),
      n_candidate_specific_bumpy_no_hk_artifact =
        sum(candidate_specific_bumpy_no_hk_artifact == TRUE, na.rm = TRUE),
      n_candidate_specific_frame_shift_no_hk_artifact =
        sum(candidate_specific_frame_shift_no_hk_artifact == TRUE,
            na.rm = TRUE),
      n_has_ddit3_fair_control =
        sum(has_ddit3_fair_control == TRUE, na.rm = TRUE),
      n_outpaces_ddit3_by_z_ge_1 =
        sum(target_outpaces_ddit3_by_z_ge_1 == TRUE, na.rm = TRUE)
    ),
    by = target_gene
  ]
} else {
  data.table()
}

myc_rrm1_infected_library_context_summary <- if (
  nrow(myc_rrm1_infected_library_metrics) > 0
) {
  context_dt <- copy(myc_rrm1_infected_library_metrics)
  for (col in c("viral_study", "viral_AUTHOR", "viral_CELL_LINE",
                "viral_TISSUE", "viral_case_condition",
                "viral_control_conditions")) {
    if (!col %in% names(context_dt)) context_dt[, (col) := NA_character_]
  }
  for (col in c("study", "BioProject", "AUTHOR", "CELL_LINE", "TISSUE",
                "CONDITION")) {
    if (!col %in% names(context_dt)) context_dt[, (col) := NA_character_]
  }
  context_dt[
    ,
    context_study := fifelse(
      !is.na(viral_study) & nzchar(viral_study),
      viral_study,
      fifelse("study" %in% names(context_dt) &
                !is.na(study) & nzchar(study),
              study,
              fifelse("BioProject" %in% names(context_dt) &
                        !is.na(BioProject) & nzchar(BioProject),
                      BioProject, "unknown_study"))
    )
  ]
  context_dt[
    ,
    context_author := fifelse(!is.na(viral_AUTHOR) & nzchar(viral_AUTHOR),
                              viral_AUTHOR,
                              fifelse(!is.na(AUTHOR) & nzchar(AUTHOR),
                                      AUTHOR, "unknown_author"))
  ]
  context_dt[
    ,
    context_cell_line := fifelse(
      !is.na(viral_CELL_LINE) & nzchar(viral_CELL_LINE),
      viral_CELL_LINE,
      fifelse(!is.na(CELL_LINE) & nzchar(CELL_LINE),
              CELL_LINE, "unknown_cell_line")
    )
  ]
  context_dt[
    ,
    context_tissue := fifelse(!is.na(viral_TISSUE) & nzchar(viral_TISSUE),
                              viral_TISSUE,
                              fifelse(!is.na(TISSUE) & nzchar(TISSUE),
                                      TISSUE, "unknown_tissue"))
  ]
  context_dt[
    ,
    context_case_condition := fifelse(
      !is.na(viral_case_condition) & nzchar(viral_case_condition),
      viral_case_condition,
      fifelse(!is.na(CONDITION) & nzchar(CONDITION),
              CONDITION, "unknown_condition")
    )
  ]
  out <- context_dt[
    ,
    .(
      n_infected_libraries = .N,
      n_candidate_specific_bumpy =
        sum(candidate_specific_bumpy_no_hk_artifact == TRUE, na.rm = TRUE),
      n_candidate_specific_frame_shift =
        sum(candidate_specific_frame_shift_no_hk_artifact == TRUE,
            na.rm = TRUE),
      n_technical_frame_bias =
        sum(technical_frame_bias_possible == TRUE, na.rm = TRUE),
      n_no_strong_signal =
        sum(infected_library_qc_class ==
              "no_strong_target_specific_signal", na.rm = TRUE),
      fraction_candidate_specific_bumpy =
        sum(candidate_specific_bumpy_no_hk_artifact == TRUE, na.rm = TRUE) /
        .N,
      fraction_candidate_specific_frame_shift =
        sum(candidate_specific_frame_shift_no_hk_artifact == TRUE,
            na.rm = TRUE) / .N,
      fraction_technical_frame_bias =
        sum(technical_frame_bias_possible == TRUE, na.rm = TRUE) / .N,
      median_target_counts = safe_median(target_cds_total_counts),
      median_target_canonical_frame_fraction =
        safe_median(target_canonical_frame_fraction),
      median_housekeeping_canonical_frame_fraction =
        safe_median(housekeeping_canonical_frame_fraction),
      median_target_top1pct_bumpiness_z_vs_depth =
        safe_median(target_top1pct_bumpiness_z_vs_depth),
      median_housekeeping_top1pct_bumpiness_z_vs_depth =
        safe_median(housekeeping_median_top1pct_bumpiness_z_vs_depth),
      median_target_minus_housekeeping_bumpiness_z =
        safe_median(target_minus_housekeeping_top1pct_bumpiness_z_vs_depth),
      n_has_ddit3_fair_control =
        sum(has_ddit3_fair_control == TRUE, na.rm = TRUE),
      n_outpaces_ddit3_by_z_ge_1 =
        sum(target_outpaces_ddit3_by_z_ge_1 == TRUE, na.rm = TRUE),
      representative_candidate_runs =
        collapse_text(Run[candidate_specific_bumpy_no_hk_artifact == TRUE],
                      max_n = 5L),
      representative_technical_runs =
        collapse_text(Run[technical_frame_bias_possible == TRUE], max_n = 5L)
    ),
    by = .(target_gene, context_study, context_author, context_cell_line,
           context_tissue, context_case_condition,
           viral_control_conditions)
  ]
  out[
    ,
    context_signal_class := fcase(
      n_infected_libraries < 2,
      "single_library_context",
      fraction_technical_frame_bias >= 0.50,
      "context_technical_bias_possible",
      n_candidate_specific_bumpy >= 2 &
        fraction_candidate_specific_bumpy >= 0.50,
      "context_enriched_candidate_bumpy",
      n_candidate_specific_bumpy >= 2,
      "context_has_repeated_candidate_bumpy",
      n_candidate_specific_frame_shift >= 2,
      "context_has_repeated_frame_shift",
      default = "context_no_strong_candidate_signal"
    )
  ]
  setorder(out, target_gene, -n_candidate_specific_bumpy,
           -fraction_candidate_specific_bumpy,
           -median_target_minus_housekeeping_bumpiness_z)
  out[]
} else {
  data.table()
}

summarize_design_library_role <- function(dt) {
  if (nrow(dt) == 0) {
    return(list(
      n_evaluable_libraries = 0L,
      n_candidate_specific_bumpy = 0L,
      n_candidate_specific_frame_shift = 0L,
      n_technical_frame_bias = 0L,
      fraction_candidate_specific_bumpy = NA_real_,
      fraction_technical_frame_bias = NA_real_,
      median_target_counts = NA_real_,
      median_target_canonical_frame_fraction = NA_real_,
      median_housekeeping_canonical_frame_fraction = NA_real_,
      median_target_minus_housekeeping_bumpiness_z = NA_real_,
      median_target_bumpiness_z = NA_real_,
      median_housekeeping_bumpiness_z = NA_real_,
      representative_candidate_runs = NA_character_
    ))
  }
  list(
    n_evaluable_libraries = uniqueN(dt$Run),
    n_candidate_specific_bumpy =
      sum(dt$candidate_specific_bumpy_no_hk_artifact == TRUE, na.rm = TRUE),
    n_candidate_specific_frame_shift =
      sum(dt$candidate_specific_frame_shift_no_hk_artifact == TRUE,
          na.rm = TRUE),
    n_technical_frame_bias =
      sum(dt$technical_frame_bias_possible == TRUE, na.rm = TRUE),
    fraction_candidate_specific_bumpy =
      sum(dt$candidate_specific_bumpy_no_hk_artifact == TRUE, na.rm = TRUE) /
      uniqueN(dt$Run),
    fraction_technical_frame_bias =
      sum(dt$technical_frame_bias_possible == TRUE, na.rm = TRUE) /
      uniqueN(dt$Run),
    median_target_counts = safe_median(dt$target_cds_total_counts),
    median_target_canonical_frame_fraction =
      safe_median(dt$target_canonical_frame_fraction),
    median_housekeeping_canonical_frame_fraction =
      safe_median(dt$housekeeping_canonical_frame_fraction),
    median_target_minus_housekeeping_bumpiness_z =
      safe_median(dt$target_minus_housekeeping_top1pct_bumpiness_z_vs_depth),
    median_target_bumpiness_z =
      safe_median(dt$target_top1pct_bumpiness_z_vs_depth),
    median_housekeeping_bumpiness_z =
      safe_median(dt$housekeeping_median_top1pct_bumpiness_z_vs_depth),
    representative_candidate_runs =
      collapse_text(dt$Run[dt$candidate_specific_bumpy_no_hk_artifact == TRUE],
                    max_n = 5L)
  )
}

myc_rrm1_viral_design_housekeeping_contrast_summary <- if (
  nrow(myc_rrm1_all_library_metrics) > 0 &&
    nrow(viral_clean[gene_symbol %chin% c("MYC", "RRM1")]) > 0
) {
  design_rows <- viral_clean[gene_symbol %chin% c("MYC", "RRM1")]
  out <- rbindlist(lapply(seq_len(nrow(design_rows)), function(i) {
    row <- design_rows[i]
    case_runs <- parse_run_signature(row$case_run_signature)
    control_runs <- parse_run_signature(row$control_run_signature)
    case_dt <- myc_rrm1_all_library_metrics[
      target_gene == row$gene_symbol & Run %chin% case_runs
    ]
    control_dt <- myc_rrm1_all_library_metrics[
      target_gene == row$gene_symbol & Run %chin% control_runs
    ]
    case_summary <- as.data.table(
      summarize_design_library_role(case_dt)
    )
    control_summary <- as.data.table(
      summarize_design_library_role(control_dt)
    )
    names(case_summary) <- paste0("case_", names(case_summary))
    names(control_summary) <- paste0("control_", names(control_summary))
    data.table(
      target_gene = row$gene_symbol,
      study = row$study,
      AUTHOR = row$AUTHOR,
      CELL_LINE = row$CELL_LINE,
      TISSUE = row$TISSUE,
      case_condition = row$case_condition,
      control_conditions = row$control_conditions,
      n_case_runs_in_design = length(case_runs),
      n_control_runs_in_design = length(control_runs)
    )[
      ,
      cbind(.SD, case_summary, control_summary)
    ][
      ,
      `:=`(
        case_minus_control_fraction_candidate_specific_bumpy =
          case_fraction_candidate_specific_bumpy -
          control_fraction_candidate_specific_bumpy,
        case_minus_control_median_target_minus_housekeeping_bumpiness_z =
          case_median_target_minus_housekeeping_bumpiness_z -
          control_median_target_minus_housekeeping_bumpiness_z,
        case_minus_control_median_target_canonical_frame_fraction =
          case_median_target_canonical_frame_fraction -
          control_median_target_canonical_frame_fraction,
        case_vs_control_target_minus_housekeeping_bumpiness_wilcox_p =
          safe_wilcox_p(
            case_dt$target_minus_housekeeping_top1pct_bumpiness_z_vs_depth,
            control_dt$target_minus_housekeeping_top1pct_bumpiness_z_vs_depth
          )
      )
    ]
  }), fill = TRUE)
  out[
    ,
    case_vs_control_target_minus_housekeeping_bumpiness_wilcox_q :=
      p.adjust(case_vs_control_target_minus_housekeeping_bumpiness_wilcox_p,
               "BH")
  ]
  out[
    ,
    design_housekeeping_qc_class := fcase(
      case_n_evaluable_libraries < 2 | control_n_evaluable_libraries < 2,
      "insufficient_matched_evaluable_libraries",
      case_fraction_technical_frame_bias >= 0.50 |
        control_fraction_technical_frame_bias >= 0.50,
      "technical_frame_bias_context",
      case_n_candidate_specific_bumpy >= 2 &
        case_minus_control_fraction_candidate_specific_bumpy > 0 &
        case_minus_control_median_target_minus_housekeeping_bumpiness_z > 0,
      "case_enriched_candidate_bumpy_vs_control",
      control_n_candidate_specific_bumpy >= 2 &
        case_minus_control_fraction_candidate_specific_bumpy < 0,
      "control_enriched_candidate_bumpy",
      default = "no_clear_case_specific_bumpiness"
    )
  ]
  setorder(out, target_gene,
           -case_minus_control_fraction_candidate_specific_bumpy,
           -case_minus_control_median_target_minus_housekeeping_bumpiness_z)
  out[]
} else {
  data.table()
}

myc_rrm1_hotspot_group_profiles <- data.table()
myc_rrm1_infected_library_top_hotspots <- data.table()
myc_rrm1_recurrent_hotspots <- data.table()
myc_rrm1_hotspot_recurrence_summary <- data.table()

if (nrow(myc_rrm1_infected_library_metrics) > 0) {
  hotspot_profile_list <- list()
  per_run_hotspot_list <- list()
  recurrent_hotspot_list <- list()

  for (gene in hotspot_genes) {
    mat <- hotspot_mats[[gene]]
    cds_positions <- hotspot_cds_positions[[gene]]
    if (is.null(mat) || is.null(cds_positions)) next
    gene_metrics <- myc_rrm1_infected_library_metrics[target_gene == gene]
    if (nrow(gene_metrics) == 0) next
    group_runs <- list(
      target_specific_bumpy =
        gene_metrics[candidate_specific_bumpy_no_hk_artifact == TRUE, Run],
      target_specific_frame_shift =
        gene_metrics[candidate_specific_frame_shift_no_hk_artifact == TRUE,
                     Run],
      technical_frame_bias =
        gene_metrics[technical_frame_bias_possible == TRUE, Run],
      no_strong_signal =
        gene_metrics[infected_library_qc_class ==
                       "no_strong_target_specific_signal", Run],
      all_infected = gene_metrics$Run,
      viral_control = viral_control_union
    )
    group_runs <- lapply(group_runs, unique)

    hotspot_profile_list[[gene]] <- rbindlist(lapply(
      names(group_runs),
      function(group) {
        cds_hotspot_profile(mat, group_runs[[group]], cds_positions,
                            gene, group)
      }
    ), fill = TRUE)

    group_lookup <- setNames(gene_metrics$infected_library_qc_class,
                             gene_metrics$Run)
    per_run <- per_run_cds_hotspots(
      mat = mat,
      runs = gene_metrics$Run,
      cds_positions = cds_positions,
      gene = gene,
      group_lookup = group_lookup
    )
    if (nrow(per_run) > 0) {
      context_cols <- intersect(
        c("Run", "target_gene", "target_cds_total_counts",
          "target_canonical_frame_fraction",
          "target_top1pct_bumpiness_z_vs_depth",
          "target_minus_housekeeping_top1pct_bumpiness_z_vs_depth",
          "technical_frame_bias_possible",
          "metadata_canonical_frame_fraction", "metadata_dominant_frame",
          "metadata_possible_misshift",
          "candidate_specific_bumpy_no_hk_artifact",
          "candidate_specific_frame_shift_no_hk_artifact",
          "has_ddit3_fair_control",
          "target_minus_ddit3_top1pct_bumpiness_z_vs_depth",
          "CELL_LINE", "TISSUE", "CONDITION", "GENE", "AUTHOR"),
        names(gene_metrics)
      )
      per_run <- merge(
        per_run,
        unique(gene_metrics[, ..context_cols], by = "Run"),
        by = "Run",
        all.x = TRUE,
        sort = FALSE
      )
    }
    per_run_hotspot_list[[gene]] <- per_run

    recurrent <- rbindlist(lapply(
      c("target_specific_bumpy", "target_specific_frame_shift",
        "technical_frame_bias", "no_strong_signal"),
      function(group) {
        summarize_recurrent_hotspots(per_run, group, group_runs[[group]])
      }
    ), fill = TRUE)
    if (nrow(recurrent) > 0) {
      recurrent_wide <- dcast(
        recurrent,
        gene_symbol + cds_relative_position + codon_index +
          transcript_position + cds_frame ~ hotspot_group,
        value.var = "fraction_runs_with_hotspot",
        fill = 0
      )
      names(recurrent_wide) <- sub("^target_specific_bumpy$",
                                   "target_specific_bumpy_fraction",
                                   names(recurrent_wide))
      names(recurrent_wide) <- sub("^target_specific_frame_shift$",
                                   "target_specific_frame_shift_fraction",
                                   names(recurrent_wide))
      names(recurrent_wide) <- sub("^technical_frame_bias$",
                                   "technical_frame_bias_fraction",
                                   names(recurrent_wide))
      names(recurrent_wide) <- sub("^no_strong_signal$",
                                   "no_strong_signal_fraction",
                                   names(recurrent_wide))
      recurrent <- merge(
        recurrent,
        recurrent_wide,
        by = c("gene_symbol", "cds_relative_position", "codon_index",
               "transcript_position", "cds_frame"),
        all.x = TRUE,
        sort = FALSE
      )
      if (!"target_specific_bumpy_fraction" %in% names(recurrent)) {
        recurrent[, target_specific_bumpy_fraction := 0]
      }
      if (!"no_strong_signal_fraction" %in% names(recurrent)) {
        recurrent[, no_strong_signal_fraction := 0]
      }
      if (!"technical_frame_bias_fraction" %in% names(recurrent)) {
        recurrent[, technical_frame_bias_fraction := 0]
      }
      recurrent[
        ,
        target_specific_minus_no_signal_fraction :=
          target_specific_bumpy_fraction - no_strong_signal_fraction
      ]
      recurrent[
        ,
        target_specific_minus_technical_fraction :=
          target_specific_bumpy_fraction - technical_frame_bias_fraction
      ]
      setorder(recurrent, gene_symbol, hotspot_group,
               -fraction_runs_with_hotspot, -total_position_counts)
    }
    recurrent_hotspot_list[[gene]] <- recurrent
  }

  myc_rrm1_hotspot_group_profiles <- rbindlist(hotspot_profile_list,
                                               fill = TRUE)
  myc_rrm1_infected_library_top_hotspots <- rbindlist(per_run_hotspot_list,
                                                      fill = TRUE)
  myc_rrm1_recurrent_hotspots <- rbindlist(recurrent_hotspot_list,
                                           fill = TRUE)
  if (nrow(myc_rrm1_recurrent_hotspots) > 0) {
    target_recurrent <- myc_rrm1_recurrent_hotspots[
      hotspot_group == "target_specific_bumpy"
    ]
    myc_rrm1_hotspot_recurrence_summary <- target_recurrent[
      ,
      .(
        target_specific_n_libraries = max(group_n_runs, na.rm = TRUE),
        max_fraction_runs_with_hotspot =
          max(fraction_runs_with_hotspot, na.rm = TRUE),
        n_positions_recurrent_ge_20pct =
          sum(fraction_runs_with_hotspot >= 0.20, na.rm = TRUE),
        n_positions_enriched_vs_no_signal_ge_10pct =
          sum(target_specific_minus_no_signal_fraction >= 0.10,
              na.rm = TRUE),
        n_positions_enriched_vs_technical_ge_10pct =
          sum(target_specific_minus_technical_fraction >= 0.10,
              na.rm = TRUE),
        top_cds_position =
          cds_relative_position[which.max(fraction_runs_with_hotspot)],
        top_codon_index = codon_index[which.max(fraction_runs_with_hotspot)],
        top_cds_frame = cds_frame[which.max(fraction_runs_with_hotspot)],
        top_n_runs_with_hotspot =
          n_runs_with_hotspot[which.max(fraction_runs_with_hotspot)],
        top_no_strong_signal_fraction =
          no_strong_signal_fraction[which.max(fraction_runs_with_hotspot)],
        top_technical_frame_bias_fraction =
          technical_frame_bias_fraction[which.max(fraction_runs_with_hotspot)]
      ),
      by = gene_symbol
    ]
    myc_rrm1_hotspot_recurrence_summary[
      target_specific_n_libraries < 3,
      recurrence_interpretation :=
        "insufficient_target_specific_libraries"
    ]
    myc_rrm1_hotspot_recurrence_summary[
      target_specific_n_libraries >= 3 &
        max_fraction_runs_with_hotspot >= 0.50,
      recurrence_interpretation := "strong_recurrent_hotspot"
    ]
    myc_rrm1_hotspot_recurrence_summary[
      target_specific_n_libraries >= 3 &
        max_fraction_runs_with_hotspot < 0.50 &
        max_fraction_runs_with_hotspot >= 0.20,
      recurrence_interpretation := "weak_or_moderate_recurrence"
    ]
    myc_rrm1_hotspot_recurrence_summary[
      target_specific_n_libraries >= 3 &
        max_fraction_runs_with_hotspot < 0.20,
      recurrence_interpretation := "no_recurrent_hotspot"
    ]
  }
}

target_profiles <- rbindlist(target_profiles, fill = TRUE)

fwrite(housekeeping_gene_summary,
       file.path(output_dir, "housekeeping_baseline_gene_summary.csv"))
fwrite(housekeeping_qc,
       file.path(output_dir, "housekeeping_baseline_run_qc.csv"))
fwrite(run_metrics, file.path(output_dir, "cds_frame_bumpiness_run_metrics.csv"))
fwrite(global_qc, file.path(output_dir, "run_global_cds_frame_qc.csv"))
fwrite(target_run_metrics,
       file.path(output_dir, "target_gene_frame_bumpiness_run_metrics.csv"))
fwrite(viral_design_summaries,
       file.path(output_dir, "target_gene_viral_design_frame_bumpiness_summary.csv"))
fwrite(target_summary,
       file.path(output_dir, "target_gene_frame_bumpiness_summary.csv"))
fwrite(shared_viral_libraries,
       file.path(output_dir, "shared_viral_library_candidate_vs_ddit3_qc.csv"))
fwrite(shared_ddit3_controls,
       file.path(output_dir,
                 "shared_viral_candidate_vs_ddit3_run_metrics.csv"))
fwrite(shared_ddit3_summary,
       file.path(output_dir, "shared_viral_candidate_vs_ddit3_summary.csv"))
fwrite(target_housekeeping_controls,
       file.path(output_dir,
                 "viral_target_vs_housekeeping_run_metrics.csv"))
fwrite(target_housekeeping_summary,
       file.path(output_dir, "viral_target_vs_housekeeping_summary.csv"))
fwrite(myc_rrm1_infected_library_metrics,
       file.path(output_dir,
                 "myc_rrm1_infected_library_frame_bumpiness_metrics.csv"))
fwrite(myc_rrm1_infected_library_summary,
       file.path(output_dir,
                 "myc_rrm1_infected_library_frame_bumpiness_summary.csv"))
fwrite(myc_rrm1_infected_library_context_summary,
       file.path(output_dir,
                 "myc_rrm1_infected_library_context_summary.csv"))
fwrite(myc_rrm1_viral_design_housekeeping_contrast_summary,
       file.path(output_dir,
                 "myc_rrm1_viral_design_housekeeping_contrast_summary.csv"))
fwrite(myc_rrm1_hotspot_group_profiles,
       file.path(output_dir, "myc_rrm1_cds_hotspot_group_profiles.csv"))
fwrite(myc_rrm1_infected_library_top_hotspots,
       file.path(output_dir, "myc_rrm1_infected_library_top_hotspots.csv"))
fwrite(myc_rrm1_recurrent_hotspots,
       file.path(output_dir, "myc_rrm1_recurrent_cds_hotspots.csv"))
fwrite(myc_rrm1_hotspot_recurrence_summary,
       file.path(output_dir, "myc_rrm1_hotspot_recurrence_summary.csv"))
fwrite(target_profiles,
       file.path(output_dir, "target_gene_cds_aggregate_profiles.csv"))

if (nrow(target_run_metrics) > 0) {
  plot_dt <- target_run_metrics[
    cds_total_counts >= 30 &
      is.finite(cds_top1pct_share) &
      is.finite(global_canonical_frame_fraction)
  ]
  plot_dt[, viral_union_role := "other"]
  plot_dt[Run %chin% viral_control_union, viral_union_role := "viral_control"]
  plot_dt[Run %chin% viral_case_union, viral_union_role := "viral_case"]
  plot_dt[, viral_union_role := factor(
    viral_union_role,
    levels = c("viral_case", "viral_control", "other")
  )]

  p_depth <- ggplot(
    plot_dt,
    aes(x = log10_cds_total_counts,
        y = cds_top1pct_share,
        color = viral_union_role)
  ) +
    geom_point(alpha = 0.45, size = 1.1) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.6,
                color = "#20242a") +
    facet_wrap(~ gene_symbol, scales = "free_x") +
    scale_color_manual(values = c(
      viral_case = "#9b1c31",
      viral_control = "#1f5aa6",
      other = "#9097a3"
    ), drop = FALSE) +
    labs(
      title = "Target-gene CDS bumpiness versus coverage depth",
      subtitle = "Top-1% position share should fall as total CDS coverage rises if peaks are mainly sampling noise",
      x = "log10 middle-CDS counts + 1",
      y = "Top-1% CDS positions share",
      color = "Run group"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir, "target_gene_bumpiness_vs_depth.png"),
         p_depth, width = 10.5, height = 5.8, dpi = 180)
  ggsave(file.path(figure_dir, "target_gene_bumpiness_vs_depth.pdf"),
         p_depth, width = 10.5, height = 5.8)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_depth,
      file.path(figure_dir, "target_gene_bumpiness_vs_depth.html"),
      title = "Target-gene CDS bumpiness versus coverage depth"
    )
  }

  p_frame <- ggplot(
    plot_dt,
    aes(x = global_canonical_frame_fraction,
        y = canonical_frame_fraction,
        color = viral_union_role)
  ) +
    geom_abline(slope = 1, intercept = 0, color = "#7b8494",
                linewidth = 0.35) +
    geom_hline(yintercept = 0.45, linetype = 2, color = "#7b8494",
               linewidth = 0.35) +
    geom_vline(xintercept = 0.45, linetype = 2, color = "#7b8494",
               linewidth = 0.35) +
    geom_point(alpha = 0.55, size = 1.1) +
    facet_wrap(~ gene_symbol) +
    scale_color_manual(values = c(
      viral_case = "#9b1c31",
      viral_control = "#1f5aa6",
      other = "#9097a3"
    ), drop = FALSE) +
    labs(
      title = "Gene-specific CDS frame versus global CDS-frame QC",
      subtitle = "If bumpy off-frame peaks are from library-level p-shift errors, points should also show low global canonical frame",
      x = "Run-level global canonical CDS-frame fraction",
      y = "Gene middle-CDS canonical frame fraction",
      color = "Run group"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir, "target_gene_frame_vs_global_frame.png"),
         p_frame, width = 10.5, height = 5.8, dpi = 180)
  ggsave(file.path(figure_dir, "target_gene_frame_vs_global_frame.pdf"),
         p_frame, width = 10.5, height = 5.8)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_frame,
      file.path(figure_dir, "target_gene_frame_vs_global_frame.html"),
      title = "Target-gene frame versus global frame"
    )
  }
}

if (nrow(housekeeping_run_metrics) > 0) {
  hk_gene_plot_dt <- housekeeping_run_metrics[
    is.finite(top1pct_bumpiness_z_vs_depth)
  ]
  if (nrow(hk_gene_plot_dt) > 0) {
    hk_order <- housekeeping_gene_summary[
      selected_for_housekeeping_baseline == TRUE,
      gene_symbol[order(median_top1pct_bumpiness_z_vs_depth)]
    ]
    hk_gene_plot_dt[
      ,
      gene_symbol := factor(gene_symbol, levels = hk_order)
    ]
    p_hk_genes <- ggplot(
      hk_gene_plot_dt,
      aes(x = gene_symbol, y = top1pct_bumpiness_z_vs_depth)
    ) +
      geom_hline(yintercept = 0, color = "#7b8494", linewidth = 0.3) +
      geom_hline(yintercept = 1, linetype = 2, color = "#9b1c31",
                 linewidth = 0.3) +
      geom_boxplot(fill = "#d9e6ee", color = "#36424f",
                   outlier.alpha = 0.12, outlier.size = 0.5,
                   linewidth = 0.35) +
      coord_flip() +
      labs(
        title = "Housekeeping-gene bumpiness baseline",
        subtitle = "Matched-depth top-1% bumpiness z across runs",
        x = NULL,
        y = "Matched-depth top-1% bumpiness z"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())
    ggsave(file.path(figure_dir,
                     "housekeeping_baseline_gene_bumpiness_distribution.png"),
           p_hk_genes, width = 7.5, height = 5.8, dpi = 180)
    ggsave(file.path(figure_dir,
                     "housekeeping_baseline_gene_bumpiness_distribution.pdf"),
           p_hk_genes, width = 7.5, height = 5.8)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_hk_genes,
        file.path(figure_dir,
                  "housekeeping_baseline_gene_bumpiness_distribution.html"),
        title = "Housekeeping bumpiness baseline"
      )
    }
  }
}

if (nrow(target_housekeeping_controls) > 0) {
  hk_compare_plot_dt <- target_housekeeping_controls[
    is.finite(target_top1pct_bumpiness_z_vs_depth) &
      is.finite(housekeeping_median_top1pct_bumpiness_z_vs_depth)
  ]
  if (nrow(hk_compare_plot_dt) > 0) {
    p_hk_bump <- ggplot(
      hk_compare_plot_dt,
      aes(
        x = housekeeping_median_top1pct_bumpiness_z_vs_depth,
        y = target_top1pct_bumpiness_z_vs_depth,
        color = housekeeping_possible_misshift,
        shape = global_possible_misshift
      )
    ) +
      geom_abline(slope = 1, intercept = 0, color = "#7b8494",
                  linewidth = 0.35) +
      geom_vline(xintercept = 1, linetype = 2, color = "#7b8494",
                 linewidth = 0.35) +
      geom_hline(yintercept = 1, linetype = 2, color = "#7b8494",
                 linewidth = 0.35) +
      geom_point(alpha = 0.8, size = 2) +
      facet_wrap(~ target_gene, scales = "free") +
      scale_color_manual(
        values = c(`TRUE` = "#9b1c31", `FALSE` = "#2f6f8f"),
        na.value = "#9097a3",
        drop = FALSE
      ) +
      scale_shape_manual(
        values = c(`TRUE` = 17, `FALSE` = 16),
        na.value = 16,
        drop = FALSE
      ) +
      labs(
        title = "Target CDS bumpiness versus housekeeping baseline in viral-case runs",
        subtitle = "Points above the diagonal are bumpier than always-on genes from the same run",
        x = "Housekeeping median matched-depth top-1% bumpiness z",
        y = "Target matched-depth top-1% bumpiness z",
        color = "Housekeeping\nframe misshift",
        shape = "Global frame\nmisshift"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir,
                     "viral_target_vs_housekeeping_bumpiness.png"),
           p_hk_bump, width = 10, height = 5.2, dpi = 180)
    ggsave(file.path(figure_dir,
                     "viral_target_vs_housekeeping_bumpiness.pdf"),
           p_hk_bump, width = 10, height = 5.2)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_hk_bump,
        file.path(figure_dir,
                  "viral_target_vs_housekeeping_bumpiness.html"),
        title = "Target versus housekeeping bumpiness"
      )
    }
  }

  hk_frame_plot_dt <- target_housekeeping_controls[
    is.finite(target_canonical_frame_fraction) &
      is.finite(housekeeping_canonical_frame_fraction)
  ]
  if (nrow(hk_frame_plot_dt) > 0) {
    p_hk_frame <- ggplot(
      hk_frame_plot_dt,
      aes(
        x = housekeeping_canonical_frame_fraction,
        y = target_canonical_frame_fraction,
        color = housekeeping_possible_misshift,
        shape = global_possible_misshift
      )
    ) +
      geom_abline(slope = 1, intercept = 0, color = "#7b8494",
                  linewidth = 0.35) +
      geom_vline(xintercept = 0.45, linetype = 2, color = "#7b8494",
                 linewidth = 0.35) +
      geom_hline(yintercept = 0.45, linetype = 2, color = "#7b8494",
                 linewidth = 0.35) +
      geom_point(alpha = 0.8, size = 2) +
      facet_wrap(~ target_gene) +
      scale_color_manual(
        values = c(`TRUE` = "#9b1c31", `FALSE` = "#2f6f8f"),
        na.value = "#9097a3",
        drop = FALSE
      ) +
      scale_shape_manual(
        values = c(`TRUE` = 17, `FALSE` = 16),
        na.value = 16,
        drop = FALSE
      ) +
      labs(
        title = "Target CDS frame versus housekeeping baseline in viral-case runs",
        subtitle = "A technical p-shift should move housekeeping and target frame fractions together",
        x = "Housekeeping aggregate canonical frame fraction",
        y = "Target canonical frame fraction",
        color = "Housekeeping\nframe misshift",
        shape = "Global frame\nmisshift"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir, "viral_target_vs_housekeeping_frame.png"),
           p_hk_frame, width = 10, height = 5.2, dpi = 180)
    ggsave(file.path(figure_dir, "viral_target_vs_housekeeping_frame.pdf"),
           p_hk_frame, width = 10, height = 5.2)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_hk_frame,
        file.path(figure_dir, "viral_target_vs_housekeeping_frame.html"),
        title = "Target versus housekeeping frame"
      )
    }
  }
}

if (nrow(myc_rrm1_infected_library_metrics) > 0) {
  myc_rrm1_plot_dt <- myc_rrm1_infected_library_metrics[
    is.finite(target_top1pct_bumpiness_z_vs_depth) &
      is.finite(target_minus_housekeeping_top1pct_bumpiness_z_vs_depth)
  ]
  if (nrow(myc_rrm1_plot_dt) > 0) {
    p_myc_rrm1 <- ggplot(
      myc_rrm1_plot_dt,
      aes(
        x = target_minus_housekeeping_top1pct_bumpiness_z_vs_depth,
        y = target_canonical_frame_fraction,
        color = infected_library_qc_class,
        size = log10(target_cds_total_counts + 1)
      )
    ) +
      geom_hline(yintercept = 0.45, linetype = 2, color = "#7b8494",
                 linewidth = 0.35) +
      geom_vline(xintercept = 1, linetype = 2, color = "#7b8494",
                 linewidth = 0.35) +
      geom_point(alpha = 0.8) +
      facet_wrap(~ target_gene) +
      scale_color_manual(
        values = c(
          candidate_bumpy_hk_and_ddit3_supported = "#7a0d19",
          candidate_bumpy_vs_housekeeping = "#c6462f",
          candidate_frame_shift_vs_housekeeping = "#c57b00",
          technical_frame_bias_possible = "#6f7785",
          no_strong_target_specific_signal = "#2f6f8f"
        ),
        labels = c(
          candidate_bumpy_hk_and_ddit3_supported = "Bumpy vs HK + DDIT3",
          candidate_bumpy_vs_housekeeping = "Bumpy vs HK",
          candidate_frame_shift_vs_housekeeping = "Frame shift vs HK",
          technical_frame_bias_possible = "Technical frame bias",
          no_strong_target_specific_signal = "No strong signal"
        ),
        drop = FALSE
      ) +
      labs(
        title = "MYC/RRM1 infected-library frame and bumpiness metrics",
        subtitle = "X-axis is target bumpiness above same-run housekeeping baseline; y-axis is target canonical-frame fraction",
        x = "Target minus housekeeping top-1% bumpiness z",
        y = "Target canonical frame fraction",
        color = "Library QC class",
        size = "log10 target\nCDS counts + 1"
      ) +
      theme_minimal(base_size = 11) +
      guides(color = guide_legend(nrow = 3, byrow = TRUE),
             size = guide_legend(nrow = 2)) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir,
                     "myc_rrm1_infected_library_frame_bumpiness_metrics.png"),
           p_myc_rrm1, width = 11.5, height = 6.1, dpi = 180)
    ggsave(file.path(figure_dir,
                     "myc_rrm1_infected_library_frame_bumpiness_metrics.pdf"),
           p_myc_rrm1, width = 11.5, height = 6.1)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_myc_rrm1,
        file.path(figure_dir,
                  "myc_rrm1_infected_library_frame_bumpiness_metrics.html"),
        title = "MYC/RRM1 infected-library frame and bumpiness metrics"
      )
    }
  }
}

if (nrow(myc_rrm1_infected_library_context_summary) > 0) {
  context_plot_dt <- myc_rrm1_infected_library_context_summary[
    n_infected_libraries >= 2
  ]
  if (nrow(context_plot_dt) > 0) {
    context_plot_dt[
      ,
      context_label := paste0(
        context_cell_line, " / ", context_case_condition,
        "\n", sub("-homo_sapiens$", "", context_study)
      )
    ]
    context_plot_dt[
      ,
      context_label := reorder(
        context_label,
        fraction_candidate_specific_bumpy +
          0.01 * n_candidate_specific_bumpy
      )
    ]
    p_context <- ggplot(
      context_plot_dt,
      aes(x = context_label,
          y = fraction_candidate_specific_bumpy,
          fill = context_signal_class)
    ) +
      geom_col(width = 0.78, color = "#28313b", linewidth = 0.18) +
      geom_text(
        aes(label = paste0(n_candidate_specific_bumpy, "/",
                           n_infected_libraries)),
        hjust = -0.08,
        size = 2.8,
        color = "#20242a"
      ) +
      coord_flip(ylim = c(0, 1.05), clip = "off") +
      facet_wrap(~ target_gene, scales = "free_y") +
      scale_fill_manual(
        values = c(
          context_enriched_candidate_bumpy = "#9b1c31",
          context_has_repeated_candidate_bumpy = "#c6462f",
          context_has_repeated_frame_shift = "#c57b00",
          context_technical_bias_possible = "#6f7785",
          context_no_strong_candidate_signal = "#2f6f8f",
          single_library_context = "#b7bec7"
        ),
        labels = c(
          context_enriched_candidate_bumpy = "Enriched bumpy",
          context_has_repeated_candidate_bumpy = "Repeated bumpy",
          context_has_repeated_frame_shift = "Repeated frame shift",
          context_technical_bias_possible = "Technical bias",
          context_no_strong_candidate_signal = "No strong signal",
          single_library_context = "Single library"
        ),
        drop = FALSE
      ) +
      labs(
        title = "MYC/RRM1 bumpy infected libraries by viral context",
        subtitle = "Labels show candidate-specific bumpy libraries / evaluable infected libraries after same-run housekeeping QC",
        x = NULL,
        y = "Fraction candidate-specific bumpy",
        fill = "Context class"
      ) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      theme_minimal(base_size = 10.5) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(5.5, 28, 5.5, 5.5)
      )
    ggsave(file.path(figure_dir,
                     "myc_rrm1_infected_library_context_summary.png"),
           p_context, width = 11.5, height = 7.2, dpi = 180)
    ggsave(file.path(figure_dir,
                     "myc_rrm1_infected_library_context_summary.pdf"),
           p_context, width = 11.5, height = 7.2)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_context,
        file.path(figure_dir,
                  "myc_rrm1_infected_library_context_summary.html"),
        title = "MYC/RRM1 infected-library context summary"
      )
    }
  }
}

if (nrow(myc_rrm1_viral_design_housekeeping_contrast_summary) > 0) {
  design_plot_dt <- myc_rrm1_viral_design_housekeeping_contrast_summary[
    case_n_evaluable_libraries >= 2 &
      control_n_evaluable_libraries >= 2 &
      is.finite(case_minus_control_median_target_minus_housekeeping_bumpiness_z)
  ]
  if (nrow(design_plot_dt) > 0) {
    design_plot_dt[
      ,
      design_label := paste0(
        CELL_LINE, " / ", case_condition, " vs ", control_conditions,
        "\n", sub("-homo_sapiens$", "", study)
      )
    ]
    design_plot_dt[
      ,
      design_label := reorder(
        design_label,
        case_minus_control_median_target_minus_housekeeping_bumpiness_z
      )
    ]
    p_design <- ggplot(
      design_plot_dt,
      aes(
        x = design_label,
        y = case_minus_control_median_target_minus_housekeeping_bumpiness_z,
        fill = design_housekeeping_qc_class
      )
    ) +
      geom_hline(yintercept = 0, color = "#7b8494", linewidth = 0.35) +
      geom_col(width = 0.78, color = "#28313b", linewidth = 0.18) +
      coord_flip() +
      facet_wrap(~ target_gene, scales = "free_y") +
      scale_fill_manual(
        values = c(
          case_enriched_candidate_bumpy_vs_control = "#9b1c31",
          control_enriched_candidate_bumpy = "#1f5aa6",
          technical_frame_bias_context = "#6f7785",
          insufficient_matched_evaluable_libraries = "#b7bec7",
          no_clear_case_specific_bumpiness = "#2f6f8f"
        ),
        labels = c(
          case_enriched_candidate_bumpy_vs_control = "Case-enriched",
          control_enriched_candidate_bumpy = "Control-enriched",
          technical_frame_bias_context = "Technical bias",
          insufficient_matched_evaluable_libraries = "Insufficient matched",
          no_clear_case_specific_bumpiness = "No clear case signal"
        ),
        drop = FALSE
      ) +
      labs(
        title = "Matched viral case-control bumpiness after housekeeping QC",
        subtitle = "Positive values mean infected/case libraries are bumpier than matched controls after subtracting same-run housekeeping baseline",
        x = NULL,
        y = "Case - control median target-minus-housekeeping bumpiness z",
        fill = "Design class"
      ) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      theme_minimal(base_size = 10.5) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir,
                     "myc_rrm1_viral_design_housekeeping_contrast.png"),
           p_design, width = 11.5, height = 7.2, dpi = 180)
    ggsave(file.path(figure_dir,
                     "myc_rrm1_viral_design_housekeeping_contrast.pdf"),
           p_design, width = 11.5, height = 7.2)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_design,
        file.path(figure_dir,
                  "myc_rrm1_viral_design_housekeeping_contrast.html"),
        title = "MYC/RRM1 matched viral housekeeping contrast"
      )
    }
  }
}

if (nrow(myc_rrm1_recurrent_hotspots) > 0) {
  recurrent_plot_dt <- myc_rrm1_recurrent_hotspots[
    hotspot_group == "target_specific_bumpy" &
      group_n_runs >= 3 &
      fraction_runs_with_hotspot > 0
  ]
  if (nrow(recurrent_plot_dt) > 0) {
    recurrent_plot_dt <- recurrent_plot_dt[
      order(gene_symbol, -fraction_runs_with_hotspot,
            -total_position_counts)
    ][
      ,
      head(.SD, 15L),
      by = gene_symbol
    ]
    recurrent_plot_dt[
      ,
      hotspot_label := paste0("CDS ", cds_relative_position,
                              "\n", cds_frame)
    ]
    recurrent_plot_dt[
      ,
      hotspot_label := reorder(hotspot_label, fraction_runs_with_hotspot)
    ]
    p_recurrent <- ggplot(
      recurrent_plot_dt,
      aes(x = hotspot_label, y = fraction_runs_with_hotspot,
          color = cds_frame, size = total_position_counts)
    ) +
      geom_segment(
        aes(x = hotspot_label, xend = hotspot_label,
            y = 0, yend = fraction_runs_with_hotspot,
            color = cds_frame),
        linewidth = 0.35,
        alpha = 0.8,
        inherit.aes = FALSE
      ) +
      geom_point(alpha = 0.9) +
      coord_flip() +
      facet_wrap(~ gene_symbol, scales = "free_y") +
      scale_color_manual(values = c(
        frame0 = "#1f5aa6",
        frame1 = "#c57b00",
        frame2 = "#9b1c31"
      )) +
      labs(
        title = "Recurrent MYC/RRM1 CDS hotspots in target-specific infected libraries",
        subtitle = "Only genes with at least 3 target-specific bumpy libraries are plotted",
        x = NULL,
        y = "Fraction of target-specific bumpy libraries",
        color = "CDS frame",
        size = "Total hotspot\ncounts"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir, "myc_rrm1_recurrent_cds_hotspots.png"),
           p_recurrent, width = 10.5, height = 6.2, dpi = 180)
    ggsave(file.path(figure_dir, "myc_rrm1_recurrent_cds_hotspots.pdf"),
           p_recurrent, width = 10.5, height = 6.2)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_recurrent,
        file.path(figure_dir, "myc_rrm1_recurrent_cds_hotspots.html"),
        title = "MYC/RRM1 recurrent CDS hotspots"
      )
    }
  }
}

if (nrow(myc_rrm1_hotspot_group_profiles) > 0) {
  profile_plot_dt <- myc_rrm1_hotspot_group_profiles[
    hotspot_group %chin% c("target_specific_bumpy",
                           "target_specific_frame_shift",
                           "technical_frame_bias",
                           "no_strong_signal",
                           "viral_control")
  ]
  profile_plot_dt[
    ,
    hotspot_group := factor(
      hotspot_group,
      levels = c("target_specific_bumpy",
                 "target_specific_frame_shift",
                 "technical_frame_bias",
                 "no_strong_signal",
                 "viral_control"),
      labels = c("Target-specific bumpy",
                 "Target-specific frame shift",
                 "Technical frame bias",
                 "No strong signal",
                 "Viral controls")
    )
  ]
  if (nrow(profile_plot_dt) > 0) {
    p_hotspot_profile <- ggplot(
      profile_plot_dt,
      aes(x = cds_relative_position, y = normalized_density,
          color = cds_frame)
    ) +
      geom_line(linewidth = 0.32, alpha = 0.95) +
      facet_grid(gene_symbol ~ hotspot_group, scales = "free_x") +
      scale_color_manual(values = c(
        frame0 = "#1f5aa6",
        frame1 = "#c57b00",
        frame2 = "#9b1c31"
      )) +
      labs(
        title = "MYC/RRM1 CDS coverage profiles by infected-library QC class",
        subtitle = "Profiles are normalized within gene and group; recurrent group-specific peaks are stronger candidate biology than isolated spikes",
        x = "CDS-relative position",
        y = "Normalized density per 1000 CDS counts",
        color = "CDS frame"
      ) +
      theme_minimal(base_size = 10.5) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    ggsave(file.path(figure_dir, "myc_rrm1_cds_hotspot_group_profiles.png"),
           p_hotspot_profile, width = 13.5, height = 6.8, dpi = 180)
    ggsave(file.path(figure_dir, "myc_rrm1_cds_hotspot_group_profiles.pdf"),
           p_hotspot_profile, width = 13.5, height = 6.8)
    if (exists("dominant_save_ggplotly")) {
      dominant_save_ggplotly(
        p_hotspot_profile,
        file.path(figure_dir, "myc_rrm1_cds_hotspot_group_profiles.html"),
        title = "MYC/RRM1 CDS hotspot group profiles"
      )
    }
  }
}

if (nrow(shared_ddit3_controls) > 0) {
  p_shared <- ggplot(
    shared_ddit3_controls[
      is.finite(candidate_top1pct_bumpiness_z_vs_depth) &
        is.finite(ddit3_top1pct_bumpiness_z_vs_depth)
    ],
    aes(
      x = ddit3_top1pct_bumpiness_z_vs_depth,
      y = candidate_top1pct_bumpiness_z_vs_depth,
      color = global_possible_misshift
    )
  ) +
    geom_abline(slope = 1, intercept = 0, color = "#7b8494",
                linewidth = 0.35) +
    geom_vline(xintercept = 1, linetype = 2, color = "#7b8494",
               linewidth = 0.35) +
    geom_hline(yintercept = 1, linetype = 2, color = "#7b8494",
               linewidth = 0.35) +
    geom_point(alpha = 0.8, size = 2) +
    facet_wrap(~ candidate_gene, scales = "free") +
    scale_color_manual(
      values = c(`TRUE` = "#9b1c31", `FALSE` = "#2f6f8f"),
      na.value = "#9097a3",
      drop = FALSE
    ) +
    labs(
      title = "Candidate CDS bumpiness versus DDIT3 in the same viral-case runs",
      subtitle = "Only runs with at least 30 middle-CDS counts for both DDIT3 and the candidate are shown",
      x = "DDIT3 matched-depth top-1% bumpiness z",
      y = "Candidate matched-depth top-1% bumpiness z",
      color = "Global frame\nmisshift"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir, "shared_viral_candidate_vs_ddit3_bumpiness.png"),
         p_shared, width = 9.5, height = 4.8, dpi = 180)
  ggsave(file.path(figure_dir, "shared_viral_candidate_vs_ddit3_bumpiness.pdf"),
         p_shared, width = 9.5, height = 4.8)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_shared,
      file.path(figure_dir,
                "shared_viral_candidate_vs_ddit3_bumpiness.html"),
      title = "Candidate versus DDIT3 bumpiness"
    )
  }
}

if (nrow(target_profiles) > 0) {
  profile_dt <- copy(target_profiles)
  profile_dt[, profile_group := factor(
    profile_group,
    levels = c("viral_case_union", "viral_control_union", "all_3857_runs")
  )]
  p_profile <- ggplot(
    profile_dt[
      profile_group %chin% c("viral_case_union", "viral_control_union",
                             "all_3857_runs")
    ],
    aes(x = cds_relative_position, y = normalized_density,
        color = cds_frame)
  ) +
    geom_line(linewidth = 0.35, alpha = 0.95) +
    facet_grid(gene_symbol ~ profile_group, scales = "free_x") +
    scale_color_manual(values = c(
      frame0 = "#1f5aa6",
      frame1 = "#c57b00",
      frame2 = "#9b1c31"
    )) +
    labs(
      title = "CDS coverage profiles by CDS-relative frame",
      subtitle = "Normalized within each gene/group; all-runs profile tests whether apparent spikes are depth-related",
      x = "CDS-relative position",
      y = "Normalized density per 1000 CDS counts",
      color = "CDS frame"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figure_dir, "target_gene_cds_frame_profiles.png"),
         p_profile, width = 12.5, height = 7.8, dpi = 180)
  ggsave(file.path(figure_dir, "target_gene_cds_frame_profiles.pdf"),
         p_profile, width = 12.5, height = 7.8)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p_profile,
      file.path(figure_dir, "target_gene_cds_frame_profiles.html"),
      title = "Target-gene CDS frame profiles"
    )
  }
}

message("Saved run metrics: ",
        file.path(output_dir, "cds_frame_bumpiness_run_metrics.csv"))
message("Saved target summaries: ",
        file.path(output_dir, "target_gene_frame_bumpiness_summary.csv"))
message("Saved viral design summaries: ",
        file.path(output_dir,
                  "target_gene_viral_design_frame_bumpiness_summary.csv"))
