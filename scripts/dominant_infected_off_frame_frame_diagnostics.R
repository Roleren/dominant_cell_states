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
strict_review_dir <- file.path(analysis_dir, "dominant_infected_off_frame_review")
strict_hotspot_dir <- file.path(analysis_dir, "dominant_infected_off_frame_hotspots")
frame_risk_dir <- file.path(analysis_dir, "dominant_infected_off_frame_frame_risk")
output_dir <- file.path(analysis_dir,
                        "dominant_infected_off_frame_frame_diagnostics")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)
coverage_cache_file <- file.path(analysis_dir, "scripts", "dominant_coverage_cache.R")
if (!file.exists(coverage_cache_file)) {
  stop("Missing dominant coverage cache helper: ", coverage_cache_file,
       call. = FALSE)
}
source(coverage_cache_file)

message("Infected off-frame hotspot frame diagnostics")
message("  1. Combine strict and frame-risk hotspot queues")
message("  2. Measure start, stop, and hotspot-window frame behavior")
message("  3. Score simple-p-shift-like versus localized/mixed behavior")

required_files <- c(
  file.path(strict_review_dir, "infected_off_frame_hotspot_review_queue.csv"),
  file.path(strict_hotspot_dir, "infected_off_frame_hotspot_run_positions.csv"),
  file.path(analysis_dir, "dominant_infected_off_frame_cds_scan",
            "infected_off_frame_cds_top_candidate_runs.csv"),
  file.path(frame_risk_dir, "infected_off_frame_frame_risk_hotspots.csv"),
  file.path(frame_risk_dir, "infected_off_frame_frame_risk_run_positions.csv"),
  file.path(analysis_dir,
            "human_dominant_cell_states_clean_cds_gene_diagnostics.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "),
       call. = FALSE)
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_text <- function(x, max_n = 8L) {
  x <- unique(safe_chr(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  suffix <- if (length(x) > max_n) ";..." else ""
  paste0(paste(head(x, max_n), collapse = ";"), suffix)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_frac <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
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

window_positions <- function(center, before, after, max_position) {
  positions <- seq.int(center - before, center + after)
  positions[positions >= 1L & positions <= max_position]
}

frame_label <- function(pos, cds_start) {
  paste0("frame", as.integer((pos - cds_start) %% 3L))
}

dominant_frame_from_counts <- function(counts, positions, cds_start) {
  if (!length(counts) || !any(is.finite(counts)) || sum(counts, na.rm = TRUE) <= 0) {
    return(NA_character_)
  }
  frames <- frame_label(positions, cds_start)
  sums <- tapply(counts, frames, sum, na.rm = TRUE)
  if (!length(sums)) return(NA_character_)
  names(sums)[which.max(sums)]
}

count_at <- function(mat, pos, run) {
  if (is.na(pos) || pos < 1L || pos > nrow(mat) || !run %in% colnames(mat)) {
    return(NA_real_)
  }
  as.numeric(mat[pos, run, drop = TRUE])
}

measure_run_hotspot <- function(mat, run, cds_positions, hotspot_position,
                                cds_relative_position, codon_index,
                                hotspot_frame) {
  if (!run %in% colnames(mat)) return(data.table())
  cds_positions <- sort(unique(as.integer(cds_positions)))
  cds_start <- min(cds_positions)
  cds_end <- max(cds_positions)
  max_position <- nrow(mat)

  hotspot_position <- as.integer(hotspot_position)
  hotspot_count <- count_at(mat, hotspot_position, run)
  codon_start_rel <- (as.integer(codon_index) - 1L) * 3L + 1L
  codon_frame0_position <- if (codon_start_rel >= 1L &&
                               codon_start_rel <= length(cds_positions)) {
    cds_positions[[codon_start_rel]]
  } else {
    NA_integer_
  }
  codon_frame0_count <- count_at(mat, codon_frame0_position, run)

  local_positions <- window_positions(hotspot_position, 15L, 15L,
                                      max_position)
  local_counts <- as.numeric(mat[local_positions, run, drop = TRUE])
  local_total <- sum(local_counts, na.rm = TRUE)
  local_frame <- frame_label(local_positions, cds_start)
  local_frame_counts <- tapply(local_counts, local_frame, sum, na.rm = TRUE)
  for (frame in paste0("frame", 0:2)) {
    if (!frame %in% names(local_frame_counts)) local_frame_counts[[frame]] <- 0
  }
  local_frame_counts <- local_frame_counts[paste0("frame", 0:2)]
  local_median <- stats::median(local_counts[is.finite(local_counts)],
                                na.rm = TRUE)
  if (!is.finite(local_median)) local_median <- NA_real_
  local_rank <- rank(-local_counts, ties.method = "min")
  local_hit_index <- match(hotspot_position, local_positions)
  local_hotspot_rank <- if (!is.na(local_hit_index)) {
    as.integer(local_rank[[local_hit_index]])
  } else {
    NA_integer_
  }

  start_positions <- window_positions(cds_start, 15L, 45L, max_position)
  start_counts <- as.numeric(mat[start_positions, run, drop = TRUE])
  start_top_index <- if (length(start_counts) &&
                         any(is.finite(start_counts) & start_counts > 0)) {
    which.max(start_counts)
  } else {
    NA_integer_
  }
  start_top_position <- if (!is.na(start_top_index)) {
    start_positions[[start_top_index]]
  } else {
    NA_integer_
  }
  start_top_count <- if (!is.na(start_top_index)) {
    start_counts[[start_top_index]]
  } else {
    NA_real_
  }
  start_total <- sum(start_counts, na.rm = TRUE)
  start_dominant_frame <- dominant_frame_from_counts(start_counts,
                                                     start_positions,
                                                     cds_start)

  end_positions <- seq.int(max(cds_start, cds_end - 60L), cds_end)
  end_counts <- as.numeric(mat[end_positions, run, drop = TRUE])
  end_top_index <- if (length(end_counts) &&
                       any(is.finite(end_counts) & end_counts > 0)) {
    which.max(end_counts)
  } else {
    NA_integer_
  }
  end_top_position <- if (!is.na(end_top_index)) {
    end_positions[[end_top_index]]
  } else {
    NA_integer_
  }
  end_top_count <- if (!is.na(end_top_index)) {
    end_counts[[end_top_index]]
  } else {
    NA_real_
  }
  end_total <- sum(end_counts, na.rm = TRUE)
  end_dominant_frame <- dominant_frame_from_counts(end_counts,
                                                   end_positions,
                                                   cds_start)

  data.table(
    Run = run,
    cds_start_transcript_position = cds_start,
    cds_end_transcript_position = cds_end,
    hotspot_position_count = hotspot_count,
    hotspot_window_total_count = local_total,
    hotspot_window_frame0_count = as.numeric(local_frame_counts[["frame0"]]),
    hotspot_window_frame1_count = as.numeric(local_frame_counts[["frame1"]]),
    hotspot_window_frame2_count = as.numeric(local_frame_counts[["frame2"]]),
    hotspot_window_dominant_frame =
      dominant_frame_from_counts(local_counts, local_positions, cds_start),
    hotspot_local_rank_31nt = local_hotspot_rank,
    hotspot_peak_share_31nt = if (is.finite(local_total) &&
                                  local_total > 0) {
      hotspot_count / local_total
    } else {
      NA_real_
    },
    hotspot_peak_over_local_median = if (is.finite(local_median) &&
                                         local_median > 0) {
      hotspot_count / local_median
    } else {
      NA_real_
    },
    codon_frame0_transcript_position = codon_frame0_position,
    codon_frame0_count = codon_frame0_count,
    hotspot_to_codon_frame0_ratio = if (is.finite(codon_frame0_count) &&
                                        codon_frame0_count > 0) {
      hotspot_count / codon_frame0_count
    } else {
      NA_real_
    },
    start_window_total_count = start_total,
    start_window_dominant_frame = start_dominant_frame,
    start_top_transcript_position = start_top_position,
    start_top_offset_from_cds_start = start_top_position - cds_start,
    start_top_frame = if (!is.na(start_top_position)) {
      frame_label(start_top_position, cds_start)
    } else {
      NA_character_
    },
    start_top_count = start_top_count,
    start_top_share = if (is.finite(start_total) && start_total > 0) {
      start_top_count / start_total
    } else {
      NA_real_
    },
    end_window_total_count = end_total,
    end_window_dominant_frame = end_dominant_frame,
    end_top_transcript_position = end_top_position,
    end_top_offset_from_cds_end = end_top_position - cds_end,
    end_top_frame = if (!is.na(end_top_position)) {
      frame_label(end_top_position, cds_start)
    } else {
      NA_character_
    },
    end_top_count = end_top_count,
    end_top_share = if (is.finite(end_total) && end_total > 0) {
      end_top_count / end_total
    } else {
      NA_real_
    },
    start_end_top_same_frame = !is.na(start_top_position) &&
      !is.na(end_top_position) &&
      frame_label(start_top_position, cds_start) ==
      frame_label(end_top_position, cds_start),
    start_end_dominant_same_frame = !is.na(start_dominant_frame) &&
      !is.na(end_dominant_frame) &&
      start_dominant_frame == end_dominant_frame,
    hotspot_matches_start_top_frame = !is.na(hotspot_frame) &&
      !is.na(start_top_position) &&
      hotspot_frame == frame_label(start_top_position, cds_start),
    hotspot_matches_end_top_frame = !is.na(hotspot_frame) &&
      !is.na(end_top_position) &&
      hotspot_frame == frame_label(end_top_position, cds_start)
  )
}

make_window_profile <- function(mat, run, cds_positions, center_position,
                                center_label, before, after) {
  if (!run %in% colnames(mat)) return(data.table())
  cds_positions <- sort(unique(as.integer(cds_positions)))
  cds_start <- min(cds_positions)
  positions <- window_positions(center_position, before, after, nrow(mat))
  if (!length(positions)) return(data.table())
  data.table(
    Run = run,
    window = center_label,
    transcript_position = positions,
    offset = positions - center_position,
    cds_relative_position = match(positions, cds_positions),
    cds_frame = frame_label(positions, cds_start),
    count = as.numeric(mat[positions, run, drop = TRUE])
  )
}

strict_queue <- fread(
  file.path(strict_review_dir, "infected_off_frame_hotspot_review_queue.csv"),
  showProgress = FALSE
)
strict_run_positions <- fread(
  file.path(strict_hotspot_dir, "infected_off_frame_hotspot_run_positions.csv"),
  showProgress = FALSE
)
strict_run_metadata <- fread(
  file.path(analysis_dir, "dominant_infected_off_frame_cds_scan",
            "infected_off_frame_cds_top_candidate_runs.csv"),
  showProgress = FALSE
)
frame_risk_hotspots <- fread(
  file.path(frame_risk_dir, "infected_off_frame_frame_risk_hotspots.csv"),
  showProgress = FALSE
)
frame_risk_run_positions <- fread(
  file.path(frame_risk_dir,
            "infected_off_frame_frame_risk_run_positions.csv"),
  showProgress = FALSE
)

strict_queue[, gene_symbol := toupper(gene_symbol)]
strict_run_positions[, gene_symbol := toupper(gene_symbol)]
strict_run_metadata[, gene_symbol := toupper(gene_symbol)]
frame_risk_hotspots[, gene_symbol := toupper(gene_symbol)]
frame_risk_run_positions[, gene_symbol := toupper(gene_symbol)]

strict_def <- strict_queue[
  ,
  .(
    evidence_layer = "strict_no_frame_risk",
    gene_symbol,
    tx_id,
    codon_index,
    cds_relative_position,
    cds_percent,
    transcript_position,
    cds_frame,
    source_score = hotspot_review_priority_score,
    source_class = hotspot_review_class,
    source_tier = review_tier,
    source_run_hits = n_runs_with_hotspot_candidate_infected,
    source_group_runs = group_n_runs_candidate_infected,
    source_contexts = contexts_candidate_infected,
    source_run_examples = run_examples_candidate_infected
  )
]

frame_risk_def <- frame_risk_hotspots[
  ,
  .(
    evidence_layer = "frame_risk",
    gene_symbol,
    tx_id,
    codon_index,
    cds_relative_position,
    cds_percent,
    transcript_position,
    cds_frame,
    source_score = frame_risk_hotspot_score,
    source_class = "frame_risk_hotspot",
    source_tier = fifelse(n_frame_risk_runs_with_hotspot >= 2,
                          "recurrent_frame_risk",
                          "single_run_frame_risk"),
    source_run_hits = n_frame_risk_runs_with_hotspot,
    source_group_runs = group_n_frame_risk_runs,
    source_contexts = contexts,
    source_run_examples = run_examples
  )
]

hotspot_def <- rbindlist(list(strict_def, frame_risk_def), fill = TRUE)
hotspot_def <- unique(hotspot_def,
                      by = c("evidence_layer", "gene_symbol", "tx_id",
                             "transcript_position", "cds_frame"))
hotspot_def[, hotspot_id := sprintf(
  "%s|%s|%s|txpos:%s|codon:%s|%s",
  evidence_layer, gene_symbol, tx_id, transcript_position, codon_index,
  cds_frame
)]

strict_hits <- strict_run_positions[
  hotspot_run_group == "candidate_infected",
  .(
    evidence_layer = "strict_no_frame_risk",
    gene_symbol,
    tx_id,
    Run,
    transcript_position,
    cds_frame,
    position_count,
    position_count_share,
    hotspot_rank_in_run,
    run_cds_total_counts,
    run_group = hotspot_run_group
  )
]
strict_meta_cols <- intersect(
  c(
    "gene_symbol", "tx_id", "Run", "BioProject", "study", "AUTHOR",
    "CELL_LINE", "TISSUE", "CONDITION", "INHIBITOR", "TIMEPOINT",
    "pathogen", "growth_context", "sample_title", "sample_source",
    "cds_total_counts", "canonical_frame_fraction", "off_frame_fraction",
    "top1pct_off_frame_fraction", "top1pct_off_frame_enrichment",
    "target_minus_housekeeping_bumpiness_z",
    "run_off_frame_peak_score",
    "housekeeping_canonical_frame_fraction",
    "housekeeping_possible_misshift",
    "metadata_canonical_frame_fraction",
    "metadata_dominant_frame",
    "metadata_possible_misshift",
    "global_canonical_frame_fraction",
    "global_possible_misshift"
  ),
  names(strict_run_metadata)
)
strict_hits <- merge(
  strict_hits,
  strict_run_metadata[, ..strict_meta_cols],
  by = c("gene_symbol", "tx_id", "Run"),
  all.x = TRUE,
  sort = FALSE
)
frame_risk_hits <- frame_risk_run_positions[
  is_off_frame == TRUE,
  .(
    evidence_layer = "frame_risk",
    gene_symbol,
    tx_id,
    Run,
    transcript_position,
    cds_frame,
    position_count,
    position_count_share,
    hotspot_rank_in_run,
    run_cds_total_counts,
    run_group = "frame_risk_infected",
    BioProject,
    study,
    AUTHOR,
    CELL_LINE,
    TISSUE,
    CONDITION,
    INHIBITOR,
    TIMEPOINT,
    pathogen,
    growth_context,
    sample_title,
    sample_source,
    metadata_canonical_frame_fraction,
    metadata_dominant_frame,
    metadata_possible_misshift,
    global_canonical_frame_fraction,
    global_possible_misshift,
    housekeeping_possible_misshift
  )
]

hit_dt <- rbindlist(list(strict_hits, frame_risk_hits), fill = TRUE)
hit_dt <- merge(
  hit_dt,
  hotspot_def[
    ,
    .(
      hotspot_id, evidence_layer, gene_symbol, tx_id, transcript_position,
      cds_frame, codon_index, cds_relative_position, cds_percent
    )
  ],
  by = c("evidence_layer", "gene_symbol", "tx_id", "transcript_position",
         "cds_frame"),
  all.x = FALSE,
  all.y = FALSE,
  sort = FALSE
)
hit_dt <- unique(hit_dt, by = c("hotspot_id", "Run"))

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
run_order <- runIDs(df)
# collection_dir_from_exp() returns RiboCrypt's generic "collection_tables"
# dir; this project's own indexed coverage-page cache (see
# dominant_cell_state_pack_fst_pages.R's documented page_source_dir) is one
# level further, at "<that>_indexed".
fst_index <- file.path(paste0(collection_dir_from_exp(df), "_indexed"), "coverage_index.fst")
if (!file.exists(fst_index)) stop("Missing FST coverage index: ", fst_index)

diagnostics <- fread(
  file.path(analysis_dir,
            "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"),
  showProgress = FALSE
)
diagnostics[, gene_symbol := toupper(gene_symbol)]
selected_tx <- unique(hotspot_def[, .(gene_symbol, tx_id)])
selected_tx <- selected_tx[!is.na(tx_id) & nzchar(tx_id)]

cds_all <- loadRegion(df, part = "cds", names.keep = selected_tx$tx_id)
leader_all <- loadRegion(df, part = "leaders", names.keep = selected_tx$tx_id)

run_diagnostics <- list()
window_profiles <- list()

for (i in seq_len(nrow(selected_tx))) {
  gene <- selected_tx$gene_symbol[[i]]
  tx <- selected_tx$tx_id[[i]]
  gene_hits <- hit_dt[gene_symbol == gene & tx_id == tx]
  if (!nrow(gene_hits) || !tx %in% names(cds_all)) next

  message("Measuring frame diagnostics for ", gene, " (", tx, ")")
  cds_tx <- cds_all[tx]
  leader_tx <- if (tx %in% names(leader_all)) {
    leader_all[tx]
  } else {
    GRangesList()
  }
  display_region <- make_display_region(leader_tx, cds_tx, tx)
  cds_positions <- safe_pmap_positions(cds_tx, display_region)
  if (!length(cds_positions)) next

  coverage <- dominant_cached_coverage_by_transcript(
    display_region = display_region,
    fst_index = fst_index,
    gene_symbol = gene,
    tx_id = tx,
    analysis_dir = analysis_dir
  )
  mat <- as.matrix(coverage)
  storage.mode(mat) <- "numeric"
  if (!identical(colnames(mat), run_order)) {
    mat <- mat[, run_order, drop = FALSE]
  }

  measured <- rbindlist(lapply(seq_len(nrow(gene_hits)), function(j) {
    row <- gene_hits[j]
    out <- measure_run_hotspot(
      mat = mat,
      run = row$Run,
      cds_positions = cds_positions,
      hotspot_position = row$transcript_position,
      cds_relative_position = row$cds_relative_position,
      codon_index = row$codon_index,
      hotspot_frame = row$cds_frame
    )
    if (!nrow(out)) return(data.table())
    cbind(row, out[, setdiff(names(out), "Run"), with = FALSE])
  }), fill = TRUE)
  run_diagnostics[[paste(gene, tx, sep = "|")]] <- measured

  profile_hits <- gene_hits[
    evidence_layer == "frame_risk" &
      gene_symbol == "B2M" &
      codon_index == 55 &
      cds_frame == "frame1"
  ]
  if (nrow(profile_hits)) {
    cds_start <- min(cds_positions)
    cds_end <- max(cds_positions)
    profiles <- rbindlist(lapply(seq_len(nrow(profile_hits)), function(j) {
      row <- profile_hits[j]
      rbindlist(list(
        make_window_profile(mat, row$Run, cds_positions, cds_start,
                            "CDS start", 15L, 45L),
        make_window_profile(mat, row$Run, cds_positions,
                            row$transcript_position,
                            "B2M codon55 hotspot", 30L, 30L),
        make_window_profile(mat, row$Run, cds_positions, cds_end,
                            "CDS end", 60L, 0L)
      ), fill = TRUE)[
        ,
        `:=`(
          hotspot_id = row$hotspot_id,
          gene_symbol = row$gene_symbol,
          tx_id = row$tx_id,
          codon_index = row$codon_index,
          hotspot_cds_frame = row$cds_frame
        )
      ]
    }), fill = TRUE)
    window_profiles[[paste(gene, tx, sep = "|")]] <- profiles
  }
}

run_diagnostics <- rbindlist(run_diagnostics, fill = TRUE)
window_profiles <- rbindlist(window_profiles, fill = TRUE)

if (!nrow(run_diagnostics)) {
  fwrite(data.table(), file.path(output_dir,
                                 "infected_off_frame_hotspot_frame_diagnostics_runs.csv"))
  fwrite(data.table(), file.path(output_dir,
                                 "infected_off_frame_hotspot_frame_diagnostics_hotspots.csv"))
  quit(save = "no")
}

run_diagnostics[, hotspot_frame_matches_metadata_dominant :=
                  !is.na(metadata_dominant_frame) &
                  cds_frame == metadata_dominant_frame]
run_diagnostics[, start_or_end_in_frame :=
                  start_top_frame == "frame0" | end_top_frame == "frame0" |
                  start_window_dominant_frame == "frame0" |
                  end_window_dominant_frame == "frame0"]
run_diagnostics[, start_end_mixed :=
                  (!is.na(start_top_frame) & !is.na(end_top_frame) &
                     start_top_frame != end_top_frame) |
                  (!is.na(start_window_dominant_frame) &
                     !is.na(end_window_dominant_frame) &
                     start_window_dominant_frame != end_window_dominant_frame)]
run_diagnostics[, simple_pshift_like_run :=
                  hotspot_frame_matches_metadata_dominant &
                  start_end_top_same_frame &
                  start_end_dominant_same_frame &
                  start_top_frame != "frame0" &
                  end_top_frame != "frame0"]
run_diagnostics[, localized_mixed_like_run :=
                  is.finite(hotspot_peak_over_local_median) &
                  hotspot_peak_over_local_median >= 5 &
                  (start_end_mixed | end_top_frame == "frame0" |
                     end_window_dominant_frame == "frame0")]

hotspot_summary <- merge(
  hotspot_def,
  run_diagnostics[
    ,
    .(
      diagnostic_n_runs = uniqueN(Run),
      median_hotspot_position_count = safe_median(hotspot_position_count),
      median_hotspot_peak_share_31nt = safe_median(hotspot_peak_share_31nt),
      median_hotspot_peak_over_local_median =
        safe_median(hotspot_peak_over_local_median),
      median_hotspot_to_codon_frame0_ratio =
        safe_median(hotspot_to_codon_frame0_ratio),
      fraction_local_rank1 = safe_frac(hotspot_local_rank_31nt <= 1),
      fraction_hotspot_frame_matches_metadata =
        safe_frac(hotspot_frame_matches_metadata_dominant),
      fraction_start_top_frame0 = safe_frac(start_top_frame == "frame0"),
      fraction_end_top_frame0 = safe_frac(end_top_frame == "frame0"),
      fraction_start_end_top_same_frame =
        safe_frac(start_end_top_same_frame),
      fraction_start_end_dominant_same_frame =
        safe_frac(start_end_dominant_same_frame),
      fraction_start_end_mixed = safe_frac(start_end_mixed),
      fraction_simple_pshift_like_runs = safe_frac(simple_pshift_like_run),
      fraction_localized_mixed_like_runs =
        safe_frac(localized_mixed_like_run),
      start_top_offsets = collapse_text(start_top_offset_from_cds_start,
                                        max_n = 8L),
      start_top_frames = collapse_text(start_top_frame, max_n = 4L),
      end_top_offsets = collapse_text(end_top_offset_from_cds_end,
                                      max_n = 8L),
      end_top_frames = collapse_text(end_top_frame, max_n = 4L),
      metadata_dominant_frames = collapse_text(metadata_dominant_frame,
                                               max_n = 4L),
      run_examples_measured = collapse_text(Run, max_n = 8L)
    ),
    by = hotspot_id
  ],
  by = "hotspot_id",
  all.x = TRUE,
  sort = FALSE
)

hotspot_summary[, simple_pshift_support_score := rowMeans(
  cbind(
    pmin(1, pmax(0, fraction_hotspot_frame_matches_metadata)),
    pmin(1, pmax(0, fraction_start_end_top_same_frame)),
    pmin(1, pmax(0, fraction_start_end_dominant_same_frame)),
    1 - pmin(1, pmax(0, fraction_end_top_frame0))
  ),
  na.rm = TRUE
)]
hotspot_summary[, localized_mixed_support_score := rowMeans(
  cbind(
    pmin(1, pmax(0, fraction_local_rank1)),
    pmin(1, pmax(0, log1p(median_hotspot_peak_over_local_median) /
                   log1p(20))),
    pmin(1, pmax(0, fraction_start_end_mixed)),
    pmin(1, pmax(0, fraction_end_top_frame0))
  ),
  na.rm = TRUE
)]
hotspot_summary[, frame_diagnostic_class := fifelse(
  simple_pshift_support_score >= 0.65 &
    simple_pshift_support_score >
    localized_mixed_support_score + 0.10,
  "simple_global_frame_shift_like",
  fifelse(
    localized_mixed_support_score >= 0.55 &
      localized_mixed_support_score >= simple_pshift_support_score,
    "localized_or_mixed_frame_behavior",
    "ambiguous_needs_read_length_decomposition"
  )
)]
hotspot_summary[, browser_priority_note := fifelse(
  evidence_layer == "frame_risk" &
    gene_symbol == "B2M" &
    codon_index == 55 &
    cds_frame == "frame1",
  "User-validated B2M codon55 hotspot; use as frame-risk biology-vs-p-shift stress test.",
  fifelse(
    frame_diagnostic_class == "localized_or_mixed_frame_behavior",
    "Localized peak with mixed start/end-frame behavior; browser review before denoising.",
    fifelse(
      frame_diagnostic_class == "simple_global_frame_shift_like",
      "Consistent with simple frame-shift-like behavior; deprioritize unless browser contradicts.",
      "Needs read-length-specific decomposition or matched controls."
    )
  )
)]
setorder(hotspot_summary, -localized_mixed_support_score,
         simple_pshift_support_score, -source_score)

class_summary <- hotspot_summary[
  ,
  .(
    n_hotspots = .N,
    n_genes = uniqueN(gene_symbol),
    median_simple_pshift_support_score =
      safe_median(simple_pshift_support_score),
    median_localized_mixed_support_score =
      safe_median(localized_mixed_support_score),
    top_examples = collapse_text(
      paste0(gene_symbol, " codon ", codon_index, " ", cds_frame),
      max_n = 10L
    )
  ),
  by = .(evidence_layer, frame_diagnostic_class)
][order(evidence_layer, frame_diagnostic_class)]

summary_metrics <- data.table(
  metric = c(
    "n_hotspots_diagnosed",
    "n_run_hotspot_diagnostics",
    "n_strict_hotspots_diagnosed",
    "n_frame_risk_hotspots_diagnosed",
    "n_localized_or_mixed_hotspots",
    "n_simple_global_frame_shift_like_hotspots",
    "n_ambiguous_hotspots"
  ),
  value = as.character(c(
    nrow(hotspot_summary[!is.na(diagnostic_n_runs)]),
    nrow(run_diagnostics),
    hotspot_summary[evidence_layer == "strict_no_frame_risk" &
                      !is.na(diagnostic_n_runs), .N],
    hotspot_summary[evidence_layer == "frame_risk" &
                      !is.na(diagnostic_n_runs), .N],
    hotspot_summary[
      frame_diagnostic_class == "localized_or_mixed_frame_behavior", .N
    ],
    hotspot_summary[
      frame_diagnostic_class == "simple_global_frame_shift_like", .N
    ],
    hotspot_summary[
      frame_diagnostic_class ==
        "ambiguous_needs_read_length_decomposition", .N
    ]
  ))
)

fwrite(run_diagnostics, file.path(
  output_dir, "infected_off_frame_hotspot_frame_diagnostics_runs.csv"
))
fwrite(hotspot_summary, file.path(
  output_dir, "infected_off_frame_hotspot_frame_diagnostics_hotspots.csv"
))
fwrite(class_summary, file.path(
  output_dir, "infected_off_frame_hotspot_frame_diagnostics_class_summary.csv"
))
fwrite(summary_metrics, file.path(
  output_dir, "infected_off_frame_hotspot_frame_diagnostics_summary_metrics.csv"
))
if (nrow(window_profiles)) {
  fwrite(window_profiles, file.path(
    output_dir, "infected_off_frame_b2m_codon55_window_profiles.csv"
  ))
}

plot_dt <- hotspot_summary[!is.na(diagnostic_n_runs)]
if (nrow(plot_dt)) {
  plot_dt[, point_label := fifelse(
    gene_symbol == "B2M" & codon_index == 55,
    paste0(gene_symbol, " codon ", codon_index),
    fifelse(
      rank(-localized_mixed_support_score, ties.method = "first") <= 12,
      paste0(gene_symbol, " c", codon_index),
      ""
    )
  )]
  p <- ggplot(
    plot_dt,
    aes(
      x = simple_pshift_support_score,
      y = localized_mixed_support_score,
      size = pmax(source_run_hits, 1),
      color = evidence_layer,
      shape = frame_diagnostic_class,
      text = paste0(
        gene_symbol, " codon ", codon_index, " ", cds_frame,
        "<br>layer=", evidence_layer,
        "<br>class=", frame_diagnostic_class,
        "<br>simple score=", round(simple_pshift_support_score, 3),
        "<br>localized/mixed score=",
        round(localized_mixed_support_score, 3),
        "<br>runs=", diagnostic_n_runs,
        "<br>start frames=", start_top_frames,
        "<br>end frames=", end_top_frames,
        "<br>metadata frames=", metadata_dominant_frames
      )
    )
  ) +
    geom_hline(yintercept = 0.55, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0.65, linetype = "dashed", color = "grey70") +
    geom_point(alpha = 0.82) +
    geom_text(
      aes(label = point_label),
      hjust = -0.05,
      vjust = 0.5,
      size = 2.6,
      show.legend = FALSE,
      check_overlap = TRUE
    ) +
    scale_color_manual(values = c(
      strict_no_frame_risk = "#2166ac",
      frame_risk = "#8b0000"
    )) +
    scale_size_continuous(range = c(2.2, 8)) +
    coord_cartesian(xlim = c(0, 1.04), ylim = c(0, 1.04),
                    clip = "off") +
    labs(
      title = "Off-Frame Hotspot Frame Diagnostics",
      subtitle = paste(
        "Simple-shift score uses metadata/start/end frame agreement;",
        "localized score rewards narrow recurrent peaks with mixed or in-frame end behavior."
      ),
      x = "Simple global frame-shift-like support",
      y = "Localized or mixed behavior support",
      color = "Evidence layer",
      shape = "Diagnostic class",
      size = "Run hits"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right",
          plot.margin = margin(8, 32, 8, 8))
  ggsave(file.path(figure_dir,
                   "infected_off_frame_hotspot_frame_diagnostic_quadrants.png"),
         p, width = 10.5, height = 7, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir,
                "infected_off_frame_hotspot_frame_diagnostic_quadrants.html"),
      title = "infected_off_frame_hotspot_frame_diagnostic_quadrants"
    )
  }
}

if (nrow(window_profiles)) {
  b2m_runs <- unique(window_profiles$Run)
  keep_runs <- head(b2m_runs, 6L)
  b2m_plot <- window_profiles[Run %chin% keep_runs]
  b2m_plot[, window := factor(window,
                              levels = c("CDS start",
                                         "B2M codon55 hotspot",
                                         "CDS end"))]
  q <- ggplot(
    b2m_plot,
    aes(x = offset, y = count, fill = cds_frame,
        text = paste0(
          "Run=", Run,
          "<br>window=", window,
          "<br>offset=", offset,
          "<br>tx pos=", transcript_position,
          "<br>frame=", cds_frame,
          "<br>count=", count
        ))
  ) +
    geom_col(width = 0.9, color = NA) +
    facet_grid(Run ~ window, scales = "free_x") +
    scale_fill_manual(values = c(
      frame0 = "#2166ac",
      frame1 = "#b2182b",
      frame2 = "#1b9e77"
    )) +
    labs(
      title = "B2M Codon 55 Frame-Risk Diagnostic Windows",
      subtitle = "Aggregate p-shifted FST coverage; start, hotspot, and CDS-end windows are shown separately.",
      x = "Offset from window anchor",
      y = "Coverage count",
      fill = "CDS frame"
    ) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(figure_dir,
                   "infected_off_frame_b2m_codon55_diagnostic_windows.png"),
         q, width = 12, height = 7, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      q,
      file.path(figure_dir,
                "infected_off_frame_b2m_codon55_diagnostic_windows.html"),
      title = "infected_off_frame_b2m_codon55_diagnostic_windows"
    )
  }
}

message("Saved frame diagnostic outputs to: ", output_dir)
message("Diagnostic class summary:")
print(class_summary)
message("Top localized/mixed hotspots:")
print(head(hotspot_summary[
  ,
  .(
    evidence_layer, gene_symbol, codon_index, cds_frame,
    frame_diagnostic_class, diagnostic_n_runs,
    simple_pshift_support_score, localized_mixed_support_score,
    start_top_frames, end_top_frames, metadata_dominant_frames,
    browser_priority_note
  )
], 25L))
