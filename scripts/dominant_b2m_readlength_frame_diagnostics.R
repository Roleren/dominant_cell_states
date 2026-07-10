#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
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
input_dir <- file.path(analysis_dir, "readlength_inputs")
output_dir <- file.path(analysis_dir,
                        "dominant_infected_off_frame_readlength_diagnostics")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

htmlwidget_helper <- file.path(analysis_dir, "scripts", "dominant_htmlwidgets.R")
if (file.exists(htmlwidget_helper)) source(htmlwidget_helper)

message("B2M read-length frame diagnostics")
message("  1. Read collapsed BAMs with ORFik::readBam() score weights")
message("  2. Apply study/run-specific p-shift tables by read length")
message("  3. Measure B2M CDS-start, codon55, and CDS-end frame behavior")

target_runs <- data.table(
  Run = c("SRR2052924", "SRR27450688"),
  BioProject = c("PRJNA285961", "PRJNA1062064"),
  bam_file = file.path(input_dir, c("SRR2052924.bam", "SRR27450688.bam")),
  shift_file = file.path(input_dir, c(
    "PRJNA285961_shifting_tables.rds",
    "PRJNA1062064_shifting_tables.rds"
  ))
)

summary_out <- file.path(
  output_dir, "b2m_readlength_frame_diagnostics_summary_metrics.csv"
)
missing_inputs <- target_runs[
  !file.exists(bam_file) | !file.exists(shift_file)
]
if (nrow(missing_inputs)) {
  fwrite(missing_inputs, file.path(output_dir,
                                   "b2m_readlength_missing_inputs.csv"))
  fwrite(data.table(
    metric = c("status", "n_missing_inputs"),
    value = c("missing_readlength_inputs", as.character(nrow(missing_inputs)))
  ), summary_out)
  message("Missing read-length inputs; wrote missing-input summary.")
  quit(save = "no", status = 0)
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

collapse_text <- function(x, max_n = 8L) {
  x <- unique(as.character(x[!is.na(x)]))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  suffix <- if (length(x) > max_n) ";..." else ""
  paste0(paste(head(x, max_n), collapse = ";"), suffix)
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

load_shift_table_for_run <- function(shift_file, run) {
  shift_list <- readRDS(shift_file)
  if (!is.list(shift_list)) {
    stop("Shift file is not a list: ", shift_file, call. = FALSE)
  }
  key <- names(shift_list)[basename(names(shift_list)) == paste0(run, ".ofst")]
  if (length(key) != 1L) {
    stop("Could not find exactly one shift-table entry for ", run,
         " in ", shift_file, ". Found: ", length(key), call. = FALSE)
  }
  shifts <- as.data.table(shift_list[[key]])
  if (!all(c("fraction", "offsets_start") %in% names(shifts))) {
    stop("Shift table for ", run,
         " must contain fraction and offsets_start.", call. = FALSE)
  }
  shifts[, `:=`(
    Run = run,
    source_shift_entry = key,
    fraction = as.integer(fraction),
    offsets_start = as.integer(offsets_start)
  )]
  shifts[]
}

make_bam_which <- function(display_region, flank = 120L) {
  gr <- unlist(display_region, use.names = FALSE)
  GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(pmax(1L, start(gr) - flank), end(gr) + flank),
    strand = "*"
  ) |>
    reduce()
}

ensure_bam_index <- function(bam_file) {
  index_file <- paste0(bam_file, ".bai")
  if (file.exists(index_file)) return(index_file)
  message("Indexing BAM: ", bam_file)
  out <- tryCatch(
    Rsamtools::indexBam(bam_file),
    error = function(e) {
      warning("Could not index BAM ", bam_file, ": ", conditionMessage(e))
      NA_character_
    }
  )
  as.character(out)
}

read_collapsed_bam_region <- function(bam_file, which_gr) {
  ensure_bam_index(bam_file)
  param <- ScanBamParam(
    which = which_gr,
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )
  reads <- ORFik::readBam(bam_file, param = param)
  if (!"score" %in% colnames(mcols(reads))) {
    mcols(reads)$score <- rep(1L, length(reads))
  }
  reads
}

coverage_by_read_length <- function(display_region, shifted_reads) {
  lengths <- sort(unique(as.integer(shifted_reads$size)))
  lengths <- lengths[is.finite(lengths)]
  rbindlist(lapply(lengths, function(len) {
    len_reads <- shifted_reads[shifted_reads$size == len]
    if (!length(len_reads)) return(data.table())
    ORFik::coveragePerTiling(
      display_region,
      len_reads,
      as.data.table = TRUE,
      keep.names = FALSE,
      weight = "score",
      drop.zero.dt = FALSE,
      fraction = len
    )
  }), fill = TRUE)
}

label_windows <- function(profile, cds_start, cds_end, hotspot_position = 194L) {
  window_defs <- data.table(
    window = c("CDS start", "B2M codon55 hotspot", "CDS end"),
    anchor = c(cds_start, hotspot_position, cds_end),
    left = c(15L, 30L, 60L),
    right = c(45L, 30L, 0L)
  )
  rbindlist(lapply(seq_len(nrow(window_defs)), function(i) {
    row <- window_defs[i]
    profile[
      position >= row$anchor - row$left &
        position <= row$anchor + row$right
    ][
      ,
      `:=`(
        window = row$window,
        window_anchor = row$anchor,
        window_offset = position - row$anchor
      )
    ]
  }), fill = TRUE)
}

shift_sensitivity <- function(display_region, shifted_reads, run, cds_start,
                              hotspot_position = 194L, deltas = -6L:6L) {
  rbindlist(lapply(deltas, function(delta) {
    shifted_delta <- GenomicRanges::shift(shifted_reads, shift = delta)
    cov <- coverage_by_read_length(display_region, shifted_delta)
    if (!nrow(cov)) return(data.table())
    cov[, `:=`(
      Run = run,
      psite_delta = delta,
      cds_frame = paste0("frame", as.integer((position - cds_start) %% 3L))
    )]
    codon_positions <- (hotspot_position - 1L):(hotspot_position + 1L)
    local_positions <- (hotspot_position - 30L):(hotspot_position + 30L)
    codon <- cov[position %in% codon_positions,
                 .(count = sum(count, na.rm = TRUE)),
                 by = .(Run, psite_delta, position, cds_frame)]
    codon_wide <- dcast(
      codon,
      Run + psite_delta ~ cds_frame,
      value.var = "count",
      fill = 0
    )
    for (frame in paste0("frame", 0:2)) {
      if (!frame %in% names(codon_wide)) codon_wide[, (frame) := 0]
    }
    local <- cov[position %in% local_positions,
                 .(count = sum(count, na.rm = TRUE)),
                 by = .(Run, psite_delta, cds_frame)]
    local_wide <- dcast(
      local,
      Run + psite_delta ~ cds_frame,
      value.var = "count",
      fill = 0
    )
    for (frame in paste0("frame", 0:2)) {
      if (!frame %in% names(local_wide)) local_wide[, (frame) := 0]
    }
    local_wide[, local_total := frame0 + frame1 + frame2]
    local_wide[, `:=`(
      local_frame0_fraction = frame0 / pmax(local_total, 1),
      local_frame1_fraction = frame1 / pmax(local_total, 1),
      local_frame2_fraction = frame2 / pmax(local_total, 1)
    )]
    setnames(codon_wide, paste0("frame", 0:2),
             paste0("codon55_", paste0("frame", 0:2), "_count"))
    setnames(local_wide, paste0("frame", 0:2),
             paste0("local_", paste0("frame", 0:2), "_count"))
    pos194 <- cov[position == hotspot_position,
                  .(position194_count = sum(count, na.rm = TRUE)),
                  by = .(Run, psite_delta)]
    merge(merge(codon_wide, local_wide, by = c("Run", "psite_delta"),
                all = TRUE),
          pos194, by = c("Run", "psite_delta"), all = TRUE)
  }), fill = TRUE)
}

df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
tx_id <- "ENST00000648006"
gene_symbol <- "B2M"
cds_tx <- loadRegion(df, part = "cds", names.keep = tx_id)
leader_tx <- loadRegion(df, part = "leaders", names.keep = tx_id)
display_region <- make_display_region(leader_tx, cds_tx, tx_id)
cds_positions <- safe_pmap_positions(cds_tx, display_region)
if (!length(cds_positions)) {
  stop("Could not map B2M CDS to transcript coordinates.", call. = FALSE)
}
cds_start <- min(cds_positions)
cds_end <- max(cds_positions)
hotspot_position <- 194L
codon55_positions <- (hotspot_position - 1L):(hotspot_position + 1L)
which_gr <- make_bam_which(display_region)

all_profiles <- list()
all_window_profiles <- list()
all_window_summary <- list()
all_sensitivity <- list()
all_read_summary <- list()
shift_tables_used <- list()

for (i in seq_len(nrow(target_runs))) {
  row <- target_runs[i]
  run <- row$Run
  message("Processing ", run)
  shifts <- load_shift_table_for_run(row$shift_file, run)
  shift_tables_used[[run]] <- shifts

  reads <- read_collapsed_bam_region(row$bam_file, which_gr)
  read_widths <- ORFik::readWidths(reads, along.reference = TRUE)
  total_score <- sum(mcols(reads)$score, na.rm = TRUE)
  message("  regional collapsed alignments: ", length(reads),
          "; score-sum=", total_score)

  shifted <- ORFik::shiftFootprints(reads, shifts, sort = FALSE)
  profile <- coverage_by_read_length(display_region, shifted)
  if (!nrow(profile)) next
  profile[, `:=`(
    Run = run,
    BioProject = row$BioProject,
    gene_symbol = gene_symbol,
    tx_id = tx_id,
    read_length = as.integer(fraction),
    cds_frame = paste0("frame", as.integer((position - cds_start) %% 3L)),
    is_cds = position %between% c(cds_start, cds_end),
    codon_index = fifelse(
      position %between% c(cds_start, cds_end),
      floor((position - cds_start) / 3L) + 1L,
      NA_integer_
    )
  )]
  profile <- profile[
    ,
    .(Run, BioProject, gene_symbol, tx_id, read_length, position,
      count, cds_frame, is_cds, codon_index)
  ]
  all_profiles[[run]] <- profile

  window_profile <- label_windows(profile, cds_start, cds_end,
                                  hotspot_position)
  all_window_profiles[[run]] <- window_profile

  window_summary <- window_profile[
    ,
    .(count = sum(count, na.rm = TRUE)),
    by = .(Run, BioProject, window, read_length, cds_frame)
  ]
  window_summary[, window_read_length_total := sum(count),
                 by = .(Run, window, read_length)]
  window_summary[, frame_fraction :=
                   count / pmax(window_read_length_total, 1)]
  all_window_summary[[run]] <- window_summary

  sens <- shift_sensitivity(
    display_region = display_region,
    shifted_reads = shifted,
    run = run,
    cds_start = cds_start,
    hotspot_position = hotspot_position
  )
  sens[, BioProject := row$BioProject]
  all_sensitivity[[run]] <- sens

  read_summary <- data.table(
    Run = run,
    BioProject = row$BioProject,
    bam_file = row$bam_file,
    bam_size_bytes = file.info(row$bam_file)$size,
    n_regional_collapsed_alignments = length(reads),
    regional_collapsed_score_sum = total_score,
    n_shifted_psites = length(shifted),
    shifted_score_sum = sum(mcols(shifted)$score, na.rm = TRUE),
    raw_read_lengths = collapse_text(names(table(read_widths))),
    shifted_read_lengths = collapse_text(sort(unique(shifted$size))),
    cds_start_transcript_position = cds_start,
    cds_end_transcript_position = cds_end,
    b2m_codon55_frame0_transcript_position = codon55_positions[[1]],
    b2m_codon55_frame1_transcript_position = codon55_positions[[2]],
    b2m_codon55_frame2_transcript_position = codon55_positions[[3]]
  )
  all_read_summary[[run]] <- read_summary
}

profiles <- rbindlist(all_profiles, fill = TRUE)
window_profiles <- rbindlist(all_window_profiles, fill = TRUE)
window_summary <- rbindlist(all_window_summary, fill = TRUE)
offset_sensitivity <- rbindlist(all_sensitivity, fill = TRUE)
read_summary <- rbindlist(all_read_summary, fill = TRUE)
shift_tables_used <- rbindlist(shift_tables_used, fill = TRUE)

codon55_by_length <- profiles[
  position %in% codon55_positions,
  .(count = sum(count, na.rm = TRUE)),
  by = .(Run, BioProject, read_length, position, cds_frame)
]
codon55_summary <- codon55_by_length[
  ,
  .(
    total_codon55_count = sum(count, na.rm = TRUE),
    frame0_count = sum(count[cds_frame == "frame0"], na.rm = TRUE),
    frame1_count = sum(count[cds_frame == "frame1"], na.rm = TRUE),
    frame2_count = sum(count[cds_frame == "frame2"], na.rm = TRUE),
    dominant_codon55_read_length =
      read_length[which.max(count)][1],
    dominant_codon55_position = position[which.max(count)][1],
    dominant_codon55_frame = cds_frame[which.max(count)][1],
    dominant_codon55_count = max(count, na.rm = TRUE)
  ),
  by = .(Run, BioProject)
]
codon55_summary[, `:=`(
  frame0_fraction = frame0_count / pmax(total_codon55_count, 1),
  frame1_fraction = frame1_count / pmax(total_codon55_count, 1),
  frame2_fraction = frame2_count / pmax(total_codon55_count, 1),
  frame1_over_frame0 = frame1_count / pmax(frame0_count, 1)
)]

bam_codon55_position_counts <- profiles[
  position %in% codon55_positions,
  .(
    bam_default_pshift_count = sum(count, na.rm = TRUE),
    bam_default_dominant_read_length =
      read_length[which.max(count)][1]
  ),
  by = .(Run, BioProject, transcript_position = position, cds_frame)
]
fst_run_positions_file <- file.path(
  analysis_dir,
  "dominant_infected_off_frame_frame_risk",
  "infected_off_frame_frame_risk_run_positions.csv"
)
fst_vs_bam <- data.table()
if (file.exists(fst_run_positions_file)) {
  fst_run_positions <- fread(fst_run_positions_file, showProgress = FALSE)
  fst_run_positions[, gene_symbol := toupper(gene_symbol)]
  fst_codon55 <- fst_run_positions[
    gene_symbol == "B2M" &
      Run %chin% target_runs$Run &
      codon_index == 55 &
      transcript_position %in% codon55_positions,
    .(
      fst_aggregate_position_count = sum(position_count, na.rm = TRUE),
      fst_hotspot_rank_in_run = min(hotspot_rank_in_run, na.rm = TRUE)
    ),
    by = .(BioProject, gene_symbol, tx_id, Run, transcript_position,
           cds_frame)
  ]
  fst_vs_bam <- merge(
    fst_codon55,
    bam_codon55_position_counts,
    by = c("Run", "BioProject", "transcript_position", "cds_frame"),
    all = TRUE,
    sort = FALSE
  )
  fst_vs_bam[is.na(gene_symbol), `:=`(
    gene_symbol = "B2M",
    tx_id = "ENST00000648006"
  )]
  fst_vs_bam[, count_delta_bam_minus_fst :=
               bam_default_pshift_count - fst_aggregate_position_count]
  setorder(fst_vs_bam, Run, transcript_position, cds_frame)
}

summary_metrics <- rbind(
  data.table(
    metric = c(
      "status",
      "n_runs_processed",
      "n_window_profile_rows",
      "n_offset_sensitivity_rows"
    ),
    value = c(
      "ok",
      as.character(uniqueN(profiles$Run)),
      as.character(nrow(window_profiles)),
      as.character(nrow(offset_sensitivity))
    )
  ),
  codon55_summary[
    ,
    .(
      metric = paste0(Run, "_codon55_", c(
        "total_count", "frame0_fraction", "frame1_fraction",
        "frame2_fraction", "dominant_read_length", "dominant_position",
        "dominant_frame", "frame1_over_frame0"
      )),
      value = as.character(c(
        total_codon55_count, round(frame0_fraction, 4),
        round(frame1_fraction, 4), round(frame2_fraction, 4),
        dominant_codon55_read_length, dominant_codon55_position,
        dominant_codon55_frame, round(frame1_over_frame0, 4)
      ))
    ),
    by = Run
  ][, .(metric, value)]
)

fwrite(shift_tables_used, file.path(output_dir,
                                    "b2m_readlength_shift_tables_used.csv"))
fwrite(read_summary, file.path(output_dir,
                               "b2m_readlength_run_summary.csv"))
fwrite(profiles, file.path(output_dir,
                           "b2m_readlength_psite_profile.csv"))
fwrite(window_profiles, file.path(output_dir,
                                  "b2m_readlength_window_profiles.csv"))
fwrite(window_summary, file.path(output_dir,
                                 "b2m_readlength_window_frame_summary.csv"))
fwrite(codon55_by_length, file.path(output_dir,
                                    "b2m_codon55_readlength_counts.csv"))
fwrite(codon55_summary, file.path(output_dir,
                                  "b2m_codon55_readlength_summary.csv"))
fwrite(offset_sensitivity, file.path(output_dir,
                                     "b2m_codon55_psite_offset_sensitivity.csv"))
if (nrow(fst_vs_bam)) {
  fwrite(fst_vs_bam, file.path(output_dir,
                               "b2m_codon55_fst_vs_bam_pshift_counts.csv"))
}
fwrite(summary_metrics, summary_out)

if (nrow(window_profiles)) {
  plot_window <- window_profiles[
    window == "B2M codon55 hotspot" & count > 0
  ]
  p <- ggplot(
    plot_window,
    aes(
      x = window_offset,
      y = factor(read_length),
      fill = count,
      text = paste0(
        "Run=", Run,
        "<br>read length=", read_length,
        "<br>offset=", window_offset,
        "<br>tx position=", position,
        "<br>frame=", cds_frame,
        "<br>count=", count
      )
    )
  ) +
    geom_tile() +
    geom_vline(xintercept = 0, color = "white", linewidth = 0.45) +
    facet_grid(Run ~ cds_frame, scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "white", high = "#8b0000") +
    labs(
      title = "B2M Codon 55 Read-Length P-Site Signal",
      subtitle = "Collapsed BAM reads are weighted by ORFik::readBam() score and shifted by the study/run shift table.",
      x = "Offset from transcript position 194",
      y = "Read length",
      fill = "P-site count"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(figure_dir, "b2m_codon55_readlength_heatmap.png"),
         p, width = 11.5, height = 6.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "b2m_codon55_readlength_heatmap.html"),
      title = "b2m_codon55_readlength_heatmap"
    )
  }

  top_lengths <- codon55_by_length[
    ,
    .(codon55_count = sum(count)),
    by = .(Run, read_length)
  ][order(Run, -codon55_count)][
    ,
    head(.SD, 8L),
    by = Run
  ]$read_length
  frame_plot <- window_summary[read_length %in% top_lengths]
  q <- ggplot(
    frame_plot,
    aes(
      x = factor(read_length),
      y = frame_fraction,
      fill = cds_frame,
      text = paste0(
        "Run=", Run,
        "<br>window=", window,
        "<br>read length=", read_length,
        "<br>frame=", cds_frame,
        "<br>fraction=", round(frame_fraction, 3),
        "<br>count=", count
      )
    )
  ) +
    geom_col(position = "stack") +
    facet_grid(Run ~ window) +
    scale_fill_manual(values = c(
      frame0 = "#2166ac",
      frame1 = "#b2182b",
      frame2 = "#1b9e77"
    )) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "B2M Read-Length Frame Fractions by Window",
      subtitle = "CDS-start, codon55, and CDS-end windows after default p-shifting.",
      x = "Read length",
      y = "Frame fraction",
      fill = "CDS frame"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figure_dir, "b2m_readlength_window_frame_fractions.png"),
         q, width = 12, height = 6.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      q,
      file.path(figure_dir,
                "b2m_readlength_window_frame_fractions.html"),
      title = "b2m_readlength_window_frame_fractions"
    )
  }
}

if (nrow(offset_sensitivity)) {
  offset_long <- melt(
    offset_sensitivity,
    id.vars = c("Run", "BioProject", "psite_delta"),
    measure.vars = c("codon55_frame0_count", "codon55_frame1_count",
                     "codon55_frame2_count"),
    variable.name = "codon55_frame",
    value.name = "count"
  )
  offset_long[, codon55_frame := sub("^codon55_", "", codon55_frame)]
  offset_long[, codon55_frame := sub("_count$", "", codon55_frame)]
  r <- ggplot(
    offset_long,
    aes(
      x = psite_delta,
      y = count,
      color = codon55_frame,
      group = codon55_frame,
      text = paste0(
        "Run=", Run,
        "<br>delta=", psite_delta,
        "<br>frame=", codon55_frame,
        "<br>count=", count
      )
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ Run, scales = "free_y") +
    scale_color_manual(values = c(
      frame0 = "#2166ac",
      frame1 = "#b2182b",
      frame2 = "#1b9e77"
    )) +
    labs(
      title = "B2M Codon 55 P-Site Offset Sensitivity",
      subtitle = "Delta shifts the default p-site coordinate before recounting codon55 frame positions.",
      x = "P-site coordinate delta from default",
      y = "Codon55 count",
      color = "Codon55 frame"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(figure_dir, "b2m_codon55_psite_offset_sensitivity.png"),
         r, width = 9.5, height = 5.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      r,
      file.path(figure_dir,
                "b2m_codon55_psite_offset_sensitivity.html"),
      title = "b2m_codon55_psite_offset_sensitivity"
    )
  }
}

message("Saved B2M read-length diagnostics to: ", output_dir)
message("Codon55 summary:")
print(codon55_summary)
