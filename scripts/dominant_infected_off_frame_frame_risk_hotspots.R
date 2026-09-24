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
input_dir <- file.path(analysis_dir, "dominant_infected_off_frame_cds_scan")
output_dir <- file.path(analysis_dir, "dominant_infected_off_frame_frame_risk")
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

message("Infected off-frame frame-risk hotspot review")
message("  1. Keep strong off-frame peaks removed by all-CDS frame QC")
message("  2. Localize exact hotspot positions for browser review")
message("  3. Keep these separate from strict no-frame-bias candidates")

frame_risk_file <- file.path(
  input_dir, "infected_off_frame_cds_frame_risk_runs.csv"
)
diagnostics_file <- file.path(
  analysis_dir, "human_dominant_cell_states_clean_cds_gene_diagnostics.csv"
)
for (file in c(frame_risk_file, diagnostics_file)) {
  if (!file.exists(file)) stop("Missing required file: ", file, call. = FALSE)
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

middle_cds_positions <- function(cds_positions, min_edge = 45L,
                                 fraction_edge = 0.1) {
  cds_positions <- sort(unique(as.integer(cds_positions)))
  n <- length(cds_positions)
  if (n < 120L) return(cds_positions)
  edge <- max(as.integer(min_edge), floor(n * fraction_edge))
  if (n <= 2L * edge + 30L) {
    edge <- max(15L, floor(n * 0.05))
  }
  cds_positions[(edge + 1L):(n - edge)]
}

per_run_top_positions <- function(mat, runs, cds_positions, use_positions,
                                  gene, tx_id, top_fraction = 0.01,
                                  max_positions_per_run = 30L) {
  runs <- intersect(runs, colnames(mat))
  if (!length(runs)) return(data.table())
  cds_positions <- sort(unique(as.integer(cds_positions)))
  use_positions <- sort(unique(as.integer(use_positions)))
  use_positions <- use_positions[
    use_positions >= 1L & use_positions <= nrow(mat)
  ]
  if (!length(use_positions)) return(data.table())
  full_rel <- match(use_positions, cds_positions)
  frame <- as.integer((use_positions - min(cds_positions)) %% 3L)
  rbindlist(lapply(runs, function(run) {
    counts <- as.numeric(mat[use_positions, run, drop = TRUE])
    total <- sum(counts, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) return(data.table())
    q99 <- as.numeric(stats::quantile(counts, probs = 1 - top_fraction,
                                      na.rm = TRUE, names = FALSE))
    idx <- which(is.finite(counts) & counts > 0 & counts >= q99)
    if (!length(idx)) return(data.table())
    idx <- idx[order(counts[idx], decreasing = TRUE)]
    idx <- head(idx, max_positions_per_run)
    data.table(
      gene_symbol = gene,
      tx_id = tx_id,
      Run = run,
      cds_relative_position = full_rel[idx],
      cds_percent = full_rel[idx] / length(cds_positions),
      codon_index = floor((full_rel[idx] - 1L) / 3L) + 1L,
      transcript_position = use_positions[idx],
      cds_frame = paste0("frame", frame[idx]),
      is_off_frame = frame[idx] != 0L,
      position_count = counts[idx],
      position_count_share = counts[idx] / total,
      hotspot_rank_in_run = seq_along(idx),
      run_cds_total_counts = total,
      top_threshold_count = q99
    )
  }), fill = TRUE)
}

aggregate_profile <- function(mat, runs, cds_positions, use_positions,
                              gene, tx_id) {
  runs <- intersect(runs, colnames(mat))
  if (!length(runs)) return(data.table())
  cds_positions <- sort(unique(as.integer(cds_positions)))
  use_positions <- sort(unique(as.integer(use_positions)))
  use_positions <- use_positions[
    use_positions >= 1L & use_positions <= nrow(mat)
  ]
  if (!length(use_positions)) return(data.table())
  counts <- as.numeric(rowSums(mat[use_positions, runs, drop = FALSE],
                               na.rm = TRUE))
  total <- sum(counts, na.rm = TRUE)
  rel <- match(use_positions, cds_positions)
  data.table(
    gene_symbol = gene,
    tx_id = tx_id,
    n_runs = length(runs),
    cds_relative_position = rel,
    cds_percent = rel / length(cds_positions),
    codon_index = floor((rel - 1L) / 3L) + 1L,
    transcript_position = use_positions,
    cds_frame = paste0("frame", (use_positions - min(cds_positions)) %% 3L),
    counts = counts,
    normalized_density = if (is.finite(total) && total > 0) {
      counts / total * 1000
    } else {
      NA_real_
    }
  )
}

frame_risk_runs <- fread(frame_risk_file, showProgress = FALSE)
if (!nrow(frame_risk_runs)) {
  fwrite(data.table(), file.path(output_dir,
                                 "infected_off_frame_frame_risk_hotspots.csv"))
  fwrite(data.table(), file.path(output_dir,
                                 "infected_off_frame_frame_risk_gene_summary.csv"))
  quit(save = "no")
}
frame_risk_runs[, gene_symbol := toupper(gene_symbol)]
fwrite(frame_risk_runs, file.path(output_dir,
                                  "infected_off_frame_frame_risk_runs.csv"))

candidate_genes <- unique(frame_risk_runs$gene_symbol)
diagnostics <- fread(diagnostics_file, showProgress = FALSE)
diagnostics[, gene_symbol := toupper(gene_symbol)]
diagnostics <- diagnostics[
  expression_source == "fst_clean_cds" &
    gene_symbol %chin% candidate_genes &
    !is.na(tx_id) & nzchar(tx_id)
]
selected_tx <- unique(diagnostics[, .(gene_symbol, tx_id)])
selected_tx <- selected_tx[match(candidate_genes, gene_symbol), ]
selected_tx <- selected_tx[!is.na(gene_symbol)]

df <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
run_order <- runIDs(df)
# collection_dir_from_exp() returns RiboCrypt's generic "collection_tables"
# dir; this project's own indexed coverage-page cache (see
# dominant_cell_state_pack_fst_pages.R's documented page_source_dir) is one
# level further, at "<that>_indexed".
fst_index <- file.path(paste0(collection_dir_from_exp(df), "_indexed"), "coverage_index.fst")
if (!file.exists(fst_index)) stop("Missing FST coverage index: ", fst_index)

cds_all <- loadRegion(df, part = "cds", names.keep = selected_tx$tx_id)
leader_all <- loadRegion(df, part = "leaders", names.keep = selected_tx$tx_id)

all_hotspots <- list()
all_profiles <- list()

for (i in seq_len(nrow(selected_tx))) {
  gene <- selected_tx$gene_symbol[[i]]
  tx_id <- selected_tx$tx_id[[i]]
  runs <- frame_risk_runs[gene_symbol == gene, unique(Run)]
  if (!length(runs) || !tx_id %in% names(cds_all)) next
  message("Localizing frame-risk hotspots for ", gene, " (", tx_id, ")")
  cds_tx <- cds_all[tx_id]
  leader_tx <- if (tx_id %in% names(leader_all)) {
    leader_all[tx_id]
  } else {
    GRangesList()
  }
  display_region <- make_display_region(leader_tx, cds_tx, tx_id)
  cds_positions <- safe_pmap_positions(cds_tx, display_region)
  use_positions <- middle_cds_positions(cds_positions)
  if (!length(use_positions)) next

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

  all_hotspots[[gene]] <- per_run_top_positions(
    mat = mat,
    runs = runs,
    cds_positions = cds_positions,
    use_positions = use_positions,
    gene = gene,
    tx_id = tx_id
  )
  all_profiles[[gene]] <- aggregate_profile(
    mat = mat,
    runs = runs,
    cds_positions = cds_positions,
    use_positions = use_positions,
    gene = gene,
    tx_id = tx_id
  )
}

hotspots <- rbindlist(all_hotspots, fill = TRUE)
profiles <- rbindlist(all_profiles, fill = TRUE)
if (!nrow(hotspots)) {
  fwrite(data.table(), file.path(output_dir,
                                 "infected_off_frame_frame_risk_hotspots.csv"))
  fwrite(data.table(), file.path(output_dir,
                                 "infected_off_frame_frame_risk_gene_summary.csv"))
  quit(save = "no")
}

hotspots <- merge(
  hotspots,
  frame_risk_runs,
  by = c("gene_symbol", "tx_id", "Run"),
  all.x = TRUE,
  sort = FALSE
)

off_hotspots <- hotspots[is_off_frame == TRUE]
hotspot_summary <- off_hotspots[
  ,
  .(
    n_frame_risk_runs_with_hotspot = uniqueN(Run),
    total_position_counts = sum(position_count, na.rm = TRUE),
    median_position_count_share = safe_median(position_count_share),
    max_position_count_share = safe_max(position_count_share),
    best_rank_in_run = min(hotspot_rank_in_run, na.rm = TRUE),
    run_examples = collapse_text(Run[order(hotspot_rank_in_run)],
                                 max_n = 8L),
    contexts = collapse_text(paste(BioProject, CELL_LINE, CONDITION,
                                   TIMEPOINT, sep = " | "),
                             max_n = 6L),
    metadata_dominant_frames = collapse_text(metadata_dominant_frame,
                                             max_n = 4L),
    median_metadata_canonical_frame_fraction =
      safe_median(metadata_canonical_frame_fraction),
    median_global_canonical_frame_fraction =
      safe_median(global_canonical_frame_fraction)
  ),
  by = .(gene_symbol, tx_id, cds_relative_position, cds_percent,
         codon_index, transcript_position, cds_frame)
]
group_sizes <- frame_risk_runs[, .(group_n_frame_risk_runs = uniqueN(Run)),
                               by = gene_symbol]
hotspot_summary <- merge(hotspot_summary, group_sizes, by = "gene_symbol",
                         all.x = TRUE, sort = FALSE)
hotspot_summary[, fraction_frame_risk_runs_with_hotspot :=
                  n_frame_risk_runs_with_hotspot /
                  pmax(group_n_frame_risk_runs, 1)]
hotspot_summary[, frame_risk_hotspot_score :=
                  log1p(n_frame_risk_runs_with_hotspot) +
                  pmax(0, fraction_frame_risk_runs_with_hotspot) +
                  pmax(0, median_position_count_share)]
setorder(hotspot_summary, -frame_risk_hotspot_score,
         -n_frame_risk_runs_with_hotspot, gene_symbol)

gene_summary <- hotspot_summary[
  ,
  .(
    n_frame_risk_hotspots = .N,
    n_recurrent_frame_risk_hotspots =
      sum(n_frame_risk_runs_with_hotspot >= 2),
    best_frame_risk_hotspot_score = max(frame_risk_hotspot_score,
                                        na.rm = TRUE),
    best_codon_index = codon_index[which.max(frame_risk_hotspot_score)],
    best_cds_frame = cds_frame[which.max(frame_risk_hotspot_score)],
    best_transcript_position =
      transcript_position[which.max(frame_risk_hotspot_score)],
    best_fraction_frame_risk_runs_with_hotspot =
      fraction_frame_risk_runs_with_hotspot[which.max(frame_risk_hotspot_score)],
    best_run_examples = run_examples[which.max(frame_risk_hotspot_score)],
    metadata_dominant_frames = collapse_text(metadata_dominant_frames,
                                             max_n = 4L)
  ),
  by = .(gene_symbol, tx_id)
]
run_summary <- frame_risk_runs[
  ,
  .(
    n_frame_risk_runs = uniqueN(Run),
    median_top1pct_off_frame_fraction =
      safe_median(top1pct_off_frame_fraction),
    median_off_frame_enrichment = safe_median(top1pct_off_frame_enrichment),
    median_metadata_canonical_frame_fraction =
      safe_median(metadata_canonical_frame_fraction),
    run_metadata_dominant_frames = collapse_text(metadata_dominant_frame,
                                                 max_n = 4L),
    run_examples = collapse_text(Run[order(-run_off_frame_peak_score)],
                                 max_n = 8L)
  ),
  by = .(gene_symbol, tx_id)
]
gene_summary <- merge(gene_summary, run_summary,
                      by = c("gene_symbol", "tx_id"), all.x = TRUE,
                      sort = FALSE)
setorder(gene_summary, -n_recurrent_frame_risk_hotspots,
         -best_frame_risk_hotspot_score, gene_symbol)

summary_metrics <- data.table(
  metric = c(
    "n_frame_risk_genes",
    "n_frame_risk_runs",
    "n_frame_risk_hotspots",
    "n_recurrent_frame_risk_hotspots"
  ),
  value = as.character(c(
    uniqueN(frame_risk_runs$gene_symbol),
    uniqueN(frame_risk_runs$Run),
    nrow(hotspot_summary),
    hotspot_summary[n_frame_risk_runs_with_hotspot >= 2, .N]
  ))
)

fwrite(hotspots, file.path(output_dir,
                           "infected_off_frame_frame_risk_run_positions.csv"))
fwrite(profiles, file.path(output_dir,
                           "infected_off_frame_frame_risk_profiles.csv"))
fwrite(hotspot_summary, file.path(output_dir,
                                  "infected_off_frame_frame_risk_hotspots.csv"))
fwrite(gene_summary, file.path(output_dir,
                               "infected_off_frame_frame_risk_gene_summary.csv"))
fwrite(summary_metrics, file.path(output_dir,
                                  "infected_off_frame_frame_risk_summary_metrics.csv"))

top_hotspots <- head(hotspot_summary, 50L)
if (nrow(top_hotspots) > 0) {
  gene_order <- gene_summary$gene_symbol
  top_hotspots[, gene_label := factor(gene_symbol, levels = rev(gene_order))]
  p <- ggplot(
    top_hotspots,
    aes(x = cds_percent, y = gene_label,
        size = n_frame_risk_runs_with_hotspot,
        fill = fraction_frame_risk_runs_with_hotspot,
        text = paste0(
          gene_symbol, "<br>codon=", codon_index,
          "<br>tx pos=", transcript_position,
          "<br>frame=", cds_frame,
          "<br>runs=", n_frame_risk_runs_with_hotspot, " / ",
          group_n_frame_risk_runs,
          "<br>metadata frames=", metadata_dominant_frames,
          "<br>median metadata frame0=",
          round(median_metadata_canonical_frame_fraction, 3),
          "<br>runs=", run_examples
        ))
  ) +
    geom_point(shape = 21, color = "#303030", alpha = 0.92) +
    scale_fill_gradient(low = "white", high = "#8b0000",
                        na.value = "grey90") +
    scale_size_continuous(range = c(2.5, 9)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Frame-Risk Infected Off-Frame CDS Hotspots",
      subtitle = "Strong off-frame peaks removed from the strict queue because all-CDS/housekeeping/global frame QC is abnormal.",
      x = "Relative CDS position",
      y = NULL,
      fill = "Fraction of\nframe-risk runs",
      size = "Run hits"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
  ggsave(file.path(figure_dir, "infected_off_frame_frame_risk_hotspots.png"),
         p, width = 11, height = 7.5, dpi = 150)
  if (exists("dominant_save_ggplotly")) {
    dominant_save_ggplotly(
      p,
      file.path(figure_dir, "infected_off_frame_frame_risk_hotspots.html"),
      title = "infected_off_frame_frame_risk_hotspots"
    )
  }
}

message("Saved frame-risk hotspot outputs to: ", output_dir)
message("Top frame-risk genes:")
print(head(gene_summary[
  ,
  .(gene_symbol, n_frame_risk_runs, n_recurrent_frame_risk_hotspots,
    best_codon_index, best_cds_frame, best_transcript_position,
    best_fraction_frame_risk_runs_with_hotspot, metadata_dominant_frames,
    run_metadata_dominant_frames)
], 25L))
