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
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(GenomicRanges)
  library(IRanges)
})

helper_candidates <- c(
  file.path("scripts", "dominant_coverage_cache.R"),
  file.path("dominant_cell_states", "scripts", "dominant_coverage_cache.R")
)
helper_file <- helper_candidates[file.exists(helper_candidates)][1]
if (!is.na(helper_file)) source(helper_file)

htmlwidget_helper <- c(
  file.path("scripts", "dominant_htmlwidgets.R"),
  "dominant_htmlwidgets.R",
  file.path("dominant_cell_states", "scripts", "dominant_htmlwidgets.R")
)
htmlwidget_helper <- htmlwidget_helper[file.exists(htmlwidget_helper)][1]
if (!is.na(htmlwidget_helper)) source(htmlwidget_helper)

if (!exists("dominant_safe_name", mode = "function")) {
  dominant_safe_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
    x <- gsub("^_+|_+$", "", x)
    ifelse(nzchar(x), x, "unknown")
  }
}

analysis_dir <- if (file.exists("dominant_ouorf_cds_coupling")) {
  "."
} else if (file.exists(file.path("dominant_cell_states",
                                 "dominant_ouorf_cds_coupling"))) {
  "dominant_cell_states"
} else {
  "."
}

coupling_review_file <- file.path(
  analysis_dir,
  "dominant_ouorf_cds_coupling",
  "dominant_ouorf_cds_coupling_candidate_review.csv"
)
output_dir <- file.path(analysis_dir, "dominant_downstream_tis_rescue")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

all_candidates_output <- file.path(
  output_dir,
  "downstream_tis_rescue_candidates.csv"
)
best_candidates_output <- file.path(
  output_dir,
  "downstream_tis_rescue_best_candidates.csv"
)
feature_summary_output <- file.path(
  output_dir,
  "downstream_tis_rescue_feature_summary.csv"
)
summary_metrics_output <- file.path(
  output_dir,
  "downstream_tis_rescue_summary_metrics.csv"
)
summary_text_output <- file.path(
  output_dir,
  "downstream_tis_rescue_summary.txt"
)

message("Downstream-TIS rescue audit:")
message("  1. Read residual-positive ouORF/CDS rescue-review candidates")
message("  2. Reuse cached transcript coverage and sequence")
message("  3. Scan downstream clean-CDS segments for ATG/near-cognate starts")
message("  4. Add sequence geometry and local coverage-step evidence")

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}

safe_div <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}

bounded_seq <- function(from, to, upper) {
  from <- max(as.integer(from), 1L, na.rm = TRUE)
  to <- min(as.integer(to), as.integer(upper), na.rm = TRUE)
  if (!is.finite(from) || !is.finite(to) || from > to) return(integer())
  seq.int(from, to)
}

kozak_strength <- function(seq_chars, tx_start) {
  minus3 <- if (tx_start > 3L) seq_chars[[tx_start - 3L]] else NA_character_
  plus4 <- if (tx_start + 3L <= length(seq_chars)) {
    seq_chars[[tx_start + 3L]]
  } else {
    NA_character_
  }
  strong_minus <- !is.na(minus3) && minus3 %chin% c("A", "G")
  strong_plus <- !is.na(plus4) && plus4 == "G"
  strength <- if (strong_minus && strong_plus) {
    "strong"
  } else if (strong_minus || strong_plus) {
    "moderate"
  } else {
    "weak"
  }
  list(minus3 = minus3, plus4 = plus4, strength = strength)
}

first_in_frame_stop <- function(display_seq, tx_start, scan_limit) {
  stop_codons <- c("TAA", "TAG", "TGA")
  max_start <- min(nchar(display_seq) - 2L, scan_limit - 2L)
  if (!is.finite(max_start) || tx_start + 3L > max_start) {
    return(list(stop_tx_start = NA_integer_, stop_codon = NA_character_))
  }
  starts <- seq.int(tx_start + 3L, max_start, by = 3L)
  codons <- substring(display_seq, starts, starts + 2L)
  hit <- which(codons %chin% stop_codons)[1]
  if (is.na(hit)) {
    return(list(stop_tx_start = NA_integer_, stop_codon = NA_character_))
  }
  list(stop_tx_start = as.integer(starts[[hit]]), stop_codon = codons[[hit]])
}

precompute_stop_index <- function(display_seq, scan_limit) {
  stop_codons <- c("TAA", "TAG", "TGA")
  max_start <- min(nchar(display_seq) - 2L, scan_limit - 2L)
  if (!is.finite(max_start) || max_start < 1L) {
    return(vector("list", 3L))
  }
  starts <- seq_len(max_start)
  codons <- substring(display_seq, starts, starts + 2L)
  hits <- starts[codons %chin% stop_codons]
  hit_codons <- codons[codons %chin% stop_codons]
  out <- vector("list", 3L)
  for (frame in 0:2) {
    keep <- ((hits - 1L) %% 3L) == frame
    out[[frame + 1L]] <- data.table(
      stop_tx_start = as.integer(hits[keep]),
      stop_codon = hit_codons[keep]
    )
  }
  out
}

next_in_frame_stop <- function(stop_index, tx_start) {
  frame <- as.integer((tx_start - 1L) %% 3L)
  stops <- stop_index[[frame + 1L]]
  if (is.null(stops) || nrow(stops) == 0L) {
    return(list(stop_tx_start = NA_integer_, stop_codon = NA_character_))
  }
  idx <- findInterval(tx_start, stops$stop_tx_start) + 1L
  if (idx > nrow(stops)) {
    return(list(stop_tx_start = NA_integer_, stop_codon = NA_character_))
  }
  list(
    stop_tx_start = as.integer(stops$stop_tx_start[[idx]]),
    stop_codon = stops$stop_codon[[idx]]
  )
}

coverage_window_metrics <- function(mat, downstream_start, clean_end, tis_start) {
  n_pos <- nrow(mat)
  pre_pos <- bounded_seq(max(downstream_start, tis_start - 90L),
                         tis_start - 1L, n_pos)
  start_pos <- bounded_seq(tis_start, min(tis_start + 29L, clean_end), n_pos)
  post_pos <- bounded_seq(tis_start, min(tis_start + 179L, clean_end), n_pos)
  before_pos <- bounded_seq(downstream_start, tis_start - 1L, n_pos)
  downstream_pos <- bounded_seq(downstream_start, clean_end, n_pos)

  total_count <- function(pos) {
    if (!length(pos)) return(NA_real_)
    sum(mat[pos, , drop = FALSE], na.rm = TRUE)
  }
  density <- function(pos) {
    if (!length(pos)) return(NA_real_)
    total_count(pos) / length(pos)
  }
  col_sum <- function(pos) {
    if (!length(pos)) return(rep(NA_real_, ncol(mat)))
    colSums(mat[pos, , drop = FALSE], na.rm = TRUE)
  }

  pre_count <- total_count(pre_pos)
  start_count <- total_count(start_pos)
  post_count <- total_count(post_pos)
  before_count <- total_count(before_pos)
  downstream_count <- total_count(downstream_pos)
  pre_density <- density(pre_pos)
  start_density <- density(start_pos)
  post_density <- density(post_pos)
  before_density <- density(before_pos)
  downstream_density <- density(downstream_pos)

  pre_sums <- col_sum(pre_pos)
  post_sums <- col_sum(post_pos)
  pre_run_density <- pre_sums / max(length(pre_pos), 1L)
  post_run_density <- post_sums / max(length(post_pos), 1L)
  run_ok <- is.finite(pre_sums) & is.finite(post_sums) &
    (pre_sums + post_sums) >= 10
  run_log2 <- rep(NA_real_, length(pre_sums))
  if (any(run_ok)) {
    run_log2[run_ok] <- log2((post_run_density[run_ok] + 0.01) /
                               (pre_run_density[run_ok] + 0.01))
  }

  data.table(
    pre_window_bases = length(pre_pos),
    start_window_bases = length(start_pos),
    post_window_bases = length(post_pos),
    before_tis_segment_bases = length(before_pos),
    downstream_segment_bases = length(downstream_pos),
    total_pre_window_counts = pre_count,
    total_start_window_counts = start_count,
    total_post_window_counts = post_count,
    total_before_tis_segment_counts = before_count,
    total_downstream_segment_counts = downstream_count,
    total_pre_window_density = pre_density,
    total_start_window_density = start_density,
    total_post_window_density = post_density,
    total_before_tis_segment_density = before_density,
    total_downstream_segment_density = downstream_density,
    total_post_vs_pre_density_log2 = if (is.finite(pre_density) &&
                                        is.finite(post_density)) {
      log2((post_density + 1e-6) / (pre_density + 1e-6))
    } else {
      NA_real_
    },
    total_start_vs_pre_density_log2 = if (is.finite(pre_density) &&
                                         is.finite(start_density)) {
      log2((start_density + 1e-6) / (pre_density + 1e-6))
    } else {
      NA_real_
    },
    n_runs_shape_evaluable = sum(run_ok),
    fraction_runs_post_gt_pre_1_25 = if (any(run_ok)) {
      mean(post_run_density[run_ok] > (pre_run_density[run_ok] * 1.25))
    } else {
      NA_real_
    },
    median_run_post_vs_pre_log2 = if (any(run_ok)) {
      median(run_log2[run_ok], na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

coverage_window_metrics_from_totals <- function(row_totals, downstream_start,
                                                clean_end, tis_start) {
  n_pos <- length(row_totals)
  pre_pos <- bounded_seq(max(downstream_start, tis_start - 90L),
                         tis_start - 1L, n_pos)
  start_pos <- bounded_seq(tis_start, min(tis_start + 29L, clean_end), n_pos)
  post_pos <- bounded_seq(tis_start, min(tis_start + 179L, clean_end), n_pos)
  before_pos <- bounded_seq(downstream_start, tis_start - 1L, n_pos)
  downstream_pos <- bounded_seq(downstream_start, clean_end, n_pos)

  total_count <- function(pos) {
    if (!length(pos)) return(NA_real_)
    sum(row_totals[pos], na.rm = TRUE)
  }
  density <- function(pos) {
    if (!length(pos)) return(NA_real_)
    total_count(pos) / length(pos)
  }

  pre_count <- total_count(pre_pos)
  start_count <- total_count(start_pos)
  post_count <- total_count(post_pos)
  before_count <- total_count(before_pos)
  downstream_count <- total_count(downstream_pos)
  pre_density <- density(pre_pos)
  start_density <- density(start_pos)
  post_density <- density(post_pos)
  before_density <- density(before_pos)
  downstream_density <- density(downstream_pos)

  data.table(
    pre_window_bases = length(pre_pos),
    start_window_bases = length(start_pos),
    post_window_bases = length(post_pos),
    before_tis_segment_bases = length(before_pos),
    downstream_segment_bases = length(downstream_pos),
    total_pre_window_counts = pre_count,
    total_start_window_counts = start_count,
    total_post_window_counts = post_count,
    total_before_tis_segment_counts = before_count,
    total_downstream_segment_counts = downstream_count,
    total_pre_window_density = pre_density,
    total_start_window_density = start_density,
    total_post_window_density = post_density,
    total_before_tis_segment_density = before_density,
    total_downstream_segment_density = downstream_density,
    total_post_vs_pre_density_log2 = if (is.finite(pre_density) &&
                                        is.finite(post_density)) {
      log2((post_density + 1e-6) / (pre_density + 1e-6))
    } else {
      NA_real_
    },
    total_start_vs_pre_density_log2 = if (is.finite(pre_density) &&
                                         is.finite(start_density)) {
      log2((start_density + 1e-6) / (pre_density + 1e-6))
    } else {
      NA_real_
    },
    n_runs_shape_evaluable = NA_integer_,
    fraction_runs_post_gt_pre_1_25 = NA_real_,
    median_run_post_vs_pre_log2 = NA_real_
  )
}

empty_coverage_metrics <- function() {
  data.table(
    pre_window_bases = NA_integer_,
    start_window_bases = NA_integer_,
    post_window_bases = NA_integer_,
    before_tis_segment_bases = NA_integer_,
    downstream_segment_bases = NA_integer_,
    total_pre_window_counts = NA_real_,
    total_start_window_counts = NA_real_,
    total_post_window_counts = NA_real_,
    total_before_tis_segment_counts = NA_real_,
    total_downstream_segment_counts = NA_real_,
    total_pre_window_density = NA_real_,
    total_start_window_density = NA_real_,
    total_post_window_density = NA_real_,
    total_before_tis_segment_density = NA_real_,
    total_downstream_segment_density = NA_real_,
    total_post_vs_pre_density_log2 = NA_real_,
    total_start_vs_pre_density_log2 = NA_real_,
    n_runs_shape_evaluable = NA_integer_,
    fraction_runs_post_gt_pre_1_25 = NA_real_,
    median_run_post_vs_pre_log2 = NA_real_
  )
}

cache_file_for_gene_tx <- function(gene_symbol, tx_id) {
  cache_dir <- file.path(analysis_dir, ".coverage_cache")
  if (!dir.exists(cache_dir)) return(NA_character_)
  pattern <- paste0(
    "^",
    dominant_safe_name(gene_symbol),
    "__",
    dominant_safe_name(tx_id),
    "__.*\\.rds$"
  )
  files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  files[order(info$mtime, decreasing = TRUE, na.last = TRUE)][1]
}

display_from_cache <- function(cache_obj) {
  sig <- cache_obj$signature
  gr <- GenomicRanges::GRanges(
    seqnames = sig$seqnames,
    ranges = IRanges::IRanges(sig$start, sig$end),
    strand = sig$strand
  )
  GenomicRanges::GRangesList(
    stats::setNames(list(gr), sig$transcript_names[[1]])
  )
}

sequence_from_cache <- function(cache_obj, df) {
  tryCatch(
    {
      display <- display_from_cache(cache_obj)
      as.character(ORFik::txSeqsFromFa(display, df, TRUE)[[1]])
    },
    error = function(e) NA_character_
  )
}

scan_downstream_starts <- function(row, cache_obj, df) {
  display_seq <- sequence_from_cache(cache_obj, df)
  if (is.na(display_seq) || !nzchar(display_seq)) {
    return(data.table())
  }

  coverage <- cache_obj$coverage
  mat <- as.matrix(coverage)
  storage.mode(mat) <- "numeric"
  row_totals <- rowSums(mat, na.rm = TRUE)
  seq_chars <- strsplit(toupper(display_seq), "", fixed = TRUE)[[1]]
  display_seq <- paste(seq_chars, collapse = "")
  n_pos <- min(length(seq_chars), nrow(mat))

  canonical_cds_tx_start <- if ("ouorf_start_relative_to_cds_start" %in%
                                names(row) &&
                                is.finite(row$ouorf_start_relative_to_cds_start)) {
    as.integer(row$tx_start - row$ouorf_start_relative_to_cds_start)
  } else {
    as.integer(row$clean_cds_tx_start)
  }
  downstream_start <- max(as.integer(row$tx_end + 1L),
                          as.integer(row$clean_cds_tx_start))
  clean_end <- min(as.integer(row$clean_cds_tx_end), n_pos)
  scan_from <- downstream_start
  scan_to <- clean_end - 2L
  if (!is.finite(scan_from) || !is.finite(scan_to) || scan_from > scan_to) {
    return(data.table())
  }

  start_codons <- c("ATG", "CTG", "GTG", "TTG", "ACG", "ATT", "ATC", "ATA")
  near_cognate <- setdiff(start_codons, "ATG")
  starts <- seq.int(scan_from, scan_to)
  codons <- substring(display_seq, starts, starts + 2L)
  keep <- codons %chin% start_codons
  if (!any(keep)) return(data.table())

  starts <- starts[keep]
  codons <- codons[keep]
  stop_index <- precompute_stop_index(display_seq, clean_end + 3L)
  cache_path <- cache_file_for_gene_tx(row$gene_symbol, row$tx_id)
  out <- rbindlist(lapply(seq_along(starts), function(i) {
    tis <- as.integer(starts[[i]])
    codon <- codons[[i]]
    kz <- kozak_strength(seq_chars, tis)
    stop <- next_in_frame_stop(stop_index, tis)
    frame_relative_to_canonical_cds <-
      as.integer((tis - canonical_cds_tx_start) %% 3L)
    frame_relative_to_ouorf_start <- as.integer((tis - row$tx_start) %% 3L)
    predicted_stop <- if (is.finite(stop$stop_tx_start)) {
      stop$stop_tx_start + 2L
    } else {
      NA_integer_
    }
    predicted_bases <- if (is.finite(predicted_stop)) {
      predicted_stop - tis + 1L
    } else {
      clean_end - tis + 1L
    }
    data.table(
      gene_symbol = row$gene_symbol,
      tx_id = row$tx_id,
      feature_id = row$feature_id,
      ouorf_tx_start = as.integer(row$tx_start),
      ouorf_tx_end = as.integer(row$tx_end),
      clean_cds_tx_start = as.integer(row$clean_cds_tx_start),
      clean_cds_tx_end = as.integer(row$clean_cds_tx_end),
      canonical_cds_tx_start = canonical_cds_tx_start,
      downstream_scan_start = downstream_start,
      candidate_tis_tx_start = tis,
      candidate_start_codon = codon,
      candidate_start_codon_class = if (codon == "ATG") {
        "ATG"
      } else if (codon %chin% near_cognate) {
        "near_cognate"
      } else {
        "other"
      },
      kozak_minus3 = kz$minus3,
      kozak_plus4 = kz$plus4,
      candidate_kozak_strength = kz$strength,
      distance_from_ouorf_stop = as.integer(tis - row$tx_end - 1L),
      frame_relative_to_canonical_cds =
        frame_relative_to_canonical_cds,
      frame_preserves_canonical_cds =
        frame_relative_to_canonical_cds == 0L,
      frame_relative_to_ouorf_start = frame_relative_to_ouorf_start,
      first_inframe_stop_tx_start = as.integer(stop$stop_tx_start),
      first_inframe_stop_codon = stop$stop_codon,
      predicted_truncated_orf_nt = as.integer(predicted_bases),
      predicted_truncated_orf_aa = as.integer(floor(predicted_bases / 3L)),
      predicted_orf_reaches_clean_cds_end =
        is.na(stop$stop_tx_start) || stop$stop_tx_start >= clean_end - 2L,
      cache_file = cache_path,
      sequence_length = length(seq_chars),
      coverage_rows = nrow(mat),
      original_ouorf_start_codon = row$start_codon,
      original_ouorf_start_codon_class = row$start_codon_class,
      original_ouorf_frame_relative_to_cds =
        row$frame_relative_to_cds,
      original_ouorf_kozak_strength = row$kozak_strength,
      raw_log2_fpkm_spearman_rho =
        row$raw_log2_fpkm_spearman_rho,
      total_route_residual_spearman_rho =
        row$total_route_residual_spearman_rho,
      relative_route_use_spearman_rho =
        row$relative_route_use_spearman_rho,
      review_priority_score = row$review_priority_score
    )
  }), fill = TRUE)

  out[, `:=`(
    codon_score = fifelse(candidate_start_codon == "ATG", 3,
                          fifelse(candidate_start_codon_class == "near_cognate",
                                  1, 0)),
    kozak_score = fcase(
      candidate_kozak_strength == "strong", 2,
      candidate_kozak_strength == "moderate", 1,
      default = 0
    ),
    frame_score = fifelse(frame_preserves_canonical_cds, 3, 0),
    length_score = fifelse(predicted_truncated_orf_aa >= 100, 2,
                           fifelse(predicted_truncated_orf_aa >= 50, 1, 0)),
    distance_score = fcase(
      distance_from_ouorf_stop >= 10 & distance_from_ouorf_stop <= 300, 2,
      distance_from_ouorf_stop > 300 & distance_from_ouorf_stop <= 900, 1,
      distance_from_ouorf_stop >= 0 & distance_from_ouorf_stop < 10, 0,
      default = 0
    )
  )]
  out[, sequence_geometry_supported :=
        codon_score > 0 &
        frame_preserves_canonical_cds &
        candidate_kozak_strength %chin% c("strong", "moderate") &
        predicted_truncated_orf_aa >= 50 &
        predicted_orf_reaches_clean_cds_end]

  empty_cov <- empty_coverage_metrics()
  for (col in names(empty_cov)) {
    out[, (col) := empty_cov[[col]][1]]
  }
  out[, coverage_evaluated := FALSE]
  out[, sequence_precoverage_score :=
        codon_score + kozak_score + frame_score + length_score +
        distance_score + fifelse(sequence_geometry_supported, 2, 0)]
  max_coverage_candidates <- 160L
  eval_idx <- unique(c(
    out[sequence_geometry_supported == TRUE, .I],
    out[order(-sequence_precoverage_score,
              -as.integer(frame_preserves_canonical_cds),
              candidate_tis_tx_start),
        .I[seq_len(min(.N, max_coverage_candidates))]]
  ))
  if (length(eval_idx)) {
    coverage_tables <- rbindlist(lapply(eval_idx, function(idx) {
      cbind(
        data.table(row_index = idx),
        coverage_window_metrics_from_totals(
          row_totals = row_totals,
          downstream_start = downstream_start,
          clean_end = clean_end,
          tis_start = out$candidate_tis_tx_start[[idx]]
        )
      )
    }), fill = TRUE)
    for (col in setdiff(names(coverage_tables), "row_index")) {
      set(out, i = coverage_tables$row_index, j = col,
          value = coverage_tables[[col]])
    }
    set(out, i = eval_idx, j = "coverage_evaluated", value = TRUE)
  }

  out[, coverage_shape_score := fcase(
    is.finite(total_post_vs_pre_density_log2) &
      total_post_vs_pre_density_log2 >= log2(2), 2,
    is.finite(total_post_vs_pre_density_log2) &
      total_post_vs_pre_density_log2 >= log2(1.25), 1,
    default = 0
  )]
  out[, coverage_shape_supported :=
        coverage_evaluated == TRUE &
        pre_window_bases >= 15L &
        post_window_bases >= 30L &
        is.finite(total_post_vs_pre_density_log2) &
        total_post_vs_pre_density_log2 >= log2(1.25) &
        is.finite(total_pre_window_counts) &
        is.finite(total_post_window_counts) &
        (total_pre_window_counts + total_post_window_counts) >= 100]
  out[, rescue_review_class := fcase(
    sequence_geometry_supported & coverage_shape_supported,
    "sequence_and_shape_downstream_tis_review",
    sequence_geometry_supported,
    "sequence_only_downstream_tis_review",
    coverage_shape_supported,
    "shape_only_without_good_start_review",
    default = "weak_or_no_downstream_tis_support"
  )]
  out[, downstream_tis_evidence_score :=
        codon_score + kozak_score + frame_score + length_score +
        distance_score + coverage_shape_score +
        fifelse(sequence_geometry_supported, 2, 0) +
        fifelse(coverage_shape_supported, 2, 0)]
  out[]
}

write_plot_pair <- function(plot, name, width = 9, height = 6) {
  png_file <- file.path(figure_dir, paste0(name, ".png"))
  pdf_file <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 220,
         bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
}

write_plotly <- function(plot, name, title = name) {
  if (exists("dominant_save_ggplotly", mode = "function")) {
    dominant_save_ggplotly(
      plot,
      file.path(figure_dir, paste0(name, ".html")),
      title = title
    )
  }
}

wrap_text <- function(x, width = 110L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

stop_if_missing(coupling_review_file)
review <- fread(coupling_review_file)
review <- review[truncated_cds_rescue_review == TRUE]

if (nrow(review) == 0L) {
  empty <- data.table()
  fwrite(empty, all_candidates_output)
  fwrite(empty, best_candidates_output)
  fwrite(empty, feature_summary_output)
  fwrite(data.table(metric = "truncated_cds_rescue_review_rows", value = 0),
         summary_metrics_output)
  writeLines("No truncated-CDS rescue-review rows were available.",
             summary_text_output)
  quit(save = "no", status = 0)
}

df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)

candidate_tables <- vector("list", nrow(review))
missing_cache <- data.table()
for (i in seq_len(nrow(review))) {
  row <- review[i]
  message("  scanning ", row$gene_symbol, " ", row$feature_id,
          " (", row$tx_id, ")")
  cache_file <- cache_file_for_gene_tx(row$gene_symbol, row$tx_id)
  if (is.na(cache_file) || !file.exists(cache_file)) {
    missing_cache <- rbind(
      missing_cache,
      data.table(
        gene_symbol = row$gene_symbol,
        tx_id = row$tx_id,
        feature_id = row$feature_id,
        missing_reason = "coverage_cache_file_missing"
      )
    )
    next
  }
  cache_obj <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  if (is.null(cache_obj) || is.null(cache_obj$coverage)) {
    missing_cache <- rbind(
      missing_cache,
      data.table(
        gene_symbol = row$gene_symbol,
        tx_id = row$tx_id,
        feature_id = row$feature_id,
        missing_reason = "coverage_cache_unreadable"
      )
    )
    next
  }
  candidate_tables[[i]] <- scan_downstream_starts(row, cache_obj, df)
  message("    downstream start candidates: ", nrow(candidate_tables[[i]]))
}

all_candidates <- rbindlist(candidate_tables, fill = TRUE)
if (nrow(all_candidates) > 0L) {
  setorder(
    all_candidates,
    gene_symbol,
    feature_id,
    -downstream_tis_evidence_score,
    -sequence_geometry_supported,
    -coverage_shape_supported,
    candidate_tis_tx_start
  )
  all_candidates[, candidate_rank_for_feature := seq_len(.N),
                 by = .(gene_symbol, tx_id, feature_id)]
}

if (nrow(all_candidates) > 0L) {
  best_candidates <- all_candidates[candidate_rank_for_feature == 1L]
} else {
  best_candidates <- data.table()
}

feature_summary <- review[, .(
  gene_symbol,
  tx_id,
  feature_id,
  ouorf_tx_start = tx_start,
  ouorf_tx_end = tx_end,
  clean_cds_tx_start,
  clean_cds_tx_end,
  original_ouorf_start_codon = start_codon,
  original_ouorf_start_codon_class = start_codon_class,
  original_ouorf_frame_relative_to_cds = frame_relative_to_cds,
  original_ouorf_kozak_strength = kozak_strength,
  raw_log2_fpkm_spearman_rho,
  total_route_residual_spearman_rho,
  relative_route_use_spearman_rho,
  review_priority_score
)]
if (nrow(all_candidates) > 0L) {
  candidate_counts <- all_candidates[, .(
    n_downstream_start_codons = .N,
    n_frame_preserving_start_codons = sum(frame_preserves_canonical_cds,
                                          na.rm = TRUE),
    n_sequence_geometry_supported = sum(sequence_geometry_supported,
                                        na.rm = TRUE),
    n_coverage_shape_supported = sum(coverage_shape_supported,
                                     na.rm = TRUE),
    n_sequence_and_shape_supported =
      sum(sequence_geometry_supported & coverage_shape_supported,
          na.rm = TRUE)
  ), by = .(gene_symbol, tx_id, feature_id)]
  best_small <- best_candidates[, .(
    gene_symbol,
    tx_id,
    feature_id,
    best_candidate_tis_tx_start = candidate_tis_tx_start,
    best_candidate_start_codon = candidate_start_codon,
    best_candidate_start_codon_class = candidate_start_codon_class,
    best_candidate_kozak_strength = candidate_kozak_strength,
    best_distance_from_ouorf_stop = distance_from_ouorf_stop,
    best_frame_relative_to_canonical_cds =
      frame_relative_to_canonical_cds,
    best_predicted_truncated_orf_aa = predicted_truncated_orf_aa,
    best_total_post_vs_pre_density_log2 =
      total_post_vs_pre_density_log2,
    best_fraction_runs_post_gt_pre_1_25 =
      fraction_runs_post_gt_pre_1_25,
    best_sequence_geometry_supported = sequence_geometry_supported,
    best_coverage_shape_supported = coverage_shape_supported,
    best_rescue_review_class = rescue_review_class,
    best_downstream_tis_evidence_score = downstream_tis_evidence_score
  )]
  feature_summary <- merge(
    feature_summary,
    candidate_counts,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
  feature_summary <- merge(
    feature_summary,
    best_small,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
if (nrow(missing_cache) > 0L) {
  feature_summary <- merge(
    feature_summary,
    missing_cache,
    by = c("gene_symbol", "tx_id", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
}
for (col in c(
  "n_downstream_start_codons",
  "n_frame_preserving_start_codons",
  "n_sequence_geometry_supported",
  "n_coverage_shape_supported",
  "n_sequence_and_shape_supported"
)) {
  if (!col %in% names(feature_summary)) feature_summary[, (col) := 0L]
  feature_summary[is.na(get(col)), (col) := 0L]
}
if (!"missing_reason" %in% names(feature_summary)) {
  feature_summary[, missing_reason := NA_character_]
}

summary_metrics <- data.table(
  metric = c(
    "truncated_cds_rescue_review_rows",
    "rows_with_cache_missing_or_unreadable",
    "downstream_tis_candidate_rows",
    "features_with_downstream_start_codons",
    "features_with_frame_preserving_start_codons",
    "features_with_sequence_geometry_supported_start",
    "features_with_coverage_shape_supported_start",
    "features_with_sequence_and_shape_supported_start"
  ),
  value = c(
    nrow(review),
    nrow(missing_cache),
    nrow(all_candidates),
    if (nrow(feature_summary)) {
      sum(feature_summary$n_downstream_start_codons > 0L, na.rm = TRUE)
    } else 0,
    if (nrow(feature_summary)) {
      sum(feature_summary$n_frame_preserving_start_codons > 0L, na.rm = TRUE)
    } else 0,
    if (nrow(feature_summary)) {
      sum(feature_summary$n_sequence_geometry_supported > 0L, na.rm = TRUE)
    } else 0,
    if (nrow(feature_summary)) {
      sum(feature_summary$n_coverage_shape_supported > 0L, na.rm = TRUE)
    } else 0,
    if (nrow(feature_summary)) {
      sum(feature_summary$n_sequence_and_shape_supported > 0L, na.rm = TRUE)
    } else 0
  )
)

fwrite(all_candidates, all_candidates_output)
fwrite(best_candidates, best_candidates_output)
fwrite(feature_summary, feature_summary_output)
fwrite(summary_metrics, summary_metrics_output)

if (nrow(best_candidates) > 0L) {
  best_plot <- copy(best_candidates)
  best_plot[, label := paste(gene_symbol, feature_id, sep = " ")]
  best_plot[, plot_y := total_post_vs_pre_density_log2]
  best_plot[!is.finite(plot_y), plot_y := 0]
  class_levels <- c(
    "Seq + shape",
    "Seq only",
    "Shape only",
    "Weak / none"
  )
  best_plot[, review_label := fcase(
    rescue_review_class == "sequence_and_shape_downstream_tis_review",
    "Seq + shape",
    rescue_review_class == "sequence_only_downstream_tis_review",
    "Seq only",
    rescue_review_class == "shape_only_without_good_start_review",
    "Shape only",
    default = "Weak / none"
  )]
  best_plot[, review_label := factor(review_label, levels = class_levels)]
  p_best <- ggplot(
    best_plot,
    aes(
      x = distance_from_ouorf_stop,
      y = plot_y,
      color = review_label,
      shape = candidate_start_codon_class,
      size = pmax(downstream_tis_evidence_score, 1),
      text = paste0(
        gene_symbol, " ", feature_id,
        "<br>TIS: ", candidate_start_codon, " at tx ", candidate_tis_tx_start,
        "<br>Kozak: ", candidate_kozak_strength,
        "<br>Frame to CDS: ", frame_relative_to_canonical_cds,
        "<br>Predicted aa: ", predicted_truncated_orf_aa,
        "<br>log2 post/pre: ",
        round(total_post_vs_pre_density_log2, 3),
        "<br>Runs post>pre*1.25: ",
        round(fraction_runs_post_gt_pre_1_25, 3)
      )
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey65") +
    geom_hline(yintercept = log2(1.25), linewidth = 0.35,
               color = "grey45", linetype = "dashed") +
    geom_point(alpha = 0.88) +
    geom_text(aes(label = label), size = 3.1, vjust = -0.8,
              check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(
      values = c(
        "Seq + shape" = "#8B0000",
        "Seq only" = "#D95F02",
        "Shape only" = "#1F78B4",
        "Weak / none" = "#6B7280"
      ),
      drop = TRUE
    ) +
    scale_size_continuous(range = c(2.3, 6), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Downstream-start audit for residual-positive ouORF/CDS rows",
      subtitle = wrap_text(
        "Each point is the best downstream ATG/near-cognate candidate after one overlapping uORF. Dashed line marks 1.25x post/pre local coverage density."
      ),
      x = "Distance from ouORF stop to candidate TIS (nt)",
      y = "Local coverage step after TIS, log2(post / pre)",
      color = "Review class",
      shape = "Start class"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(12, 24, 12, 20)
    )
  write_plot_pair(p_best, "downstream_tis_rescue_best_candidates",
                  width = 10.8, height = 6.8)
  write_plotly(p_best, "downstream_tis_rescue_best_candidates",
               "Downstream-TIS rescue best candidates")

  rects <- rbindlist(list(
    best_candidates[, .(
      gene_symbol,
      feature_id,
      region = "overlapping uORF",
      xmin = ouorf_tx_start,
      xmax = ouorf_tx_end,
      ymin = 0.54,
      ymax = 0.74
    )],
    best_candidates[, .(
      gene_symbol,
      feature_id,
      region = "downstream clean CDS",
      xmin = clean_cds_tx_start,
      xmax = clean_cds_tx_end,
      ymin = 0.26,
      ymax = 0.46
    )]
  ), fill = TRUE)
  rects[, region := factor(region, levels = c("overlapping uORF",
                                              "downstream clean CDS"))]
  p_map <- ggplot() +
    geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
      color = NA,
      alpha = 0.9
    ) +
    geom_segment(
      data = best_candidates,
      aes(x = candidate_tis_tx_start, xend = candidate_tis_tx_start,
          y = 0.20, yend = 0.82),
      linewidth = 0.35,
      color = "#111827"
    ) +
    geom_point(
      data = best_candidates,
      aes(x = candidate_tis_tx_start, y = 0.84,
          color = factor(fcase(
            rescue_review_class == "sequence_and_shape_downstream_tis_review",
            "Seq + shape",
            rescue_review_class == "sequence_only_downstream_tis_review",
            "Seq only",
            rescue_review_class == "shape_only_without_good_start_review",
            "Shape only",
            default = "Weak / none"
          ), levels = class_levels),
          shape = candidate_start_codon_class),
      size = 3.2
    ) +
    geom_text(
      data = best_candidates,
      aes(x = candidate_tis_tx_start, y = 0.92,
          label = paste0(candidate_start_codon, "\n",
                         candidate_kozak_strength)),
      size = 2.65,
      lineheight = 0.9,
      check_overlap = TRUE
    ) +
    facet_wrap(~ gene_symbol, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("overlapping uORF" = "#7C3AED",
                                 "downstream clean CDS" = "#0F766E")) +
    scale_color_manual(
      values = c(
        "Seq + shape" = "#8B0000",
        "Seq only" = "#D95F02",
        "Shape only" = "#1F78B4",
        "Weak / none" = "#6B7280"
      ),
      drop = TRUE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.13))) +
    coord_cartesian(ylim = c(0.15, 1.02), expand = TRUE, clip = "off") +
    labs(
      title = "Transcript geometry for best downstream-start candidates",
      subtitle = wrap_text(
        "Region boxes have no outlines. Vertical lines mark the best candidate downstream start in the clean-CDS segment."
      ),
      x = "Transcript coordinate",
      y = NULL,
      fill = "Region",
      color = "Review class",
      shape = "Start class"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(12, 24, 12, 20)
    )
  write_plot_pair(p_map, "downstream_tis_rescue_gene_maps",
                  width = 11.5, height = 9)
  write_plotly(p_map, "downstream_tis_rescue_gene_maps",
               "Downstream-TIS rescue gene maps")
}

top_text <- if (nrow(best_candidates) > 0L) {
  best_candidates[
    order(-downstream_tis_evidence_score,
          -sequence_geometry_supported,
          -coverage_shape_supported)
  ][
    seq_len(min(.N, 10L)),
    paste0(
      gene_symbol, " ", feature_id, ": ",
      candidate_start_codon, " at tx ", candidate_tis_tx_start,
      ", frame=", frame_relative_to_canonical_cds,
      ", Kozak=", candidate_kozak_strength,
      ", class=", rescue_review_class,
      ", score=", downstream_tis_evidence_score
    )
  ]
} else {
  "No downstream ATG/near-cognate candidates were detected."
}
writeLines(
  c(
    "Downstream-TIS rescue audit summary",
    "",
    paste0("Input truncated-CDS rescue-review rows: ", nrow(review)),
    paste0("Downstream ATG/near-cognate candidate rows: ",
           nrow(all_candidates)),
    paste0("Features with sequence-geometry supported starts: ",
           summary_metrics[metric ==
                             "features_with_sequence_geometry_supported_start",
                           value]),
    paste0("Features with coverage-shape supported starts: ",
           summary_metrics[metric ==
                             "features_with_coverage_shape_supported_start",
                           value]),
    paste0("Features with both sequence and shape support: ",
           summary_metrics[metric ==
                             "features_with_sequence_and_shape_supported_start",
                           value]),
    "",
    "Top best-candidate rows:",
    top_text,
    "",
    "Interpretation: this is a review queue, not proof of downstream initiation. Total positional coverage can still be shared by canonical CDS elongation, transcript abundance, or study effects. A true rescue call still needs p-site/frame/start-peak evidence and browser review."
  ),
  summary_text_output
)

message("Saved downstream-TIS rescue audit outputs in ", output_dir)
