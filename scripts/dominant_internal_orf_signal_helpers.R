dominant_internal_signal_position_sums <- function(mat, positions) {
  positions <- sort(unique(as.integer(positions)))
  positions <- positions[
    is.finite(positions) & positions >= 1L & positions <= nrow(mat)
  ]
  if (length(positions) == 0L) {
    return(list(
      positions = integer(),
      bases = 0L,
      counts = rep(0, ncol(mat))
    ))
  }
  list(
    positions = positions,
    bases = length(positions),
    counts = as.numeric(colSums(mat[positions, , drop = FALSE], na.rm = TRUE))
  )
}

dominant_internal_signal_density_ratio <- function(numerator_counts,
                                                   numerator_bases,
                                                   denominator_counts,
                                                   denominator_bases) {
  out <- rep(NA_real_, length(numerator_counts))
  numerator_bases <- rep_len(as.numeric(numerator_bases), length(out))
  denominator_bases <- rep_len(as.numeric(denominator_bases), length(out))
  numerator_density <- rep(NA_real_, length(out))
  denominator_density <- rep(NA_real_, length(out))
  numerator_ok <- is.finite(numerator_bases) & numerator_bases > 0L
  denominator_ok <- is.finite(denominator_bases) & denominator_bases > 0L
  numerator_density[numerator_ok] <- numerator_counts[numerator_ok] /
    numerator_bases[numerator_ok]
  denominator_density[denominator_ok] <- denominator_counts[denominator_ok] /
    denominator_bases[denominator_ok]
  ok <- numerator_ok & denominator_ok &
    is.finite(numerator_density) & is.finite(denominator_density) &
    denominator_density > 0
  out[ok] <- numerator_density[ok] / denominator_density[ok]
  out
}

dominant_internal_orf_empty_signal_metrics <- function(runs) {
  n <- length(runs)
  data.table::data.table(
    Run = as.character(runs),
    internal_frame_relative_to_cds = rep(NA_integer_, n),
    internal_phase_total_raw_counts = rep(NA_real_, n),
    internal_feature_phase0_raw_counts = rep(NA_real_, n),
    internal_feature_phase1_raw_counts = rep(NA_real_, n),
    internal_feature_phase2_raw_counts = rep(NA_real_, n),
    internal_canonical_phase_raw_counts = rep(NA_real_, n),
    internal_iorf_frame_fraction = rep(NA_real_, n),
    internal_canonical_frame_fraction = rep(NA_real_, n),
    internal_iorf_vs_canonical_frame_delta = rep(NA_real_, n),
    internal_start_window_bases = rep(NA_integer_, n),
    internal_feature_core_bases = rep(NA_integer_, n),
    internal_upstream_window_bases = rep(NA_integer_, n),
    internal_start_window_raw_counts = rep(NA_real_, n),
    internal_feature_core_raw_counts = rep(NA_real_, n),
    internal_upstream_window_raw_counts = rep(NA_real_, n),
    internal_start_vs_body_density_ratio = rep(NA_real_, n),
    internal_start_vs_upstream_density_ratio = rep(NA_real_, n)
  )
}

dominant_internal_orf_signal_metrics <- function(mat, feature_positions,
                                                 cds_positions, feature_type,
                                                 category,
                                                 start_window_bases = 12L,
                                                 upstream_window_bases = 30L) {
  runs <- colnames(mat)
  if (is.null(runs)) runs <- as.character(seq_len(ncol(mat)))
  out <- dominant_internal_orf_empty_signal_metrics(runs)
  internal_feature <- identical(as.character(feature_type)[1], "internal_orf") ||
    identical(as.character(category)[1], "internal")
  if (!internal_feature || nrow(mat) == 0L || ncol(mat) == 0L) return(out)

  feature_signal <- dominant_internal_signal_position_sums(mat, feature_positions)
  cds_signal <- dominant_internal_signal_position_sums(mat, cds_positions)
  positions <- feature_signal$positions
  cds_positions <- cds_signal$positions
  if (length(positions) == 0L || length(cds_positions) == 0L) return(out)

  feature_phase <- (seq_along(positions) - 1L) %% 3L
  phase_signals <- lapply(0:2, function(phase) {
    dominant_internal_signal_position_sums(
      mat,
      positions[feature_phase == phase]
    )
  })
  phase_counts <- do.call(cbind, lapply(phase_signals, `[[`, "counts"))
  phase_total <- rowSums(phase_counts, na.rm = TRUE)
  frame_relative_to_cds <- as.integer((min(positions) - min(cds_positions)) %% 3L)
  canonical_phase <- as.integer((3L - frame_relative_to_cds) %% 3L)
  canonical_counts <- phase_counts[, canonical_phase + 1L]

  start_positions <- head(
    positions,
    min(as.integer(start_window_bases), length(positions))
  )
  core_positions <- setdiff(positions, start_positions)
  upstream_candidates <- cds_positions[cds_positions < min(positions)]
  upstream_positions <- tail(
    upstream_candidates,
    min(as.integer(upstream_window_bases), length(upstream_candidates))
  )
  start_signal <- dominant_internal_signal_position_sums(mat, start_positions)
  core_signal <- dominant_internal_signal_position_sums(mat, core_positions)
  upstream_signal <- dominant_internal_signal_position_sums(mat, upstream_positions)

  out[, `:=`(
    internal_frame_relative_to_cds = frame_relative_to_cds,
    internal_phase_total_raw_counts = phase_total,
    internal_feature_phase0_raw_counts = phase_counts[, 1L],
    internal_feature_phase1_raw_counts = phase_counts[, 2L],
    internal_feature_phase2_raw_counts = phase_counts[, 3L],
    internal_canonical_phase_raw_counts = canonical_counts,
    internal_iorf_frame_fraction = data.table::fifelse(
      phase_total > 0,
      phase_counts[, 1L] / phase_total,
      NA_real_
    ),
    internal_canonical_frame_fraction = data.table::fifelse(
      phase_total > 0,
      canonical_counts / phase_total,
      NA_real_
    ),
    internal_start_window_bases = start_signal$bases,
    internal_feature_core_bases = core_signal$bases,
    internal_upstream_window_bases = upstream_signal$bases,
    internal_start_window_raw_counts = start_signal$counts,
    internal_feature_core_raw_counts = core_signal$counts,
    internal_upstream_window_raw_counts = upstream_signal$counts
  )]
  out[, internal_iorf_vs_canonical_frame_delta :=
        internal_iorf_frame_fraction - internal_canonical_frame_fraction]
  out[, internal_start_vs_body_density_ratio :=
        dominant_internal_signal_density_ratio(
          internal_start_window_raw_counts,
          internal_start_window_bases,
          internal_feature_core_raw_counts,
          internal_feature_core_bases
        )]
  out[, internal_start_vs_upstream_density_ratio :=
        dominant_internal_signal_density_ratio(
          internal_start_window_raw_counts,
          internal_start_window_bases,
          internal_upstream_window_raw_counts,
          internal_upstream_window_bases
        )]
  out[]
}
