test_that("internal ORF signal metrics capture frame and start-window evidence", {
  root <- dcs_project_root()
  source(file.path(root, "scripts", "dominant_internal_orf_signal_helpers.R"))

  mat <- matrix(1, nrow = 80L, ncol = 2L)
  colnames(mat) <- c("RUN1", "RUN2")
  cds_positions <- 10:70
  feature_positions <- 20:49
  feature_phase <- (seq_along(feature_positions) - 1L) %% 3L
  mat[feature_positions[feature_phase == 0L], ] <- 12
  mat[head(feature_positions, 12L), ] <- mat[head(feature_positions, 12L), ] + 12

  metrics <- dominant_internal_orf_signal_metrics(
    mat,
    feature_positions = feature_positions,
    cds_positions = cds_positions,
    feature_type = "internal_orf",
    category = "internal"
  )
  expect_equal(nrow(metrics), 2L)
  expect_true(all(metrics$internal_frame_relative_to_cds == 1L))
  expect_true(all(metrics$internal_iorf_frame_fraction >
                    metrics$internal_canonical_frame_fraction))
  expect_true(all(metrics$internal_start_vs_body_density_ratio > 1))
  expect_true(all(metrics$internal_start_vs_upstream_density_ratio > 1))

  non_internal <- dominant_internal_orf_signal_metrics(
    mat,
    feature_positions = feature_positions,
    cds_positions = cds_positions,
    feature_type = "leader_uorf",
    category = "uORF"
  )
  expect_true(all(is.na(non_internal$internal_phase_total_raw_counts)))

  short_internal <- dominant_internal_orf_signal_metrics(
    mat,
    feature_positions = 20:28,
    cds_positions = cds_positions,
    feature_type = "internal_orf",
    category = "internal"
  )
  expect_true(all(short_internal$internal_feature_core_bases == 0L))
  expect_true(all(is.na(short_internal$internal_start_vs_body_density_ratio)))
})
