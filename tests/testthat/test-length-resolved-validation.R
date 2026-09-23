test_that("read-length strata distinguish UV-short and canonical footprints", {
  lengths <- c(14, 15, 20, 23, 24, 28, 32, 33, 35, NA)

  expect_equal(
    dcs_read_length_stratum(lengths),
    c(
      NA, "other_15_34nt", "short_20_23nt", "short_20_23nt",
      "other_15_34nt", "canonical_28_32nt", "canonical_28_32nt",
      "other_15_34nt", NA, NA
    )
  )
})

test_that("P-site calibration QC uses ORFik offset convention", {
  relative_five_prime <- c(
    rep(c(-12L, -11L, -10L), c(80L, 10L, 10L)),
    rep(c(-12L, -11L, -10L), c(80L, 10L, 10L))
  )
  read_length <- c(rep(28L, 100L), rep(29L, 100L))
  shifts <- data.frame(
    fraction = c(28L, 29L), offsets_start = c(-12L, -12L)
  )

  qc <- dcs_psite_calibration_qc(
    relative_five_prime, read_length, shifts, downstream = 20L
  )

  expect_equal(qc$phase0_fraction[qc$read_length == 28L], 0.8)
  expect_true(qc$periodicity_gate[qc$read_length == 28L])
  expect_equal(qc$phase0_fraction[qc$read_length == 29L], 0.8)
  expect_true(qc$periodicity_gate[qc$read_length == 29L])
})

test_that("P-site calibration rejects ambiguous shift tables", {
  expect_error(
    dcs_psite_calibration_qc(
      0L, 28L,
      data.frame(fraction = c(28L, 28L), offsets_start = c(-12L, -13L))
    ),
    "duplicated"
  )
  expect_error(
    dcs_psite_calibration_qc(
      0L, 28L, data.frame(fraction = 28L, offsets_start = 12L)
    ),
    "non-positive"
  )
  expect_error(
    dcs_psite_calibration_qc(
      0L, 16L, data.frame(fraction = 16L, offsets_start = -18L)
    ),
    "smaller than the footprint length"
  )
})

test_that("position leaveout audit distinguishes broad and peak-only effects", {
  broad <- expand.grid(
    category = "broad", condition = c("control", "UV"),
    region = c("uorf", "clean_CDS"), position = 1:6,
    stringsAsFactors = FALSE
  )
  broad$counts <- with(broad, ifelse(
    condition == "UV" & region == "uorf", 10,
    ifelse(condition == "UV" & region == "clean_CDS", 5,
           ifelse(region == "uorf", 2, 12))
  ))
  broad$position[broad$region == "clean_CDS"] <-
    broad$position[broad$region == "clean_CDS"] + 20L

  peak <- data.frame(
    category = "peak",
    condition = rep(c("control", "UV"), each = 4L),
    region = rep(c("uorf", "clean_CDS"), times = 4L),
    position = rep(c(1L, 21L), times = 4L),
    counts = c(1, 10, 1, 10, 20, 10, 1, 10)
  )
  audit <- dcs_position_leaveout_audit(
    rbind(broad, peak), case_label = "UV", control_label = "control"
  )$summary

  expect_true(audit$single_position_direction_stable[
    audit$category == "broad"
  ])
  expect_true(audit$greedy_direction_stable[audit$category == "broad"])
  expect_false(audit$single_position_direction_stable[
    audit$category == "peak"
  ])
})
