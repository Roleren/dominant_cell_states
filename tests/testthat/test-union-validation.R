test_that("one-row region metrics combine without an implicit merge key", {
  uorf <- data.frame(uorf_counts = 42, uorf_density = 1.5)
  cds <- data.frame(clean_cds_counts = 150, clean_cds_density = 2.5)

  combined <- dcs_bind_one_row_metrics(uorf, cds)

  expect_equal(nrow(combined), 1L)
  expect_named(combined, c(
    "uorf_counts", "uorf_density", "clean_cds_counts", "clean_cds_density"
  ))
  expect_equal(combined$clean_cds_counts, 150)
})

test_that("one-row region metric binding rejects ambiguous inputs", {
  expect_error(
    dcs_bind_one_row_metrics(data.frame(value = 1), data.frame(value = 2)),
    "must be unique"
  )
  expect_error(
    dcs_bind_one_row_metrics(data.frame(value = 1:2), data.frame(other = 3)),
    "one-row"
  )
})

test_that("manual context identity excludes selector and ignored lane axes", {
  fields <- c("study", "CONDITION", "FRACTION", "CELL_LINE")

  expect_equal(
    dcs_manual_identity_fields(fields, "CONDITION", "FRACTION"),
    c("study", "CELL_LINE")
  )
  expect_error(
    dcs_manual_identity_fields(fields, ""),
    "one non-empty string"
  )
})

test_that("explicit replicate labels collapse technical lanes", {
  metadata <- data.frame(
    Run = paste0("run", 1:4),
    Experiment = paste0("experiment", 1:4),
    BioSample = paste0("sample", 1:4),
    SampleName = paste0("name", 1:4),
    REPLICATE = c(1, 1, 2, 2)
  )

  units <- dcs_biological_unit_ids(metadata, "case")

  expect_equal(units, c(
    "case::replicate_1", "case::replicate_1",
    "case::replicate_2", "case::replicate_2"
  ))
})

test_that("biological-unit IDs fall back when replicate labels are absent", {
  metadata <- data.frame(
    Run = c("run1", "run2", "run3"),
    Experiment = c("experiment1", "experiment1", ""),
    BioSample = c("sample1", "sample1", "sample3"),
    SampleName = c("name1", "name1", "name3"),
    REPLICATE = c("MISSING", "", NA)
  )

  units <- dcs_biological_unit_ids(metadata, "control")

  expect_equal(units, c(
    "control::experiment1", "control::experiment1", "control::sample3"
  ))
})

test_that("uORF phase positions retain only phase-agreeing overlaps", {
  annotation <- data.frame(
    tx_start = c(2L, 5L, 6L),
    tx_end = c(7L, 10L, 8L)
  )

  phase <- dcs_unambiguous_uorf_phase_positions(annotation, 12L)

  expect_equal(phase$phase0, c(2L, 5L))
  expect_equal(phase$phase1, c(3L, 9L))
  expect_equal(phase$phase2, c(4L, 10L))
  expect_equal(phase$informative, c(2L, 3L, 4L, 5L, 9L, 10L))
  expect_false(any(6:8 %in% phase$informative))
})
