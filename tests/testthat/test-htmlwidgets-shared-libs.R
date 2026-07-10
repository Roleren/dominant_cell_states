test_that("shared HTML widget libraries are not recopied when unchanged", {
  root <- dcs_project_root()
  source(file.path(root, "scripts", "dominant_htmlwidgets.R"))

  tmp <- tempfile("dominant-htmlwidgets-test-")
  private <- file.path(tmp, "private")
  shared <- file.path(tmp, "shared")
  dir.create(file.path(private, "dep-one"), recursive = TRUE)
  dir.create(file.path(shared, "dep-one"), recursive = TRUE)

  writeLines("same-content", file.path(private, "dep-one", "dep.js"))
  writeLines("same-content", file.path(shared, "dep-one", "dep.js"))
  old_time <- as.POSIXct("2001-01-01 00:00:00", tz = "UTC")
  Sys.setFileTime(file.path(shared, "dep-one", "dep.js"), old_time)

  dominant_copy_widget_libs(private, shared)
  mtime_after_skip <- file.info(file.path(shared, "dep-one", "dep.js"))$mtime
  expect_lt(abs(as.numeric(difftime(mtime_after_skip, old_time, units = "secs"))), 2)
  expect_true(file.exists(file.path(shared, "dep-one.dominant_widget_signature")))

  writeLines("new-content", file.path(private, "dep-one", "dep.js"))
  dominant_copy_widget_libs(private, shared)
  target_content <- readLines(file.path(shared, "dep-one", "dep.js"), warn = FALSE)
  expect_identical(target_content, "new-content")
})
