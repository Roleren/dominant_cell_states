test_that("legacy result links are removed, not created", {
  root <- tempfile("dcs-root-")
  dir.create(root)
  dir.create(file.path(root, "scripts"))
  dir.create(file.path(root, "inst", "shiny", "rdg_atlas_browser"),
             recursive = TRUE)
  dir.create(file.path(root, "results", "dominant_example"), recursive = TRUE)
  dir.create(file.path(root, "md"))
  writeLines("note", file.path(root, "md",
                               "dominant_cell_state_scientific_note.md"))

  legacy <- file.path(root, "dominant_example")
  expect_true(file.symlink("results/dominant_example", legacy))

  removed <- dcs_remove_legacy_result_links(root = root)

  expect_true(legacy %in% removed)
  expect_true(is.na(Sys.readlink(legacy)) || !nzchar(Sys.readlink(legacy)))
})

test_that("manual translon legacy links are removed", {
  root <- tempfile("dcs-root-")
  dir.create(file.path(root, "manual_translons"), recursive = TRUE)
  dir.create(file.path(root, "scripts"))
  dir.create(file.path(root, "inst", "shiny", "rdg_atlas_browser"),
             recursive = TRUE)
  dir.create(file.path(root, "results", "manual_translons"), recursive = TRUE)
  dir.create(file.path(root, "md"))
  writeLines("note", file.path(root, "md",
                               "dominant_cell_state_scientific_note.md"))

  link <- file.path(root, "manual_translons", "manual_translons_manifest.csv")
  expect_true(file.symlink(
    "../results/manual_translons/manual_translons_manifest.csv",
    link
  ))

  removed <- dcs_remove_legacy_result_links(root = root)

  expect_true(link %in% removed)
  expect_true(is.na(Sys.readlink(link)) || !nzchar(Sys.readlink(link)))
})
