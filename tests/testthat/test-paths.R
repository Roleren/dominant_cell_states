test_that("project and result paths resolve inside the standalone repo", {
  root <- dcs_project_root()
  expect_true(file.exists(file.path(root, "DESCRIPTION")))
  expect_true(dir.exists(file.path(root, "R")))
  expect_identical(dcs_source_path("scripts", root = root),
                   file.path(root, "scripts"))
  expect_identical(dcs_results_dir(root = root),
                   normalizePath(file.path(root, "results"),
                                 mustWork = FALSE))
  expect_identical(dcs_result_path("x.csv", root = root),
                   file.path(dcs_results_dir(root = root), "x.csv"))
})

test_that("RiboCrypt discovery refuses the analysis package DESCRIPTION", {
  old_repo <- Sys.getenv("RIBOCRYPT_REPO", unset = NA_character_)
  on.exit({
    if (is.na(old_repo)) {
      Sys.unsetenv("RIBOCRYPT_REPO")
    } else {
      Sys.setenv(RIBOCRYPT_REPO = old_repo)
    }
  }, add = TRUE)

  Sys.unsetenv("RIBOCRYPT_REPO")
  expect_identical(dcs_ribocrypt_repo(required = FALSE), "")
})
