legacy_env <- function() {
  env <- new.env()
  sys.source(test_path("fixtures", "legacy_condition_families.R"), envir = env)
  env
}

vocabulary <- c(
  "WT", "wild type", "Control", "ctrl", "mock", "DMSO", "untreated", "NTC",
  "scrambled", "shCtrl", "empty vector", "EV", "vector", "none", NA, "",
  "SARS-CoV-2 infected", "EBV", "IFNg 6h", "poly(I:C)", "arsenite", "Tunicamycin",
  "no BCAA", "noGln", "hypoxia 1%", "H2O2", "heat shock", "HS recovery", "NaCl",
  "UV", "etoposide", "senescent", "tumour", "normal adjacent", "KO", "METTL3 KD",
  "siRNA", "CRISPR", "OE", "rescue", "mutant", "SNP", "Torin1", "rapamycin",
  "treated", "LCL", "macrophage", "A549", "polysome", "cyto", "membrane",
  "timecourse 2h", "day 3", "Ribo-seq", "harringtonine", "LPS", "ver-155008"
)

test_that("condition vocabulary matches the former script helpers exactly", {
  old <- legacy_env()
  set.seed(8)
  n <- 3000
  pick <- function() sample(vocabulary, n, TRUE)
  cond <- pick(); frac <- pick(); inh <- pick(); cell <- pick()
  tissue <- pick(); gene <- pick(); cancer <- pick()
  expect_identical(dcs_clean_level(cond), old$dominant_clean_level(cond))
  expect_identical(dcs_normalise_label(cond), old$dominant_normalise_label(cond))
  expect_identical(dcs_is_control_condition(cond), old$dominant_is_control_condition(cond))
  fam_new <- dcs_classify_design_family(cond, frac, inh, cell, tissue, gene, cancer)
  fam_old <- old$dominant_classify_design_family(cond, frac, inh, cell, tissue, gene, cancer)
  expect_identical(fam_new, fam_old)
  expect_identical(dcs_classify_design_family(cond), old$dominant_classify_design_family(cond))
  expect_identical(dcs_classify_perturbation_class(cond, fam_new),
                   old$dominant_classify_perturbation_class(cond, fam_old))
  expect_identical(dcs_classify_perturbation_class(cond),
                   old$dominant_classify_perturbation_class(cond))
  expect_identical(dcs_classify_perturbation_class(design_family = c(NA, "control")),
                   old$dominant_classify_perturbation_class(design_family = c(NA, "control")))
})

test_that("the script shim keeps the dominant_* names", {
  env <- new.env()
  withr::local_envvar(DOMINANT_CELL_STATES_DIR = dcs_project_root())
  sys.source(file.path(dcs_project_root(), "scripts", "dominant_condition_families.R"),
             envir = env)
  expect_identical(env$dominant_classify_design_family("KO"), "loss_of_function")
  expect_identical(env$dominant_is_control_condition(c("WT", "KO")), c(TRUE, FALSE))
})

test_that("stratum ids paste present fields and can mark missing values", {
  dt <- data.frame(study = c("S1", "S1", NA), CELL_LINE = c("HeLa", NA, "HEK293"))
  expect_identical(dcs_stratum_id(dt, c("study", "CELL_LINE", "absent")),
                   c("S1 | HeLa", "S1 | NA", "NA | HEK293"))
  expect_identical(dcs_stratum_id(dt, c("study", "CELL_LINE"), missing = "MISSING"),
                   c("S1 | HeLa", "S1 | MISSING", "MISSING | HEK293"))
})

test_that("tight controls take precedence over loose ones within a group", {
  cond <- c("WT", "scrambled", "KD", "WT", "KO", "parental", " NTC ")
  study <- c("A", "A", "A", "B", "B", "C", "C")
  out <- dcs_control_tiers(cond, study)
  expect_identical(out$is_control, c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(out$group_has_tight, c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  two_keys <- dcs_control_tiers(cond, list(study, c("x", "y", "x", "x", "x", "x", "x")))
  expect_identical(two_keys$is_control, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
})

test_that("contrast units need a control in the same stratum", {
  dt <- data.frame(
    Run = paste0("r", 1:9),
    stratum_id = c("s1", "s1", "s1", "s1", "s2", "s2", "s3", "s3", "s3"),
    CONDITION = c("WT", "WT", "KO", "KO", "KO", "OE", "WT", "KD", NA),
    is_control_condition = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    design_family = c("control", "control", "loss_of_function", "loss_of_function",
                      "loss_of_function", "gain_of_function", "control",
                      "loss_of_function", "other_non_control"),
    stringsAsFactors = FALSE
  )
  units <- dcs_contrast_units(dt, carry = "design_family")
  expect_identical(units$stratum_id, c("s1", "s3"))
  expect_identical(units$case_condition, c("KO", "KD"))
  expect_identical(units$n_case_samples, c(2L, 1L))
  expect_identical(units$n_control_samples, c(2L, 1L))
  expect_identical(units$case_runs, c("r3;r4", "r8"))
  expect_identical(units$control_runs, c("r1;r2", "r7"))
  expect_identical(units$design_family, c("loss_of_function", "loss_of_function"))
  expect_identical(units$sample_set_signature[[1]], "s1 || KO || r3;r4 || r1;r2")
  rows <- dcs_contrast_rows(dt$stratum_id, dt$CONDITION, dt$is_control_condition,
                            "s3", "KD")
  expect_identical(rows, 7:8)
  empty <- dcs_contrast_units(dt[dt$stratum_id == "s2", ])
  expect_identical(nrow(empty), 0L)
})

test_that("contrast rows treat NA comparisons as no match", {
  # Same rule as data.table's i = stratum_id == s & (is_control == TRUE |
  # CONDITION == k), which the scripts used: NA selects nothing.
  set.seed(2)
  n <- 400
  stratum <- sample(c("a", "b", NA), n, TRUE)
  cond <- sample(c("WT", "KO", "KD", NA), n, TRUE)
  ctrl <- sample(c(TRUE, FALSE, NA), n, TRUE)
  for (s in c("a", "b")) for (k in c("KO", "KD")) {
    ref <- which(stratum == s & (ctrl == TRUE | cond == k))
    expect_identical(dcs_contrast_rows(stratum, cond, ctrl, s, k), ref)
  }
})
