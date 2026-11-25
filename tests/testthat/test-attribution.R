test_that("split band is the quantile of balanced split differences", {
  x <- c(1, 2, 3, 4)
  # splits of 4 into 2+2: {1,2}|{3,4} -> -2, {1,3}|{2,4} -> -1, {1,4}|{2,3} -> 0, and mirrors
  d <- abs(c(-2, -1, 0, 0, 1, 2))
  expect_equal(dcs_split_band(x, level = 1), 2)
  expect_equal(dcs_split_band(x), unname(quantile(d, 0.95)))
  expect_true(is.na(dcs_split_band(c(1, NA))))
  expect_equal(dcs_split_band(c(1, 3)), 2)
  set.seed(3); y <- rnorm(20)
  expect_identical(dcs_split_band(y, max_splits = 50, seed = 2), dcs_split_band(y, max_splits = 50, seed = 2))
})

test_that("band claim needs the whole interval outside the band", {
  expect_equal(dcs_band_claim(c(0.3, 0.1, -0.9, -0.5), c(0.5, 0.5, -0.4, 0.5), 0.2),
               c(TRUE, FALSE, TRUE, FALSE))
  expect_true(is.na(dcs_band_claim(1, 2, NA)))
})

test_that("band-aware gene test is calibrated under its own null", {
  set.seed(4)
  n <- 5000; v <- runif(n, 0.001, 0.01); tau <- 0.1
  d <- rnorm(n, 0, sqrt(v + tau^2))
  g <- dcs_gene_deviation(d, v, tau)
  expect_lt(mean(g$p < 0.05), 0.07)
  expect_lte(sum(g$call), 2L)
  # ignoring tau (count-only) calls many
  expect_gt(mean(dcs_gene_deviation(d, v, 0)$call), 0.1)
  d2 <- d; d2[1:100] <- d2[1:100] + 1
  expect_gt(sum(dcs_gene_deviation(d2, v, tau)$call[1:100]), 80)
})

test_that("signature match finds a planted similarity and not a shuffled one", {
  set.seed(5)
  sig <- setNames(rnorm(300), paste0("g", 1:300))
  x <- sig * 0.5 + rnorm(300, 0, 0.5)
  m <- dcs_signature_match(x, sig, n_perm = 200)
  expect_gt(m$r, 0.5); expect_lt(m$p_perm, 0.01); expect_lt(m$p_z, 1e-6)
  y <- setNames(rnorm(300), names(sig))
  expect_gt(dcs_signature_match(y, sig, n_perm = 200)$p_perm, 0.01)
  # aligned by names when both are named
  m2 <- dcs_signature_match(rev(x), sig, n_perm = 50)
  expect_equal(m2$r, m$r)
  expect_true(is.na(dcs_signature_match(1:3, 3:1)$r))
})

test_that("evidence rung and wording", {
  r <- dcs_evidence_rung(c(TRUE, TRUE, TRUE, FALSE, TRUE), feature = c(FALSE, TRUE, TRUE, TRUE, FALSE),
                         term = c(FALSE, FALSE, TRUE, TRUE, FALSE), signature = FALSE,
                         mechanism = c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(r, c("A0", "A1", "A2", "none", "A4"))
  expect_equal(dcs_claim_wording("A4", regulator = "eIF5A"), "because of eIF5A")
  expect_false(grepl("because", dcs_claim_wording("A3", regulator = "eIF5A")))
  w <- dcs_claim_wording(c("A0", "A3", "A4", "none"), feature = "x", regulator = c("a", "b", "c", "d"))
  expect_equal(w, c("deviates", "matches the signature of b", "because of c", "no deviation beyond the band"))
})
