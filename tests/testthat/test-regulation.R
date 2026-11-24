test_that("held-out draw is seeded, stratified, honours exclusions and guarantees", {
  st <- sprintf("S%03d", 1:100)
  strat <- rep(c("a", "b", "c", "d"), times = c(50, 30, 15, 5))
  tags <- data.frame(study = c("S001", "S002", "S051", "S099", "S100"),
                     tag = c("t1", "t1", "t2", "t3", "t3"))
  args <- list(study = st, stratum = strat, tags = tags,
               min_carriers = c(t1 = 2L, t2 = 2L, t3 = 2L),
               exclude = c("S001", "S099"), fraction = 0.2, seed = 7L)
  x <- do.call(dcs_draw_heldout_studies, args)
  y <- do.call(dcs_draw_heldout_studies, args)
  expect_identical(x, y)
  expect_false(any(x$heldout[x$study %in% c("S001", "S099")]))
  # t1 has one eligible carrier (S002) -> below min_carriers, not guaranteed
  expect_false(any(grepl("t1", x$reason)))
  # t3 has one eligible carrier too; t2 has one -> none guaranteed either
  expect_equal(sum(grepl("guarantee", x$reason)), 0L)
  args$min_carriers <- c(t1 = 1L, t3 = 1L)
  z <- do.call(dcs_draw_heldout_studies, args)
  expect_true(z$heldout[z$study == "S002"])
  expect_true(z$heldout[z$study == "S100"])
  share <- tapply(z$heldout, z$stratum, mean)
  expect_true(all(abs(share - 0.2) <= 0.21))
  expect_equal(sum(z$heldout[z$stratum == "a"]), 10L)
})

test_that("leave-study-out deviation matches the brute-force mean", {
  set.seed(1)
  M <- matrix(rnorm(40), 5, 8)
  st <- c("a", "a", "b", "c", "c", "c", "d", "e")
  D <- dcs_leave_study_out_deviation(M, st)
  brute <- sapply(seq_len(ncol(M)), function(i) M[, i] - rowMeans(M[, st != st[i], drop = FALSE]))
  expect_equal(unname(D), unname(brute))
})

test_that("poolability separates a replicated group from a random one", {
  set.seed(2)
  G <- 400; n <- 40
  study <- sprintf("s%02d", seq_len(n))
  group <- rep(c("real", "fake", NA), times = c(6, 6, 28))
  sig <- rnorm(G)
  M <- matrix(rnorm(G * n), G, n)
  M[, group %in% "real"] <- M[, group %in% "real"] + 0.8 * sig
  cls <- rep(c("chx", "none"), length.out = n)
  res <- dcs_poolability(M, group, study, unit_class = cls, n_perm = 300, seed = 3)
  real <- res[res$group == "real", ]
  fake <- res[res$group == "fake", ]
  expect_gt(real$replication, 0.2)
  expect_lt(real$p_perm, 0.01)
  expect_gt(fake$p_perm, 0.01)
  expect_lt(abs(fake$replication), 0.1)
  # one-unit-per-study rule and the minimum size
  res2 <- dcs_poolability(M[, 1:7], c(rep("g", 2), rep(NA, 5)), study[1:7], n_perm = 10)
  expect_null(res2)
})

test_that("window and gene contrasts are zero for identical arms", {
  a <- c(300, 500, 250, 1000, 220)
  gene <- c("g1", "g1", "g1", "g2", "g2")
  w <- dcs_window_contrast(a, a, gene, min_reads = 200)
  expect_equal(w$d, rep(0, 5))
  expect_equal(w$v, 2 / a)
  g <- dcs_gene_contrast(c(2000, 5000, 100), c(4000, 10000, 300), min_reads = 1000)
  expect_equal(g$i, 1:2)
  expect_equal(g$d, c(0, 0))
  al <- dcs_gene_contrast(c(2000, 1500), c(1000, 3000), c(4000, 3000), c(2000, 6000), min_reads = 1000)
  expect_equal(al$d, c(0, 0))
})

test_that("placebo band recovers a known excess spread", {
  set.seed(4)
  n <- 20000
  lam <- runif(n, 300, 3000)
  ex <- rnorm(n, 0, 0.15)
  a <- rpois(n, lam * exp(ex)); b <- rpois(n, lam)
  d <- log(a / lam) - log(b / lam)
  bd <- dcs_placebo_band(d, 1 / a + 1 / b, robust = FALSE)
  expect_equal(unname(bd["sd_excess"]), 0.15, tolerance = 0.05)
  bd0 <- dcs_placebo_band(log(b / lam) - log(rpois(n, lam) / lam), 2 / lam)
  expect_lt(unname(bd0["sd_excess"]), 0.03)
  expect_lt(abs(unname(bd0["false_rate_count"]) - 0.05), 0.01)
})

test_that("motif amplitude recovers planted effects and matches glm", {
  set.seed(5)
  ng <- 60; L <- 80
  gene <- rep(seq_len(ng), each = L)
  x1 <- rbinom(ng * L, 1, 0.1); x2 <- rbinom(ng * L, 1, 0.3)
  mu_g <- rep(exp(rnorm(ng, 1, 0.5)), each = L)
  y <- rpois(ng * L, mu_g * exp(0.7 * x1 - 0.3 * x2))
  X <- cbind(pp = x1, other = x2)
  fit <- dcs_motif_amplitude(y, gene, X)
  expect_true(attr(fit, "converged"))
  ref <- stats::glm(y ~ factor(gene) + x1 + x2, family = poisson())
  expect_equal(fit$beta, unname(coef(ref)[c("x1", "x2")]), tolerance = 1e-6)
  expect_equal(fit$se_model, unname(sqrt(diag(vcov(ref)))[c("x1", "x2")]), tolerance = 1e-4)
  expect_true(all(abs(fit$beta - c(0.7, -0.3)) < 4 * fit$se))
})

test_that("arm contrast gives a Welch interval", {
  r <- dcs_arm_contrast(c(1.0, 1.2, 0.2, 0.4, 0.3), rep(0.05, 5),
                        c("case", "case", "control", "control", "control"))
  expect_equal(r$delta, 1.1 - 0.3)
  tt <- stats::t.test(c(1.0, 1.2), c(0.2, 0.4, 0.3))
  expect_equal(c(r$lower, r$upper), unname(tt$conf.int[1:2]), tolerance = 1e-8)
  r1 <- dcs_arm_contrast(c(1, 0.2, 0.4), c(0.1, 0.1, 0.1), c("case", "control", "control"))
  expect_true(is.infinite(r1$df))
  expect_equal(r1$se, sqrt(0.01 / 1 + var(c(0.2, 0.4)) / 2))
})

test_that("fill pass keeps the requested carriers of a tag outside the held-out set", {
  st <- sprintf("S%02d", 1:10)
  tags <- data.frame(study = c("S01", "S02"), tag = "rare")
  for (seed in 1:20) {
    x <- dcs_draw_heldout_studies(st, rep("a", 10), tags = tags,
                                  min_carriers = c(rare = 2L),
                                  keep_free = c(rare = 1L), fraction = 0.5,
                                  seed = seed)
    expect_equal(sum(x$heldout[x$study %in% c("S01", "S02")]), 1L)
    expect_equal(sum(x$heldout), 5L)
  }
})
