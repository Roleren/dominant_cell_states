test_that("dcs_wilcox_p equals wilcox.test(exact = FALSE) exactly", {
  set.seed(11)
  for (i in seq_len(300)) {
    nx <- sample(1:12, 1); ny <- sample(1:12, 1)
    tied <- i %% 3 == 0
    x <- if (tied) sample(0:3, nx, TRUE) / 7 else rexp(nx)
    y <- if (tied) sample(0:4, ny, TRUE) / 7 else rexp(ny) * runif(1, 0.5, 2)
    if (i %% 10 == 0) x <- c(x, NA, Inf)
    if (length(unique(c(x[is.finite(x)], y))) < 2) next
    expected <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
    expect_identical(dcs_wilcox_p(x, y), expected)
  }
})

test_that("dcs_wilcox_p errors like wilcox.test on empty samples", {
  expect_error(dcs_wilcox_p(numeric(), 1:3))
  expect_error(dcs_wilcox_p(c(NA, 1)[1], 1:3))
})
