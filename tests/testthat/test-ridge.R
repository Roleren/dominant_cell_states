ridge_fixture <- function(seed = 7, n = 1500, groups = 25) {
  set.seed(seed)
  g <- sample(groups, n, TRUE)
  d <- data.frame(g = factor(g), f = factor(sample(letters[1:6], n, TRUE)),
                  a = rnorm(n), b = rexp(n))
  y <- 0.3 * d$a + (as.integer(d$f) %% 3) * 0.2 + rnorm(groups)[g] + rnorm(n)
  # Group dummies are constant (zero) on the training rows of their own fold.
  x <- Matrix::sparse.model.matrix(~ g + f + a + b + f:a, d)[, -1]
  list(x = x, y = as.numeric(y), block = paste0("block_", g))
}

tight <- list(thresh = 1e-14)

test_that("closed-form ridge reproduces converged glmnet coefficients", {
  skip_if_not_installed("glmnet")
  fx <- ridge_fixture()
  lambda <- c(2, 0.3, 0.01)
  for (l in lambda) {
    ref <- as.matrix(stats::coef(glmnet::glmnet(fx$x, fx$y, alpha = 0, lambda = l,
                                                control = tight), s = l))
    expect_lt(max(abs(dcs_ridge_fit(fx$x, fx$y, l) - ref)), 1e-5)
  }
})

test_that("closed-form cross-validation matches cv.glmnet lambda choice", {
  skip_if_not_installed("glmnet")
  fx <- ridge_fixture()
  foldid <- match(fx$block, sort(unique(fx$block)))
  ref <- glmnet::cv.glmnet(fx$x, fx$y, alpha = 0, foldid = foldid,
                           nlambda = 80, control = tight)
  fast <- dcs_cv_ridge(fx$x, fx$y, foldid, nlambda = 80)
  expect_equal(fast$lambda, ref$lambda, tolerance = 1e-12)
  expect_lt(max(abs(fast$cvm / ref$cvm - 1)), 1e-4)
  expect_lt(max(abs(fast$cvsd / ref$cvsd - 1)), 1e-3)
  # Same grid point; the grid itself agrees to the last bits.
  expect_identical(match(fast$lambda.min, fast$lambda), match(ref$lambda.min, ref$lambda))
  expect_identical(match(fast$lambda.1se, fast$lambda), match(ref$lambda.1se, ref$lambda))
})

test_that("blocked ridge engines agree on predictions and coefficients", {
  skip_if_not_installed("glmnet")
  fx <- ridge_fixture(seed = 3, n = 900, groups = 12)
  fast <- dcs_blocked_ridge(fx$x, fx$y, fx$block, nlambda = 80,
                            min_train = 100, engine = "fast")
  slow <- dcs_blocked_ridge(fx$x, fx$y, fx$block, nlambda = 80,
                            min_train = 100, engine = "glmnet")
  expect_equal(fast$lambda.min, slow$lambda.min, tolerance = 1e-12)
  expect_equal(fast$lambda.1se, slow$lambda.1se, tolerance = 1e-12)
  pf <- fast$predict(fast$lambda.min)
  ps <- slow$predict(slow$lambda.min)
  expect_identical(pf$train_rows, ps$train_rows)
  # glmnet's default convergence threshold leaves ~1e-4 differences.
  expect_lt(max(abs(pf$prediction - ps$prediction)), 1e-3 * stats::sd(fx$y))
  expect_lt(max(abs(fast$coef(fast$lambda.min) - slow$coef(slow$lambda.min))), 1e-3)
  expect_identical(rownames(fast$coef(1)), rownames(slow$coef(1)))
})

test_that("blocks below min_train are left unpredicted", {
  fx <- ridge_fixture(seed = 5, n = 300, groups = 4)
  fit <- dcs_blocked_ridge(fx$x, fx$y, fx$block, nlambda = 20,
                           min_train = 240, engine = "fast")
  p <- fit$predict(fit$lambda.min)
  train_n <- length(fx$y) - table(fx$block)[fx$block]
  expect_identical(is.na(p$prediction), as.vector(train_n < 240))
  expect_identical(is.na(p$train_rows), as.vector(train_n < 240))
})

test_that("auto engine keeps glmnet for wide designs", {
  skip_if_not_installed("glmnet")
  fx <- ridge_fixture(seed = 9, n = 400, groups = 8)
  expect_identical(dcs_blocked_ridge(fx$x, fx$y, fx$block, nlambda = 10,
                                     engine = "auto")$engine, "fast")
  expect_identical(dcs_blocked_ridge(fx$x, fx$y, fx$block, nlambda = 10,
                                     engine = "auto",
                                     max_fast_features = 5)$engine, "glmnet")
})
