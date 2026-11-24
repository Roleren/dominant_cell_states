# Reference copies of the script-local estimators replaced by dcs_random_effects()
# (dominant_rdg_model_agreement_loo.R / _transportability.R and
# dominant_rdg_branch_usage_model.R as of 2026-11-24).
safe_num <- function(x) suppressWarnings(as.numeric(x))
legacy_re_floor <- function(effect, se) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  keep <- is.finite(effect) & is.finite(se) & se > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L, effect = NA_real_, se = NA_real_, lower = NA_real_,
      upper = NA_real_, p = NA_real_, tau2 = NA_real_, i2 = NA_real_,
      direction_fraction = NA_real_, ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  vi <- pmax(se[keep], 1e-6)^2
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  q_stat <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q_stat - df) / c_val) else 0
  wi_re <- 1 / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  direction_fraction <- if (sign(mu) == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign(mu))
  }
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    p = 2 * pnorm(-abs(mu / se_mu)),
    tau2 = tau2,
    i2 = if (q_stat > 0 && df > 0) max(0, (q_stat - df) / q_stat) else 0,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

legacy_re_weighted <- function(effect, se, weight = NULL) {
  effect <- suppressWarnings(as.numeric(effect))
  se <- suppressWarnings(as.numeric(se))
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- suppressWarnings(as.numeric(weight))
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L,
      effect = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      z = NA_real_,
      p = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  vi <- se[keep]^2
  external_weight <- weight[keep] / mean(weight[keep])
  wi <- external_weight / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q - df) / c_val) else 0
  wi_re <- external_weight / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  z <- mu / se_mu
  p <- 2 * pnorm(-abs(z))
  sign_ref <- sign(mu)
  direction_fraction <- if (sign_ref == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign_ref)
  }
  i2 <- if (q > 0 && df > 0) max(0, (q - df) / q) else 0
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    z = z,
    p = p,
    tau2 = tau2,
    i2 = i2,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

test_that("dcs_random_effects reproduces the script estimators exactly", {
  skip_if_not_installed("data.table")
  library(data.table)
  set.seed(5)
  for (i in seq_len(400)) {
    k <- sample(0:9, 1)
    effect <- rnorm(k, sd = runif(1, 0, 0.3))
    se <- rexp(k) * 0.1
    if (k && i %% 7 == 0) se[1] <- 1e-8
    if (k && i %% 5 == 0) effect[sample(k, 1)] <- NA
    if (k && i %% 11 == 0) effect[] <- 0
    w <- runif(k, 0.1, 3)
    a <- as.list(legacy_re_floor(effect, se))
    b <- dcs_random_effects(effect, se, se_floor = 1e-6)
    expect_identical(b[names(a)], a)
    a <- as.list(legacy_re_weighted(effect, se))
    expect_identical(dcs_random_effects(effect, se)[names(a)], a)
    a <- as.list(legacy_re_weighted(effect, se, w))
    expect_identical(dcs_random_effects(effect, se, w)[names(a)], a)
  }
})

legacy_re_floor_weighted <- function(effect, se, weight = NULL) {
  effect <- safe_num(effect)
  se <- safe_num(se)
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- safe_num(weight)
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L,
      effect = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      p = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      q_stat = NA_real_,
      direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  sei <- pmax(se[keep], 1e-6)
  wi_external <- weight[keep] / mean(weight[keep])
  vi <- sei^2
  wi <- wi_external / vi
  fixed <- sum(wi * yi) / sum(wi)
  q_stat <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q_stat - df) / c_val) else 0
  wi_re <- wi_external / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  z <- mu / se_mu
  p <- 2 * pnorm(-abs(z))
  sign_ref <- sign(mu)
  direction_fraction <- if (sign_ref == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign_ref)
  }
  i2 <- if (q_stat > 0 && df > 0) max(0, (q_stat - df) / q_stat) else 0
  data.table(
    n = length(yi),
    effect = mu,
    se = se_mu,
    lower = lower,
    upper = upper,
    p = p,
    tau2 = tau2,
    i2 = i2,
    q_stat = q_stat,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}
legacy_re_atlas <- function(effect, se, weight = NULL) {
  effect <- suppressWarnings(as.numeric(effect))
  se <- suppressWarnings(as.numeric(se))
  if (is.null(weight)) weight <- rep(1, length(effect))
  weight <- suppressWarnings(as.numeric(weight))
  keep <- is.finite(effect) & is.finite(se) & se > 0 &
    is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(data.table(
      n = 0L,
      effect = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      q = NA_real_,
      direction_fraction = NA_real_,
      ci_excludes_zero = FALSE
    ))
  }
  yi <- effect[keep]
  vi <- se[keep]^2
  external_weight <- weight[keep] / mean(weight[keep])
  wi <- external_weight / vi
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q - df) / c_val) else 0
  wi_re <- external_weight / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  sign_ref <- sign(mu)
  direction_fraction <- if (sign_ref == 0) {
    mean(yi == 0)
  } else {
    mean(sign(yi) == sign_ref)
  }
  i2 <- if (q > 0 && df > 0) max(0, (q - df) / q) else 0
  data.table(
    n = length(yi),
    effect = mu,
    lower = lower,
    upper = upper,
    tau2 = tau2,
    i2 = i2,
    q = q,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}

test_that("dcs_random_effects reproduces the floored weighted and atlas variants", {
  library(data.table)
  set.seed(6)
  for (i in seq_len(300)) {
    k <- sample(0:9, 1)
    effect <- rnorm(k, sd = 0.2)
    se <- rexp(k) * 0.1
    if (k && i %% 7 == 0) se[1] <- 1e-8
    if (k && i %% 5 == 0) effect[sample(k, 1)] <- NA
    w <- if (i %% 2) runif(k, 0.1, 3) else NULL
    a <- as.list(legacy_re_floor_weighted(effect, se, w))
    expect_identical(dcs_random_effects(effect, se, w, se_floor = 1e-6)[names(a)], a)
    a <- as.list(legacy_re_atlas(effect, se, w))
    b <- dcs_random_effects(effect, se, w)
    b$q <- b$q_stat
    expect_identical(b[names(a)], a)
  }
})
