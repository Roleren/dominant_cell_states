# Closed-form DerSimonian-Laird random-effects pooling.
#
# Nine RDG scripts carry their own copy of this estimator (random_effects_summary)
# with small variations (external weights, a floor on the standard error, an
# extra z column) and return a one-row data.table, whose construction costs
# more than the estimate itself in leave-one-out and transport loops. This
# function covers the variants and returns a plain list; the arithmetic is
# written in the same order as the script copies so results are identical.

#' DerSimonian-Laird random-effects summary
#'
#' @param effect Per-study effect estimates.
#' @param se Their standard errors.
#' @param weight Optional external weights (rescaled to mean 1), multiplying
#'   the inverse-variance weights.
#' @param se_floor Lower bound applied to `se` after filtering (0 for none).
#' @return List with n, effect, se, lower, upper, z, p, tau2, i2,
#'   direction_fraction and ci_excludes_zero.
#' @export
dcs_random_effects <- function(effect, se, weight = NULL, se_floor = 0) {
  effect <- suppressWarnings(as.numeric(effect))
  se <- suppressWarnings(as.numeric(se))
  keep <- is.finite(effect) & is.finite(se) & se > 0
  if (!is.null(weight)) {
    weight <- suppressWarnings(as.numeric(weight))
    keep <- keep & is.finite(weight) & weight > 0
  }
  if (!any(keep)) {
    return(list(n = 0L, effect = NA_real_, se = NA_real_, lower = NA_real_,
                upper = NA_real_, z = NA_real_, p = NA_real_, tau2 = NA_real_,
                i2 = NA_real_, direction_fraction = NA_real_,
                ci_excludes_zero = FALSE))
  }
  yi <- effect[keep]
  vi <- if (se_floor > 0) pmax(se[keep], se_floor)^2 else se[keep]^2
  if (is.null(weight)) {
    wi <- 1 / vi
  } else {
    external_weight <- weight[keep] / mean(weight[keep])
    wi <- external_weight / vi
  }
  fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - fixed)^2)
  df <- length(yi) - 1
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (df > 0 && c_val > 0) max(0, (q - df) / c_val) else 0
  wi_re <- if (is.null(weight)) 1 / (vi + tau2) else external_weight / (vi + tau2)
  mu <- sum(wi_re * yi) / sum(wi_re)
  se_mu <- sqrt(1 / sum(wi_re))
  lower <- mu - 1.96 * se_mu
  upper <- mu + 1.96 * se_mu
  z <- mu / se_mu
  direction_fraction <- if (sign(mu) == 0) mean(yi == 0) else mean(sign(yi) == sign(mu))
  list(
    n = length(yi), effect = mu, se = se_mu, lower = lower, upper = upper,
    z = z, p = 2 * stats::pnorm(-abs(z)), tau2 = tau2,
    i2 = if (q > 0 && df > 0) max(0, (q - df) / q) else 0,
    direction_fraction = direction_fraction,
    ci_excludes_zero = isTRUE(lower > 0 | upper < 0)
  )
}
