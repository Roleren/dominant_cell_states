# Small statistics helpers shared by pipeline scripts.

# Base R >= 4.4 provides %||%; define it for older R.
if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

#' Two-sided Wilcoxon rank-sum p-value (normal approximation)
#'
#' Same value as `wilcox.test(x, y, exact = FALSE)$p.value` (continuity
#' correction, tie-corrected variance, missing values dropped), computed
#' with the identical sequence of floating point operations but without the
#' method dispatch, argument checks and `table()` call that dominate when the
#' test runs on thousands of small groups.
#'
#' @param x,y Numeric samples.
#' @return p-value; an error, like wilcox.test, when either sample is empty.
#' @export
dcs_wilcox_p <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 1L) stop("not enough (non-missing) 'x' observations")
  if (length(y) < 1L) stop("not enough 'y' observations")
  r <- rank(c(x, y))
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))
  statistic <- sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2
  nties <- tabulate(match(r, unique(r)))
  z <- statistic - n.x * n.y / 2
  sigma <- sqrt((n.x * n.y / 12) * ((n.x + n.y + 1) -
    sum(nties^3 - nties) / ((n.x + n.y) * (n.x + n.y - 1))))
  correction <- sign(z) * 0.5
  z <- (z - correction) / sigma
  2 * min(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
}
