# Deviation claims and attribution helpers for group-level contrasts:
# split bands from control-against-control contrasts, band-aware tests,
# signature similarity against a shuffled-signature null and the evidence
# ladder A0-A4 with its allowed wording.
#
# Plain vectors in, plain values out; reading data lives in the scripts.

#' Noise band from balanced splits of one arm
#'
#' For per-sample values of one arm (e.g. the control arm of a designed
#' contrast), the band is the `level` quantile of |mean(A) - mean(B)| over
#' balanced splits A / B of the samples (all splits, or `max_splits` drawn
#' at random when there are more).
#' @param x Numeric per-sample values (non-finite values are dropped).
#' @param level Quantile of the absolute split differences.
#' @param max_splits Largest number of splits evaluated.
#' @param seed Seed used when splits are subsampled.
#' @return Single number, NA with fewer than two finite values.
#' @export
dcs_split_band <- function(x, level = 0.95, max_splits = 200L, seed = 1L) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)
  sp <- utils::combn(n, floor(n / 2), simplify = FALSE)
  if (length(sp) > max_splits) {
    set.seed(seed)
    sp <- sp[sample.int(length(sp), max_splits)]
  }
  d <- vapply(sp, function(s) mean(x[s]) - mean(x[-s]), numeric(1))
  unname(stats::quantile(abs(d), level))
}

#' Does an interval exclude a noise band?
#'
#' A claim is an interval of the deviation that lies entirely above `band`
#' or entirely below `-band`.
#' @param lower,upper Interval bounds.
#' @param band Non-negative band (NA gives NA).
#' @return Logical vector.
#' @export
dcs_band_claim <- function(lower, upper, band) {
  out <- lower > band | upper < -band
  out[!is.finite(band)] <- NA
  out
}

#' Band-aware per-feature deviation test
#'
#' Each contrast d has variance v (counting noise) plus tau^2 (the excess
#' between-arm variance measured on placebo contrasts). Two-sided normal
#' p-values, Benjamini-Hochberg q-values across the features of one contrast.
#' @param d Contrasts (e.g. per gene).
#' @param v Their counting variances.
#' @param tau Excess SD (scalar or per feature).
#' @param alpha FDR level of a call.
#' @return data.frame with z, p, q and call.
#' @export
dcs_gene_deviation <- function(d, v, tau, alpha = 0.05) {
  stopifnot(length(d) == length(v))
  z <- d / sqrt(v + tau^2)
  p <- 2 * stats::pnorm(-abs(z))
  q <- stats::p.adjust(p, "BH")
  data.frame(z = z, p = p, q = q, call = q < alpha)
}

#' Similarity of a deviation profile to a signature, with a shuffled null
#'
#' Pearson correlation of `x` and `sig` over entries finite in both. The
#' null shuffles the signature values over the same entries `n_perm` times;
#' `p_perm` counts shuffles at least as similar (+1 smoothing) and `p_z` is
#' the upper normal tail of r standardised by the mean and SD of the shuffles
#' (resolves p-values below 1 / n_perm).
#' @param x Deviation profile (named or aligned with `sig`).
#' @param sig Signature profile.
#' @param n_perm Number of shuffles.
#' @param seed Random seed.
#' @param min_n Smallest number of shared entries.
#' @return List: n, r, p_perm, z, p_z, null_q95.
#' @export
dcs_signature_match <- function(x, sig, n_perm = 1000L, seed = 1L, min_n = 5L) {
  if (!is.null(names(x)) && !is.null(names(sig))) {
    k <- intersect(names(x), names(sig))
    x <- x[k]; sig <- sig[k]
  }
  stopifnot(length(x) == length(sig))
  ok <- is.finite(x) & is.finite(sig)
  x <- unname(x[ok]); sig <- unname(sig[ok])
  n <- length(x)
  if (n < min_n || stats::sd(x) == 0 || stats::sd(sig) == 0) {
    return(list(n = n, r = NA_real_, p_perm = NA_real_, z = NA_real_, p_z = NA_real_, null_q95 = NA_real_))
  }
  r <- stats::cor(x, sig)
  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(i) stats::cor(x, sig[sample.int(n)]), numeric(1))
  z <- (r - mean(null)) / stats::sd(null)
  list(n = n, r = r, p_perm = (1 + sum(null >= r)) / (1 + n_perm),
       z = z, p_z = stats::pnorm(-z), null_q95 = unname(stats::quantile(null, 0.95)))
}

#' Evidence ladder rung of a deviation
#'
#' A0 deviation exists; A1 associated with a feature against matched
#' controls; A2 reproduced by one term amplitude; A3 matches the signature of
#' a known regulator seen in 2+ studies; A4 the group is that perturbation,
#' the signature fits the mechanism and replicates across studies. Rungs are
#' cumulative only from A0: a higher rung needs a deviation, the others are
#' checked independently and the highest true one is returned.
#' @param deviates,feature,term,signature,mechanism Logical vectors.
#' @return Character: "none", "A0", ..., "A4".
#' @export
dcs_evidence_rung <- function(deviates, feature = FALSE, term = FALSE,
                              signature = FALSE, mechanism = FALSE) {
  n <- max(length(deviates), length(feature), length(term), length(signature), length(mechanism))
  f <- function(x) rep_len(x %in% TRUE, n)
  dv <- f(deviates)
  out <- ifelse(dv, "A0", "none")
  out[dv & f(feature)] <- "A1"
  out[dv & f(term)] <- "A2"
  out[dv & f(signature)] <- "A3"
  out[dv & f(mechanism)] <- "A4"
  out
}

#' Wording allowed at an evidence rung
#'
#' "because of" is reserved for A4.
#' @param rung Output of [dcs_evidence_rung()].
#' @param feature,regulator Names inserted in the phrase.
#' @return Character phrase.
#' @export
dcs_claim_wording <- function(rung, feature = "feature", regulator = "R") {
  unname(c(none = "no deviation beyond the band",
           A0 = "deviates",
           A1 = paste0("deviation is associated with ", feature),
           A2 = paste0("explained by term ", feature),
           A3 = paste0("matches the signature of ", regulator),
           A4 = paste0("because of ", regulator))[rung])
}
