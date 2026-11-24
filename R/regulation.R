# Group-level helpers: held-out study sets, poolability of a group across
# independent studies, control-against-control (placebo) contrasts and
# within-gene motif amplitudes with arm contrasts.
#
# All functions work on plain vectors / matrices so they can be tested on
# small synthetic data; reading Bio_data lives in the calling scripts.

# ---------------------------------------------------------------------------
# Held-out study sets
# ---------------------------------------------------------------------------

#' Draw a stratified held-out set of studies
#'
#' Two passes, both seeded. (1) Guarantee: for every tag listed in
#' `min_carriers` (smallest tags first), if the tag has at least that many
#' eligible studies and none of them is held out yet, one carrier is drawn,
#' from the `prefer`-red carriers when there are any. (2) Fill: in every
#' stratum, studies are added at random until the stratum holds
#' `fraction` of its studies (stochastic rounding of the remainder), counting
#' the studies already held out by the guarantee pass. A fill draw is skipped
#' when it would leave fewer than `keep_free[tag]` carriers of a tag outside
#' the held-out set (excluded studies count as outside).
#'
#' @param study Character vector of unique study ids.
#' @param stratum Character vector, the stratum of each study.
#' @param tags Optional data.frame with columns `study` and `tag` (a study can
#'   carry several tags).
#' @param min_carriers Named integer vector: tag -> smallest number of
#'   eligible carriers for which one held-out study is guaranteed.
#' @param prefer Optional logical vector (per study), preferred carriers for
#'   the guarantee pass (e.g. replicated, controlled designs).
#' @param keep_free Named integer vector: tag -> carriers that must stay
#'   outside the held-out set after the fill pass (e.g. 1 to keep one study
#'   of a class for development).
#' @param exclude Studies that may not be drawn (e.g. already analysed
#'   development data).
#' @param fraction Target held-out share per stratum.
#' @param seed Random seed.
#' @return data.frame with study, stratum, eligible, heldout and reason.
#' @export
dcs_draw_heldout_studies <- function(study, stratum, tags = NULL,
                                     min_carriers = integer(), prefer = NULL,
                                     keep_free = integer(),
                                     exclude = character(), fraction = 0.2,
                                     seed = 1L) {
  study <- as.character(study)
  stratum <- as.character(stratum)
  stopifnot(length(study) == length(stratum), !anyDuplicated(study),
            fraction >= 0, fraction <= 1)
  if (is.null(prefer)) prefer <- rep(FALSE, length(study))
  stopifnot(length(prefer) == length(study))
  prefer[is.na(prefer)] <- FALSE
  pick_one <- function(x) x[sample.int(length(x), 1L)]
  set.seed(seed)
  eligible <- !study %in% exclude
  held <- rep(FALSE, length(study))
  reason <- rep("", length(study))
  all_carriers <- list()
  if (!is.null(tags)) {
    tags <- as.data.frame(tags, stringsAsFactors = FALSE)
    tags <- tags[tags$study %in% study, , drop = FALSE]
    all_carriers <- lapply(split(match(tags$study, study), tags$tag), unique)
  }
  keep_free <- keep_free[names(keep_free) %in% names(all_carriers)]
  breaches <- function(k) {
    for (tg in names(keep_free)) {
      ix <- all_carriers[[tg]]
      if (k %in% ix && sum(!held[ix]) - 1L < keep_free[[tg]]) return(TRUE)
    }
    FALSE
  }
  if (length(min_carriers) && !is.null(tags)) {
    tags <- tags[tags$tag %in% names(min_carriers), , drop = FALSE]
    carriers <- split(match(tags$study, study), tags$tag)
    carriers <- lapply(carriers, function(i) sort(unique(i[eligible[i]])))
    n_c <- vapply(carriers, length, integer(1))
    for (tg in names(carriers)[order(n_c, names(carriers))]) {
      ix <- carriers[[tg]]
      if (length(ix) < min_carriers[[tg]] || any(held[ix])) next
      pool <- if (any(prefer[ix])) ix[prefer[ix]] else ix
      k <- pick_one(pool)
      held[k] <- TRUE
      reason[k] <- paste0("guarantee: ", tg)
    }
  }
  for (s in sort(unique(stratum))) {
    ix <- which(stratum == s)
    want <- fraction * length(ix)
    target <- floor(want) + (stats::runif(1) < want - floor(want))
    need <- target - sum(held[ix])
    free <- ix[eligible[ix] & !held[ix]]
    if (need > 0 && length(free)) {
      for (k in free[sample.int(length(free))]) {
        if (need <= 0) break
        if (breaches(k)) next
        held[k] <- TRUE
        reason[k] <- "stratified"
        need <- need - 1L
      }
    }
  }
  data.frame(study = study, stratum = stratum, eligible = eligible,
             heldout = held, reason = reason, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Poolability: does a group's deviation repeat across its studies?
# ---------------------------------------------------------------------------

#' Deviation of every unit from the units of all other studies
#'
#' @param M Numeric matrix, features x units (e.g. genes x study units of
#'   log CPM).
#' @param study Study of each unit (column).
#' @return Matrix of the same shape: each column minus the mean of the
#'   columns that belong to other studies.
#' @export
dcs_leave_study_out_deviation <- function(M, study) {
  study <- as.character(study)
  stopifnot(ncol(M) == length(study), length(unique(study)) >= 2)
  S <- rowSums(M)
  St <- rowsum(t(M), study, reorder = FALSE)          # studies x features
  n_st <- as.vector(table(factor(study, levels = rownames(St))))
  idx <- match(study, rownames(St))
  others <- (S - t(St[idx, , drop = FALSE])) /
    rep(ncol(M) - n_st[idx], each = nrow(M))
  M - others
}

# mean leave-one-out replication from a covariance matrix of deviations
.dcs_lso_replication <- function(V) {
  k <- nrow(V)
  if (k < 2) return(NA_real_)
  tot <- sum(V)
  r <- vapply(seq_len(k), function(u) {
    cov_u <- sum(V[u, -u])
    var_o <- tot - 2 * sum(V[u, ]) + V[u, u]
    if (var_o <= 0 || V[u, u] <= 0) return(NA_real_)
    cov_u / sqrt(V[u, u] * var_o)
  }, numeric(1))
  mean(r, na.rm = TRUE)
}

#' Poolability of groups of study units
#'
#' For each group with at least `min_studies` units (one unit per study and
#' group), the statistic is the mean over the group's units of the
#' correlation between the unit's leave-study-out deviation and the mean
#' deviation of the group's other units (leave-study-out replication). The
#' null draws the same number of units at random with the same composition
#' of `unit_class` (e.g. protocol class) and distinct studies.
#'
#' @param M Numeric matrix, features x units.
#' @param group Group of each unit (NA = in the universe but no group).
#' @param study Study of each unit.
#' @param unit_class Optional class of each unit; the null matches the
#'   group's class counts.
#' @param n_perm Null draws per group.
#' @param min_studies Smallest group size tested.
#' @param seed Random seed.
#' @return data.frame: group, studies, replication, mean_pair_r,
#'   positive_share (units whose replication r > 0), null_mean, null_q95,
#'   z (against the null), p_perm.
#' @export
dcs_poolability <- function(M, group, study, unit_class = NULL,
                            n_perm = 2000L, min_studies = 3L, seed = 1L) {
  group <- as.character(group)
  study <- as.character(study)
  if (is.null(unit_class)) unit_class <- rep("all", ncol(M))
  unit_class <- as.character(unit_class)
  stopifnot(ncol(M) == length(group), ncol(M) == length(study),
            ncol(M) == length(unit_class))
  D <- dcs_leave_study_out_deviation(M, study)
  V <- stats::cov(D)
  sdv <- sqrt(diag(V))
  C <- V / outer(sdv, sdv)
  set.seed(seed)
  per_unit <- function(ix) {
    k <- length(ix); tot <- sum(V[ix, ix])
    vapply(seq_len(k), function(j) {
      u <- ix[j]; o <- ix[-j]
      var_o <- sum(V[o, o])
      if (var_o <= 0) return(NA_real_)
      sum(V[u, o]) / sqrt(V[u, u] * var_o)
    }, numeric(1))
  }
  draw_null <- function(cls_n) {
    pick <- integer(); used <- character()
    for (cl in names(cls_n)) {
      cand <- which(unit_class == cl & !study %in% used)
      cand <- cand[sample.int(length(cand))]
      cand <- cand[!duplicated(study[cand])]
      if (length(cand) < cls_n[[cl]]) return(NULL)
      take <- cand[seq_len(cls_n[[cl]])]
      pick <- c(pick, take); used <- c(used, study[take])
    }
    pick
  }
  groups <- sort(unique(group[!is.na(group)]))
  rows <- lapply(groups, function(g) {
    ix <- which(!is.na(group) & group == g)
    ix <- ix[!duplicated(study[ix])]
    if (length(ix) < min_studies) return(NULL)
    r_u <- per_unit(ix)
    obs <- mean(r_u, na.rm = TRUE)
    cc <- C[ix, ix]
    cls_n <- table(unit_class[ix])
    null <- vapply(seq_len(n_perm), function(b) {
      p <- draw_null(cls_n)
      if (is.null(p)) NA_real_ else .dcs_lso_replication(V[p, p, drop = FALSE])
    }, numeric(1))
    null <- null[is.finite(null)]
    data.frame(group = g, studies = length(ix), replication = obs,
               mean_pair_r = mean(cc[upper.tri(cc)]),
               positive_share = mean(r_u > 0, na.rm = TRUE),
               null_mean = if (length(null)) mean(null) else NA_real_,
               null_q95 = if (length(null)) unname(stats::quantile(null, 0.95)) else NA_real_,
               z = if (length(null) > 1) (obs - mean(null)) / stats::sd(null) else NA_real_,
               p_perm = (sum(null >= obs) + 1) / (length(null) + 1),
               n_null = length(null), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Placebo (control against control) contrasts
# ---------------------------------------------------------------------------

#' Gene-centred window contrast between two arms
#'
#' d = log(a / A_gene) - log(b / B_gene) for windows with at least
#' `min_reads` in both arms, where A_gene and B_gene are the arms' totals over
#' the gene's windows (all windows, before filtering). `v` is the Poisson
#' count variance 1/a + 1/b.
#' @param a,b Window counts of the two arms (same window order).
#' @param gene Gene (or transcript) of each window.
#' @param min_reads Smallest count per window and arm.
#' @return data.frame with window index `i`, gene, d, v, a, b.
#' @export
dcs_window_contrast <- function(a, b, gene, min_reads = 200) {
  stopifnot(length(a) == length(b), length(a) == length(gene))
  A <- ave(as.numeric(a), gene, FUN = sum)
  Bt <- ave(as.numeric(b), gene, FUN = sum)
  keep <- which(a >= min_reads & b >= min_reads)
  data.frame(i = keep, gene = gene[keep],
             d = log(a[keep] / A[keep]) - log(b[keep] / Bt[keep]),
             v = 1 / a[keep] + 1 / b[keep], a = a[keep], b = b[keep],
             stringsAsFactors = FALSE)
}

#' Gene-level contrast between two arms
#'
#' Level: d = log share of the gene in arm a minus in arm b, median-centred.
#' With `a_den` / `b_den` (e.g. CDS counts) the contrast is an allocation
#' ratio instead: d = log(a / a_den) - log(b / b_den), median-centred.
#' @param a,b Gene counts of the two arms.
#' @param a_den,b_den Optional denominators per gene.
#' @param min_reads Smallest count per gene, arm (and denominator).
#' @return data.frame with i, d, v.
#' @export
dcs_gene_contrast <- function(a, b, a_den = NULL, b_den = NULL,
                              min_reads = 1000) {
  stopifnot(length(a) == length(b))
  if (is.null(a_den)) {
    keep <- which(a >= min_reads & b >= min_reads)
    d <- log(a[keep] / sum(a)) - log(b[keep] / sum(b))
    v <- 1 / a[keep] + 1 / b[keep]
  } else {
    keep <- which(a >= min_reads & b >= min_reads &
                    a_den >= min_reads & b_den >= min_reads)
    d <- log(a[keep] / a_den[keep]) - log(b[keep] / b_den[keep])
    v <- 1 / a[keep] + 1 / b[keep] + 1 / a_den[keep] + 1 / b_den[keep]
  }
  data.frame(i = keep, d = d - stats::median(d), v = v)
}

#' Noise band of a set of placebo contrasts
#'
#' Splits the observed spread of placebo contrasts into the Poisson count
#' part and an excess (between-arm) part. The observed variance is robust
#' (MAD^2) by default; the count part is the mean of `v`.
#' @param d Contrasts.
#' @param v Their count variances.
#' @param robust Use MAD instead of SD.
#' @param level Two-sided coverage of the band quantile.
#' @return Named numeric: n, sd_obs, sd_count, sd_excess (sqrt of the
#'   variance above counting noise, 0 if none), band (quantile of |d| at
#'   `level`), false_rate_count (share of |d| / sqrt(v) above 1.96: the
#'   share a count-only test would call).
#' @export
dcs_placebo_band <- function(d, v, robust = TRUE, level = 0.95) {
  ok <- is.finite(d) & is.finite(v)
  d <- d[ok]; v <- v[ok]
  if (length(d) < 3) {
    return(c(n = as.numeric(length(d)), sd_obs = NA_real_, sd_count = NA_real_,
             sd_excess = NA_real_, band = NA_real_, false_rate_count = NA_real_))
  }
  sd_obs <- if (robust) stats::mad(d) else stats::sd(d)
  var_c <- mean(v)
  c(n = as.numeric(length(d)), sd_obs = sd_obs, sd_count = sqrt(var_c),
    sd_excess = sqrt(max(0, sd_obs^2 - var_c)),
    band = unname(stats::quantile(abs(d - stats::median(d)), level)),
    false_rate_count = mean(abs(d - stats::median(d)) / sqrt(v) > 1.96))
}

# ---------------------------------------------------------------------------
# Within-gene motif amplitude
# ---------------------------------------------------------------------------

#' Within-gene log-linear amplitude of sequence features
#'
#' Poisson model y_i ~ mu_g(i) exp(X_i beta) with one free rate per gene,
#' profiled out exactly (multinomial within gene). Fisher scoring on beta;
#' standard errors are gene-clustered sandwich estimates, so neighbouring
#' positions of one gene do not count as independent.
#'
#' @param y Counts per position.
#' @param gene Gene of each position.
#' @param X Numeric (or sparse) design matrix, positions x features, no
#'   intercept. Features absent in every kept gene are dropped.
#' @param maxit,tol Fisher scoring controls.
#' @param ridge Small ridge added to the information (numerical safety).
#' @return data.frame: feature, beta, se (sandwich), se_model, n_pos (positions
#'   with the feature), plus attributes iterations and converged.
#' @export
dcs_motif_amplitude <- function(y, gene, X, maxit = 50L, tol = 1e-8,
                                ridge = 1e-8) {
  X <- as.matrix(X)
  stopifnot(length(y) == nrow(X), length(gene) == nrow(X))
  g <- as.integer(factor(gene))
  yt <- rowsum(y, g, reorder = TRUE)[, 1]
  keep_g <- yt > 0
  keep <- keep_g[g]
  y <- y[keep]; X <- X[keep, , drop = FALSE]; g <- as.integer(factor(g[keep]))
  n_pos <- colSums(X != 0)
  X <- X[, n_pos > 0, drop = FALSE]
  n_pos <- n_pos[n_pos > 0]
  yt <- rowsum(y, g, reorder = TRUE)[, 1]
  beta <- rep(0, ncol(X)); converged <- FALSE
  fit_mu <- function(beta) {
    eta <- drop(X %*% beta)
    w <- exp(eta - stats::ave(eta, g, FUN = max))
    w / rowsum(w, g, reorder = TRUE)[g, 1] * yt[g]
  }
  for (it in seq_len(maxit)) {
    mu <- fit_mu(beta)
    score <- drop(crossprod(X, y - mu))
    Xbar <- rowsum(X * mu, g, reorder = TRUE) / yt          # genes x p
    info <- crossprod(X * sqrt(mu)) - crossprod(Xbar * sqrt(yt))
    diag(info) <- diag(info) + ridge
    step <- drop(solve(info, score))
    beta <- beta + step
    if (max(abs(step)) < tol) { converged <- TRUE; break }
  }
  mu <- fit_mu(beta)
  Xbar <- rowsum(X * mu, g, reorder = TRUE) / yt
  info <- crossprod(X * sqrt(mu)) - crossprod(Xbar * sqrt(yt))
  diag(info) <- diag(info) + ridge
  inv <- solve(info)
  U <- rowsum(X * (y - mu), g, reorder = TRUE)               # gene scores
  sand <- inv %*% crossprod(U) %*% inv
  out <- data.frame(feature = colnames(X) %||% paste0("x", seq_len(ncol(X))),
                    beta = beta, se = sqrt(pmax(diag(sand), 0)),
                    se_model = sqrt(pmax(diag(inv), 0)), n_pos = unname(n_pos),
                    stringsAsFactors = FALSE)
  attr(out, "iterations") <- it
  attr(out, "converged") <- converged
  out
}

#' Contrast of per-library estimates between two arms
#'
#' Difference of arm means of a per-library estimate (e.g. a motif
#' amplitude) with a Welch interval from the replicate spread. With a single
#' library in an arm, that arm's variance falls back to the mean of the
#' supplied per-library sampling variances (`se^2`), and the interval uses
#' a normal quantile.
#' @param est,se Per-library estimate and its standard error.
#' @param arm Arm of each library ("case" / "control" or any two levels).
#' @param case Level treated as the case arm.
#' @param level Interval coverage.
#' @return One-row data.frame: n_case, n_control, mean_case, mean_control,
#'   delta, se, lower, upper, df.
#' @export
dcs_arm_contrast <- function(est, se, arm, case = "case", level = 0.95) {
  arm <- as.character(arm)
  a <- arm == case
  b <- !a
  stopifnot(sum(a) >= 1, sum(b) >= 1)
  arm_var <- function(i) {
    if (sum(i) >= 2) stats::var(est[i]) else mean(se[i]^2)
  }
  va <- arm_var(a) / sum(a); vb <- arm_var(b) / sum(b)
  delta <- mean(est[a]) - mean(est[b])
  s <- sqrt(va + vb)
  df <- if (sum(a) >= 2 && sum(b) >= 2) {
    (va + vb)^2 / (va^2 / (sum(a) - 1) + vb^2 / (sum(b) - 1))
  } else Inf
  q <- if (is.finite(df)) stats::qt(1 - (1 - level) / 2, df) else
    stats::qnorm(1 - (1 - level) / 2)
  data.frame(n_case = sum(a), n_control = sum(b), mean_case = mean(est[a]),
             mean_control = mean(est[b]), delta = delta, se = s,
             lower = delta - q * s, upper = delta + q * s, df = df)
}

# ---------------------------------------------------------------------------
# Protocol multipliers and matched null sets
# ---------------------------------------------------------------------------

#' Protocol multiplier per feature, empirical-Bayes shrunk
#'
#' For a features x studies matrix of log profiles and a two-level study
#' label (e.g. no drug vs cycloheximide), the raw multiplier is the
#' difference of the two group means per feature; its sampling variance
#' comes from the between-study spread inside each group. The shrunk
#' multiplier is raw * tau2 / (tau2 + se^2), with tau2 the robust variance
#' (MAD^2) of the raw multipliers above their median sampling variance, so a
#' few very noisy features do not hide a shared signal.
#' @param P Numeric matrix, features x studies.
#' @param group Logical (TRUE = the level whose effect is estimated) per study.
#' @return list(raw, se, shrunk, tau2).
#' @export
dcs_protocol_multiplier <- function(P, group) {
  group <- as.logical(group)
  stopifnot(ncol(P) == length(group), sum(group) >= 2, sum(!group) >= 2)
  a <- P[, group, drop = FALSE]; b <- P[, !group, drop = FALSE]
  ma <- rowMeans(a); mb <- rowMeans(b)
  va <- rowSums((a - ma)^2) / (ncol(a) - 1); vb <- rowSums((b - mb)^2) / (ncol(b) - 1)
  raw <- ma - mb
  se2 <- va / ncol(a) + vb / ncol(b)
  tau2 <- max(0, stats::mad(raw)^2 - stats::median(se2))
  list(raw = raw, se = sqrt(se2), shrunk = raw * tau2 / (tau2 + se2), tau2 = tau2)
}

#' Matched random sets of features
#'
#' For each target, draws a non-target feature from the same stratum (e.g.
#' a decile x decile bin of expression and leader reads); repeated `n_sets`
#' times. Strata with no candidate fall back to all candidates.
#' @param target Indices (or names) of the target features within `stratum`.
#' @param stratum Stratum label per feature.
#' @param candidate Logical per feature: allowed in matched sets.
#' @param n_sets Number of sets.
#' @param seed Random seed.
#' @return Integer matrix, n_sets x length(target), of feature indices.
#' @export
dcs_matched_sets <- function(target, stratum, candidate, n_sets = 200L, seed = 1L) {
  if (is.character(target)) target <- match(target, names(stratum))
  set.seed(seed)
  cand <- which(candidate)
  by_s <- split(cand, stratum[cand])
  sapply(target, function(t) {
    pool <- by_s[[as.character(stratum[t])]]
    if (!length(pool)) pool <- cand
    pool[sample.int(length(pool), n_sets, replace = TRUE)]
  })
}
