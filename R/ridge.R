# Closed-form ridge regression compatible with glmnet(alpha = 0).
#
# The RDG residual models (hierarchical context, branch covariates, graph
# geometry, mechanistic motifs) fit a Gaussian ridge model with blocked
# (leave-one-study / leave-one-gene) cross-validation. With glmnet this costs
# one coordinate-descent path per fold inside cv.glmnet() plus one more fit per
# fold and candidate lambda for the held-out predictions. For a Gaussian ridge
# every one of those fits has a closed form: one eigendecomposition of the
# standardized training Gram matrix per fold gives the whole lambda path.
#
# The functions below reproduce glmnet's conventions so that lambda selection
# and predictions agree with glmnet to its own convergence tolerance:
#   * columns standardized with the 1/n standard deviation, constant training
#     columns get coefficient 0;
#   * penalty (lambda / sd(y)) / 2 * ||beta_std||^2 on the standardized scale;
#   * lambda path: lambda_max = max_j |cov(x_j, y)| / sd(x_j) / 1e-3, a
#     geometric sequence down to lambda.min.ratio (0.01 if n < p, else 1e-4),
#     all-zero coefficients at the first lambda, and glmnet's early stop
#     (after 5 lambdas, stop when the deviance ratio gains < 1e-5 relative or
#     exceeds 0.999);
#   * cv.glmnet: each fold computes its own path and predicts at the full-data
#     lambdas by linear interpolation between path points (alignment
#     "lambda"); grouped MSE; lambda.min / lambda.1se as getOptcv.glmnet.
#
# Pure R + LAPACK so that pipeline scripts can source() this file.

.dcs_ridge_ctl <- list(fdev = 1e-5, devmax = 0.999, mnlam = 5L,
                       alpha_floor = 1e-3, const_tol = 1e-10)

# Sufficient statistics of a dense, globally centred design.
dcs_ridge_prepare <- function(x, y) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.numeric(y)
  if (nrow(x) != length(y)) stop("nrow(x) must equal length(y).", call. = FALSE)
  if (anyNA(x) || anyNA(y)) stop("x and y must not contain NA.", call. = FALSE)
  center <- colMeans(x)
  y_center <- mean(y)
  xc <- sweep(x, 2L, center, check.margin = FALSE)
  yc <- y - y_center
  list(
    xc = xc, yc = yc, center = center, y_center = y_center,
    n = nrow(x), p = ncol(x), colnames = colnames(x),
    gram = crossprod(xc), xty = drop(crossprod(xc, yc))
  )
}

# Eigen-decomposed standardized training problem for a set of training rows.
dcs_ridge_train <- function(prep, train_rows) {
  nt <- length(train_rows)
  test_rows <- setdiff(seq_len(prep$n), train_rows)
  if (length(test_rows) < length(train_rows) && length(test_rows)) {
    xb <- prep$xc[test_rows, , drop = FALSE]
    yb <- prep$yc[test_rows]
    gram <- prep$gram - crossprod(xb)
    xty <- prep$xty - drop(crossprod(xb, yb))
    sx <- -colSums(xb)
    sy <- -sum(yb)
  } else {
    xt <- prep$xc[train_rows, , drop = FALSE]
    yt <- prep$yc[train_rows]
    gram <- crossprod(xt)
    xty <- drop(crossprod(xt, yt))
    sx <- colSums(xt)
    sy <- sum(yt)
  }
  mx <- sx / nt
  my <- sy / nt
  cov_xx <- gram / nt - tcrossprod(mx)
  cov_xy <- xty / nt - mx * my
  var_x <- pmax(diag(cov_xx), 0)
  y_t <- prep$yc[train_rows]
  y_sd <- sqrt(mean((y_t - my)^2))
  # A column constant on the training rows has a (numerically) zero variance
  # relative to its mean square; glmnet excludes such columns.
  keep <- var_x > .dcs_ridge_ctl$const_tol * pmax(diag(gram) / nt, .Machine$double.xmin)
  xs <- sqrt(var_x[keep])
  r <- cov_xx[keep, keep, drop = FALSE] / tcrossprod(xs)
  b <- cov_xy[keep] / xs / y_sd
  eig <- eigen(r, symmetric = TRUE)
  d <- pmax(eig$values, 0)
  list(
    n = nt, p = prep$p, keep = keep, xs = xs, mx = mx, my = my,
    y_sd = y_sd, vectors = eig$vectors, d = d,
    z = drop(crossprod(eig$vectors, b)), b = b
  )
}

# Standardized coefficients (columns: lambdas) and deviance ratios.
dcs_ridge_solve <- function(tr, lambda) {
  pen <- lambda / tr$y_sd
  denom <- outer(tr$d, pen, "+")
  zz <- tr$z / denom
  beta_std <- tr$vectors %*% zz
  rsq <- colSums(2 * tr$z * zz) - colSums(tr$d * zz^2)
  list(beta_std = beta_std, rsq = rsq)
}

dcs_ridge_lambda_max <- function(tr) {
  if (!any(tr$keep)) return(NA_real_)
  max(abs(tr$b)) * tr$y_sd / .dcs_ridge_ctl$alpha_floor
}

# glmnet's default lambda path for one training set (with early stopping).
dcs_ridge_path <- function(tr, nlambda = 100L, lambda_min_ratio = NULL) {
  if (is.null(lambda_min_ratio)) {
    lambda_min_ratio <- if (tr$n < tr$p) 0.01 else 1e-4
  }
  lmax <- dcs_ridge_lambda_max(tr)
  # glmnet multiplies successively (alm <- alm * alf) and reports the first
  # lambda extrapolated on the log scale from the second and third.
  alf <- lambda_min_ratio^(1 / (nlambda - 1))
  lambda <- numeric(nlambda)
  lambda[1L] <- lmax
  if (nlambda >= 2L) lambda[2L] <- lmax * alf
  for (m in seq_len(nlambda)[-(1:2)]) lambda[m] <- lambda[m - 1L] * alf
  sol <- dcs_ridge_solve(tr, lambda)
  sol$beta_std[, 1L] <- 0
  sol$rsq[[1L]] <- 0
  ctl <- .dcs_ridge_ctl
  lmu <- nlambda
  for (m in seq_len(nlambda)) {
    if (m < ctl$mnlam) next
    if (sol$rsq[[m]] - sol$rsq[[m - 1L]] < ctl$fdev * sol$rsq[[m]] ||
        sol$rsq[[m]] > ctl$devmax) {
      lmu <- m
      break
    }
  }
  idx <- seq_len(lmu)
  lambda <- lambda[idx]
  if (lmu >= 3L) lambda[1L] <- exp(2 * log(lambda[2L]) - log(lambda[3L]))
  list(lambda = lambda, beta_std = sol$beta_std[, idx, drop = FALSE],
       rsq = sol$rsq[idx])
}

# Base R >= 4.4 provides %||%; define it for older R.
if (!exists("%||%", envir = baseenv(), inherits = FALSE)) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# Predictions for rows `rows` from standardized coefficient columns.
dcs_ridge_predict_rows <- function(prep, tr, beta_std, rows) {
  beta_std <- as.matrix(beta_std)
  beta <- beta_std * tr$y_sd / tr$xs
  xt <- sweep(prep$xc[rows, tr$keep, drop = FALSE], 2L, tr$mx[tr$keep],
              check.margin = FALSE)
  pred <- xt %*% beta
  sweep(pred, 2L, rep(prep$y_center + tr$my, ncol(beta)), "+",
        check.margin = FALSE)
}

#' Closed-form ridge coefficients
#'
#' Coefficients (intercept first, one column per lambda) of
#' `glmnet(x, y, alpha = 0, lambda = lambda)` computed exactly.
#'
#' @param x,y Design matrix (without intercept) and response.
#' @param lambda Penalties on glmnet's scale.
#' @param prep Optional result of `dcs_ridge_prepare()` to reuse.
#' @export
dcs_ridge_fit <- function(x, y, lambda, prep = NULL) {
  if (is.null(prep)) prep <- dcs_ridge_prepare(x, y)
  tr <- dcs_ridge_train(prep, seq_len(prep$n))
  sol <- dcs_ridge_solve(tr, lambda)
  beta <- matrix(0, prep$p, length(lambda))
  beta[tr$keep, ] <- sol$beta_std * tr$y_sd / tr$xs
  x_mean <- prep$center + tr$mx
  intercept <- prep$y_center + tr$my - drop(crossprod(x_mean, beta))
  out <- rbind(intercept, beta)
  rownames(out) <- c("(Intercept)", prep$colnames %||% paste0("V", seq_len(prep$p)))
  colnames(out) <- paste0("s", seq_along(lambda) - 1L)
  out
}

.dcs_lambda_interp <- function(lambda, s) {
  if (length(lambda) == 1L) {
    n <- length(s)
    return(list(left = rep(1L, n), right = rep(1L, n), frac = rep(1, n)))
  }
  k <- length(lambda)
  sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
  lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
  sfrac[sfrac < min(lambda)] <- min(lambda)
  sfrac[sfrac > max(lambda)] <- max(lambda)
  coord <- stats::approx(lambda, seq_along(lambda), sfrac)$y
  left <- floor(coord)
  right <- ceiling(coord)
  sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
  sfrac[left == right] <- 1
  sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
  list(left = left, right = right, frac = sfrac)
}

# Many small eigendecompositions are much slower with a multi-threaded BLAS
# on a shared server (measured 9x for p = 464: thread start-up and
# contention); run them single-threaded and restore the setting afterwards.
.dcs_single_thread_blas <- function() {
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) return(function() NULL)
  old <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1L)
  function() RhpcBLASctl::blas_set_num_threads(old)
}

#' Blocked cross-validated ridge regression (glmnet-compatible, closed form)
#'
#' Equivalent to `glmnet::cv.glmnet(x, y, alpha = 0, foldid = foldid,
#' nlambda = nlambda, standardize = TRUE, intercept = TRUE)` followed by
#' per-fold refits at selected lambdas, computed from one eigendecomposition
#' per fold.
#'
#' @param x Numeric (dense or sparse) design matrix without intercept column.
#' @param y Numeric response.
#' @param foldid Integer fold labels 1..K.
#' @param nlambda Length of the lambda path.
#' @param min_train Folds with fewer training rows get `NA` from
#'   `fold_prediction()` (the scripts' historical `length(train_idx) < 100`
#'   rule); the CV curve itself uses every fold, as cv.glmnet does.
#' @param lambda_min_ratio Smallest lambda as a fraction of lambda_max
#'   (NULL: glmnet's default).
#' @return List with `lambda`, `cvm`, `cvsd`, `lambda.min`, `lambda.1se`,
#'   `fold_prediction(lambda_values)` (held-out predictions with an exact
#'   refit per fold, as `glmnet(lambda = s)` per fold) and `fit(lambda_values)`
#'   (full-data coefficients).
#' @export
dcs_cv_ridge <- function(x, y, foldid, nlambda = 100L, min_train = 0L,
                         lambda_min_ratio = NULL) {
  restore <- .dcs_single_thread_blas()
  on.exit(restore(), add = TRUE)
  prep <- dcs_ridge_prepare(x, y)
  foldid <- as.integer(foldid)
  nfolds <- max(foldid)
  full <- dcs_ridge_train(prep, seq_len(prep$n))
  full_path <- dcs_ridge_path(full, nlambda = nlambda,
                             lambda_min_ratio = lambda_min_ratio)
  lambda <- full_path$lambda
  folds <- lapply(seq_len(nfolds), function(k) {
    test <- which(foldid == k)
    train <- which(foldid != k)
    tr <- dcs_ridge_train(prep, train)
    list(test = test, train = train, tr = tr)
  })
  predmat <- matrix(NA_real_, prep$n, length(lambda))
  for (k in seq_len(nfolds)) {
    f <- folds[[k]]
    if (!length(f$test)) next
    path <- dcs_ridge_path(f$tr, nlambda = nlambda,
                           lambda_min_ratio = lambda_min_ratio)
    pred_path <- dcs_ridge_predict_rows(prep, f$tr, path$beta_std, f$test)
    ip <- .dcs_lambda_interp(path$lambda, lambda)
    predmat[f$test, ] <- sweep(pred_path[, ip$left, drop = FALSE], 2L,
                               ip$frac, "*") +
      sweep(pred_path[, ip$right, drop = FALSE], 2L, 1 - ip$frac, "*")
  }
  y_num <- as.numeric(y)
  cvraw <- (y_num - predmat)^2
  fold_sizes <- tabulate(foldid, nfolds)
  fold_mean <- t(vapply(seq_len(nfolds), function(k) {
    colMeans(cvraw[foldid == k, , drop = FALSE])
  }, numeric(length(lambda))))
  if (length(lambda) == 1L) fold_mean <- matrix(fold_mean, ncol = 1L)
  ok <- fold_sizes > 0
  w <- fold_sizes[ok]
  fm <- fold_mean[ok, , drop = FALSE]
  cvm <- colSums(fm * w) / sum(w)
  cvsd <- sqrt(colSums(sweep(fm, 2L, cvm)^2 * w) / sum(w) / (sum(ok) - 1))
  cvmin <- min(cvm, na.rm = TRUE)
  lambda_min <- max(lambda[cvm <= cvmin], na.rm = TRUE)
  idmin <- match(lambda_min, lambda)
  lambda_1se <- max(lambda[cvm <= (cvm + cvsd)[idmin]], na.rm = TRUE)

  fold_prediction <- function(lambda_values) {
    restore <- .dcs_single_thread_blas()
    on.exit(restore(), add = TRUE)
    out <- matrix(NA_real_, prep$n, length(lambda_values))
    for (f in folds) {
      if (!length(f$test) || length(f$train) < min_train) next
      sol <- dcs_ridge_solve(f$tr, lambda_values)
      out[f$test, ] <- dcs_ridge_predict_rows(prep, f$tr, sol$beta_std, f$test)
    }
    out
  }
  list(
    lambda = lambda, cvm = cvm, cvsd = cvsd,
    lambda.min = lambda_min, lambda.1se = lambda_1se,
    fold_prediction = fold_prediction,
    fit = function(lambda_values) dcs_ridge_fit(prep = prep, lambda = lambda_values),
    n_folds = nfolds
  )
}

#' Blocked ridge model shared by the RDG residual scripts
#'
#' Fits a Gaussian ridge (alpha = 0) with blocked cross-validation for lambda
#' and returns closures for held-out predictions and full-data coefficients.
#' `engine = "fast"` uses the closed-form path above, `"glmnet"` the original
#' cv.glmnet + per-block glmnet refits. `"auto"` (default, or the
#' `DOMINANT_RIDGE_ENGINE` environment variable) uses the closed form up to
#' `max_fast_features` columns, where one eigendecomposition per block is
#' cheap; wider sparse designs keep glmnet.
#'
#' @param x Design matrix (sparse or dense) without intercept.
#' @param y Numeric response.
#' @param blocks Block label per row; one fold per sorted unique label.
#' @param nlambda Lambda path length.
#' @param min_train Blocks with fewer training rows are not predicted.
#' @param lambda_min_ratio Smallest lambda as a fraction of lambda_max
#'   (default as glmnet: 1e-4, or 0.01 when n < p).
#' @param extend_grid When the CV optimum is the smallest lambda of the path,
#'   divide `lambda_min_ratio` by 100 and refit until the optimum is interior
#'   or the ratio would fall below `min_lambda_min_ratio`.
#' @param min_lambda_min_ratio Lower limit for the extension.
#' @return List with `engine`, `lambda.min`, `lambda.1se`, `n_blocks`,
#'   `lambda_min_ratio` (used), `lambda_min_interior`,
#'   `predict(lambda)` returning list(prediction, train_rows) with `NA` for
#'   skipped blocks, and `coef(lambda)` returning a one-column matrix with
#'   the intercept first.
#' @export
dcs_blocked_ridge <- function(x, y, blocks, nlambda = 80L, min_train = 100L,
                              engine = Sys.getenv("DOMINANT_RIDGE_ENGINE", "auto"),
                              max_fast_features = 2000L,
                              lambda_min_ratio = NULL, extend_grid = TRUE,
                              min_lambda_min_ratio = 1e-10) {
  engine <- match.arg(engine, c("auto", "fast", "glmnet"))
  if (identical(engine, "auto")) {
    engine <- if (ncol(x) <= max_fast_features) "fast" else "glmnet"
  }
  blocks <- as.character(blocks)
  levels_b <- sort(unique(blocks))
  foldid <- match(blocks, levels_b)
  n <- length(y)
  feature_names <- colnames(x)
  default_ratio <- if (n < ncol(x)) 0.01 else 1e-4
  ratio <- lambda_min_ratio %||% default_ratio

  cv_fit <- function(ratio) {
    if (identical(engine, "fast")) {
      dcs_cv_ridge(x, y, foldid, nlambda = nlambda, min_train = min_train,
                   lambda_min_ratio = ratio)
    } else {
      glmnet::cv.glmnet(x, y, family = "gaussian", alpha = 0,
                        foldid = foldid, standardize = TRUE,
                        intercept = TRUE, nlambda = nlambda,
                        lambda.min.ratio = ratio)
    }
  }
  # A CV optimum on the smallest lambda of the path means the grid, not the
  # data, chose the penalty: widen the path 100-fold until it is interior.
  cv <- cv_fit(ratio)
  while (isTRUE(extend_grid) && cv$lambda.min <= min(cv$lambda) &&
         ratio / 100 >= min_lambda_min_ratio) {
    ratio <- ratio / 100
    cv <- cv_fit(ratio)
  }
  grid <- list(lambda_min_ratio = ratio,
               lambda_min_interior = cv$lambda.min > min(cv$lambda))

  if (identical(engine, "fast")) {
    train_n <- n - tabulate(foldid, length(levels_b))
    predict_fun <- function(lambda) {
      pred <- drop(cv$fold_prediction(lambda))
      tr <- train_n[foldid]
      tr[tr < min_train] <- NA_integer_
      list(prediction = pred, train_rows = as.integer(tr))
    }
    coef_fun <- function(lambda) {
      out <- cv$fit(lambda)
      rownames(out) <- c("(Intercept)", feature_names %||% rownames(out)[-1])
      out
    }
  } else {
    predict_fun <- function(lambda) {
      pred <- rep(NA_real_, n)
      train_rows <- rep(NA_integer_, n)
      for (block_id in levels_b) {
        test_idx <- which(blocks == block_id)
        train_idx <- which(blocks != block_id)
        if (!length(test_idx) || length(train_idx) < min_train) next
        fit <- glmnet::glmnet(x[train_idx, , drop = FALSE], y[train_idx],
                              family = "gaussian", alpha = 0, lambda = lambda,
                              standardize = TRUE, intercept = TRUE)
        pred[test_idx] <- as.numeric(stats::predict(
          fit, newx = x[test_idx, , drop = FALSE], s = lambda))
        train_rows[test_idx] <- length(train_idx)
      }
      list(prediction = pred, train_rows = train_rows)
    }
    coef_fun <- function(lambda) {
      fit <- glmnet::glmnet(x, y, family = "gaussian", alpha = 0,
                            lambda = lambda, standardize = TRUE,
                            intercept = TRUE)
      as.matrix(stats::coef(fit, s = lambda))
    }
  }
  c(list(engine = engine, lambda.min = cv$lambda.min,
         lambda.1se = cv$lambda.1se, lambda = cv$lambda, cvm = cv$cvm,
         n_blocks = length(levels_b), predict = predict_fun, coef = coef_fun),
    grid)
}
