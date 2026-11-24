# Input checking, starting values and construction of the fitted object.

.is_prepare <- function(X, y, family, alpha, intercept, weights, offset, unpenalized, control) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("'X' must be a numeric matrix", call. = FALSE)
  if (intercept) X <- cbind(1, X)
  nobs <- nrow(X)
  nvars <- ncol(X)
  if (nvars == 0L) stop("the model matrix is empty", call. = FALSE)
  if (anyNA(X) || any(!is.finite(X))) stop("'X' contains missing or infinite values", call. = FALSE)
  nms <- colnames(X)
  if (is.null(nms)) nms <- paste0("X", seq_len(nvars) - intercept)
  if (intercept) nms[1L] <- "(Intercept)"
  colnames(X) <- nms
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1)
    stop("'alpha' must be a number in [0, 1] (0 = ridge, 1 = lasso)", call. = FALSE)

  # unpenalized: logical vector over the columns of X (intercept always unpenalized)
  unp <- logical(nvars)
  if (intercept) unp[1L] <- TRUE
  if (!is.null(unpenalized) && length(unpenalized)) {
    if (is.list(unpenalized) || is.factor(unpenalized))
      stop("'unpenalized' must be a vector of names, indices or a logical vector", call. = FALSE)
    cand <- if (intercept) nms[-1L] else nms
    if (is.logical(unpenalized)) {
      if (length(unpenalized) != length(cand))
        stop("a logical 'unpenalized' must have one entry per column of 'X'", call. = FALSE)
      id <- which(unpenalized)
    } else if (is.character(unpenalized)) {
      id <- match(unpenalized, cand)
      if (anyNA(id)) id[is.na(id)] <- pmatch(unpenalized[is.na(id)], cand)
      if (anyNA(id))
        stop(sprintf("the following names are not columns of the model matrix: %s",
                     paste(unpenalized[is.na(id)], collapse = ", ")), call. = FALSE)
    } else {
      id <- unpenalized
      if (any(abs(id - round(id)) > sqrt(.Machine$double.eps)))
        stop("some elements of 'unpenalized' are not integers", call. = FALSE)
      if (any(id < 1 | id > length(cand)))
        stop(sprintf("the elements of 'unpenalized' must be between 1 and %d", length(cand)), call. = FALSE)
    }
    unp[id + intercept] <- TRUE
  }

  if (is.null(offset)) offset <- numeric(nobs)
  offset <- as.double(offset)
  if (length(offset) != nobs)
    stop(sprintf("number of offsets is %d, should equal %d (number of observations)", length(offset), nobs),
         call. = FALSE)
  if (is.null(weights)) {
    weights <- rep(1, nobs)
  } else {
    if (!is.numeric(weights) || length(weights) != nobs)
      stop(sprintf("'weights' must be a numeric vector of length %d", nobs), call. = FALSE)
    weights[is.nan(weights)] <- 0
    if (any(weights < 0)) stop("negative weights are not allowed", call. = FALSE)
    if (all(weights == 0)) stop("all the weights are equal to zero", call. = FALSE)
  }

  family <- .is_family(family)
  fam <- family$family

  setting <- control
  if (setting$nfolds < 3) stop("'nfolds' must be at least 3", call. = FALSE)
  if (setting$tol <= 0) stop("'tol' must be positive", call. = FALSE)
  if (setting$itmax < 1) stop("'itmax' must be a positive integer", call. = FALSE)
  setting$estpai <- FALSE
  if (any(setting$c > 1)) stop("the mixture weight 'c' must be in [0, 1], or negative to be estimated",
                               call. = FALSE)
  if (setting$c[1L] < 0) {
    setting$estpai <- TRUE
    setting$c <- 0.5
  }
  setting$c <- rep_len(setting$c, nvars)

  # response
  n <- rep.int(1, nobs)
  if (is.factor(y) || is.character(y) || NCOL(y) == 2L) {
    if (!fam %in% c("binomial", "quasibinomial"))
      stop(sprintf("the %s family does not accept a factor or a two-column response", sQuote(fam)),
           call. = FALSE)
  }
  if (fam %in% c("binomial", "quasibinomial")) {
    if (NCOL(y) == 1L) {
      if (is.character(y)) y <- factor(y)
      if (is.factor(y)) y <- as.numeric(y != levels(y)[1L])
      y <- as.double(y)
      y[weights == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1", call. = FALSE)
      if (fam == "binomial" && any(abs(weights * y - round(weights * y)) > 0.001))
        warning("non-integer #successes in a binomial glm!", call. = FALSE)
    } else if (NCOL(y) == 2L) {
      if (fam == "binomial" && any(abs(y - round(y)) > 0.001))
        warning("non-integer counts in a binomial glm!", call. = FALSE)
      n <- y[, 1L] + y[, 2L]
      y <- ifelse(n == 0, 0, y[, 1L] / n)
      weights <- weights * n
    } else {
      stop("for the binomial family, y must be a vector of 0 and 1's or a two-column matrix of successes and failures",
           call. = FALSE)
    }
  } else {
    y <- as.double(drop(y))
  }
  if (anyNA(y)) stop("the response contains missing values", call. = FALSE)
  if (length(y) != nobs) stop("incompatible dimensions of 'X' and 'y'", call. = FALSE)
  if (fam %in% c("poisson", "quasipoisson") && any(y < 0))
    stop("negative values not allowed for the Poisson family", call. = FALSE)
  if (fam == "Gamma" && any(y <= 0))
    stop("non-positive values not allowed for the Gamma family", call. = FALSE)
  if (fam == "inverse.gaussian" && any(y <= 0))
    stop("non-positive values not allowed for the inverse.gaussian family", call. = FALSE)

  list(X = X, y = y, family = family, alpha = alpha, intercept = intercept, weights = weights,
       offset = offset, unpenalized = unp, nobs = nobs, nvars = nvars, n = n, nms = nms,
       setting = setting)
}

# run code with a temporary seed, restoring the RNG state afterwards
.is_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  force(expr)
}

# family argument for glmnet
.is_glmnet_family <- function(family) {
  f <- family$family
  if (f %in% c("gaussian", "binomial", "poisson") && identical(family$link, stats::make.link(
    switch(f, gaussian = "identity", binomial = "logit", poisson = "log"))$name)) return(f)
  if (f == "quasibinomial" && family$link == "logit") return("binomial")
  if (f == "quasipoisson" && family$link == "log") return("poisson")
  family
}

# Starting values (standardised scale) and, if 'lambda' is missing, its selection by
# cross-validation with glmnet.  lambda (islasso scale) = lambda (glmnet scale) * sum(weights).
.is_start <- function(prep, Xs, lambda) {
  X <- Xs
  y <- prep$y
  w <- prep$weights
  family <- prep$family
  intercept <- prep$intercept
  setting <- prep$setting
  p <- ncol(X)
  sw <- sum(w)
  pf <- as.numeric(!prep$unpenalized)
  interval <- NULL
  if (!is.null(setting$b0)) {
    if (missing(lambda) || is.null(lambda)) stop("a value for 'lambda' must be given with 'b0'", call. = FALSE)
    b0 <- as.double(setting$b0)
    if (length(b0) != p) stop(sprintf("'b0' must have length %d", p), call. = FALSE)
    # b0 is on the original scale: map to the standardised one
    A <- .is_backmap(prep$xm, prep$xs, intercept)
    est <- drop(solve(A, b0))
    return(list(beta = est, lambda = lambda, interval = NULL))
  }
  xx <- if (intercept) X[, -1L, drop = FALSE] else X
  pfx <- if (intercept) pf[-1L] else pf
  if (ncol(xx) < 2L) {
    if (missing(lambda) || is.null(lambda))
      stop("with less than two covariates a value for 'lambda' must be given", call. = FALSE)
    est <- rep(0, p)
    if (intercept) est[1L] <- family$linkfun(min(max(sum(w * y) / sw, 1e-4), if (family$family %in% c("binomial", "quasibinomial")) 1 - 1e-4 else Inf))
    return(list(beta = est, lambda = lambda, interval = NULL))
  }
  gfam <- .is_glmnet_family(family)
  if (is.character(gfam) && gfam == "binomial") {
    y2 <- cbind(1 - y, y) * w
    w2 <- rep(1, length(y))
  } else {
    y2 <- y
    w2 <- w
  }
  if (sum(pfx) == 0) pfx[] <- 1     # glmnet needs at least one penalised column
  args <- list(x = xx, y = y2, family = gfam, alpha = prep$alpha,
               weights = w2, offset = prep$offset, standardize = FALSE, intercept = intercept,
               penalty.factor = pfx)
  fit_est <- function(obj, s) {
    est <- as.vector(stats::coef(obj, s = s))
    if (!intercept) est <- est[-1L]
    est
  }
  out <- tryCatch({
    if (missing(lambda) || is.null(lambda)) {
      foldid <- .is_with_seed(setting$seed, sample(rep_len(seq_len(setting$nfolds), length(y))))
      obj <- suppressWarnings(do.call(glmnet::cv.glmnet, c(args, list(foldid = foldid))))
      list(beta = fit_est(obj, "lambda.min"), lambda = obj$lambda.min * sw,
           interval = range(obj$lambda) * sw)
    } else {
      obj <- suppressWarnings(do.call(glmnet::glmnet, args))
      s <- min(max(lambda / sw, min(obj$lambda)), max(obj$lambda))
      list(beta = fit_est(obj, s), lambda = lambda, interval = range(obj$lambda) * sw)
    }
  }, error = function(e) {
    if (missing(lambda) || is.null(lambda))
      stop("cross-validation with glmnet failed (", conditionMessage(e), "); please supply 'lambda'",
           call. = FALSE)
    NULL
  })
  if (is.null(out)) {
    est <- rep(0, p)
    if (intercept) est[1L] <- family$linkfun(sum(w * y) / sw)
    out <- list(beta = est, lambda = lambda, interval = NULL)
  }
  out
}

# Build the fitted 'islasso' object from a fit on the standardised scale.
.is_finalize <- function(prep, fit, Xs, lambda, start) {
  X <- prep$X
  y <- prep$y
  family <- prep$family
  weights <- prep$weights
  offset <- prep$offset
  intercept <- prep$intercept
  nobs <- prep$nobs
  nvars <- prep$nvars
  nms <- prep$nms
  setting <- prep$setting
  xm <- prep$xm
  xs <- prep$xs
  A <- .is_backmap(xm, xs, intercept)
  beta <- drop(A %*% fit$theta)

  # covariance: correlation of the sandwich at the solution, variances of the
  # (damped) fixed-point iteration, as in islasso <= 1.6.2
  Vs <- .is_vcov_core(fit$mode, fit$G, fit$Xt, fit$theta, fit$se, fit$cm, fit$lam, prep$alpha)$V
  d <- sqrt(pmax(diag(Vs), 0))
  R <- Vs / outer(d, d)
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  sv <- sqrt(fit$v)
  Vstd <- fit$dispersion * R * outer(sv, sv)
  V <- A %*% Vstd %*% t(A)
  V <- (V + t(V)) / 2
  se <- sqrt(pmax(diag(V), 0))
  dimnames(V) <- list(nms, nms)

  eta <- drop(X %*% beta) + offset
  mu <- family$linkinv(eta)
  me <- family$mu.eta(eta)
  residuals <- (y - mu) / me
  w <- weights * me^2 / family$variance(mu)
  w[!is.finite(w) | weights == 0] <- 0
  dev <- sum(family$dev.resids(y, mu, weights))
  rank <- fit$edf
  resdf <- sum(weights > 0) - rank

  aicv <- NA_real_
  if (nobs > nvars && is.function(family$aic)) {
    aicv <- tryCatch(family$aic(y, prep$n, mu, weights, dev), error = function(e) NA_real_)
  }
  aic.model <- if (is.finite(aicv)) aicv + 2 * rank else dev + 2 * rank

  wtdmu <- if (intercept) sum(weights * y) / sum(weights) else family$linkinv(offset)
  nulldev <- sum(family$dev.resids(y, wtdmu, weights))
  nulldf <- sum(weights > 0) - as.integer(intercept)

  pz <- .is_penalty_R(fit$theta, fit$se, fit$cm, fit$lam, prep$alpha)
  gstd <- -drop(crossprod(Xs, w * residuals)) + pz$gradient
  gradient <- drop(solve(t(A), gstd))   # theta = A^{-1} beta, so dF/dbeta = A^{-T} dF/dtheta
  pen <- prep$alpha * pz$P + (1 - prep$alpha) * fit$theta
  pen[prep$unpenalized] <- 0

  unp <- prep$unpenalized
  names(beta) <- names(se) <- names(gradient) <- names(pen) <- names(unp) <- nms
  hi <- fit$h
  names(hi) <- nms

  internal <- list(
    n = nobs, p = nvars, lambda.seq = fit$lam, vcov = V, vcov.std = Vstd,
    gradient = gradient, hi = hi, intercept = intercept, unpenalized = unp,
    nms = nms, estc = setting$estpai, lmbd.interval = start$interval,
    theta = fit$theta, se.std = fit$se, v = fit$v, c = fit$cm, lam = fit$lam, xm = xm, xs = xs,
    aic.ok = is.finite(aicv), mode = fit$mode, inner.iter = fit$inner_iter, tol = fit$tol, conv = fit$conv, pen = pen,
    n.trials = prep$n
  )
  structure(list(
    coefficients = beta, se = se, dispersion = fit$dispersion,
    residuals = residuals, fitted.values = mu, rank = rank,
    family = family, linear.predictors = eta, deviance = dev, aic = aic.model,
    null.deviance = nulldev, iter = fit$iter, weights = w, prior.weights = weights,
    df.residual = resdf, df.null = nulldf, y = y, converged = fit$conv == 0L,
    model = NULL, call = NULL, formula = NULL, terms = NULL, data = NULL, offset = offset,
    contrasts = NULL, control = setting, internal = internal, xlevels = NULL,
    lambda = lambda, alpha = prep$alpha
  ), class = "islasso")
}
