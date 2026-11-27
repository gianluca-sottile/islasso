#' Induced Smoothed Lasso Regularization Path
#'
#' Fits the induced smoothed lasso (or elastic net) over a grid of values of the tuning
#' parameter, with warm starts.
#'
#' @aliases islasso.path
#' @aliases print.islasso.path
#' @aliases islasso.path.fit
#' @aliases coef.islasso.path
#' @aliases deviance.islasso.path
#' @aliases fitted.islasso.path
#' @aliases residuals.islasso.path
#' @aliases logLik.islasso.path
#' @aliases print.logLik.islasso.path
#' @aliases model.matrix.islasso.path
#'
#' @inheritParams islasso
#' @param lambda optional vector of values of the tuning parameter; by default a sequence
#'   of `nlambda` values, log-equispaced between `lambda.min.ratio * lambda.max` and
#'   `lambda.max`, the smallest value giving (practically) null penalized coefficients.
#' @param nlambda number of values of `lambda` (default 100).
#' @param lambda.min.ratio ratio between the smallest and the largest value of `lambda`;
#'   default `1e-2` if the number of observations is smaller than the number of covariates and
#'   `1e-4` otherwise.
#'
#' @details
#' Each model of the path is fitted as [islasso()] would do at the same value of `lambda`
#' (starting values from the \pkg{glmnet} path, standard errors started at their default
#' values), so that the path and a single fit agree; with the Gaussian family the matrix
#' \eqn{X^\top W X} is computed once for the whole path. The models are independent and can be
#' fitted in parallel with `control = is.control(ncores = )` (socket cluster of the 'parallel'
#' package; worthwhile for large problems, since starting the workers takes about a second). Every model
#' is kept in the output, together with its convergence flag (column `conv` of `Info`), except
#' those for which the algorithm failed (a warning is given).
#' The information criteria stored in `GoF` are
#' \eqn{AIC = -2\ell + 2\,df}, \eqn{BIC = -2\ell + \log(n)\,df},
#' \eqn{AICc = AIC + 2\,df(df+1)/(n-df-1)}, \eqn{eBIC = -2\ell + \{\log(n) + 2g\log(p)\}\,df},
#' \eqn{GCV = n\,D/(n-df)^2} (with \eqn{D} the deviance) and \eqn{GIC = -2\ell + \log\log(n)\log(p)\,df},
#' where \eqn{df} is the effective number of parameters (sum of the pseudo hat values) and
#' \eqn{-2\ell} is replaced by the deviance when \eqn{p \ge n} or the family has no likelihood.
#'
#' @return An object of class `"islasso.path"`, a list with components
#'   \item{Coef, SE, Weight, Gradient}{matrices (one row per value of `lambda`, in increasing
#'     order) of estimates, standard errors, mixture weights and gradients of the smoothed
#'     estimating equation (a convergence diagnostic, numerically close to zero);}
#'   \item{Score}{matrix of standardized scores
#'     \eqn{\{x_j^\top W(z - X\hat\beta)/\lambda_j - (1-\alpha)\hat\beta_j\}/\alpha} (standardized
#'     covariates), equal at the solution to the smoothed sign \eqn{P_j \in [-1, 1]}: the
#'     induced-smoothing analogue of the lasso / elastic-net optimality conditions (0 for unpenalized
#'     columns, `NA` for the ridge penalty, `alpha = 0`);}
#'   \item{Info}{matrix with columns `lambda`, `df`, `phi`, `deviance`, `logLik`, `iter`, `err`, `conv`;}
#'   \item{GoF}{matrix of information criteria (AIC, BIC, AICc, eBIC, GCV, GIC);}
#'   \item{Linear.predictors, Fitted.values, Residuals}{matrices with one row per value of `lambda`;}
#'   \item{Input}{list with `n`, `p`, `intercept`, `unpenalized`;}
#'   \item{call, family, prior.weights, y, offset, alpha, control, formula, model, terms, data,
#'     xlevels, contrasts}{the inputs.}
#'
#' @references
#' Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2020). The induced smoothed lasso: a
#' practical framework for hypothesis testing in high dimensional regression.
#' *Statistical Methods in Medical Research*, **29**(3), 765--777. \doi{10.1177/0962280219842890}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#'
#' @seealso [plot.islasso.path()], [predict.islasso.path()], [summary.islasso.path()],
#'   [GoF.islasso.path()], [islasso()], [is.control()]
#'
#' @examples
#' n <- 100; p <- 30; p1 <- 10
#' beta.veri <- sort(round(c(seq(.5, 3, length.out = p1/2),
#'                          seq(-1, -2, length.out = p1/2)), 2))
#' beta <- c(beta.veri, rep(0, p - p1))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso.path(y ~ ., data = sim1$data, family = gaussian(), nlambda = 30L)
#' o
#' summary(o, lambda = 10, pval = 0.05)
#' coef(o, lambda = 10)
#' GoF.islasso.path(o, plot = FALSE)$lambda.min
#'
#' \donttest{
#' plot(o, yvar = "coefficients", gof = "AIC")
#' # standardized scores: +-1 for the relevant covariates, inside (-1, 1) for the others
#' plot(o, yvar = "score", gof = "AIC", label = TRUE)
#'
#' # Poisson
#' beta <- c(1, 1, 1, rep(0, p - 3))
#' sim3 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1, family = poisson())
#' o3 <- islasso.path(y ~ ., data = sim3$data, family = poisson(), nlambda = 30L)
#' temp <- GoF.islasso.path(o3, plot = FALSE)
#' summary(o3, pval = 0.05, lambda = temp$lambda.min["BIC"])
#' }
#'
#' @export
islasso.path <- function(formula, family = gaussian(), lambda = NULL, nlambda = 100,
                         lambda.min.ratio = ifelse(nobs < nvars, 1e-2, 1e-4),
                         alpha = 1, data, weights, subset, offset, contrasts = NULL,
                         unpenalized, control = is.control()) {
  this.call <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) stop("the model matrix is empty", call. = FALSE)
  y <- model.response(mf, "any")
  x <- model.matrix(mt, mf, contrasts)
  nobs <- nrow(x)
  nvars <- ncol(x)
  force(lambda.min.ratio)
  offset <- as.vector(model.offset(mf))
  weights <- as.vector(model.weights(mf))
  intercept <- attr(mt, "intercept") != 0
  if (intercept) x <- x[, -1L, drop = FALSE]
  if (missing(unpenalized)) unpenalized <- NULL
  fit <- islasso.path.fit(X = x, y = y, family = family, lambda = lambda, nlambda = nlambda,
                          lambda.min.ratio = lambda.min.ratio, alpha = alpha, intercept = intercept,
                          weights = weights, offset = offset, unpenalized = unpenalized,
                          control = control)
  fit$call <- this.call
  fit$formula <- formula
  fit$model <- mf
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- attr(model.matrix(mt, mf, contrasts), "contrasts")
  fit$xlevels <- .getXlevels(mt, mf)
  fit$na.action <- attr(mf, "na.action")
  fit
}

#' @rdname islasso.path
#' @export
islasso.path.fit <- function(X, y, family = gaussian(), lambda = NULL, nlambda = 100,
                             lambda.min.ratio = ifelse(NROW(X) < NCOL(X), 1e-2, 1e-4), alpha = 1,
                             intercept = FALSE, weights = NULL, offset = NULL, unpenalized = NULL,
                             control = is.control()) {
  control <- do.call(is.control, control)
  prep <- .is_prepare(X, y, family, alpha, intercept, weights, offset, unpenalized, control)
  st <- .is_standardize(prep$X, prep$weights, intercept, isTRUE(control$stand))
  prep$xm <- st$xm
  prep$xs <- st$xs
  Xs <- st$X
  family <- prep$family
  nobs <- prep$nobs
  nvars <- prep$nvars
  setting <- prep$setting
  unp <- prep$unpenalized

  null <- .is_null_fit(Xs, prep, unp)
  if (is.null(lambda)) {
    lmax <- .is_lambda_max(Xs, prep, null$mu, unp)
    if (!is.numeric(lambda.min.ratio) || lambda.min.ratio <= 0 || lambda.min.ratio >= 1)
      stop("'lambda.min.ratio' must be in (0, 1)", call. = FALSE)
    lambda <- exp(seq(log(lmax), log(lmax * lambda.min.ratio), length.out = nlambda))
  } else {
    if (!is.numeric(lambda) || any(lambda < 0)) stop("'lambda' must be non-negative", call. = FALSE)
    lambda <- sort(unique(lambda), decreasing = TRUE)
  }
  nlambda <- length(lambda)

  mode <- .is_mode(nobs, nvars, setting)
  G <- bz <- NULL
  if (.is_gaussian_identity(family) && mode == 0L) {
    sw <- sqrt(prep$weights)
    Xt <- Xs * sw
    G <- crossprod(Xt)
    bz <- drop(crossprod(Xt, sw * (prep$y - prep$offset)))
  }

  # starting values: glmnet at each lambda, as in islasso() (otherwise the unpenalized fit)
  starts <- .is_start_path(prep, Xs, lambda)
  theta_null <- null$theta
  v_init <- if (is.null(setting$V0)) rep(0.01, nvars) else rep_len(diag(as.matrix(setting$V0)), nvars)

  ll_ok <- nobs > nvars && is.function(family$aic)
  nok <- sum(prep$weights > 0)
  disp <- family$family %in% c("gaussian", "Gamma", "inverse.gaussian")
  trace <- setting$trace
  setting_in <- setting
  setting_in$trace <- 0L

  # one model of the path; the standard errors (and mixture weights) always start from their
  # initial values, so that every model coincides with islasso() at the same lambda (the IS fixed
  # point of the null coefficients depends on the starting standard errors)
  fit_one <- function(k, theta = NULL) {
    lam <- lambda[k] * as.numeric(!unp)
    if (is.null(theta)) theta <- if (!is.null(starts)) starts[, k] else theta_null
    theta[theta == 0] <- 1e-5
    fit <- tryCatch(
      .is_fit_single(Xs, prep$y, family, lam, prep$alpha, prep$weights, prep$offset, intercept,
                     setting_in, theta, v_init, setting$c, G = G, bz = bz, mode = mode),
      error = function(e) conditionMessage(e))
    if (is.character(fit)) return(list(error = fit))
    theta <- fit$theta
    cm_k <- fit$cm
    beta <- .is_theta2beta(theta, st$xm, st$xs, intercept)
    se <- sqrt(fit$dispersion * fit$v) / st$xs
    if (intercept && nvars > 1L) {
      # the intercept on the original scale combines all the standardised coefficients:
      # variance a' Vstd a with Vstd = phi * diag(sqrt(v / vn)) S diag(sqrt(v / vn)), S = sandwich
      a <- .is_intercept_row(st$xm, st$xs)
      b <- a * ifelse(fit$vn > 0, sqrt(fit$v / fit$vn), 0)
      q <- .is_quad(fit$mode, fit$G, fit$Xt, theta, fit$se.used, cm_k, fit$lam, prep$alpha, b)
      se[1L] <- sqrt(fit$dispersion * max(q, 0))
    }
    eta <- drop(prep$X %*% beta) + prep$offset
    mu <- family$linkinv(eta)
    me <- family$mu.eta(eta)
    dev <- sum(family$dev.resids(prep$y, mu, prep$weights))
    df <- fit$edf
    ll <- dev
    llfam <- FALSE
    if (ll_ok) {
      a <- tryCatch(family$aic(prep$y, prep$n, mu, prep$weights, dev), error = function(e) NA_real_)
      if (is.finite(a)) { ll <- a - 2 * disp; llfam <- TRUE }  # -2 loglik (family$aic adds 2 for phi)
    }
    w <- prep$weights * me^2 / family$variance(mu)
    w[!is.finite(w)] <- 0
    pz <- .is_penalty_R(theta, fit$se.used, cm_k, fit$lam, prep$alpha)
    score <- drop(crossprod(Xs, w * (prep$y - mu) / me))
    gstd <- -score + pz$gradient
    # standardized score {x_j' W (z - X beta) / lambda_j - (1 - alpha) beta_j} / alpha, equal at the
    # solution to the smoothed sign P_j in [-1, 1] (lasso: x_j' W (z - X beta) / lambda_j); 0 for the
    # unpenalized columns, NA for the ridge penalty (alpha = 0), for which it is not defined
    zscore <- if (prep$alpha > 0) {
      ifelse(fit$lam > 0,
             (score / pmax(fit$lam, .Machine$double.xmin) - (1 - prep$alpha) * theta) / prep$alpha, 0)
    } else rep(NA_real_, length(theta))
    list(theta = theta, coef = beta, se = se, w = cm_k, grad = .is_grad_back(gstd, st$xm, st$xs, intercept),
         score = zscore,
         info = c(lambda[k], df, fit$dispersion, dev, -ll / 2, fit$iter, fit$tol, fit$conv),
         gof = .is_gof(if (llfam) ll + 2 * disp else ll, dev, df, nok, nvars, setting$g),
         eta = eta, mu = mu, res = (prep$y - mu) / me)
  }

  ncores <- max(1L, as.integer(setting$ncores))
  if (ncores > 1L && !is.null(starts)) {
    # socket cluster (portable, and safe with multi-threaded BLAS libraries, unlike forking)
    cl <- parallel::makePSOCKcluster(min(ncores, nlambda))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(lp) invisible(.libPaths(lp)), .libPaths())
    fits <- parallel::parLapply(cl, seq_len(nlambda), fit_one)
  } else {
    if (trace == 2L) cat("\n   Ind|Max   lambda       df        phi      iter   error\n")
    # sequential fits: same starting values as in parallel (glmnet at each lambda), so that the
    # results do not depend on 'ncores'; without glmnet starts, warm start of the estimates
    fits <- vector("list", nlambda)
    theta <- NULL
    for (k in seq_len(nlambda)) {
      fits[[k]] <- fit_one(k, if (is.null(starts)) theta else NULL)
      theta <- fits[[k]]$theta
      if (trace == 2L && !is.null(fits[[k]]$info)) {
        f <- fits[[k]]$info
        cat(sprintf("%4d|%4d %10.4f %8.3f %10.4f %6d %9.2e\n", k, nlambda, f[1], f[2], f[3],
                    as.integer(f[6]), f[7]))
      } else if (trace == 1L) {
        cat("\r", k, "|", nlambda, "      ")
      }
    }
    if (trace > 0L) cat("\n")
  }
  ok <- vapply(fits, function(f) is.list(f) && !is.null(f$coef), logical(1L))
  msg <- unique(unlist(lapply(fits[!ok], function(f) if (is.list(f)) f$error else as.character(f))))
  if (!any(ok)) stop("the algorithm failed for all the values of lambda: ", paste(msg, collapse = "; "),
                     call. = FALSE)
  if (any(!ok))
    warning("the algorithm failed for ", sum(!ok), " values of lambda, removed from the output (",
            paste(msg, collapse = "; "), ")", call. = FALSE)
  fits <- fits[ok]
  stack <- function(name) do.call(rbind, lapply(fits, `[[`, name))
  out_coef <- stack("coef"); out_se <- stack("se"); out_w <- stack("w"); out_grad <- stack("grad")
  out_score <- stack("score")
  out_info <- stack("info"); out_gof <- stack("gof")
  out_eta <- stack("eta"); out_mu <- stack("mu"); out_res <- stack("res")
  if (any(out_info[, 8L] == 1))
    warning("maximum number of iterations reached for some values of lambda (conv = 1 in 'Info')",
            call. = FALSE)

  keep <- rev(seq_len(nrow(out_coef)))   # increasing lambda
  nms <- prep$nms
  out_info <- out_info[keep, , drop = FALSE]
  colnames(out_info) <- c("lambda", "df", "phi", "deviance", "logLik", "iter", "err", "conv")
  out_gof <- out_gof[keep, , drop = FALSE]
  colnames(out_gof) <- c("AIC", "BIC", "AICc", "eBIC", "GCV", "GIC")
  mk <- function(M) { M <- M[keep, , drop = FALSE]; colnames(M) <- nms; M }
  structure(list(
    Coef = mk(out_coef), SE = mk(out_se), Residuals = out_res[keep, , drop = FALSE],
    Fitted.values = out_mu[keep, , drop = FALSE], Info = out_info, GoF = out_gof,
    Linear.predictors = out_eta[keep, , drop = FALSE], prior.weights = prep$weights,
    Weight = mk(out_w), Gradient = mk(out_grad), Score = mk(out_score), y = prep$y, model = NULL, call = match.call(),
    family = family, formula = NULL, terms = NULL, data = NULL, offset = prep$offset,
    control = setting, contrasts = NULL, xlevels = NULL, alpha = prep$alpha,
    Input = list(n = nobs, p = nvars, intercept = intercept, unpenalized = unp, n.trials = prep$n)
  ), class = c("islasso.path", "islasso"))
}

# glmnet starting values (standardised scale) for all the values of lambda; NULL if unavailable
.is_start_path <- function(prep, Xs, lambda) {
  intercept <- prep$intercept
  xx <- if (intercept) Xs[, -1L, drop = FALSE] else Xs
  if (ncol(xx) < 2L) return(NULL)
  w <- prep$weights
  sw <- sum(w)
  pf <- as.numeric(!prep$unpenalized)
  pfx <- if (intercept) pf[-1L] else pf
  if (sum(pfx) == 0) return(NULL)
  gfam <- .is_glmnet_family(prep$family)
  if (is.character(gfam) && gfam == "binomial") {
    y2 <- cbind(1 - prep$y, prep$y) * w
    w2 <- rep(1, length(prep$y))
  } else {
    y2 <- prep$y
    w2 <- w
  }
  obj <- tryCatch(suppressWarnings(glmnet::glmnet(
    x = xx, y = y2, family = gfam, alpha = prep$alpha, weights = w2, offset = prep$offset,
    standardize = FALSE, intercept = intercept, penalty.factor = pfx,
    lambda = sort(lambda / sw, decreasing = TRUE))), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  B <- tryCatch(as.matrix(stats::coef(obj, s = lambda / sw)), error = function(e) NULL)
  if (is.null(B) || ncol(B) != length(lambda) || anyNA(B)) return(NULL)
  if (!intercept) B <- B[-1L, , drop = FALSE]
  unname(B)
}

# fit of the unpenalized part (intercept and 'unpenalized' columns) on the standardised scale
.is_null_fit <- function(Xs, prep, unp) {
  family <- prep$family
  p <- ncol(Xs)
  theta <- numeric(p)
  if (any(unp)) {
    fit <- suppressWarnings(stats::glm.fit(Xs[, unp, drop = FALSE], prep$y, weights = prep$weights,
                                           offset = prep$offset, family = family,
                                           intercept = FALSE))
    b <- fit$coefficients
    b[is.na(b)] <- 0
    theta[unp] <- b
    mu <- fit$fitted.values
  } else {
    mu <- family$linkinv(prep$offset)
  }
  list(theta = theta, mu = mu)
}

# smallest lambda with null penalized coefficients (lasso part), on the islasso scale
.is_lambda_max <- function(Xs, prep, mu, unp) {
  family <- prep$family
  eta <- family$linkfun(mu)
  r <- prep$weights * (prep$y - mu) * family$mu.eta(eta) / family$variance(mu)
  r[!is.finite(r)] <- 0
  g <- abs(drop(crossprod(Xs, r)))
  g[unp] <- 0
  lmax <- max(g) / max(prep$alpha, 1e-3)
  if (!(lmax > 0)) lmax <- 1
  lmax
}

.is_gof <- function(ll, dev, df, n, p, g) {
  c(AIC = ll + 2 * df,
    BIC = ll + log(n) * df,
    AICc = ll + 2 * df + 2 * df * (df + 1) / max(n - df - 1, 1e-8),
    eBIC = ll + (log(n) + 2 * g * log(p)) * df,
    GCV = n * dev / max(n - df, 1e-8)^2,
    GIC = ll + log(log(n)) * log(p) * df)
}

#' Select the Tuning Parameter by Information Criteria
#'
#' Finds, along an [islasso.path()] fit, the values of `lambda` minimising the information
#' criteria AIC, BIC, AICc, eBIC, GCV and GIC (see [islasso.path()] for their definitions).
#'
#' @param object a fitted model of class `"islasso.path"`.
#' @param plot logical: should the criteria be plotted against \eqn{\log\lambda}?
#' @param ... not used.
#'
#' @details The degrees of freedom are the effective number of parameters (trace of the pseudo
#' hat matrix), which is not necessarily an integer.
#'
#' @return A list with components
#' \item{gof}{the matrix of criteria along the path;}
#' \item{minimum}{the positions of the minima;}
#' \item{lambda.min}{the values of `lambda` minimising each criterion;}
#' \item{plot}{the `ggplot` object (if `plot = TRUE`).}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#' @seealso [islasso.path()], [summary.islasso.path()], [predict.islasso.path()]
#'
#' @examples
#' n <- 100; p <- 30
#' beta <- c(runif(10, -2, 2), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' fit <- islasso.path(y ~ ., data = sim$data, family = gaussian(), nlambda = 30)
#' GoF.islasso.path(fit, plot = FALSE)$lambda.min
#'
#' @export
GoF.islasso.path <- function(object, plot = TRUE, ...) {
  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)
  gof.name <- c("AIC", "BIC", "AICc", "eBIC", "GCV", "GIC")
  gof <- object$GoF
  id.min <- apply(gof, 2, which.min)
  lambda.min <- lambda.seq[id.min]
  names(lambda.min) <- gof.name
  out <- list(gof = gof, minimum = id.min, lambda.min = lambda.min)
  if (plot) {
    df <- data.frame(lambda = log(lambda.seq), value = c(gof),
                     measure = factor(rep(gof.name, each = nlambda), levels = gof.name))
    vl <- data.frame(measure = factor(gof.name, levels = gof.name), opt = log(lambda.min))
    p <- ggplot(df, aes(x = lambda, y = value)) +
      geom_line() +
      geom_point(size = 1) +
      facet_wrap(~ measure, scales = "free_y", nrow = 2, ncol = 3) +
      theme_bw() +
      xlab(expression(log(lambda))) +
      ylab("") +
      geom_vline(data = vl, aes(xintercept = opt), linetype = "dashed", colour = "red") +
      theme(strip.background = element_rect(fill = "white"))
    out$plot <- p
    print(p)
  }
  out
}
