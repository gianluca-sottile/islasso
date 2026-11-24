#' The Induced Smoothed Lasso
#'
#' Fits regression models with a smoothed lasso (or elastic-net) penalty under the
#' induced smoothing paradigm, which yields point estimates, standard errors and Wald
#' tests for all the coefficients.
#'
#' @aliases islasso print.islasso islasso.fit vcov.islasso deviance.islasso residuals.islasso
#'   logLik.islasso model.matrix.islasso cooks.distance.islasso extractAIC.islasso
#'   family.islasso formula.islasso influence.islasso model.frame.islasso nobs.islasso
#'   rstandard.islasso rstudent.islasso variable.names.islasso weights.islasso
#'
#' @param formula a symbolic description of the model (see [stats::formula()]).
#' @param family the error distribution and link function: a family object, a family
#'   function or its name. Any family accepted by [stats::glm()] can be used, for instance
#'   `gaussian()`, `binomial()` (links logit, probit, cauchit, log, cloglog), `poisson()`,
#'   `Gamma()`, `inverse.gaussian()`, the `quasi*` families and `MASS::negative.binomial()`
#'   with known `theta`.
#' @param lambda the tuning parameter (a non-negative number). If missing, it is chosen by
#'   cross-validation with [glmnet::cv.glmnet()] (see `nfolds` and `seed` in [is.control()]).
#' @param alpha the elastic-net mixing parameter, \eqn{0 \le \alpha \le 1}: `alpha = 1` is the
#'   lasso, `alpha = 0` the ridge penalty.
#' @param data an optional data frame containing the variables in the model.
#' @param weights optional vector of prior weights.
#' @param subset an optional vector specifying a subset of observations.
#' @param offset an optional offset in the linear predictor.
#' @param unpenalized names or indices (with respect to the columns of the model matrix
#'   without intercept) of the coefficients left unpenalized. The intercept is never penalized.
#' @param contrasts an optional list of contrasts for the factors (see [stats::model.matrix()]).
#' @param control a list of control parameters, see [is.control()].
#' @param X,y for `islasso.fit()`: the model matrix (without the column of ones) and the response.
#' @param intercept for `islasso.fit()`: logical, should an (unpenalized) intercept be included?
#'
#' @details
#' For a fixed \eqn{\lambda} the lasso estimating equation
#' \eqn{-X^\top W (z - X\beta) + \lambda\,\mathrm{sign}(\beta) = 0} is not smooth. The induced
#' smoothing (IS) approach (Brown and Wang, 2005) replaces \eqn{\mathrm{sign}(\beta_j)} by its
#' expectation under the approximate sampling distribution of the estimator, modelled as the
#' mixture \eqn{c\,N(0,1) + (1-c)\,N(0,10^{-10})} of a continuous component and a component
#' concentrated at zero:
#' \deqn{P_j = c\,\{2\Phi(\beta_j/s_j) - 1\} + (1-c)\,\{2\Phi(10^5 \beta_j/s_j) - 1\},}
#' where \eqn{s_j} is the standard error of \eqn{\hat\beta_j}. The smoothed equation
#' \eqn{\tilde U(\beta) = -X^\top W(z - X\beta) + \lambda\{\alpha P + (1-\alpha)\beta\} = 0} admits
#' the derivative \eqn{\tilde U'(\beta) = X^\top W X + \lambda\,\mathrm{diag}\{\alpha P' + 1 - \alpha\}},
#' and the covariance matrix is obtained by the sandwich formula
#' \eqn{V = \phi\,\tilde U'^{-1} X^\top W X\, \tilde U'^{-1}}. Point estimates and standard errors are
#' computed alternately until the standard errors stabilise.
#'
#' At fixed standard errors, \eqn{\tilde U} is the gradient of a smooth convex function, which is
#' minimised by a monotone majorise-minimise algorithm (default) or by Newton's method with line
#' search (`inner = "newton"` in [is.control()]). When the number of columns exceeds twice the
#' number of observations, the linear algebra uses the Woodbury identity, with cost
#' \eqn{O(n^2 p)} instead of \eqn{O(p^3)}. Generalised linear models are fitted by iteratively
#' reweighted least squares, so any family with the usual components (`linkinv`, `mu.eta`,
#' `variance`, `dev.resids`) is supported.
#'
#' With `stand = TRUE` (default in [is.control()]) the covariates are standardised (weighted
#' mean and standard deviation with divisor \eqn{\sum w}, as in \pkg{glmnet}) and the penalty
#' acts on the standardised coefficients; the results are reported on the original scale. The
#' scale of `lambda` is that of the deviance: \eqn{\lambda = \sum_i w_i\,\lambda_{glmnet}}.
#'
#' Since the smoothed penalty is never exactly flat, no estimate is exactly zero: the Wald
#' statistics, rather than the point estimates, are used to judge which covariates matter.
#' For coefficients whose estimate is shrunk to (practically) zero the induced-smoothing
#' fixed point drives the standard error towards zero as well, while the Wald ratio stays
#' stable: their standard errors depend on the convergence tolerance and should not be used to
#' build confidence intervals. See `vignette("islasso-methodology", package = "islasso")`.
#'
#' @return An object of class `"islasso"`, a list with components including
#'   \item{coefficients, se}{estimates and standard errors on the original scale;}
#'   \item{dispersion}{the dispersion parameter;}
#'   \item{fitted.values, linear.predictors, residuals, weights, prior.weights}{as in [stats::glm()]
#'     (`residuals` are the working residuals and `weights` the working weights);}
#'   \item{deviance, null.deviance, aic, rank, df.residual, df.null}{deviance-based summaries; `rank`
#'     is the effective number of parameters, the sum of the pseudo hat values;}
#'   \item{converged, iter}{convergence flag and number of iterations;}
#'   \item{lambda, alpha, family, control, call, formula, terms, model, data, xlevels, contrasts, offset, y}{the inputs;}
#'   \item{internal}{other internal quantities (full covariance matrix, pseudo hat values
#'     `hi`, gradient, mixture weights, standardisation constants, ...).}
#'
#' @references
#' Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2020). The induced smoothed lasso: a
#' practical framework for hypothesis testing in high dimensional regression.
#' *Statistical Methods in Medical Research*, **29**(3), 765--777. \doi{10.1177/0962280219842890}
#'
#' Sottile G., Cilluffo G., Muggeo V.M.R. (2019). The R package islasso: estimation and
#' hypothesis testing in lasso regression. Technical report. \doi{10.13140/RG.2.2.16360.11521}
#'
#' Brown B.M., Wang Y.G. (2005). Standard errors and covariance matrices for smoothed rank
#' estimators. *Biometrika*, **92**(1), 149--158. \doi{10.1093/biomet/92.1.149}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#'
#' @seealso [summary.islasso()], [confint.islasso()], [predict.islasso()], [anova.islasso()],
#'   [plot.islasso()], [aic.islasso()], [relax.islasso()], [islasso.path()], [is.control()],
#'   [simulXy()]
#'
#' @examples
#' n <- 100; p <- 50
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso(y ~ ., data = sim1$data, family = gaussian(), lambda = 10)
#' summary(o, pval = 0.05)
#' head(coef(o))
#' logLik(o)
#'
#' \donttest{
#' # lambda chosen by 5-fold cross-validation (reproducible through 'seed')
#' o <- islasso(y ~ ., data = sim1$data, control = is.control(seed = 1))
#' summary(relax.islasso(o), pval = 0.05)
#'
#' # binomial
#' beta <- c(1, 1, 1, rep(0, p - 3))
#' sim2 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1, size = 100,
#'                 family = binomial())
#' o2 <- islasso(cbind(success, failure) ~ ., data = sim2$data, family = binomial(),
#'               lambda = 5)
#' summary(o2, pval = 0.05)
#'
#' # Poisson
#' sim3 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1, family = poisson())
#' o3 <- islasso(y ~ ., data = sim3$data, family = poisson(), lambda = 5)
#' summary(o3, pval = 0.05)
#'
#' # Gamma with log link
#' sim4 <- simulXy(n = n, p = p, beta = beta, interc = -1, seed = 1, dispersion = 0.1,
#'                 family = Gamma(link = "log"))
#' o4 <- islasso(y ~ ., data = sim4$data, family = Gamma(link = "log"), lambda = 5)
#' summary(o4, pval = 0.05)
#' }
#'
#' @export
islasso <- function(formula, family = gaussian, lambda, alpha = 1, data, weights, subset, offset,
                    unpenalized, contrasts = NULL, control = is.control()) {
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
  X <- model.matrix(mt, mf, contrasts)
  offset <- as.vector(model.offset(mf))
  weights <- as.vector(model.weights(mf))
  intercept <- attr(mt, "intercept") != 0
  if (intercept) X <- X[, -1L, drop = FALSE]
  if (!missing(lambda) && (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0))
    stop("'lambda' must be a non-negative number", call. = FALSE)
  if (missing(unpenalized)) unpenalized <- NULL
  fit <- islasso.fit(X = X, y = y, family = family, lambda = lambda, alpha = alpha,
                     intercept = intercept, weights = weights, offset = offset,
                     unpenalized = unpenalized, control = control)
  fit$model <- mf
  fit$call <- this.call
  fit$formula <- formula
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- attr(model.matrix(mt, mf, contrasts), "contrasts")
  fit$xlevels <- .getXlevels(mt, mf)
  fit
}

#' @rdname islasso
#' @export
islasso.fit <- function(X, y, family = gaussian, lambda, alpha = 1, intercept = FALSE, weights = NULL,
                        offset = NULL, unpenalized = NULL, control = is.control()) {
  control <- do.call(is.control, control)
  prep <- .is_prepare(X, y, family, alpha, intercept, weights, offset, unpenalized, control)
  st <- .is_standardize(prep$X, prep$weights, intercept, isTRUE(control$stand))
  prep$xm <- st$xm
  prep$xs <- st$xs
  start <- .is_start(prep, st$X, if (missing(lambda)) NULL else lambda)
  lambda <- start$lambda
  lam <- lambda * as.numeric(!prep$unpenalized)
  theta <- start$beta
  theta[theta == 0] <- 1e-5
  v0 <- if (is.null(control$V0)) rep(0.01, prep$nvars) else {
    V0 <- as.matrix(control$V0)
    if (length(V0) == 1L) rep(V0[1L], prep$nvars) else diag(V0)
  }
  fit <- .is_fit_single(st$X, prep$y, prep$family, lam, prep$alpha, prep$weights, prep$offset,
                        intercept, prep$setting, theta, v0, prep$setting$c)
  if (fit$conv == 1L) warning("maximum number of iterations reached", call. = FALSE)
  .is_finalize(prep, fit, st$X, lambda, start)
}

#' Control Parameters for islasso Fitting
#'
#' Collects the parameters controlling the fitting process of [islasso()] and
#' [islasso.path()].
#'
#' @param sigma2 fixed value of the dispersion parameter; if negative (default) it is
#'   estimated (Pearson statistic over the residual degrees of freedom). Ignored for the
#'   binomial, Poisson and negative binomial families, where it is 1.
#' @param tol tolerance of the convergence criteria (default `1e-5`).
#' @param itmax maximum number of iterations (default 1000).
#' @param stand logical: should the covariates be standardised (default `TRUE`)? With
#'   `stand = TRUE` the penalty acts on the standardised coefficients, as in \pkg{glmnet};
#'   estimates are always returned on the original scale.
#' @param trace integer: 0 (default) silent, 1 compact and 2 detailed printing of the iterations.
#' @param nfolds number of folds for the cross-validation used when `lambda` is missing.
#' @param seed optional integer seed for the cross-validation folds (the random number
#'   generator state of the session is restored afterwards).
#' @param adaptive logical: adaptive lasso with weights based on the pseudo hat values
#'   (experimental).
#' @param g numeric in \[0, 1\]: the parameter of the extended BIC (`g = 0` gives the BIC).
#' @param b0 optional starting values for the coefficients (original scale, intercept first).
#' @param V0 optional starting covariance matrix (or a number for its diagonal) on the
#'   standardised scale; default `0.01`.
#' @param c the mixture weight of the induced smoothing (default 0.5), a number in \[0, 1\]
#'   or a vector with one value per coefficient; a negative value means that the weights are
#'   re-estimated at each iteration.
#' @param inner the algorithm solving the smoothed estimating equation at fixed standard
#'   errors: `"mm"` (default, monotone majorise-minimise iterations, the undamped version of
#'   the fixed-point iterations of islasso 1.6) or `"newton"` (Newton's method with line search).
#' @param solver linear algebra: `"auto"` (default) uses `"woodbury"` when the number of
#'   columns exceeds twice the number of observations and `"dense"` otherwise.
#' @param ... further arguments, ignored (allows `do.call(is.control, control)`).
#'
#' @return A list with the control parameters.
#' @author Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#' @seealso [islasso()], [islasso.path()]
#' @examples
#' is.control(tol = 1e-6, seed = 1)
#' @export
is.control <- function(sigma2 = -1, tol = 1e-05, itmax = 1000, stand = TRUE, trace = 0, nfolds = 5,
                       seed = NULL, adaptive = FALSE, g = 0.5, b0 = NULL, V0 = NULL, c = 0.5,
                       inner = c("mm", "newton"), solver = c("auto", "dense", "woodbury"), ...) {
  inner <- match.arg(inner)
  solver <- match.arg(solver)
  if (!is.numeric(tol) || tol <= 0) stop("'tol' must be positive", call. = FALSE)
  if (!is.numeric(itmax) || itmax < 1) stop("'itmax' must be a positive integer", call. = FALSE)
  if (!is.numeric(g) || g < 0 || g > 1) stop("'g' must be in [0, 1]", call. = FALSE)
  list(sigma2 = sigma2, tol = tol, itmax = as.integer(itmax), trace = as.integer(trace),
       stand = as.logical(stand), nfolds = nfolds, seed = seed, adaptive = as.logical(adaptive),
       g = g, b0 = b0, V0 = V0, c = c, inner = inner, solver = solver)
}
