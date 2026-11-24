#' Optimal Tuning Parameter by Information Criteria
#'
#' Selects the tuning parameter `lambda` of an [islasso()] model by minimising an information
#' criterion (AIC, BIC/eBIC, AICc, GCV or GIC; see [islasso.path()] for the definitions), as an
#' alternative to cross-validation.
#'
#' @param object a fitted model of class `"islasso"` (optional if `y` and `X` are given).
#' @param method the criterion to minimise.
#' @param interval the search interval for `lambda`; by default the range of the
#'   \pkg{glmnet} path used to obtain the starting values of `object`.
#' @param g numeric in \[0, 1\]: the parameter of the extended BIC (`g = 0` gives the BIC).
#' @param y,X response and model matrix (without the column of ones), used if `object` is missing.
#' @param intercept,family,alpha,offset,weights,unpenalized,control as in [islasso.fit()],
#'   used if `object` is missing.
#' @param trace logical: print the progress of the optimisation?
#'
#' @details The criterion is minimised over \eqn{\log\lambda} with [stats::optimize()]; the
#' degrees of freedom are the effective number of parameters (trace of the pseudo hat matrix).
#'
#' @return The value of `lambda` minimising the criterion.
#' @author Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#' @seealso [islasso()], [islasso.path()], [GoF.islasso.path()]
#'
#' @examples
#' \donttest{
#' n <- 100; p <- 30
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso(y ~ ., data = sim1$data, family = gaussian(), lambda = 5)
#' lambda.bic <- aic.islasso(o, method = "BIC", interval = c(1, 50), trace = FALSE)
#' o2 <- update(o, lambda = lambda.bic)
#' }
#'
#' @export
aic.islasso <- function(object, method = c("AIC", "BIC", "AICc", "GCV", "GIC"), interval, g = 0,
                        y, X, intercept = FALSE, family = gaussian(), alpha = 1, offset, weights,
                        unpenalized, control = is.control(), trace = TRUE) {
  method <- match.arg(method)
  if (!is.numeric(g) || g < 0 || g > 1) stop("'g' must be in [0, 1]", call. = FALSE)
  if (missing(object)) {
    if (missing(interval)) stop("please specify an interval to search for a minimum", call. = FALSE)
    if (missing(y)) stop("the response 'y' is missing", call. = FALSE)
    if (missing(X)) stop("the model matrix 'X' is missing", call. = FALSE)
    if (missing(offset)) offset <- NULL
    if (missing(weights)) weights <- NULL
    if (missing(unpenalized)) unpenalized <- NULL
    ny <- if (is.matrix(y)) nrow(y) else length(y)
    if (!is.null(weights)) ny <- sum(weights > 0)
  } else {
    if (!inherits(object, "islasso")) stop("'object' must be of class 'islasso'", call. = FALSE)
    if (missing(interval)) interval <- object$internal$lmbd.interval
    if (is.null(interval)) stop("please specify an interval to search for a minimum", call. = FALSE)
    X <- model.matrix(object)
    intercept <- object$internal$intercept
    if (intercept) X <- X[, -1L, drop = FALSE]
    unp <- object$internal$unpenalized
    unpenalized <- which(if (intercept) unp[-1L] else unp)
    family <- object$family
    alpha <- object$alpha
    offset <- object$offset
    control <- object$control
    if (isTRUE(control$estpai)) control$c <- -1 else if (length(unique(control$c)) == 1L) control$c <- control$c[1L]
    # response and prior weights as given to the original fit
    y <- object$y
    weights <- object$prior.weights
    ny <- sum(weights > 0)
  }
  if (length(interval) != 2L || any(interval <= 0)) stop("'interval' must contain two positive values",
                                                          call. = FALSE)
  control$trace <- 0L
  p <- NCOL(X) + intercept
  crit <- function(llam) {
    lambda <- exp(llam)
    obj <- tryCatch(islasso.fit(X = X, y = y, family = family, lambda = lambda, alpha = alpha,
                                intercept = intercept, weights = weights, offset = offset,
                                unpenalized = unpenalized, control = control),
                    error = function(e) NULL)
    if (is.null(obj)) return(Inf)
    ll <- obj$aic - 2 * obj$rank
    val <- .is_gof(ll, obj$deviance, obj$rank, ny, p, g)[[switch(method, BIC = if (g == 0) "BIC" else "eBIC",
                                                                    method)]]
    if (trace) cat("lambda = ", formatC(lambda, digits = 4, width = 10, format = "f"), "  ", method,
                   " = ", formatC(val, digits = 5, width = 12, format = "f"), "\n", sep = "")
    val
  }
  if (trace) cat("\nOptimization through ", if (method == "BIC" && g != 0) "eBIC" else method, "\n\n", sep = "")
  opt <- stats::optimize(crit, interval = log(sort(interval)))
  exp(opt$minimum)
}

#' Relaxed islasso Fit
#'
#' Refits an [islasso()] model leaving unpenalized a subset of covariates, chosen by their
#' p-values or given explicitly, while the other covariates are (practically) removed by a large
#' penalty.
#'
#' @param fit a fitted model of class `"islasso"`.
#' @param id optional names or indices (with respect to the covariates, intercept excluded)
#'   of the covariates to leave unpenalized; by default those with p-value not greater than `pval`.
#' @param pval the significance level used when `id` is missing (default 0.05).
#' @param ... further arguments replacing those of the original call.
#'
#' @details The refit uses a value of `lambda` equal to 100 times the smallest value that
#' shrinks all the penalized coefficients to zero, so that the unselected covariates have
#' (practically) null estimates.
#'
#' @return An object of class `"islasso"`.
#' @seealso [islasso()], [summary.islasso()]
#' @examples
#' \donttest{
#' n <- 100; p <- 30
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso(y ~ ., data = sim1$data, lambda = 10)
#' summary(relax.islasso(o), pval = 0.05)
#' summary(relax.islasso(o, id = c("X1", "X2")), pval = 0.05)
#' }
#' @export
relax.islasso <- function(fit, id = NULL, pval = 0.05, ...) {
  if (!inherits(fit, "islasso")) stop("'fit' must be of class 'islasso'", call. = FALSE)
  intercept <- fit$internal$intercept
  nms <- fit$internal$nms
  cov_nms <- if (intercept) nms[-1L] else nms
  p <- length(cov_nms)
  unp <- fit$internal$unpenalized
  unp <- if (intercept) unp[-1L] else unp
  if (is.null(id)) {
    pv <- 2 * stats::pnorm(-abs(fit$coefficients / fit$se))
    if (intercept) pv <- pv[-1L]
    id <- which(pv <= pval)
  } else if (is.character(id)) {
    m <- match(id, cov_nms)
    if (anyNA(m)) stop("some names in 'id' are not covariates of the model", call. = FALSE)
    id <- m
  } else if (!is.numeric(id) || any(id < 1 | id > p)) {
    stop(sprintf("the values in 'id' must be between 1 and %d", p), call. = FALSE)
  }
  unp[id] <- TRUE
  if (all(unp)) stop("all the covariates are selected: the relaxed model would be unpenalized; use glm()",
                     call. = FALSE)
  if (!any(unp[id])) warning("no covariate selected", call. = FALSE)
  # large lambda: 100 times lambda.max given the unpenalized covariates
  X <- model.matrix(fit)
  prep <- .is_prepare(if (intercept) X[, -1L, drop = FALSE] else X, fit$y, fit$family, fit$alpha,
                      intercept, fit$prior.weights, fit$offset, which(unp), fit$control)
  st <- .is_standardize(prep$X, prep$weights, intercept, isTRUE(fit$control$stand))
  null <- .is_null_fit(st$X, prep, prep$unpenalized)
  lambda <- 100 * .is_lambda_max(st$X, prep, null$mu, prep$unpenalized)
  cl <- fit$call
  if (is.null(cl)) stop("'fit' has no call component", call. = FALSE)
  cl$unpenalized <- cov_nms[unp]
  cl$lambda <- lambda
  extras <- list(...)
  for (nm in names(extras)) cl[[nm]] <- extras[[nm]]
  eval(cl, parent.frame())
}
