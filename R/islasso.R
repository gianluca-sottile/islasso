#' The Induced Smoothed Lasso
#'
#' Fits regression models with a smoothed L1 penalty under the induced smoothing paradigm.
#' Supports linear, logistic, Poisson, and Gamma responses. Enables reliable standard errors
#' and Wald-based inference.
#'
#' @aliases islasso print.islasso islasso.fit vcov.islasso deviance.islasso residuals.islasso
#'          logLik.islasso model.matrix.islasso cooks.distance.islasso extractAIC.islasso
#'          family.islasso formula.islasso influence.islasso model.frame.islasso nobs.islasso
#'          rstandard.islasso rstudent.islasso variable.names.islasso weights.islasso
#'
#' @param formula A symbolic formula describing the model.
#' @param family Response distribution. Can be \code{gaussian}, \code{binomial}, \code{poisson}, or \code{Gamma}.
#' @param lambda Regularization parameter. If missing, it is estimated via \code{\link[glmnet]{cv.glmnet}}.
#' @param alpha Elastic-net mixing parameter (\eqn{0 \le \alpha \le 1}).
#' @param data A data frame or environment containing the variables in the model.
#' @param weights Observation weights. Defaults to 1.
#' @param subset Optional vector specifying a subset of rows to include.
#' @param offset Optional numeric vector of offsets in the linear predictor.
#' @param unpenalized Vector indicating variables (by name or index) to exclude from penalization.
#' @param contrasts Optional contrasts specification for factor variables.
#' @param control A list of parameters to control model fitting. See \code{\link{is.control}}.
#'
#' @details
#' The non-smooth L1 penalty is replaced by a smooth approximation, enabling inference through
#' standard errors and Wald tests. The approach controls type-I error and shows strong power
#' in various simulation settings.
#'
#' @return A list with components such as:
#'   \item{coefficients}{Estimated coefficients}
#'   \item{se}{Standard errors}
#'   \item{fitted.values}{Fitted values}
#'   \item{deviance, aic, null.deviance}{Model diagnostic metrics}
#'   \item{residuals, weights}{IWLS residuals and weights}
#'   \item{df.residual, df.null, rank}{Degrees of freedom}
#'   \item{converged}{Logical; convergence status}
#'   \item{model, call, terms, formula, data, offset}{Model objects}
#'   \item{xlevels, contrasts}{Factor handling details}
#'   \item{lambda, alpha, dispersion}{Model parameters}
#'   \item{internal}{Other internal values}
#'
#' @references
#' Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2019)
#' \emph{The Induced Smoothed Lasso: A practical framework for hypothesis testing in high dimensional regression}.
#' Statistical Methods in Medical Research. DOI: 10.1177/0962280219842890
#'
#' Sottile G., Cilluffo G., Muggeo V.M.R. (2019)
#' \emph{The R package islasso: estimation and hypothesis testing in lasso regression}.
#' Technical Report. DOI: 10.13140/RG.2.2.16360.11521
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{confint.islasso}}, \code{\link{plot.islasso}}, \code{\link{predict.islasso}},
#'          \code{\link{summary.islasso}}, \code{\link{is.control}}, \code{\link{aic.islasso}},
#'          \code{\link{anova.islasso}}, \code{\link{islasso.path}}, \code{\link{simulXy}}
#'
#' @keywords models regression
#'
#' @examples
#' n <- 100; p <- 100
#'
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso(y ~ ., data = sim1$data, family = gaussian())
#'
#' summary(o, pval = 0.05)
#' coef(o)
#' fitted(o)
#' predict(o, type="response")
#' plot(o)
#' residuals(o)
#' deviance(o)
#' AIC(o)
#' logLik(o)
#'
#' \dontrun{
#' # for the interaction
#' o <- islasso(y ~ X1 * X2, data = sim1$data, family = gaussian())
#'
#' ##### binomial ######
#' beta <- c(c(1,1,1), rep(0, p-3))
#' sim2 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1,
#'                 size = 100, family = binomial())
#' o2 <- islasso(cbind(y.success, y.failure) ~ .,
#'               data = sim2$data, family = binomial())
#' summary(o2, pval = 0.05)
#'
#' ##### poisson ######
#' beta <- c(c(1,1,1), rep(0, p-3))
#' sim3 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1,
#'                 family = poisson())
#' o3 <- islasso(y ~ ., data = sim3$data, family = poisson())
#' summary(o3, pval = 0.05)
#'
#' ##### Gamma ######
#' beta <- c(c(1,1,1), rep(0, p-3))
#' sim4 <- simulXy(n = n, p = p, beta = beta, interc = -1, seed = 1,
#'                 dispersion = 0.1, family = Gamma(link = "log"))
#' o4 <- islasso(y ~ ., data = sim4$data, family = Gamma(link = "log"))
#' summary(o4, pval = 0.05)
#' }
#'
#' @export
islasso <- function(formula, family=gaussian, lambda, alpha=1, data, weights, subset, offset,
                    unpenalized, contrasts = NULL, control = is.control()){
  this.call <- match.call()

  # Consolidate data environment check
  if(missing(data)) data <- environment(formula)

  # More efficient model frame creation
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"), names(mf), 0L)
  ioff <- m[5] != 0  # Simplified logical check
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")

  # Store offset before evaluation if it exists
  if(ioff) off <- mf$offset

  # Evaluate model frame once
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  # Check for empty model early to avoid unnecessary computation
  if(is.empty.model(mt)) stop("Model matrix is empty!")

  # Extract response and create model matrix
  y <- model.response(mf, "any")
  X <- model.matrix(mt, mf, contrasts)

  # More efficient factor/character column identification
  dataClasses <- attr(mt, "dataClasses")
  termLabels <- attr(mt, "term.labels")
  temp <- which(dataClasses[names(dataClasses) %in% termLabels] %in% c("factor", "character"))
  temp <- which(attr(X, "assign") %in% temp)

  # Optimize offset handling
  if(ioff){
    off_chars <- unlist(lapply(off, as.character))
    noff <- match(off_chars, colnames(X))
    valid_noff <- noff[!is.na(noff) & noff != 0]
    if(length(valid_noff) > 0) X <- X[, -valid_noff, drop = FALSE]
  }

  # Process offset vector
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)){
    if(length(offset) != NROW(y))
      stop(sprintf("number of offsets is %d should equal %d (number of observations)",
                   length(offset), NROW(y)))
  }

  # Consolidated input validation
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights)) {
    if(!is.numeric(weights))
      stop("'weights' must be a numeric vector")
    if(any(weights < 0))
      stop("negative weights not allowed")
  }

  if(alpha < 0 || alpha > 1)  # More efficient range check
    stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")

  # Simplified intercept handling
  intercept <- attr(mt, "intercept") != 0
  if(intercept){
    X <- X[, -1, drop = FALSE]
  }

  # Set attributes more efficiently
  attributes(X)$dataClasses <- temp

  # Consolidated lambda validation
  if(!missing(lambda)){
    if(is.character(lambda) || is.factor(lambda) || lambda < 0)
      stop("'lambda' must be a non-negative numeric value")
  }

  # Default for unpenalized if missing
  if(missing(unpenalized)) unpenalized <- NULL

  # Fit model
  fit <- islasso.fit(X=X, y=y, family=family, lambda=lambda, alpha=alpha, intercept=intercept,
                     weights=weights, offset=offset, unpenalized=unpenalized, control=control)

  # Set model attributes
  fit$model <- mf
  fit$call <- this.call
  fit$formula <- formula
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- contrasts
  fit$xlevels <- .getXlevels(mt, mf)
  class(fit) <- "islasso"

  return(fit)
}

#' @export
islasso.fit <- function(X, y, family=gaussian, lambda, alpha=1, intercept=FALSE, weights=NULL,
                        offset=NULL, unpenalized=NULL, control=is.control()){
  this.call <- match.call()

  # Call the general input checking function
  prep <- .checkinput(X, y, family, alpha, intercept, weights, offset, unpenalized, control)

  # Call the starting point function
  start <- .startpoint(prep$X, prep$y, lambda, alpha, prep$weights, prep$offset,
                       prep$mustart, prep$family, intercept, prep$setting)

  # Create Lambda vector more efficiently
  Lambda <- numeric(prep$nvars)
  Lambda[] <- start$lambda  # Vectorized assignment
  if(!is.null(prep$unpenalized) && length(prep$unpenalized) > 0) {
    Lambda[prep$unpenalized] <- 0
  }

  # Optimize family and link determination
  tempFamily <- prep$tempFamily
  tempLink <- prep$tempLink

  # Define families once
  families <- c("binomial", "poisson", "Gamma")
  quasiFamilies <- c("quasibinomial", "quasipoisson", "Gamma")

  # Check family type more efficiently
  isRegularFamily <- tempFamily %in% families
  isQuasiFamily <- tempFamily %in% quasiFamilies

  if(isRegularFamily || isQuasiFamily) {
    # Determine family index more efficiently
    if(isRegularFamily) {
      famIndex <- match(tempFamily, families)
    } else {
      famIndex <- match(tempFamily, quasiFamilies)
    }

    # Determine link function more efficiently using a list
    linkOptions <- list(
      c("logit", "probit"),
      c("log"),
      c("inverse", "log", "identity")
    )

    link <- match(tempLink, linkOptions[[famIndex]])
    fam <- famIndex
  } else {
    fam <- 0
    link <- 0
  }

  # Call the core function
  out <- .islasso(prep, start, Lambda, fam, link)

  return(out)
}

.checkinput <- function(X, y, family, alpha, intercept, weights, offset, unpenalized, control) {
  # Convert X to matrix and add intercept if needed
  X <- as.matrix(X)
  X <- if(intercept) cbind(1, X) else X
  nobs <- nrow(X)
  nvars <- ncol(X)

  # Handle column names
  nms <- colnames(X)
  if(is.null(nms) && nvars > 0) {
    nms <- paste0("X", seq_len(nvars))
  }
  if(intercept) nms[1] <- "(Intercept)"
  colnames(X) <- nms

  # Validate alpha parameter
  if(alpha < 0 || alpha > 1) {
    stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")
  }

  # Handle unpenalized variables
  if(is.null(unpenalized)) {
    unpenalized <- logical(nvars)
    if(intercept) unpenalized[1] <- TRUE
  } else {
    if(!is.vector(unpenalized)) stop("'unpenalized' is not a vector")
    if(is.list(unpenalized)) stop("'unpenalized' can not be a list")
    if(is.factor(unpenalized)) stop("'unpenalized' can not be a factor")

    if(is.character(unpenalized)) {
      temp_nms <- if(intercept) nms[-1] else nms
      unpenalized_id <- pmatch(unpenalized, temp_nms)
      if(any(is.na(unpenalized_id))) {
        stop(sprintf("the following names are not in colnames(X): %s",
                     paste(unpenalized[is.na(unpenalized_id)], collapse = ", ")))
      }
      unpenalized <- sort(unpenalized_id)
    } else {
      unpenalized <- sort(unpenalized)
    }

    # Validate unpenalized indices
    if(any(abs(unpenalized - round(unpenalized)) > .Machine$double.eps^0.5)) {
      stop("some element of 'unpenalized' is not an integer")
    }
    if(any(unpenalized < 0)) {
      stop("some element of 'unpenalized' is smaller than zero")
    }
    if(any(unpenalized > (nvars - intercept))) {
      stop("some element of 'unpenalized' is greater than the number of columns of the matrix 'X'")
    }

    temp <- logical(nvars - intercept)
    temp[unpenalized] <- TRUE
    if(intercept) temp <- c(TRUE, temp)
    unpenalized <- temp
  }

  # Handle offset and weights
  offset <- if(is.null(offset)) numeric(nobs) else offset

  if(is.null(weights)) {
    weights <- rep(1, nobs)
  } else {
    if(!is.vector(weights)) stop("argument 'weights' is not a vector")
    if(is.list(weights)) stop("argument 'weights' can not be a list")
    if(is.factor(weights)) stop("argument 'weights' can not be a factor")
    if(is.character(weights)) stop("vector 'weights' can not be a character")
    if(length(weights) != nobs) stop(sprintf("the length of the vector 'weights' is not equal to %s", nobs))

    # Handle NaN values in weights
    weights[is.nan(weights)] <- 0
    if(all(weights == 0)) stop("all the entries of the vector 'weights' are equal to zero")
  }

  # Process family
  if(is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if(is.function(family)) {
    family <- do.call(family, args = list(), envir = parent.frame())
  }
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized!")
  }

  tempFamily <- family$family
  tempLink <- family$link

  # Validate and set control parameters
  setting <- control
  if(setting$nfolds < 3) stop("'nfolds' should be greater than 3")
  if(setting$tol <= 0) stop("'tol' should be a non-negative value")
  if(setting$itmax < 0) stop("'itmax' should be a non-negative value")

  setting$estpai <- FALSE
  if(setting$c[1] > 1) {
    stop("'mixing parameter' should be fixed in (0,1) or estimated using a negative value")
  }
  if(setting$c[1] < 0) {
    setting$estpai <- TRUE
    setting$c <- 0.5
  }
  setting$c <- rep(setting$c, nvars)

  # Set sigma2 based on family
  setting$sigma2 <- switch(tempFamily,
                           binomial = 1,
                           quasibinomial = -1,
                           poisson = 1,
                           quasipoisson = -1,
                           setting$sigma2)

  # Validate link function
  okLinks <- c("identity", "logit", "log", "probit", "inverse")
  if(!(tempLink %in% okLinks)) {
    stop(sprintf("%s link not recognized", sQuote(tempLink)))
  }

  # Check compatibility between family and link
  validLinks <- switch(tempFamily,
                       gaussian = "identity",
                       poisson = , quasipoisson = "log",
                       binomial = , quasibinomial = c("logit", "probit"),
                       Gamma = c("log", "identity", "inverse"),
                       okLinks)

  if(!(tempLink %in% validLinks)) {
    stop(sprintf("The %s family does not accept the link %s.\nThe accepted link(s): %s",
                 sQuote(tempFamily), sQuote(tempLink),
                 paste(sQuote(validLinks), collapse = ", ")))
  }

  # Check response variables
  if((is.character(y) || is.factor(y) || is.matrix(y)) &&
     !tempFamily %in% c("binomial", "quasibinomial")) {
    stop(sprintf("The %s family does not accept %s as response variable",
                 sQuote(tempFamily), ifelse(is.character(y), "a character",
                                            ifelse(is.factor(y), "a factor", "a matrix"))))
  }

  # Initialize from family
  n <- rep.int(1, nobs)
  mustart <- NULL

  # Process based on family type
  if(tempFamily == "gaussian") {
    mustart <- y
  } else if(tempFamily %in% c("binomial", "quasibinomial")) {
    if(NCOL(y) == 1) {
      if(is.factor(y)) y <- y != levels(y)[1L]
      y[weights == 0] <- 0
      if(any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")

      mustart <- (weights * y + 0.5) / (weights + 1)
      m <- weights * y
      if(tempFamily == "binomial" && any(abs(m - round(m)) > 0.001)) {
        warning(sprintf("non-integer #successes in a %s glm!", "binomial"))
      }
    } else if(NCOL(y) == 2) {
      if(tempFamily == "binomial" && any(abs(y - round(y)) > 0.001)) {
        warning(sprintf("non-integer counts in a %s glm!", "binomial"))
      }

      y1 <- y[, 1L]
      n <- y1 + y[, 2L]
      y <- y1 / n
      y[n == 0] <- 0
      weights <- weights * n
      mustart <- (n * y + 0.5) / (n + 1)
    } else {
      stop(sprintf("for the '%s' family, y must be a vector of 0 and 1's or a 2 column matrix where col 1 is no. successes and col 2 is no. failures", "binomial"))
    }
  } else if(tempFamily %in% c("poisson", "quasipoisson")) {
    if(any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    mustart <- y + 0.1
  } else if(tempFamily == "Gamma") {
    if(any(y <= 0)) stop("non-positive values not allowed for the 'Gamma' family")
    mustart <- y
  }

  # Ensure response is a vector
  y <- drop(y)

  return(list(
    y = y,
    X = X,
    intercept = intercept,
    alpha = alpha,
    mustart = mustart,
    weights = weights,
    offset = offset,
    unpenalized = unpenalized,
    nobs = nobs,
    nvars = nvars,
    n = n,
    nms = nms,
    family = family,
    tempFamily = tempFamily,
    tempLink = tempLink,
    setting = setting
  ))
}

.startpoint <- function(X, y, lambda, alpha, weights, offset, mustart, family, intercept, setting) {
  nobs <- nrow(X)
  nvars <- ncol(X)

  # Determine family type
  tempFamily <- family$family
  if(tempFamily %in% c("quasibinomial", "quasipoisson")) {
    tempFamily <- sub("^quasi", "", tempFamily)
  } else if(!tempFamily %in% c("gaussian", "binomial", "poisson")) {
    tempFamily <- family
  }

  # Initialize coefficient estimates
  if(is.null(setting$b0)) {
    if((nvars - intercept) < 2) {
      if(missing(lambda)) stop("Insert a positive value for lambda")
      est <- rep(0.1, nvars)
      interval <- NULL
    } else {
      # Prepare data for glmnet
      x <- as.matrix(X)
      y2 <- y
      weights2 <- weights

      # Special handling for binomial family
      if(family$family %in% c("binomial", "quasibinomial")) {
        y2 <- weights * cbind(1-y, y)
        weights2 <- rep(1, nobs)
      }

      # Remove intercept from predictors if needed
      if(intercept) x <- x[, -1, drop = FALSE]

      # Standardize predictors if requested
      if(setting$stand) {
        x_mean <- colMeans(x)
        x_centered <- sweep(x, 2, x_mean, "-")
        x_sd <- apply(x_centered, 2, function(x) sqrt(sum(x^2) / nobs))
        x <- sweep(x_centered, 2, x_sd, "/")
      }

      # Fit model with glmnet
      if(missing(lambda)) {
        obj <- suppressWarnings(
          cv.glmnet(
            x = x,
            y = y2,
            family = tempFamily,
            nfolds = setting$nfolds,
            standardize = FALSE,
            intercept = intercept,
            offset = offset,
            weights = weights2,
            alpha = alpha
          )
        )
        lambda <- obj$lambda.min * nobs
        est <- as.vector(coef(obj, s = "lambda.min"))
      } else {
        obj <- suppressWarnings(
          glmnet(
            x = x,
            y = y2,
            family = tempFamily,
            alpha = alpha,
            weights = weights2,
            standardize = FALSE,
            intercept = intercept,
            offset = offset
          )
        )
        est <- as.vector(coef(obj, s = lambda / nobs))
      }

      # Adjust coefficients
      if(!intercept) est <- est[-1]

      # Get lambda interval for potential reuse
      interval <- range(rev(obj$lambda)) * nobs
    }
  } else {
    if(missing(lambda)) stop("Insert a positive value for lambda")
    est <- setting$b0
    interval <- NULL
  }

  # Prepare covariance and other return values
  covar <- if(is.null(setting$V0)) diag(0.01, nvars) else setting$V0
  se <- sqrt(diag(covar))

  # Calculate fitted values
  eta <- family$linkfun(mustart) + offset
  mu <- family$linkinv(eta)
  residuals <- (y - mu) / family$mu.eta(eta)

  return(list(
    fam = family$family,
    lambda = lambda,
    interval = interval,
    beta = est,
    covar = covar,
    se = se,
    eta = eta,
    mu = mu,
    residuals = residuals
  ))
}

.islasso <- function(prep, start, Lambda, fam, link) {
  # --- Estrai parametri principali ---
  X <- prep$X
  storage.mode(X) <- "double"
  y <- prep$y
  storage.mode(y) <- "double"
  alpha <- as.double(prep$alpha)
  offset <- prep$offset
  storage.mode(offset) <- "double"
  weights <- prep$weights
  storage.mode(weights) <- "double"
  intercept <- as.integer(prep$intercept)
  nobs <- as.integer(prep$nobs)
  nvars <- as.integer(prep$nvars)
  tempFamily <- prep$tempFamily
  family <- prep$family

  beta <- start$beta
  beta[beta == 0] <- 1e-5
  storage.mode(beta) <- "double"

  se <- start$se
  storage.mode(se) <- "double"
  covar <- start$covar
  storage.mode(covar) <- "double"

  # --- Parametri di controllo ---
  setting <- prep$setting
  Lambda <- Lambda
  storage.mode(Lambda) <- "double"
  setting$c <- setting$c
  storage.mode(setting$c) <- "double"
  setting$estpai <- as.integer(setting$estpai)
  sigma2 <- as.double(setting$sigma2)
  itmax <- as.integer(setting$itmax)
  tol <- as.double(setting$tol)
  trace <- as.integer(setting$trace)
  setting$adaptive <- as.integer(setting$adaptive)
  setting$stand <- as.integer(setting$stand)

  # --- Invoca funzione C++ ottimizzata ---
  # fit <- if (tempFamily == "gaussian") {
  #   islasso_cpp(X, y, Lambda, alpha, beta, se, covar, sigma2,
  #               pi = setting$c, weights, offset, setting$estpai,
  #               stand = setting$stand, intercept, itmax, itmaxse = itmax,
  #               tol, adaptive = setting$adaptive, trace = setting$trace)
  # } else {
  #   islasso_glm_cpp(X, y, Lambda, alpha, beta, se, covar, sigma2,
  #                   pi = setting$c, weights, offset, fam, link, estpi = setting$estpai,
  #                   stand = setting$stand, intercept, itmax, itmaxse = itmax, tol,
  #                   adaptive = setting$adaptive, trace = setting$trace)
  # }
  #
  # fit <- lapply(fit, drop)

  fit <- if (tempFamily == "gaussian") {
    .Fortran(C_islasso, X = X, y = y, n = nobs, p = nvars, theta = beta, se = se,
             cov = covar, lambda = Lambda, alpha = alpha, pi = setting$c, estpi = setting$estpai,
             itmax = itmax, itmaxse = itmax, tol = tol, sigma2 = sigma2, trace = trace,
             adaptive = setting$adaptive, offset = offset, conv = integer(1), stand = setting$stand,
             intercept = intercept, eta = y, mu = y, varmu = y, mu_eta_val = y, w = y, res = y, dev = double(1),
             weights = weights, hi = beta, edf = double(1), xtw = matrix(0.0, nrow = nvars, ncol = nobs),
             xtx = covar, grad = beta, hess = covar, invH = covar, pen = beta)
  } else {
    .Fortran(C_islasso_glm, X = X, y = y, n = nobs, p = nvars, theta = beta, se = se,
             cov = covar, lambda = Lambda, alpha = alpha, pi = setting$c, estpi = setting$estpai,
             itmax = itmax, itmaxse = itmax, tol = tol, sigma2 = sigma2, trace = trace,
             adaptive = setting$adaptive, offset = offset, conv = integer(1), stand = setting$stand,
             intercept = intercept, eta = y, mu = y, varmu = y, mu_eta_val = y, w = y, res = y, dev = double(1),
             weights = weights, hi = beta, edf = double(1), xtw = matrix(0.0, nrow = nvars, ncol = nobs),
             xtx = covar, grad = beta, hess = covar, invH = covar, pen = beta, fam = fam, link = link)
  }
  # Check convergence
  if(fit$conv == -1) stop("Infinite values attained, try to change lambda value!!")
  if(fit$conv == 1) warning("Maximum number of iterations attained!!")
  if(fit$conv == 2) stop("Safe exit from ISLASSO algorithm after an inversion problem, try to change lambda value!!")

  setting$c <- fit$pi

  # --- Funzioni di famiglia ---
  dev.resids <- family$dev.resids
  aic_fun <- family$aic
  linkinv <- family$linkinv

  # --- Estrazione risultati ---
  beta <- fit$theta
  se <- fit$se
  covar <- fit$cov
  s2 <- fit$sigma2

  eta <- fit$eta
  mu <- fit$mu
  residuals <- fit$res
  w <- fit$w
  dev <- sum(dev.resids(y, mu, weights))

  rank <- fit$edf
  resdf <- nobs - rank
  aic.model <- if (nobs > nvars && tempFamily %in% c("gaussian", "binomial", "poisson", "Gamma"))
    aic_fun(y, prep$n, mu, weights, dev) + 2 * rank else dev + 2 * rank

  # --- Matrici derivate ---
  XtW <- fit$xtw
  XtX <- fit$xtx
  P <- fit$pen
  if (any(prep$unpenalized)) P[prep$unpenalized] <- 0
  gradient <- fit$grad
  hi <- fit$hi
  H <- fit$hess
  invH <- fit$invH

  # --- Calcola QR e struttura R ---
  z <- eta - offset + residuals
  design <- rbind(X * sqrt(w), sqrt(Lambda) * diag(sqrt(P / beta)))
  response <- c(z * sqrt(w), rep(0, nvars))
  fit$qr <- .lm.fit(x = design, y = response)

  nr <- min(nobs, nvars)
  Rmat <- diag(nvars)
  Rmat[seq_len(nr), ] <- fit$qr$qr[seq_len(nr), seq_len(nvars)]
  Rmat[row(Rmat) > col(Rmat)] <- 0

  fit$effects <- fit$qr$effects[1:nobs]

  # --- Modello nullo ---
  wtdmu <- if (intercept) sum(weights * y) / sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  nulldf <- nobs - as.integer(intercept)

  # --- Etichette e nomi ---
  nms <- prep$nms
  names(gradient) <- names(se) <- names(beta) <- colnames(XtX) <-
    rownames(XtX) <- colnames(covar) <- rownames(covar) <- nms

  internal <- list(
    n = nobs, p = nvars, lambda.seq = fit$lambda,
    XtW = XtW, XtX = XtX, invH = invH, vcov = s2 * covar,
    gradient = gradient, hessian = H, hi = hi,
    intercept = intercept, unpenalized = prep$unpenalized,
    fam = fam, link = link, nms = nms,
    estc = setting$estpai, lmbd.interval = start$interval
  )

  out <- list(
    coefficients = beta, se = se, dispersion = s2,
    residuals = residuals, fitted.values = mu, effects = fit$effect,
    R = Rmat, rank = rank,
    qr = structure(fit$qr[c("qr", "qraux", "pivot", "tol", "rank")], class = "qr"),
    family = family, linear.predictors = eta, deviance = dev, aic = aic.model,
    null.deviance = nulldev, iter = fit$itmax, weights = w, prior.weights = weights,
    df.residual = resdf, df.null = nulldf, y = y, converged = fit$conv,
    model = NULL, call = NULL, formula = NULL, terms = NULL,
    data = NULL, offset = offset, contrasts = NULL, control = setting,
    internal = internal, xlevels = NULL, lambda = start$lambda, alpha = alpha
  )

  out
}

#' Control Settings for islasso Model Fitting
#'
#' Auxiliary function used to configure and customize the fitting process of \code{\link{islasso}} models.
#'
#' @param sigma2 Numeric. Fixed value of the dispersion parameter. If \code{-1} (default), it is estimated from data.
#' @param tol Numeric. Tolerance level to declare convergence. Default is \code{1e-5}.
#' @param itmax Integer. Maximum number of iterations. Default is \code{1000}.
#' @param stand Logical. If \code{TRUE} (default), standardizes covariates before fitting. Returned coefficients remain on the original scale.
#' @param trace Integer. Controls verbosity of the iterative procedure:
#'   \itemize{
#'     \item \code{0} - no printing,
#'     \item \code{1} - compact printing,
#'     \item \code{2} - detailed printing,
#'     \item \code{3} - compact printing with Fisher scoring info (only for GLM).
#'   }
#' @param nfolds Integer. Number of folds for CV if \code{lambda} is missing in \code{islasso}. Defaults to \code{5}.
#' @param seed Optional. Integer seed for reproducibility in cross-validation.
#' @param adaptive Logical. If \code{TRUE}, fits an adaptive LASSO. (Experimental)
#' @param g Numeric in \code{[0,1]}. Governs BIC selection: \code{g = 0} is standard BIC; \code{g = 0.5} is extended BIC.
#' @param b0 Optional. Starting values for regression coefficients. If \code{NULL}, uses \code{glmnet} estimates.
#' @param V0 Optional. Initial covariance matrix. Defaults to identity matrix if \code{NULL}.
#' @param c Numeric. Controls the weight in the induced smoothed LASSO. Default is \code{0.5}; use \code{-1} to recompute at every iteration.
#'
#' @return A list of control parameters for use in \code{\link{islasso}}.
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso}}
#'
#' @export
is.control <- function(
    sigma2 = -1,         # Initial error variance (negative means estimate from data)
    tol = 1E-05,         # Convergence tolerance
    itmax = 1E+3,        # Maximum number of iterations
    stand = TRUE,        # Standardize predictors
    trace = 0,           # Verbosity level (0=none)
    nfolds = 5,          # Number of cross-validation folds
    seed = NULL,         # Random seed for reproducibility
    adaptive = FALSE,    # Use adaptive lasso
    g = 0.5,             # Gamma parameter for eBIC
    b0 = NULL,           # Prior mean for coefficients
    V0 = NULL,           # Prior variance for coefficients
    c = 0.5              # Shrinkage parameter
) {
  # Return a list of control parameters
  list(
    sigma2 = sigma2,
    tol = tol,
    itmax = itmax,
    trace = trace,
    stand = stand,
    nfolds = nfolds,
    seed = seed,
    adaptive = adaptive,
    g = g,
    b0 = b0,
    V0 = V0,
    c = c
  )
}


#' Optimization for Lambda Selection
#'
#' Minimizes information criteria to select the optimal tuning parameter \code{lambda} for \code{\link{islasso}} models.
#' Supports AIC, BIC, AICc, GCV, and GIC.
#'
#' @param object Fitted model of class \code{"islasso"}.
#' @param method Criterion to minimize. Options are \code{"AIC"}, \code{"BIC"}, \code{"AICc"}, \code{"GCV"}, \code{"GIC"}.
#' @param interval Numeric vector (length 2) giving lower and upper bounds for \code{lambda} optimization. Optional if \code{object} includes prior cross-validation.
#' @param g Numeric in \code{[0,1]}. Governs BIC generalization: \code{g = 0} is classic BIC, \code{g = 0.5} is extended BIC.
#' @param y Response vector. Required only if \code{object} is missing.
#' @param X Design matrix. Required only if \code{object} is missing.
#' @param intercept Logical. Whether to include intercept in \code{X}. Used if \code{object} is missing.
#' @param family Error distribution. Accepted: \code{gaussian}, \code{binomial}, \code{poisson}. Uses canonical link.
#' @param alpha Elastic-net mixing parameter, \code{0 <= alpha <= 1}. Lasso: \code{alpha = 1}; Ridge: \code{alpha = 0}.
#' @param offset Optional numeric vector. Adds known linear predictor component.
#' @param weights Optional weights for observations. Defaults to 1.
#' @param unpenalized Logical vector indicating variables to exclude from penalization.
#' @param control List of control parameters. See \code{\link{is.control}}.
#' @param trace Logical. If \code{TRUE}, prints progress of optimization. Default is \code{TRUE}.
#'
#' @details
#' Instead of using cross-validation, this function selects the best \code{lambda} by minimizing criteria like AIC or BIC.
#' Degrees of freedom are computed as the trace of the hat matrix (not necessarily an integer).
#'
#' @return Optimal \code{lambda} value as numeric.
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso}}, \code{\link{islasso.fit}}, \code{\link{summary.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 100
#' beta <- c(rep(2, 20), rep(0, p - 20))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso(y ~ ., data = sim1$data, family = gaussian())
#'
#' \dontrun{
#' # Use the evaluation interval of the fit
#' lambda_aic <- aic.islasso(o, method = "AIC")
#'
#' # Overwrites the evaluation interval for lambda
#' lambda_bic <- aic.islasso(o, interval = c(0.1, 30), method = "BIC")
#'
#' # Overwrites the evaluation interval for lambda using eBIC criterion
#' lambda_ebic <- aic.islasso(o, interval = c(0.1, 30), method = "BIC", g = 0.5)
#' }
#'
#' @export
aic.islasso <- function(
    object,              # An islasso object (optional)
    method = c("AIC", "BIC", "AICc", "GCV", "GIC"),  # Information criterion
    interval,            # Lambda search interval
    g = 0,               # Parameter for extended BIC
    y,                   # Response variable
    X,                   # Model matrix
    intercept = FALSE,   # Include intercept
    family = gaussian(), # Model family
    alpha = 1,           # Elastic net mixing parameter (1=lasso, 0=ridge)
    offset,              # Model offset
    weights,             # Observation weights
    unpenalized,         # Variables to leave unpenalized
    control = is.control(), # Control parameters
    trace = TRUE         # Show optimization progress
) {
  # Input validation and parameter extraction
  if (missing(object)) {
    if (missing(interval))
      stop("Please specify an interval to search for a minimum")
    if (missing(y))
      stop("Model response is missing")
    if (missing(X))
      stop("Model matrix is missing")

    n <- NROW(X)
    p <- NCOL(X)

    if (alpha < 0 || alpha > 1)
      stop("Alpha parameter must be in [0,1]")

    # Set default values
    if (missing(offset))
      offset <- rep(0, n)
    if (missing(weights))
      weights <- rep(1, n)
    if (missing(unpenalized))
      unpenalized <- NULL
  } else {
    # Extract information from existing object
    if (missing(interval))
      interval <- object$internal$lmbd.interval
    if (is.null(interval))
      stop("Please specify an interval to search for a minimum")

    n <- object$internal$n
    p <- object$internal$p
    X <- model.matrix(object)
    y <- object$y
    intercept <- object$internal$intercept
    nms <- object$internal$nms
    unpenalized <- object$internal$unpenalized
    alpha <- object$alpha

    # Handle intercept
    if (intercept) {
      X <- X[, -1, drop = FALSE]
      nms <- nms[-1]
      unpenalized <- unpenalized[-1]
    }

    unpenalized <- nms[unpenalized]
    family <- object$family
    offset <- object$offset
    weights <- object$weights
    control <- object$control

    # Check control parameters
    if (control$estpai)
      control$c <- -1
    if (length(control$c) != 1L)
      control$c <- max(control$c)
  }

  # Disable tracing in the control object to avoid excessive output
  control$trace <- 0

  # Match the method argument
  method <- match.arg(method)

  # Define formulas for different information criteria
  k <- switch(method,
              "AIC" = "ll + 2 * df",
              "BIC" = "ll + (log(n) + 2 * g * log(p)) * df",
              "AICc" = "ll + 2 * n * df * (df + 1) / (n - df - 1) + 2 * df",
              "GCV" = "ll / ((1 - df / n)^2)",
              "GIC" = "ll + log(log(n)) * log(p) * df")

  # Rename method if using extended BIC
  if (method == "BIC" && g != 0)
    method <- "eBIC"

  # Validate g parameter
  if (g < 0 || g > 1)
    stop("Gamma parameter must be set in [0,1]")

  # Print optimization method if trace is enabled
  if (trace)
    cat(paste0("\nOptimization through ", method, "\n\n"))

  # Define the objective function for optimization
  fun <- function(lambda, X, y, alpha, family, intercept, weights, offset, unpenalized, control, k, n, p, trace) {
    # Try to fit the model
    obj <- try(islasso.fit(
      X = X,
      y = y,
      family = family,
      lambda = lambda,
      alpha = alpha,
      intercept = intercept,
      weights = weights,
      offset = offset,
      unpenalized = unpenalized,
      control = control
    ), silent = TRUE)

    # Handle errors in model fitting
    if (inherits(obj, "try-error")) {
      return(Inf)
    } else {
      # Calculate the information criterion
      ll <- obj$aic - 2 * obj$rank
      df <- obj$rank
      temp <- eval(parse(text = k))

      # Print current lambda and criterion value if trace is enabled
      if (trace) {
        cat("lambda = ", formatC(lambda, digits = 4, width = 8, format = "f"),
            method, "= ", formatC(temp, digits = 5, width = 10, format = "f"), "\n")
      }

      return(temp)
    }
  }

  # Optimize to find the best lambda
  opt_result <- optimize(
    f = fun,
    interval = interval,
    X = X,
    y = y,
    alpha = alpha,
    family = family,
    intercept = intercept,
    weights = weights,
    offset = offset,
    unpenalized = unpenalized,
    control = control,
    k = k,
    n = n,
    p = p,
    trace = trace
  )

  lambda.min <- opt_result$minimum

  return(lambda.min)
}
