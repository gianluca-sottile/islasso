#' Induced Smoothed Lasso Regularization Path
#'
#' Fits a sequence of penalized regression models using the Induced Smoothing Lasso approach over a grid of lambda values.
#' Supports elastic-net penalties and generalized linear models: Gaussian, Binomial, Poisson, and Gamma.
#'
#' @aliases islasso.path print.islasso.path islasso.path.fit coef.islasso.path
#'          deviance.islasso.path fitted.islasso.path residuals.islasso.path
#'          logLik.islasso.path print.logLik.islasso.path model.matrix.islasso.path
#'
#' @param formula Model formula of type \code{response ~ predictors}.
#' @param family Response distribution. Supported families: \code{gaussian()}, \code{binomial()}, \code{poisson()}, \code{Gamma()}.
#' @param lambda Optional numeric vector of lambda values. If not provided, a sequence is automatically generated.
#' @param nlambda Integer. Number of lambda values to generate if \code{lambda} is missing. Default is \code{100}.
#' @param lambda.min.ratio Smallest lambda as a fraction of \code{lambda.max}. Default: \code{1e-2} if \code{nobs < nvars}, else \code{1e-3}.
#' @param alpha Elastic-net mixing parameter: \code{alpha = 1} is lasso, \code{alpha = 0} is ridge.
#' @param data Data frame containing model variables.
#' @param weights Optional observation weights.
#' @param subset Optional logical or numeric vector to subset observations.
#' @param offset Optional vector of prior known components for the linear predictor.
#' @param contrasts Optional contrast settings for factor variables.
#' @param unpenalized Optional vector of variable names or indices excluded from penalization.
#' @param control A list of control parameters via \code{\link{is.control}}.
#'
#' @details
#' This function fits a regularization path of models using the induced smoothing paradigm, replacing the non-smooth L1 penalty with a differentiable surrogate.
#' Standard errors are returned for all lambda points, allowing for Wald-based hypothesis testing. The regularization path spans a range of lambda values,
#' either user-defined or automatically computed.
#'
#' @return A list with components:
#' \item{call}{Matched function call.}
#' \item{Info}{Matrix with diagnostics: lambda, deviance, degrees of freedom, dispersion, iterations, convergence status.}
#' \item{GoF}{Model goodness-of-fit metrics: AIC, BIC, AICc, GCV, GIC, eBIC.}
#' \item{Coef}{Matrix of coefficients across lambda values.}
#' \item{SE}{Matrix of standard errors.}
#' \item{Weights}{Matrix of mixing weights for the smoothed penalty.}
#' \item{Gradient}{Matrix of gradients for the smoothed penalty.}
#' \item{Linear.predictors, Fitted.values, Residuals}{Matrices of fitted quantities across the path.}
#' \item{Input}{List of input arguments and design matrix.}
#' \item{control, formula, model, terms, data, xlevels, contrasts}{Standard model components.}
#'
#' @references
#' Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2019).
#' \emph{The Induced Smoothed Lasso: A practical framework for hypothesis testing in high dimensional regression}.
#' Statistical Methods in Medical Research. DOI: 10.1177/0962280219842890
#'
#' Sottile G., Cilluffo G., Muggeo V.M.R. (2019).
#' \emph{The R package islasso: estimation and hypothesis testing in lasso regression}.
#' Technical Report. DOI: 10.13140/RG.2.2.16360.11521
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso}}, \code{\link{summary.islasso.path}}, \code{\link{coef.islasso.path}},
#'          \code{\link{predict.islasso.path}}, \code{\link{GoF.islasso.path}}
#'
#' @examples
#' n <- 100; p <- 30; p1 <- 10  # number of nonzero coefficients
#'
#' beta.veri <- sort(round(c(seq(.5, 3, length.out = p1/2),
#'                          seq(-1, -2, length.out = p1/2)), 2))
#' beta <- c(beta.veri, rep(0, p - p1))
#' sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso.path(y ~ ., data = sim1$data,
#'                   family = gaussian(), nlambda = 30L)
#' o
#'
#' summary(o, lambda = 10, pval = 0.05)
#' coef(o, lambda = 10)
#' fitted(o, lambda = 10)
#' predict(o, type = "response", lambda = 10)
#' plot(o, yvar = "coef")
#' residuals(o, lambda = 10)
#' deviance(o, lambda = 10)
#' logLik(o, lambda = 10)
#' GoF.islasso.path(o)
#'
#' \dontrun{
#' ##### binomial ######
#' beta <- c(1, 1, 1, rep(0, p - 3))
#' sim2 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1,
#'                 size = 100, family = binomial())
#' o2 <- islasso.path(cbind(y.success, y.failure) ~ ., data = sim2$data,
#'                    family = binomial(), lambda = seq(0.1, 100, l = 50L))
#' temp <- GoF.islasso.path(o2)
#' summary(o2, pval = 0.05, lambda = temp$lambda.min["BIC"])
#'
#' ##### poisson ######
#' beta <- c(1, 1, 1, rep(0, p - 3))
#' sim3 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1,
#'                 family = poisson())
#' o3 <- islasso.path(y ~ ., data = sim3$data, family = poisson(), nlambda = 30L)
#' temp <- GoF.islasso.path(o3)
#' summary(o3, pval = 0.05, lambda = temp$lambda.min["BIC"])
#'
#' ##### Gamma ######
#' beta <- c(1, 1, 1, rep(0, p - 3))
#' sim4 <- simulXy(n = n, p = p, beta = beta, interc = -1, seed = 1,
#'                 family = Gamma(link = "log"))
#' o4 <- islasso.path(y ~ ., data = sim4$data, family = Gamma(link = "log"),
#'                    nlambda = 30L)
#' temp <- GoF.islasso.path(o4)
#' summary(o4, pval = .05, lambda = temp$lambda.min["BIC"])
#' }
#'
#' @export
islasso.path <- function(formula, family = gaussian(), lambda = NULL, nlambda = 100,
                         lambda.min.ratio = ifelse(nobs < nvars, 1E-2, 1E-04),
                         alpha = 1, data, weights, subset, offset, contrasts = NULL,
                         unpenalized, control = is.control()) {
  this.call <- match.call()

  # Handle missing data argument
  if(missing(data)) data <- environment(formula)

  # Create model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"), names(mf), 0L)
  ioff <- if(m[5] != 0) TRUE else FALSE
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  if(ioff) off <- mf$offset
  mf <- eval(mf, parent.frame())

  # Extract model components
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")

  # Create model matrix
  x <- if(!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  } else {
    stop("Model matrix is empty!")
  }

  # Get dimensions early for use in lambda.min.ratio
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])

  # Factor handling
  temp <- which(attr(mt, "dataClasses")[names(attr(mt, "dataClasses")) %in%
                                          attr(mt, "term.labels")] %in% c("factor", "character"))
  temp <- which(attr(x, "assign") %in% temp)

  # Handle offset
  if(ioff){
    noff <- match(unlist(lapply(off, as.character)), colnames(x))
    if(!all(is.na(noff))) x <- x[, -noff[which(noff!=0)]]
  }
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)){
    if(length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }

  # Handle weights
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")

  # Validate alpha
  if(alpha > 1 | alpha < 0)
    stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")

  # Handle intercept
  if(attr(mt,"intercept") == 0){
    intercept <- FALSE
  } else {
    intercept <- TRUE
    x <- x[, -1, drop = FALSE]
  }

  attributes(x)$dataClasses <- temp

  # Handle unpenalized variables
  if(missing(unpenalized)) unpenalized <- NULL

  # Fit the model
  fit <- islasso.path.fit(X = x, y = y, family = family, lambda = lambda, nlambda = nlambda,
                          lambda.min.ratio = lambda.min.ratio, alpha = alpha, intercept = intercept,
                          weights = weights, offset = offset, unpenalized = unpenalized, control = control)

  # Add necessary components to fit
  fit$call <- this.call
  fit$formula <- formula
  fit$model <- mf
  fit$terms <- mt
  fit$data <- data
  fit$contrasts <- contrasts
  fit$xlevels <- .getXlevels(mt, mf)
  class(fit) <- c('islasso.path', "islasso")

  return(fit)
}

#' @export
islasso.path.fit <- function(X, y, family = gaussian(), lambda, nlambda,
                             lambda.min.ratio, alpha = 1, intercept = FALSE,
                             weights = NULL, offset = NULL, unpenalized = NULL,
                             control = is.control()){
  this.call <- match.call()

  # Check inputs - this is a core step, no optimization needed as it's called once
  prep <- checkinput.islasso.path(X, y, family, lambda, nlambda, lambda.min.ratio,
                                  alpha, intercept, weights, offset, unpenalized, control)

  # Get starting point - also a one-time call
  start <- startpoint.islasso.path(prep$X, prep$y, prep$lambda, alpha, prep$weights,
                                   prep$offset, prep$mustart, prep$family, intercept, prep$setting)

  # Determine family and link
  fam <- 0
  link <- 0
  fam_names <- c("binomial", "poisson", "Gamma")
  fam_names2 <- c("quasibinomial", "quasipoisson", "Gamma")

  if((a <- prep$tempFamily %in% fam_names) | (prep$tempFamily %in% fam_names2)){
    fam <- ifelse(a, match(prep$tempFamily, fam_names), match(prep$tempFamily, fam_names2))

    link_options <- switch(fam,
                           "1" = c("logit", "probit"),
                           "2" = c("log"),
                           "3" = c("inverse", "log", "identity"))

    link <- match(prep$tempLink, link_options)
  }

  out <- islasso.path.fit.glm(prep, start, start$lambda, fam, link)

  out$control <- prep$setting
  return(out)
}

checkinput.islasso.path <- function(X, y, family, lambda, nlambda, lambda.min.ratio, alpha, intercept,
                                    weights, offset, unpenalized, control) {
  # Helper function for lambda calculation - moved inside to avoid namespace clutter
  get_start <- function(x, y, weights, family, intercept, offset, exclude, alpha) {
    nobs <- nrow(x)
    nvars <- ncol(x)
    is.offset <- !all(offset == 0)

    # Calculate mu
    if (intercept) {
      if (is.offset) {
        suppressWarnings(tempfit <- glm(y ~ 1, family = family, weights = weights, offset = offset))
        mu <- tempfit$fitted.values
      } else {
        mu <- rep(weighted.mean(y, weights), times = nobs)
      }
    } else {
      mu <- family$linkinv(offset)
    }

    # Calculate penalty components
    nulldev <- sum(family$dev.resids(y, mu, weights))
    ju <- rep(1, nvars)
    ju[exclude] <- 0

    # More efficient calculation of g
    r <- y - mu
    eta <- family$linkfun(mu)
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)

    # Normalize weights
    w_norm <- weights / sum(weights)

    # Vectorized computation
    rv <- r/v * m.e * w_norm
    g <- abs(colSums(rv * x)) * ju

    lambda_max <- max(g) / max(alpha, 0.001)
    return(list(nulldev = nulldev, mu = mu, lambda_max = lambda_max * nobs))
  }

  # Convert X to matrix (efficient handling)
  X <- as.matrix(X)
  X <- if(intercept) cbind(1, X) else X
  nobs <- nrow(X)
  nvars <- ncol(X)

  # Create column names if needed
  nms <- colnames(X)
  if(is.null(nms) & ncol(X) != 0) nms <- paste0("X", 1:nvars)
  if(intercept) nms[1] <- "(Intercept)"
  colnames(X) <- nms

  # Check alpha value
  if(alpha > 1 | alpha < 0)
    stop("alpha must be in [0, 1] (0 for ridge penalty, 1 for lasso penalty)")

  # Handle unpenalized variables efficiently
  if(is.null(unpenalized)){
    unpenalized <- logical(nvars)
    if(intercept) unpenalized[1] <- TRUE
  } else {
    if(!is.vector(unpenalized)) stop("'unpenalized' is not a vector")
    if(is.list(unpenalized)) stop("'unpenalized' can not be a list")
    if(is.factor(unpenalized)) stop("'unpenalized' can not be a factor")

    if(is.character(unpenalized)){
      temp_nms <- if(intercept) nms[-1] else nms
      unpenalized_id <- pmatch(unpenalized, temp_nms)
      if(any(is.na(unpenalized_id)))
        stop(gettextf("the following names are not in colnames(X): %s",
                      paste(unpenalized[is.na(unpenalized_id)], collapse = ", ")))
      unpenalized_ids <- sort(unpenalized_id)

      # Create logical vector directly
      temp <- logical(nvars - intercept)
      temp[unpenalized_ids] <- TRUE
      if(intercept) temp <- c(TRUE, temp)
      unpenalized <- temp
    } else {
      # Validate numeric inputs
      unpenalized <- sort(unpenalized)
      if(any(abs(unpenalized - round(unpenalized)) > .Machine$double.eps^0.5))
        stop("some element of 'unpenalized' is not an integer")
      if(any(unpenalized <= 0))
        stop("some element of 'unpenalized' is smaller than zero")
      if(any(unpenalized > (nvars-1*intercept)))
        stop("some element of 'unpenalized' is greater than the number of columns of the matrix 'X'")

      # Create logical vector directly
      temp <- logical(nvars - intercept)
      temp[unpenalized] <- TRUE
      if(intercept) temp <- c(TRUE, temp)
      unpenalized <- temp
    }
  }

  # Handle offset and weights efficiently
  if(is.null(offset)) offset <- numeric(nobs)

  if(is.null(weights)){
    weights <- rep(1, nobs)
  } else {
    if(!is.vector(weights)) stop("argument 'weights' is not a vector")
    if(is.list(weights)) stop("argument 'weights' can not be a list")
    if(is.factor(weights)) stop("argument 'weights' can not be a factor")
    if(is.character(weights)) stop("vector 'weights' can not be a character")
    if(length(weights) != nobs) stop("the length of the vector 'weights' is not equal to ", sQuote(nobs))

    # Handle NaN weights more efficiently
    if(any(is.nan(weights))) weights[is.nan(weights)] <- 0
    if(all(weights == 0)) stop("all the entries of the vector 'weights' are equal to zero")
  }

  # Process family object
  if(is.character(family))
    family <- get(family, mode="function", envir=parent.frame())
  if(is.function(family))
    family <- do.call(family, args=list(), envir=parent.frame())
  if(is.null(family$family)){
    print(family)
    stop("'family' not recognized!")
  }

  tempFamily <- family$family
  tempLink <- family$link

  # Configure control settings
  setting <- control
  if(setting$itmax < 0) stop("'itmax' should be a non-negative value")

  # Handle mixing parameter
  if(setting$c[1] > 1) stop("'mixing parameter' should be fixed in (0,1) or estimated using a negative value")
  if(setting$c[1] < 0) {
    setting$c <- rep(1, nvars)
    setting$fix.c <- FALSE
  } else {
    if(length(setting$c) != nvars) setting$c <- rep(setting$c[1], nvars)
    setting$fix.c <- TRUE
  }

  # Validate gamma parameter
  if(setting$g < 0 | setting$g > 1) {
    warning("gamma parameter have to be set in (0,1). Default parameter is 0.5")
    setting$g <- .5
  }

  # Set sigma2 based on family
  family_sigma <- c(binomial = 1, quasibinomial = -1, poisson = 1, quasipoisson = -1)
  if(tempFamily %in% names(family_sigma)) {
    setting$sigma2 <- family_sigma[tempFamily]
  }

  # Check link compatibility more efficiently
  okLinks <- c("identity", "logit", "log", "probit", "inverse")
  if(!(tempLink %in% okLinks))
    stop(gettextf("%s link not recognized", sQuote(tempLink)), domain = NA)

  # More concise link validation
  family_link_map <- list(
    gaussian = "identity",
    binomial = c("logit", "probit"),
    quasibinomial = c("logit", "probit"),
    poisson = "log",
    quasipoisson = "log",
    Gamma = c("log", "identity", "inverse")
  )

  if(tempFamily %in% names(family_link_map)) {
    valid_links <- family_link_map[[tempFamily]]
    if(!(tempLink %in% valid_links)) {
      stop(gettextf("The %s family does not accept the link %s.\n  The accepted link(s): %s",
                    sQuote(tempFamily), sQuote(tempLink),
                    paste(sQuote(valid_links), collapse = ", ")), domain = NA)
    }
  }

  # Check response type compatibility
  if((is.character(y) || is.factor(y) || is.matrix(y)) &&
     !tempFamily %in% c("binomial", "quasibinomial")) {
    stop(gettextf("The %s family does not accept a %s as response variable",
                  sQuote(tempFamily), class(y)[1]), domain = NA)
  }

  # Initialize from family function
  etastart <- 0
  mustart <- NULL

  # Handle initialization for different families
  if(tempFamily == "gaussian"){
    n <- rep.int(1, nobs)
    mustart <- y
  } else if(tempFamily %in% c("binomial", "quasibinomial")) {
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (weights * y + 0.5)/(weights + 1)
      m <- weights * y
      if (tempFamily == "binomial" && any(abs(m - round(m)) > 0.001))
        warning(gettextf("non-integer #successes in a %s glm!", "binomial"), domain = NA)
    }
    else if (NCOL(y) == 2) {
      if (tempFamily == "binomial" && any(abs(y - round(y)) > 0.001))
        warning(gettextf("non-integer counts in a %s glm!", "binomial"), domain = NA)
      n <- (y1 <- y[, 1L]) + y[, 2L]
      y <- y1 / n
      if (any(n0 <- n == 0)) y[n0] <- 0
      weights <- weights * n
      mustart <- (n * y + 0.5)/(n + 1)
    }
    else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures",
                       "binomial"), domain = NA)
  } else if(tempFamily %in% c("poisson", "quasipoisson")) {
    if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  } else if(tempFamily == "Gamma") {
    if (any(y <= 0)) stop("non-positive values not allowed for the 'Gamma' family")
    n <- rep.int(1, nobs)
    mustart <- y
  }

  # Ensure y is a vector (not a matrix)
  y <- drop(y)

  # Lambda sequence generation
  if(is.null(lambda)) {
    lmax <- get_start(X, y, weights, family, intercept, offset, unpenalized, alpha)$lambda_max
    lmin <- lmax * lambda.min.ratio * 1.01
    lambda <- exp(seq(log(lmin), log(lmax), length.out = nlambda))
  } else {
    nlambda <- length(lambda)
  }

  # Return prepared data
  return(list(
    y = y,
    X = X,
    intercept = intercept,
    lambda = lambda,
    nlambda = nlambda,
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

startpoint.islasso.path <- function(X, y, lambda, alpha, weights, offset, mustart, family, intercept, setting) {
  nobs <- nrow(X)
  nvars <- ncol(X)

  # Determine tempFamily efficiently
  tempFamily <- family$family
  if(tempFamily %in% c("quasibinomial", "quasipoisson")) {
    tempFamily <- sub("quasi", "", tempFamily)
  }

  # Initialize beta vector
  if((nvars - intercept) < 2){
    est <- rep(0.1, nvars)
  } else {
    # Copy data to avoid modifying original
    x <- X
    y2 <- y
    weights2 <- weights

    # Special handling for binomial family
    if(family$family %in% c("binomial", "quasibinomial")) {
      y2 <- weights * cbind(1 - y, y)
      weights2 <- rep(1, nobs)
    }

    # Remove intercept for standardization
    if(intercept) x <- x[, -1, drop = FALSE]

    # Standardize predictors if needed
    if(setting$stand) {
      x_mean <- colMeans(x)
      x_centered <- sweep(x, 2, x_mean, "-")
      # More efficient variance calculation
      x_sd <- sqrt(colSums(x_centered^2) / nobs)
      x <- sweep(x_centered, 2, x_sd, "/")
    }

    # Use glmnet for initial estimates (much faster than direct calculation)
    if(tempFamily != "Gamma") {
      obj <- suppressWarnings(glmnet(x = x, y = y2, family = tempFamily, alpha = alpha,
                                     weights = weights2, standardize = FALSE,
                                     intercept = intercept, offset = offset))
    } else {
      obj <- suppressWarnings(glmnet(x = x, y = y2, family = family, alpha = alpha,
                                     weights = weights2, standardize = FALSE,
                                     intercept = intercept, offset = offset))
    }

    # Extract coefficients at appropriate lambda
    est <- as.vector(coef(obj, s = max(min(obj$lambda), min(lambda) / nobs)))

    # Remove intercept if needed
    if(!intercept) est <- est[-1]
  }

  # Calculate interval and other necessary values
  interval <- range(lambda)

  # Initialize covariance matrix
  covar <- diag(0.01, nrow = nvars, ncol = nvars)
  se <- sqrt(diag(covar))

  # Calculate linear predictor, fitted values and residuals
  eta <- family$linkfun(mustart) + offset
  mu <- family$linkinv(eta)
  residuals <- (y - mu) / family$mu.eta(eta)

  # Return list of starting values
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

islasso.path.fit.glm <- function(prep, start, lambda, fam, link) {
  # --- Estrazione e preparazione ---
  X <- prep$X
  storage.mode(X) <- "double"
  y <- prep$y
  storage.mode(y) <- "double"
  offset <- prep$offset
  storage.mode(offset) <- "double"
  weights <- prep$weights
  storage.mode(weights) <- "double"
  family <- prep$family
  unpenalized <- prep$unpenalized
  intercept <- as.integer(prep$intercept)
  nvars <- as.integer(prep$nvars)
  nobs <- as.integer(prep$nobs)
  alpha <- as.double(prep$alpha)
  tempFamily <- prep$tempFamily

  b <- start$beta
  b[b == 0] <- 1e-5
  storage.mode(b) <- "double"
  se <- start$se
  storage.mode(se) <- "double"
  cov.unscaled <- start$covar
  storage.mode(cov.unscaled) <- "double"

  setting <- prep$setting
  c <- as.double(setting$c)
  storage.mode(c) <- "double"
  est.c <- as.integer(!setting$fix.c)
  tol <- as.double(setting$tol)
  maxIter <- as.integer(setting$itmax)
  trace <- as.integer(setting$trace)
  sigma2 <- as.double(setting$sigma2)
  g <- as.double(setting$g)
  adaptive <- as.integer(setting$adaptive)
  stand <- as.integer(setting$stand)

  storage.mode(lambda) <- "double"
  nlambda <- length(lambda)

  gradient <- double(nvars)

  # --- Output containers ---
  outputInfo <- matrix(NA, nrow = nlambda, ncol = 7)
  outputGoF <- matrix(NA, nrow = nlambda, ncol = 6)
  outputCoef <- matrix(NA, nrow = nlambda, ncol = nvars)
  outputSE <- matrix(NA, nrow = nlambda, ncol = nvars)
  outputWeight <- matrix(NA, nrow = nlambda, ncol = nvars)
  outputGrad <- matrix(NA, nrow = nlambda, ncol = nvars)
  outputLinPred <- outputFitted <- outputResid <- NULL
  CONV <- integer(nlambda)

  # --- Funzioni famiglia ---
  linkinv <- family$linkinv
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic_fun <- family$aic
  mu.eta <- family$mu.eta

  # --- Tracciamento ---
  if (trace == 2L) cat("\n\n   Ind|Max  \tlambda     \tdf        \tphi      \tIter\tError")
  if (trace > 0) time0 <- Sys.time()

  for (rep in seq_len(nlambda)) {
    lmbd <- lambda[rep] * (!unpenalized)

    active <- abs(b) > 1e-8

    if (sum(active) == 0 || (intercept && all(which(active) == 1))) break

    # if (tempFamily == "gaussian") {
    #   temp.fit <- islasso_cpp(
    #     X = X[, active, drop = FALSE], y = y, lambda = lmbd[active], alpha = alpha,
    #     theta = b[active], se = se[active], covar = cov.unscaled[active, active], sigma2 = sigma2,
    #     pi = c[active], weights = weights, offset = offset,
    #     estpi = est.c, stand = stand, intercept = intercept,
    #     itmax = maxIter, itmaxse = maxIter,
    #     tol = tol, adaptive = adaptive, trace = 0L
    #   )
    # } else {
    #   temp.fit <- islasso_glm_cpp(
    #     X = X[, active, drop = FALSE], y = y, lambda = lmbd[active], alpha = alpha,
    #     theta = b[active], se = se[active], covar = cov.unscaled[active, active], sigma2 = sigma2,
    #     pi = c[active], weights = weights, offset = offset, fam = fam, link = link,
    #     estpi = est.c, stand = stand, intercept = intercept,
    #     itmax = maxIter, itmaxse = maxIter,
    #     tol = tol, adaptive = adaptive, trace = 0L
    #   )
    # }
    #
    # temp.fit <- lapply(temp.fit, drop)

    temp.fit <- if (tempFamily == "gaussian") {
      .Fortran(C_islasso, X = X[, active, drop = FALSE], y = y, n = nobs, p = sum(active), theta = b[active], se = se[active],
               cov = cov.unscaled[active, active], lambda = lmbd[active], alpha = alpha, pi = c[active], estpi = est.c,
               itmax = maxIter, itmaxse = maxIter, tol = tol, sigma2 = sigma2, trace = 0L,
               adaptive = adaptive, offset = offset, conv = integer(1), stand = stand,
               intercept = intercept, eta = y, mu = y, varmu = y, mu_eta_val = y, w = y, res = y, dev = double(1),
               weights = weights, hi = b[active], edf = double(1), xtw = matrix(0.0, nrow = sum(active), ncol = nobs),
               xtx = cov.unscaled[active, active], grad = b[active], hess = cov.unscaled[active, active], invH = cov.unscaled[active, active], pen = b[active])
    } else {
      .Fortran(C_islasso_glm, X = X[, active, drop = FALSE], y = y, n = nobs, p = sum(active), theta = b[active], se = se[active],
               cov = cov.unscaled[active, active], lambda = lmbd[active], alpha = alpha, pi = c[active], estpi = est.c,
               itmax = maxIter, itmaxse = maxIter, tol = tol, sigma2 = sigma2, trace = 0L,
               adaptive = adaptive, offset = offset, conv = integer(1), stand = stand,
               intercept = intercept, eta = y, mu = y, varmu = y, mu_eta_val = y, w = y, res = y, dev = double(1),
               weights = weights, hi = b[active], edf = double(1), xtw = matrix(0.0, nrow = sum(active), ncol = nobs),
               xtx = cov.unscaled[active, active], grad = b[active], hess = cov.unscaled[active, active], invH = cov.unscaled[active, active], pen = b[active], fam = fam, link = link)
    }

    CONV[rep] <- conv <- temp.fit$conv
    if (conv > 1) next

    iter <- temp.fit$itmax
    err <- temp.fit$tol

    b[active] <- temp.fit$theta
    se[active] <- temp.fit$se
    c[active] <- temp.fit$pi
    cov.unscaled[active, active] <- temp.fit$cov
    gradient[active] <- temp.fit$grad

    # eta <- drop(X %*% b + offset)
    eta <- temp.fit$eta
    # mu <- linkinv(eta)
    mu <- temp.fit$mu
    # mu.eta.val <- mu.eta(eta)
    mu.eta.val <- temp.fit$mu_eta_val
    # varmu <- variance(mu)
    varmu <- temp.fit$varmu
    # w <- (weights * mu.eta.val^2) / varmu
    w <-  temp.fit$w
    # residuals <- (y - mu) / mu.eta.val
    residuals <- temp.fit$res

    # XtW <- t(w * X)
    # XtW <- temp.fit$xtw
    # XtX <- XtW %*% X
    # XtX <- temp.fit$xtx

    # H <- compute_hessian(XtX, b, se, lmbd, c, alpha, 2L)
    # H <- .Fortran(C_hessian, theta = b, se = se, lambda = lmbd, xtx = XtX, pi = c, p = nvars, hess = XtX, alpha = alpha)$hess
    # H <- temp.fit$hess
    # invH <- inv_sym_matrix(H)
    # invH <- .Fortran(C_inv_lapack, n = nvars, a = H, ainv = H, info = integer(1), ipiv = integer(nvars), work = double(nvars*nvars))$ainv
    # invH <- temp.fit$invH

    df <- temp.fit$edf
    dev <- sum(dev.resids(y, mu, weights))
    ll <- dev
    if (nobs > nvars && tempFamily %in% c("gaussian", "binomial", "poisson", "Gamma"))
      ll <- aic_fun(y, prep$n, mu, weights, dev)

    # --- Misure di bontà ---
    aic <- ll + 2 * df
    bic <- ll + log(nobs) * df
    aicc <- ll + 2 * nobs * df * (df + 1) / (nobs - df - 1) + 2 * df
    ebic <- ll + (log(nobs) + 2 * g * log(nvars)) * df
    gcv <- ll / ((1 - df / nobs)^2)
    gic <- ll + log(log(nobs)) * log(nvars) * df

    # --- Tracciamento progressivo ---
    if (trace == 2L) {
      cat(sprintf("\n%4d|%4d\t%10.4f\t%10.4f\t%10.4f\t%4d\t%9.6f",
                  rep, nlambda, lambda[rep], df, temp.fit$sigma2, iter, err))
    } else if (trace == 1L) {
      cat("\r", rep, "|", nlambda, rep(" ", 20))
    }

    # --- Output singolo step ---
    outputInfo[rep, ]   <- c(lambda[rep], df, temp.fit$sigma2, dev, -ll/2, iter, err)
    outputGoF[rep, ]    <- c(aic, bic, aicc, ebic, gcv, gic)
    outputCoef[rep, ]   <- b
    outputSE[rep, ]     <- se
    outputWeight[rep, ] <- c
    outputGrad[rep, ] <- gradient

    outputLinPred <- rbind(outputLinPred, eta)
    outputFitted  <- rbind(outputFitted, mu)
    outputResid   <- rbind(outputResid, residuals)
  }

  # --- Rimuovi modelli non convergenti ---
  valid <- CONV != 1L
  outputInfo     <- outputInfo[valid, , drop = FALSE]
  outputGoF      <- outputGoF[valid, , drop = FALSE]
  outputCoef     <- outputCoef[valid, , drop = FALSE]
  outputSE       <- outputSE[valid, , drop = FALSE]
  outputWeight   <- outputWeight[valid, , drop = FALSE]
  outputGrad     <- outputGrad[valid, , drop = FALSE]
  outputLinPred  <- outputLinPred[valid, , drop = FALSE]
  outputFitted   <- outputFitted[valid, , drop = FALSE]
  outputResid    <- outputResid[valid, , drop = FALSE]

  # --- Etichette ---
  colnames(outputInfo)   <- c("lambda", "df", "phi", "deviance", "logLik", "iter", "err")
  colnames(outputGoF)    <- c("AIC", "BIC", "AICc", "eBIC", "GCV", "GIC")
  colnames(outputCoef)   <- prep$nms
  colnames(outputSE)     <- prep$nms
  colnames(outputWeight) <- prep$nms
  colnames(outputGrad)   <- prep$nms

  if (trace > 0) cat("\n\n Executed in", Sys.time() - time0, "\n")

  list(
    Coef = na.omit(outputCoef), SE = na.omit(outputSE), Residuals = outputResid,
    Fitted.values = outputFitted, Info = na.omit(outputInfo), GoF = na.omit(outputGoF),
    Linear.predictors = outputLinPred, prior.weights = weights, Weight = na.omit(outputWeight),
    Gradient = na.omit(outputGrad), y = y, model = NULL, call = match.call(), family = family, formula = formula, terms = NULL,
    data = NULL, offset = offset, control = NULL, contrasts = NULL, xlevels = NULL, alpha = alpha,
    Input = list(n = nobs, p = nvars, intercept = intercept, unpenalized = unpenalized)
  )
}

#' Select Optimal Lambda via Goodness-of-Fit Criteria
#'
#' Extracts the tuning parameter \code{lambda} minimizing multiple information criteria from a fitted \code{\link{islasso.path}} object.
#' Supported criteria include AIC, BIC, AICc, eBIC, GCV, and GIC.
#'
#' @param object A fitted model of class \code{"islasso.path"}.
#' @param plot Logical. If \code{TRUE} (default), displays plots for each criterion over the lambda path.
#' @param ... Additional arguments passed to lower-level plotting or diagnostic methods.
#'
#' @details
#' This function identifies the optimal regularization parameter \code{lambda} by minimizing various information-based selection criteria.
#' Degrees of freedom are computed as the trace of the hat matrix, which may be fractional under induced smoothing.
#' This provides a robust alternative to cross-validation, especially in high-dimensional settings.
#'
#' @return A list with components:
#' \item{gof}{Matrix of goodness-of-fit values across lambda values.}
#' \item{minimum}{Index positions of the minimum for each criterion.}
#' \item{lambda.min}{Optimal lambda values that minimize each criterion.}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso.path}}, \code{\link{summary.islasso.path}}, \code{\link{predict.islasso.path}},
#'          \code{\link{coef.islasso.path}}, \code{\link{deviance.islasso.path}}, \code{\link{logLik.islasso.path}},
#'          \code{\link{residuals.islasso.path}}, \code{\link{fitted.islasso.path}}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 30
#' beta <- c(runif(10, -2, 2), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' fit <- islasso.path(y ~ ., data = sim$data, family = gaussian())
#' GoF.islasso.path(fit)
#'
#' @export
GoF.islasso.path <- function(object, plot = TRUE, ...) {
  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)
  gof.name <- c("AIC", "BIC", "AICc", "eBIC", "GCV", "GIC")
  gof <- object$GoF

  # More efficient way to find minimums
  id.min <- apply(gof, 2, which.min)
  lambda.min <- lambda.seq[id.min]
  names(lambda.min) <- gof.name

  out <- list(gof = gof, minimum = id.min, lambda.min = lambda.min)

  if (plot) {
    # Use tidyverse piping more consistently
    p <- data.frame(
      lambda = log(lambda.seq),
      value = c(gof),
      measure = factor(rep(gof.name, each = nlambda), levels = gof.name)
    ) |>
      ggplot(aes(x = lambda, y = value)) +
      geom_line() +
      geom_point(size = 1) +
      facet_wrap(~ measure, scales = "free_y", nrow = 2, ncol = 3) +
      theme_bw() +
      xlab(expression(log(lambda))) +
      ylab("") +
      geom_vline(
        data = data.frame(
          measure = factor(gof.name, levels = gof.name),
          opt = log(lambda.min)
        ),
        aes(xintercept = opt),
        linetype = "dashed",
        colour = "red"
      ) +
      theme(strip.background = element_rect(fill = "white"))

    out$plot <- p

    # Return the plot if desired
    if (inherits(plot, "logical") && plot) {
      print(p)
    }
  }

  return(out)
}
