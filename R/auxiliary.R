#' Simulate a Model Matrix and a Response Vector
#'
#' Generates covariates (multivariate normal with AR(1) correlation) and a response from a
#' Gaussian, binomial, Poisson or Gamma regression model, for simulation studies and examples.
#'
#' @param n number of observations.
#' @param p number of covariates.
#' @param interc the intercept of the linear predictor (default 0).
#' @param beta vector of `p` regression coefficients; if missing, a proportion `prop` of
#'   coefficients is drawn uniformly in `lim.b` and the others are zero.
#' @param family a family object, function or name: `gaussian`, `binomial`, `poisson` or `Gamma`
#'   (with any of their links).
#' @param prop proportion of non-null coefficients when `beta` is missing (default 0.1).
#' @param lim.b range of the non-null coefficients when `beta` is missing.
#' @param sigma standard deviation of the Gaussian response.
#' @param size number of trials of the binomial response; with `size > 1` the response is a
#'   two-column matrix of failures and successes.
#' @param rho correlation parameter: \eqn{cor(X_j, X_k) = \rho^{|j-k|}}.
#' @param scale.data logical: should the columns of the model matrix be standardised?
#' @param seed optional seed; the state of the random number generator of the session is
#'   restored on exit.
#' @param X optional model matrix replacing the generated one.
#' @param dispersion dispersion parameter of the Gamma response (default 0.1).
#'
#' @return A list with components
#' \item{data}{a data frame with the response `y` (columns `y.failure`, `y.success` when
#'   `size > 1`) and the covariates `X1`, ..., `Xp`;}
#' \item{beta0}{the intercept;}
#' \item{beta}{the regression coefficients.}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@@unipa.it}
#' @seealso [islasso()], [islasso.path()]
#' @examples
#' n <- 100; p <- 100
#' beta <- c(runif(10, -3, 3), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1234)
#' str(sim$data[, 1:5])
#' @export
simulXy <- function(n, p, interc = 0, beta, family = gaussian(), prop = 0.1, lim.b = c(-3, 3),
                    sigma = 1, size = 1, rho = 0, scale.data = TRUE, seed = NULL, X = NULL,
                    dispersion = 0.1) {
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  fam <- family$family
  if (is.null(fam) || !(fam %in% c("gaussian", "binomial", "poisson", "Gamma")))
    stop("'family' must be one of gaussian, binomial, poisson or Gamma", call. = FALSE)
  .is_with_seed(seed, {
    if (missing(beta)) {
      if (prop < 0 || prop > 1) stop("'prop' must be in [0, 1]", call. = FALSE)
      p.true <- trunc(prop * p)
      beta <- c(stats::runif(p.true, lim.b[1L], lim.b[2L]), rep(0, p - p.true))
    }
    if (length(beta) != p) stop("'beta' must have length 'p'", call. = FALSE)
    if (is.null(X)) {
      X <- modelX(n, p, rho, scale.data)
    } else {
      X <- as.matrix(X)
      if (nrow(X) != n || ncol(X) != p) stop("'X' must be a n x p matrix", call. = FALSE)
    }
    colnames(X) <- paste0("X", seq_len(p))
    eta <- as.vector(X %*% beta + interc)
    mu <- family$linkinv(eta)
    y <- switch(fam,
                gaussian = stats::rnorm(n, mean = mu, sd = sigma),
                binomial = stats::rbinom(n, size = size, prob = mu),
                poisson = stats::rpois(n, lambda = mu),
                Gamma = stats::rgamma(n, shape = 1 / dispersion, scale = mu * dispersion))
    if (fam == "binomial" && size > 1) y <- cbind(failure = size - y, success = y)
    list(data = data.frame(y = y, X), beta0 = interc, beta = beta)
  })
}

# AR(1) covariates: X_j = rho X_{j-1} + sqrt(1 - rho^2) Z_j, i.e. Z %*% chol(Sigma) with
# Sigma_jk = rho^|j-k|, computed in O(np) time and memory
modelX <- function(n, p, rho = 0.5, scale.data = TRUE) {
  X <- matrix(stats::rnorm(n * p), n, p)
  if (rho != 0 && p > 1L) {
    s <- sqrt(1 - rho^2)
    for (j in 2:p) X[, j] <- rho * X[, j - 1L] + s * X[, j]
  }
  if (scale.data) X <- scale(X)
  X
}


predislasso <- function(object, newdata, type = c("response", "terms"),
                        terms = NULL, na.action = na.pass, ...){
  type <- match.arg(type)

  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset")))
      for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }

  n <- object$internal$n
  p <- object$internal$p
  p1 <- seq_len(p)
  piv <- p1
  beta <- object$coefficients
  predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
  if (!is.null(offset)) predictor <- predictor + offset

  if (type == "terms") {
    if (!mmDone) {
      mm <- model.matrix(object)
      mmDone <- TRUE
    }
    aa <- attr(mm, "assign")
    ll <- attr(tt, "term.labels")
    hasintercept <- attr(tt, "intercept") > 0L
    if (hasintercept)
      ll <- c("(Intercept)", ll)
    aaa <- factor(aa, labels = ll)
    asgn <- split(order(aa), aaa)
    if (hasintercept) {
      asgn$"(Intercept)" <- NULL
      avx <- colMeans(mm)
      termsconst <- sum(avx[piv] * beta[piv])
    }
    nterms <- length(asgn)
    if (nterms > 0) {
      predictor <- matrix(ncol = nterms, nrow = NROW(X))
      dimnames(predictor) <- list(rownames(X), names(asgn))
      if (hasintercept)
        X <- sweep(X, 2L, avx, check.margin = FALSE)
      unpiv <- rep.int(0L, NCOL(X))
      unpiv[piv] <- p1
      for (i in seq.int(1L, nterms, length.out = nterms)) {
        iipiv <- asgn[[i]]
        ii <- unpiv[iipiv]
        iipiv[ii == 0L] <- 0L
        predictor[, i] <- if (any(iipiv > 0L))
          X[, iipiv, drop = FALSE] %*% beta[iipiv]
        else 0
      }
      if (!is.null(terms)) predictor <- predictor[, terms, drop = FALSE]
    }
    else predictor <- ip <- matrix(0, n, 0L)
    attr(predictor, "constant") <- if (hasintercept)
      termsconst
    else 0
  }
  if (missing(newdata) && !is.null(na.act <- object$na.action)) predictor <- napredict(na.act, predictor)
  attr(predictor, "X") <- X
  return(predictor)
}

'%||%' <- function (L, R) if (is.null(L)) R else L


# observation leverages x_i' H^-1 x_i (weighted, standardised scale) from the sandwich Hessian
.is_leverage <- function(object) {
  it <- object$internal
  X <- model.matrix(object)
  Xs <- sweep(sweep(X, 2L, it$xm, "-", check.margin = FALSE), 2L, it$xs, "/", check.margin = FALSE)
  w <- object$weights
  Xt <- Xs * sqrt(w)
  mode <- if (is.null(it$mode)) 0L else it$mode
  G <- if (mode == 0L) crossprod(Xt) else NULL
  h <- .is_vcov_core(mode, G, Xt, it$theta, it$se.std, it$c, it$lam, object$alpha, hobs = TRUE)$hobs
  h <- pmin(pmax(h, 0), 1)
  names(h) <- rownames(X)
  h
}

# arrange a list of ggplot objects in a grid (no extra dependency)
.is_arrange <- function(plots, ncol = 2L) {
  nr <- ceiling(length(plots) / ncol)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, ncol)))
  for (k in seq_along(plots)) {
    print(plots[[k]], vp = grid::viewport(layout.pos.row = (k - 1L) %/% ncol + 1L,
                                          layout.pos.col = (k - 1L) %% ncol + 1L))
  }
  grid::popViewport()
  invisible(NULL)
}
