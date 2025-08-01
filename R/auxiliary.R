qqNorm <- function(x,
                   probs = seq(0.005, 0.995, length.out = 200),
                   centre = FALSE, scale = FALSE,
                   leg = TRUE, mean = 0, sd = 1,
                   dF = FALSE, ylab = NULL,
                   color = "black", ...) {
  stopifnot(is.numeric(x))

  if (centre) x <- x - mean(x)
  if (scale)  x <- (x - mean(x)) / sd(x) * sd + mean

  emp_q <- quantile(x, probs, names = FALSE)
  teor_q <- qnorm(probs, mean = mean, sd = sd)

  df <- data.frame(Theoretical = teor_q, Empirical = emp_q)

  if (dF) {
    df_cdf <- data.frame(
      x = sort(x),
      Empirical = seq_along(x) / length(x),
      Theoretical = pnorm(sort(x), mean = mean, sd = sd)
    )

    p <- ggplot(df_cdf, aes(x = x)) +
      geom_step(aes(y = Empirical), color = color, ...) +
      geom_line(aes(y = Theoretical), color = "red", linewidth = 0.9) +
      labs(x = "", y = "Distribution Function") +
      theme_minimal()
  } else {
    p <- ggplot(df, aes(x = Theoretical, y = Empirical)) +
      geom_point(color = color, ...) +
      geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
      labs(x = "Theoretical Quantiles",
           y = ylab %||% "Empirical Quantiles") +
      theme_minimal()
  }

  if (leg && !dF) {
    emp_stats <- sprintf("emp. mean = %.3f\nemp. sd = %.3f", mean(x), sd(x))
    theor_stats <- sprintf("theor. mean = %.3f\ntheor. sd = %.3f", mean, sd)
    p <- p +
      annotate("text", x = min(teor_q), y = max(emp_q),
               label = emp_stats, hjust = 0, vjust = 1, size = 3.2) +
      annotate("text", x = max(teor_q), y = min(emp_q),
               label = theor_stats, hjust = 1, vjust = 0, size = 3.2, color = "red")
  }

  return(p)
}

#' Simulate Model Matrix and Response Vector
#'
#' Generates synthetic covariates and response vector from a specified distribution for simulation studies or method validation.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Total number of covariates in the model matrix.
#' @param interc Numeric. Intercept to include in the linear predictor. Default is \code{0}.
#' @param beta Numeric vector of length \code{p}. Regression coefficients in the linear predictor.
#' @param family Distribution and link function. Allowed: \code{gaussian()}, \code{binomial()}, \code{poisson()} and , \code{Gamma()}. Can be a string, function, or family object.
#' @param prop Numeric in \code{[0,1]}. Used only if \code{beta} is missing; proportion of non-zero coefficients in \code{p}. Default is \code{0.1}.
#' @param lim.b Numeric vector of length 2. Range for coefficients if \code{beta} is missing. Default: \code{c(-3, 3)}.
#' @param sigma Standard deviation of Gaussian response. Default is \code{1}.
#' @param size Integer. Number of trials for binomial response. Default is \code{1}.
#' @param rho Numeric. Correlation coefficient for generating covariates. Used to create AR(1)-type covariance: \code{rho^|i-j|}. Default is \code{0}.
#' @param scale.data Logical. Whether to scale columns of the model matrix. Default is \code{TRUE}.
#' @param seed Optional. Integer seed for reproducibility.
#' @param X Optional. Custom model matrix. If supplied, it overrides the internally generated \code{X}.
#' @param dispersion Dispersion parameter of Gamma response. Default is \code{0.1}.
#'
#' @return A list with components:
#' \item{X}{Model matrix of dimension \code{n x p}}
#' \item{y}{Simulated response vector}
#' \item{beta}{True regression coefficients used}
#' \item{eta}{Linear predictor}
#'
#' @examples
#' n <- 100; p <- 100
#' beta <- c(runif(10, -3, 3), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1234)
#' o <- islasso(y ~ ., data = sim$data, family = gaussian())
#' summary(o, pval = 0.05)
#'
#' @export
simulXy <- function(n, p, interc = 0, beta,
                    family = gaussian(), prop = 0.1,
                    lim.b = c(-3, 3), sigma = 1, size = 1,
                    rho = 0, scale.data = TRUE,
                    seed = NULL, X = NULL, dispersion = 0.1) {

  if (!is.null(seed)) set.seed(seed)

  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()

  fam <- family$family
  if (is.null(fam) || !(fam %in% c("gaussian", "binomial", "poisson", "Gamma"))) {
    stop("Family not recognized. Use 'gaussian', 'binomial', 'poisson', or 'Gamma'.")
  }

  if (missing(beta)) {
    if (prop < 0 || prop > 1) stop("Invalid 'prop'. Must be between 0 and 1.")
    p.true <- trunc(prop * p)
    beta <- c(runif(p.true, lim.b[1], lim.b[2]), rep(0, p - p.true))
  }

  if (is.null(X)) {
    X <- modelX(n, p, rho, scale.data)
    colnames(X) <- paste0("X", seq_len(p))
  }

  eta <- as.vector(X %*% beta + interc)
  mu <- family$linkinv(eta)

  y <- switch(fam,
              gaussian = rnorm(n, mean = mu, sd = sigma),
              binomial = rbinom(n, size = size, prob = mu),
              poisson  = rpois(n, lambda = mu),
              Gamma    = rgamma(n, shape = 1 / dispersion, scale = mu * dispersion))

  if (fam == "binomial" && size > 1) {
    y <- cbind(failure = size - y, success = y)
  }

  out <- list(data = data.frame(y = y, X),
              beta0 = interc, beta = beta)
  return(out)
}


modelX <- function(n, p, rho = 0.5, scale.data = TRUE) {
  # Crea matrice di correlazione AR(1)-like
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))

  # Fattorizzazione di Cholesky
  cSigma <- chol(Sigma)

  # Generazione di dati simulati con correlazione specificata
  X <- matrix(rnorm(n * p), n, p) %*% cSigma

  # Opzionale: scala i dati
  if (scale.data) {
    X <- scale(X)
  }

  return(X)
}


lminfl <- function(mod, tol = 1e-8) {
  # Estrae la decomposizione QR
  Q <- qr.Q(mod$qr)
  n2 <- nrow(Q)
  X <- model.matrix(mod)
  n <- nrow(X)
  k <- ncol(X)
  q <- if (is.matrix(residuals(mod))) ncol(residuals(mod)) else 1
  resid <- residuals(mod)

  # Calcolo dei valori hat: diagonale della matrice degli proiettori
  hat <- rowSums(Q^2)[1:n]
  hat[hat >= (1 - tol)] <- 1.0

  # Calcolo della sigma leave-one-out per ogni osservazione e per ogni colonna (se multivariate)
  denom <- (n2 - k - 1)
  sigma <- matrix(NA, n, q)

  for (j in seq_len(q)) {
    res_j <- if (q == 1) resid else resid[, j]
    rss <- sum(res_j^2)
    for (i in seq_len(n)) {
      if (hat[i] < 1) {
        sigma[i, j] <- sqrt((rss - res_j[i]^2 / (1 - hat[i])) / denom)
      } else {
        sigma[i, j] <- sqrt(rss / denom)
      }
    }
  }

  return(list(hat = hat, sigma = drop(sigma)))
}

is.influence <- function (model, do.coef = TRUE) {
  wt.res <- weighted.residuals(model)
  e <- na.omit(wt.res)
  n <- length(wt.res)
  is.mlm <- is.matrix(e)
  if (model$rank == 0) {
    n <- length(wt.res)
    sigma <- sqrt(deviance(model)/df.residual(model))
    res <- list(hat = rep(0, n), coefficients = matrix(0, n, 0), sigma = rep(sigma, n))
  }
  else {
    e[abs(e) < 100 * .Machine$double.eps * median(abs(e))] <- 0
    mqr <- model$qr
    do.coef <- as.logical(do.coef)
    tol <- 10 * .Machine$double.eps
    # res2 <- .Call(stats:::C_influence, mqr, e, tol)
    # res <- infl(mqr, FALSE, e, tol)
    res <- lminfl(model, tol)
    res$hat <- res$hat[1:n]
    res$sigma <- res$sigma[1:n]
    if (do.coef) {
      ok <- seq_len(mqr$rank)
      Q <- qr.Q(mqr)[, ok, drop = FALSE]
      R <- qr.R(mqr)[ok, ok, drop = FALSE]
      hat <- res$hat
      invRQtt <- t(backsolve(R, t(Q)))
      k <- NCOL(Q)
      q <- NCOL(e)
      if (is.mlm) {
        cf <- array(0, c(n, k, q))
        for (j in seq_len(q)) cf[, , j] <- invRQtt[1:n, , drop = FALSE] * ifelse(hat == 1, 0, e[, j]/(1 - hat))
      }
      else cf <- invRQtt[1:n, , drop = FALSE] * ifelse(hat == 1, 0, e/(1 - hat))
      res$coefficients <- cf
    }
    drop1d <- function(a) {
      d <- dim(a)
      if (length(d) == 3L && d[[3L]] == 1L) dim(a) <- d[-3L]
      a
    }
    if (is.null(model$na.action)) {
      if (!is.mlm) {
        res$sigma <- drop(res$sigma)
        if (do.coef) res$coefficients <- drop1d(res$coefficients)
      }
    }
    else {
      hat <- naresid(model$na.action, res$hat)
      hat[is.na(hat)] <- 0
      res$hat <- hat
      if (do.coef) {
        coefficients <- naresid(model$na.action, res$coefficients)
        coefficients[is.na(coefficients)] <- 0
        res$coefficients <- if (is.mlm)
          coefficients
        else drop1d(coefficients)
      }
      sigma <- naresid(model$na.action, res$sigma)
      sigma[is.na(sigma)] <- sqrt(deviance(model)/df.residual(model))
      res$sigma <- if (is.mlm)
        sigma
      else drop(sigma)
    }
  }
  res$wt.res <- naresid(model$na.action, e)
  res$hat[res$hat > 1 - 10 * .Machine$double.eps] <- 1
  names(res$hat) <- names(res$sigma) <- names(res$wt.res)
  if (do.coef) {
    cf <- coef(model)
    if (is.mlm) {
      dnr <- dimnames(res$wt.res)
      dimnames(res$coefficients) <- list(dnr[[1L]], rownames(cf)[!apply(cf, 1L, anyNA)], dnr[[2L]])
    }
    else dimnames(res$coefficients) <- list(names(res$wt.res), names(cf)[!is.na(cf)])
  }
  res[c("hat", "coefficients", "sigma", "wt.res")]
}

islasso.diag <- function (glmfit) {
  w <- if (is.null(glmfit$prior.weights))
    rep(1, length(glmfit$residuals))
  else glmfit$prior.weights
  sd <- sqrt(glmfit$dispersion)
  dev <- residuals(glmfit, type = "deviance") / sd
  pear <- residuals(glmfit, type = "pearson") / sd
  h <- rep(0, length(w))
  h[w != 0] <- influence(glmfit)$hat
  p <- glmfit$rank
  rp <- pear/sqrt(1 - h)
  rd <- dev/sqrt(1 - h)
  cook <- (h * rp^2)/((1 - h) * p)
  res <- sign(dev) * sqrt(dev^2 + h * rp^2)
  list(res = res, rd = rd, rp = rp, cook = cook, h = h, sd = sd)
}

islasso.diag.plots <- function (glmfit, glmdiag = islasso.diag(glmfit), subset = NULL,
                                iden = FALSE, labels = NULL, ret = FALSE) {
  if (is.null(glmdiag)) glmdiag <- islasso.diag(glmfit)
  if (is.null(subset))
    subset <- seq_along(glmdiag$h)
  else if (is.logical(subset))
    subset <- seq_along(subset)[subset]
  else if (is.numeric(subset) && all(subset < 0))
    subset <- (1L:(length(subset) + length(glmdiag$h)))[subset]
  else if (is.character(subset)) {
    if (is.null(labels))
      labels <- subset
    subset <- seq_along(subset)
  }

  par(mfrow = c(2, 2))
  x1 <- predict(glmfit)
  plot(x1, glmdiag$res, xlab = "Linear predictor", ylab = "Residuals")
  pars <- vector(4L, mode = "list")
  pars[[1L]] <- par("usr")
  y2 <- glmdiag$rd
  x2 <- qnorm(ppoints(length(y2)))[rank(y2)]
  plot(x2, y2, ylab = "Quantiles of standard normal", xlab = "Ordered deviance residuals")
  abline(0, 1, lty = 2)
  pars[[2L]] <- par("usr")
  hh <- glmdiag$h/(1 - glmdiag$h)
  plot(hh, glmdiag$cook, xlab = "h/(1-h)", ylab = "Cook statistic")
  rx <- range(hh)
  ry <- range(glmdiag$cook)
  rank.fit <- glmfit$rank
  nobs <- rank.fit + glmfit$df.residual
  cooky <- 8/(nobs - 2 * rank.fit)
  hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
  if ((cooky >= ry[1L]) && (cooky <= ry[2L]))
    abline(h = cooky, lty = 2)
  if ((hy >= rx[1L]) && (hy <= rx[2L]))
    abline(v = hy, lty = 2)
  pars[[3L]] <- par("usr")
  plot(subset, glmdiag$cook, xlab = "Case", ylab = "Cook statistic")
  if ((cooky >= ry[1L]) && (cooky <= ry[2L]))
    abline(h = cooky, lty = 2)
  xx <- list(x1, x2, hh, subset)
  yy <- list(glmdiag$res, y2, glmdiag$cook, glmdiag$cook)
  pars[[4L]] <- par("usr")
  if (is.null(labels))
    labels <- names(x1)
  while (iden) {
    cat("****************************************************\n")
    cat("Please Input a screen number (1,2,3 or 4)\n")
    cat("0 will terminate the function \n")
    num <- as.numeric(readline())
    if ((length(num) > 0L) && ((num == 1) || (num == 2) ||
                               (num == 3) || (num == 4))) {
      cat(paste("Interactive Identification for screen",
                num, "\n"))
      cat("left button = Identify, center button = Exit\n")
      nm <- num + 1
      par(mfg = c(trunc(nm/2), 1 + nm%%2, 2, 2))
      par(usr = pars[[num]])
      identify(xx[[num]], yy[[num]], labels)
    }
    else iden <- FALSE
  }
  par(mfrow = c(1, 1))
  if (ret)
    glmdiag
  else invisible()
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
  piv <- if (p) object$qr$pivot[p1]
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

