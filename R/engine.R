# Numerical engine of the induced smoothed lasso.
#
# Notation (standardised scale):
#   Xt = sqrt(w) * Xs (n x p), zt = sqrt(w) * z, G = Xt'Xt, bz = Xt'zt,
#   lam = lambda * penalty factor (0 for unpenalised columns),
#   se = smoothing scales (standard errors), cm = mixture weights.
# At fixed (se, cm) the IS estimating equation is the gradient of a smooth convex
# objective, minimised by the compiled kernel C_is_inner (or by .is_inner_R).

.is_eps <- 1e-5
.is_delta <- 1e-6

.is_P <- function(b, s, cm) cm * (2 * stats::pnorm(b / s) - 1) + (1 - cm) * (2 * stats::pnorm(b / (.is_eps * s)) - 1)

# Smoothed sign P, penalty value, gradient and curvatures (see src/islasso_core.c):
# for |theta| < .is_delta the smoothed sign is linearised (as in islasso <= 1.6.2);
# 'curvature' is the exact derivative used in the sandwich formula, 'curvature.fit' the
# curvature of the (linearised) objective.
.is_penalty_R <- function(theta, se, cm, lam, alpha) {
  lin <- abs(theta) < .is_delta
  Pd <- .is_P(.is_delta, se, cm) / .is_delta
  P <- ifelse(lin, Pd * theta, .is_P(theta, se, cm))
  t <- theta / se
  dP <- 2 * (cm * stats::dnorm(t) + (1 - cm) * stats::dnorm(t / .is_eps) / .is_eps) / se
  g <- function(b, s) b * (2 * stats::pnorm(b / s) - 1) + 2 * s * stats::dnorm(b / s)
  G <- ifelse(lin, cm * g(.is_delta, se) + (1 - cm) * g(.is_delta, .is_eps * se) +
                0.5 * Pd * (theta^2 - .is_delta^2),
              cm * g(theta, se) + (1 - cm) * g(theta, .is_eps * se))
  val <- lam * (alpha * G + 0.5 * (1 - alpha) * theta^2)
  pen <- lam > 0
  list(value = sum(val[pen]),
       gradient = ifelse(pen, lam * (alpha * P + (1 - alpha) * theta), 0),
       curvature = ifelse(pen, lam * (alpha * dP + (1 - alpha)), 0),
       curvature.fit = ifelse(pen, lam * (alpha * ifelse(lin, Pd, dP) + (1 - alpha)), 0),
       P = P)
}

.is_mm_curv_R <- function(theta, se, cm, lam, alpha) {
  ab <- pmax(abs(theta), .is_delta)
  ifelse(lam > 0, lam * (alpha * .is_P(ab, se, cm) / ab + (1 - alpha)), 0)
}

# Reference implementation of the compiled kernel (dense algebra only).
.is_inner_R <- function(mode, method, G, bz, Xt, zt, theta, se, cm, lam, alpha, tol, maxit) {
  if (mode == 1L) { G <- crossprod(Xt); bz <- drop(crossprod(Xt, zt)) }
  Fq <- function(th, pv) 0.5 * sum(th * (G %*% th)) - sum(bz * th) + pv
  pz <- .is_penalty_R(theta, se, cm, lam, alpha)
  f <- Fq(theta, pz$value)
  info <- 0L
  it <- 0L
  repeat {
    it <- it + 1L
    if (it > maxit) { it <- maxit; info <- 1L; break }
    g <- drop(G %*% theta) - bz + pz$gradient
    Ds <- if (method == 0L) .is_mm_curv_R(theta, se, cm, lam, alpha) else pz$curvature.fit
    H <- G
    diag(H) <- diag(H) + Ds
    R <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(R)) { info <- 2L; break }
    d <- -backsolve(R, forwardsolve(t(R), g))
    gd <- sum(g * d)
    if (!(gd < 0)) {
      if (gd == 0) break
      d <- -g
      gd <- -sum(g^2)
    }
    sig <- if (method == 0L) 0 else 1e-4
    step <- 1
    ok <- FALSE
    for (ls in 1:60) {
      thn <- theta + step * d
      fn <- Fq(thn, .is_penalty_R(thn, se, cm, lam, alpha)$value)
      if (is.finite(fn) && fn <= f + sig * step * gd + 1e-14 * abs(f)) { ok <- TRUE; break }
      step <- step / 2
    }
    if (!ok) { info <- 3L; break }
    theta <- thn
    f <- fn
    pz <- .is_penalty_R(theta, se, cm, lam, alpha)
    if (max(abs(d)) <= tol) break
  }
  D <- .is_penalty_R(theta, se, cm, lam, alpha)$curvature
  H <- G
  diag(H) <- diag(H) + D
  Hi <- tryCatch(chol2inv(chol(H)), error = function(e) NULL)
  if (is.null(Hi)) {
    return(list(theta = theta, v = rep(NA_real_, length(theta)), h = rep(NA_real_, length(theta)),
                D = D, iter = it, info = 2L))
  }
  v <- pmax(diag(Hi) - drop((Hi^2) %*% D), 0)
  h <- 1 - D * diag(Hi)
  list(theta = theta, v = v, h = h, D = D, iter = it, info = info)
}

.is_native <- function() {
  isTRUE(getOption("islasso.engine", "C") != "R") &&
    exists("C_is_inner", envir = asNamespace("islasso"), inherits = FALSE)
}

.is_inner <- function(mode, method, G, bz, Xt, zt, theta, se, cm, lam, alpha, tol, maxit) {
  if (.is_native()) {
    .Call(C_is_inner, as.integer(mode), as.integer(method), G, bz, Xt, zt,
          as.double(theta), as.double(se), as.double(cm), as.double(lam), as.double(alpha),
          as.double(tol), as.integer(maxit))
  } else {
    .is_inner_R(mode, method, G, bz, Xt, zt, theta, se, cm, lam, alpha, tol, maxit)
  }
}

# b' H^-1 G H^-1 b without forming the full matrix
.is_quad <- function(mode, G, Xt, theta, se, cm, lam, alpha, b) {
  if (.is_native()) {
    return(.Call(C_is_quad, as.integer(mode), G, Xt, as.double(theta), as.double(se),
                 as.double(cm), as.double(lam), as.double(alpha), as.double(b)))
  }
  if (is.null(G)) G <- crossprod(Xt)
  D <- .is_penalty_R(theta, se, cm, lam, alpha)$curvature
  H <- G
  diag(H) <- diag(H) + D
  x <- solve(H, b)
  sum(x * (G %*% x))
}

# full unscaled sandwich H^-1 G H^-1 and (optionally) observation leverages
.is_vcov_core <- function(mode, G, Xt, theta, se, cm, lam, alpha, hobs = FALSE) {
  if (.is_native()) {
    out <- .Call(C_is_vcov, as.integer(mode), G, Xt, as.double(theta), as.double(se),
                 as.double(cm), as.double(lam), as.double(alpha), as.logical(hobs))
    if (out$info != 0L) stop("the sandwich matrix could not be computed", call. = FALSE)
    if (!hobs) out$hobs <- NULL
    return(out)
  }
  if (is.null(G)) G <- crossprod(Xt)
  D <- .is_penalty_R(theta, se, cm, lam, alpha)$curvature
  H <- G
  diag(H) <- diag(H) + D
  Hi <- chol2inv(chol(H))
  V <- Hi %*% G %*% Hi
  V <- (V + t(V)) / 2
  list(V = V, hobs = if (hobs) rowSums((Xt %*% Hi) * Xt) else NULL, info = 0L)
}

# ------------------------------------------------------------------------------
# standardisation (glmnet conventions: weighted, divisor sum(w))
# ------------------------------------------------------------------------------

.is_standardize <- function(X, weights, intercept, stand) {
  p <- ncol(X)
  xm <- numeric(p)
  xs <- rep(1, p)
  if (stand) {
    w <- weights / sum(weights)
    idx <- if (intercept) seq_len(p)[-1L] else seq_len(p)
    if (length(idx)) {
      Xi <- X[, idx, drop = FALSE]
      if (intercept) {
        m <- drop(crossprod(w, Xi))
        xm[idx] <- m
        xs[idx] <- sqrt(drop(crossprod(w, sweep(Xi, 2L, m, check.margin = FALSE)^2)))  # two-pass
      } else {
        xs[idx] <- sqrt(drop(crossprod(w, Xi^2)))
      }
      # constant columns (exactly, among the observations with positive weight)
      pos <- weights > 0
      const <- vapply(idx, function(j) { x <- X[pos, j]; all(x == x[1L]) }, logical(1L))
      bad <- const | !(xs[idx] > 1e-10 * pmax(1, abs(xm[idx])))
      if (intercept) xm[idx[bad]] <- 0
      xs[idx[bad]] <- 1
    }
  }
  Xs <- X
  if (stand) {
    Xs <- sweep(sweep(X, 2L, xm, "-", check.margin = FALSE), 2L, xs, "/", check.margin = FALSE)
  }
  list(X = Xs, xm = xm, xs = xs)
}

# Linear map between standardised (theta) and original (beta) coefficients, beta = A theta,
# with A = diag(1/xs) plus, with intercept, the first row (1, -xm/xs). Closed forms, O(p) or O(p^2).
.is_theta2beta <- function(theta, xm, xs, intercept) {
  beta <- theta / xs
  if (intercept && length(theta) > 1L) beta[1L] <- theta[1L] - sum(xm[-1L] * theta[-1L] / xs[-1L])
  beta
}
.is_beta2theta <- function(beta, xm, xs, intercept) {
  theta <- beta * xs
  if (intercept && length(beta) > 1L) theta[1L] <- beta[1L] + sum(xm[-1L] * beta[-1L])
  theta
}
# A V A'
.is_Vback <- function(V, xm, xs, intercept) {
  p <- length(xs)
  if (intercept && p > 1L) {
    a <- c(1, -xm[-1L] / xs[-1L])
    Va <- drop(V %*% a)
    V[1L, ] <- Va
    V[, 1L] <- Va
    V[1L, 1L] <- sum(a * Va)
  }
  V / outer(xs, xs)
}
# gradient with respect to beta: A^{-T} g (theta = A^{-1} beta, A^{-1} = diag(xs) + first row xm)
.is_grad_back <- function(g, xm, xs, intercept) {
  gb <- g * xs
  if (intercept && length(g) > 1L) gb[-1L] <- gb[-1L] + xm[-1L] * g[1L]
  gb[1L] <- if (intercept) g[1L] else gb[1L]
  gb
}
# the row of A for the intercept (for the variance of the intercept)
.is_intercept_row <- function(xm, xs) c(1, -xm[-1L] / xs[-1L])

# ------------------------------------------------------------------------------
# family helpers
# ------------------------------------------------------------------------------

.is_family <- function(family) {
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame(2L))
  if (is.function(family)) family <- family()
  if (!inherits(family, "family") || is.null(family$family))
    stop("'family' not recognized", call. = FALSE)
  for (f in c("linkfun", "linkinv", "mu.eta", "variance", "dev.resids"))
    if (!is.function(family[[f]])) stop(sprintf("the family object has no '%s' function", f), call. = FALSE)
  family
}

# dispersion fixed to one for these families (as in summary.glm)
.is_fixed_dispersion <- function(family) {
  family$family %in% c("binomial", "poisson") || grepl("^Negative Binomial", family$family)
}

.is_gaussian_identity <- function(family) family$family == "gaussian" && family$link == "identity"

# ------------------------------------------------------------------------------
# core fitting loop at a single lambda
# ------------------------------------------------------------------------------

# X: n x p model matrix (with the column of ones if intercept), already standardised.
# Returns the fit on the standardised scale.
.is_fit_single <- function(Xs, y, family, lam, alpha, weights, offset, intercept, control,
                           theta, v0, cm, G = NULL, bz = NULL, mode = NULL) {
  n <- nrow(Xs)
  p <- ncol(Xs)
  if (is.null(mode)) mode <- .is_mode(n, p, control)
  gaus <- .is_gaussian_identity(family)
  tol <- control$tol
  itmax <- control$itmax
  method <- if (identical(control$inner, "newton")) 1L else 0L
  fixphi <- .is_fixed_dispersion(family)
  sigma2 <- if (fixphi) 1 else control$sigma2
  estc <- isTRUE(control$estpai)
  adaptive <- isTRUE(control$adaptive)
  trace <- control$trace
  lam0 <- lam
  nok <- sum(weights > 0)

  se0 <- sqrt(v0)
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  variance <- family$variance
  valideta <- if (is.function(family$valideta)) family$valideta else function(eta) TRUE
  validmu <- if (is.function(family$validmu)) family$validmu else function(mu) TRUE

  irls <- function(eta) {
    mu <- linkinv(eta)
    me <- mu.eta(eta)
    vm <- variance(mu)
    w <- weights * me^2 / vm
    w[!is.finite(w) | weights == 0] <- 0
    z <- (eta - offset) + (y - mu) / me
    z[w == 0] <- 0
    list(mu = mu, me = me, vm = vm, w = w, z = z, sw = sqrt(w))
  }

  if (gaus) {
    sw <- sqrt(weights)
    Xt <- Xs * sw
    zt <- sw * (y - offset)
    if (mode == 0L && is.null(G)) {
      G <- crossprod(Xt)
      bz <- drop(crossprod(Xt, zt))
    }
    wk <- list(w = weights)
  } else {
    # first IRLS weights from 'mustart' (as in glm.fit), unless the starting values are given
    eta <- drop(Xs %*% theta) + offset
    if (is.null(control$b0) || !valideta(eta) || !validmu(linkinv(eta)) || any(!is.finite(eta))) {
      eta <- family$linkfun(.is_mustart(y, weights, family))
    }
    wk <- irls(eta)
    Xt <- Xs * wk$sw
    zt <- wk$sw * wk$z
    if (mode == 0L) { G <- crossprod(Xt); bz <- drop(crossprod(Xt, zt)) }
  }
  colss <- if (mode == 0L) pmax(diag(G), .Machine$double.eps) else pmax(colSums(Xt^2), .Machine$double.eps)

  # Gaussian family: the whole fixed-point iteration runs in compiled code
  if (gaus && trace == 0L && .is_native()) {
    out <- .Call(C_is_gauss, as.integer(mode), as.integer(method), G, bz, sum(zt^2), Xt, zt,
                 as.double(theta), as.double(v0), as.double(cm), as.double(lam), as.double(alpha),
                 as.double(sigma2), as.double(nok), as.double(tol), as.integer(itmax), estc, adaptive)
    if (out$conv == 2L)
      stop("safe exit from the IS-lasso algorithm after a numerical failure; try a different value of 'lambda'",
           call. = FALSE)
    if (out$inner_info %in% c(1L, 3L))
      warning("the inner iterations did not fully converge at the last step; consider increasing 'itmax'",
              call. = FALSE)
    eta <- drop(Xs %*% out$theta) + offset
    return(list(theta = out$theta, se = out$se, se.used = out$se.used, vn = out$vn, v = out$v,
                cm = out$cm, lam = out$lam, h = out$h, edf = sum(out$h), dispersion = out$dispersion,
                conv = out$conv, iter = out$iter, inner_iter = out$inner_iter, tol = out$tol,
                eta = eta, Xt = Xt, zt = zt, G = G, mode = mode, w = weights))
  }

  conv <- 0L
  h <- rep(1, p)
  ind <- ind2 <- NA_real_
  s2 <- if (sigma2 > 0) sigma2 else 1
  it_tot <- 0L
  inner_info <- 0L
  se_used <- se0
  vn <- v0
  i <- 0L
  if (trace == 1L) cat("\nIS-lasso algorithm (lambda = ", format(max(lam), digits = 5), ")\n", sep = "")
  while (i < itmax) {
    i <- i + 1L
    if (adaptive && i > 10L) {
      hh <- pmax(pmin(h, 1), 1e-5)
      lam <- lam0 * (1 / hh^2 - 1)
    }
    if (estc && i > 5L) cm <- stats::plogis(abs(theta / se0))
    theta_old <- theta
    se_used <- se0
    inn <- .is_inner(mode, method, G, bz, Xt, zt, theta, se0, cm, lam, alpha, tol, itmax)
    it_tot <- it_tot + inn$iter
    inner_info <- inn$info
    if (inn$info == 2L || any(!is.finite(inn$theta))) { conv <- 2L; break }
    theta <- inn$theta
    ind2 <- max(abs(theta - theta_old))
    if (gaus) {
      eta <- drop(Xs %*% theta) + offset
      mu <- eta
      res <- y - mu
      pearson <- sum(weights * res^2)
    } else {
      eta <- drop(Xs %*% theta) + offset
      mu <- linkinv(eta)
      nh <- 0L
      while ((!valideta(eta) || !validmu(mu) || any(!is.finite(mu))) && nh < 30L) {
        nh <- nh + 1L
        theta <- (theta + theta_old) / 2
        eta <- drop(Xs %*% theta) + offset
        mu <- linkinv(eta)
      }
      if (nh >= 30L) { conv <- 2L; break }
      wk <- irls(eta)
      pearson <- sum(weights * (y - wk$mu)^2 / wk$vm, na.rm = TRUE)
      Xt <- Xs * wk$sw
      zt <- wk$sw * wk$z
      if (mode == 0L) { G <- crossprod(Xt); bz <- drop(crossprod(Xt, zt)) }
    }
    h <- inn$h
    vn <- inn$v
    v <- v0 + 0.1 * (vn - v0)
    edf <- sum(h)
    redf <- nok - edf
    if (redf <= 1e-12) redf <- 1
    if (sigma2 <= 0) s2 <- pearson / redf
    se <- sqrt(pmax(s2 * v, 0))
    se <- pmax(se, 1e-8 * sqrt(s2 / colss))
    ind <- max(abs(se - se0))
    if (trace == 1L) {
      cat(sprintf("Step = %4d, DEV = %10.4f (%5.2f df), phi = %7.4f, max|dBETA| = %.2e, max|dSE| = %.2e\n",
                  i, pearson, redf, s2, ind2, ind))
    } else if (trace == 2L) {
      cat(sprintf("\n=== IS-lasso step %d ===\n  Pearson statistic = %.6f on %.2f df\n  dispersion = %.5f\n  max|BETAn - BETAo| = %.8f\n  max|SEn - SEo| = %.8f\n",
                  i, pearson, redf, s2, ind2, ind))
    }
    v0 <- v
    se0 <- se
    if (ind <= 10 * tol) break
    if (i >= itmax) { conv <- 1L; break }
  }
  if (conv == 2L) {
    stop("safe exit from the IS-lasso algorithm after a numerical failure; try a different value of 'lambda'",
         call. = FALSE)
  }
  if (trace %in% c(1L, 2L)) cat("Convergence criterion is met.\n\n")
  if (inner_info %in% c(1L, 3L))
    warning("the inner iterations did not fully converge at the last step (",
            if (inner_info == 1L) "iteration limit" else "line search", "); consider increasing 'itmax'",
            call. = FALSE)
  list(theta = theta, se = se0, se.used = se_used, vn = vn, v = v0, cm = cm, lam = lam, h = h, edf = sum(h), dispersion = s2,
       conv = conv, iter = i, inner_iter = it_tot, tol = ind, eta = eta, Xt = Xt, zt = zt, G = G,
       mode = mode, w = if (gaus) weights else wk$w)
}

.is_mode <- function(n, p, control) {
  s <- control$solver
  if (is.null(s) || identical(s, "auto")) return(if (p > 2 * n) 1L else 0L)
  if (identical(s, "woodbury")) 1L else 0L
}

# starting fitted values, as in the 'initialize' expressions of the families in 'stats'
.is_mustart <- function(y, weights, family) {
  f <- family$family
  ybar <- sum(weights * y) / sum(weights)
  mu <- if (f %in% c("binomial", "quasibinomial")) {
    (weights * y + 0.5) / (weights + 1)
  } else if (f %in% c("poisson", "quasipoisson") || grepl("^Negative Binomial", f)) {
    y + 0.1
  } else if (f %in% c("Gamma", "inverse.gaussian")) {
    y
  } else if (f == "gaussian") {
    y
  } else {
    (y + ybar) / 2
  }
  if (!is.function(family$validmu) || isTRUE(family$validmu(mu))) {
    eta <- family$linkfun(mu)
    if (all(is.finite(eta)) && (!is.function(family$valideta) || isTRUE(family$valideta(eta)))) return(mu)
  }
  rep(ybar, length(y))
}
