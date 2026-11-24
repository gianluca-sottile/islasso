test_that("penalty gradient and curvature match finite differences", {
  set.seed(1)
  p <- 12
  th <- c(rnorm(6), 1e-3 * rnorm(3), 5e-7 * c(1, -1, 0.3))
  se <- runif(p, 0.05, 0.3); cm <- rep(0.5, p); lam <- rep(3, p); lam[2] <- 0
  pz <- islasso:::.is_penalty_R(th, se, cm, lam, 0.7)
  f <- function(x) islasso:::.is_penalty_R(x, se, cm, lam, 0.7)$value
  big <- abs(th) > 1e-4
  expect_equal(pz$gradient[big], num_grad(f, th, 1e-7)[big], tolerance = 1e-5)
  gf <- function(x) islasso:::.is_penalty_R(x, se, cm, lam, 0.7)$gradient
  H <- sapply(seq_len(p), function(j) { e <- replace(numeric(p), j, 1e-7); (gf(th + e) - gf(th - e))[j] / 2e-7 })
  expect_equal(pz$curvature.fit[big], H[big], tolerance = 1e-4)
})

test_that("inner solver: MM and Newton give the same root; gradient is zero", {
  set.seed(2)
  n <- 40; p <- 10
  X <- matrix(rnorm(n * p), n); z <- X[, 1] - X[, 2] + rnorm(n)
  se <- rep(0.15, p); cm <- rep(0.5, p); lam <- rep(5, p); lam[1] <- 0
  G <- crossprod(X); bz <- drop(crossprod(X, z))
  r0 <- islasso:::.is_inner(0L, 0L, G, bz, X, z, rep(1e-3, p), se, cm, lam, 1, 1e-10, 10000)
  r1 <- islasso:::.is_inner(0L, 1L, G, bz, X, z, rep(1e-3, p), se, cm, lam, 1, 1e-10, 10000)
  expect_equal(r0$theta, r1$theta, tolerance = 1e-7)
  g <- drop(G %*% r1$theta) - bz + islasso:::.is_penalty_R(r1$theta, se, cm, lam, 1)$gradient
  expect_lt(max(abs(g)), 1e-6)
  # sandwich diagonal and hat values
  D <- islasso:::.is_penalty_R(r1$theta, se, cm, lam, 1)$curvature
  Hi <- solve(G + diag(D))
  expect_equal(r1$v, diag(Hi %*% G %*% Hi), tolerance = 1e-7)
  expect_equal(r1$h, diag(Hi %*% G), tolerance = 1e-8)
})

test_that("dense and Woodbury back ends agree (p > n)", {
  set.seed(3)
  n <- 20; p <- 50
  X <- matrix(rnorm(n * p), n); z <- X[, 1] * 2 + rnorm(n)
  se <- runif(p, 0.05, 0.2); cm <- rep(0.5, p); lam <- rep(3, p); lam[c(1, 4)] <- 0
  G <- crossprod(X); bz <- drop(crossprod(X, z))
  r0 <- islasso:::.is_inner(0L, 1L, G, bz, X, z, rep(1e-3, p), se, cm, lam, 1, 1e-10, 5000)
  r1 <- islasso:::.is_inner(1L, 1L, NULL, NULL, X, z, rep(1e-3, p), se, cm, lam, 1, 1e-10, 5000)
  expect_equal(r0$theta, r1$theta, tolerance = 1e-8)
  expect_equal(r0$h, r1$h, tolerance = 1e-8)
  expect_equal(r0$v, r1$v, tolerance = 1e-6)
  V0 <- islasso:::.is_vcov_core(0L, G, X, r0$theta, se, cm, lam, 1, hobs = TRUE)
  V1 <- islasso:::.is_vcov_core(1L, NULL, X, r0$theta, se, cm, lam, 1, hobs = TRUE)
  expect_equal(V0$V, V1$V, tolerance = 1e-8)
  expect_equal(V0$hobs, V1$hobs, tolerance = 1e-8)
  expect_equal(diag(V0$V), r0$v, tolerance = 1e-6)
})

test_that("quadratic form of the sandwich and back-transformations", {
  set.seed(6)
  n <- 30; p <- 12
  X <- matrix(rnorm(n * p), n); z <- X[, 1] + rnorm(n)
  se <- runif(p, 0.05, 0.2); cm <- rep(0.5, p); lam <- rep(3, p); lam[1] <- 0
  th <- rnorm(p) * 0.1; b <- rnorm(p)
  V <- islasso:::.is_vcov_core(0L, crossprod(X), X, th, se, cm, lam, 1)$V
  q0 <- islasso:::.is_quad(0L, crossprod(X), X, th, se, cm, lam, 1, b)
  q1 <- islasso:::.is_quad(1L, NULL, X, th, se, cm, lam, 1, b)
  expect_equal(q0, drop(t(b) %*% V %*% b), tolerance = 1e-8)
  expect_equal(q1, q0, tolerance = 1e-8)
  xm <- c(0, rnorm(p - 1)); xs <- c(1, runif(p - 1, 0.5, 2))
  A <- diag(1 / xs); A[1, -1] <- -xm[-1] / xs[-1]
  expect_equal(islasso:::.is_theta2beta(th, xm, xs, TRUE), drop(A %*% th))
  expect_equal(islasso:::.is_beta2theta(drop(A %*% th), xm, xs, TRUE), th)
  expect_equal(islasso:::.is_Vback(V, xm, xs, TRUE), A %*% V %*% t(A))
  expect_equal(islasso:::.is_grad_back(b, xm, xs, TRUE), drop(solve(t(A), b)))
})

test_that("constant covariates are detected by the standardisation", {
  set.seed(1); n <- 50
  X <- cbind(a = rnorm(n), k = 0.1, k2 = 1/3, k3 = 7.3)
  st <- islasso:::.is_standardize(cbind(1, X), rep(c(1, 2), 25), TRUE, TRUE)
  expect_equal(st$xs[3:5], c(1, 1, 1))
  expect_equal(st$xm[3:5], c(0, 0, 0))
})

test_that("compiled Gaussian loop and R reference engine give the same fit", {
  skip_if_not(islasso:::.is_native())
  d <- sim_gaus()
  o1 <- islasso(y ~ ., data = d$data, lambda = 8)
  old <- options(islasso.engine = "R")
  on.exit(options(old))
  o2 <- islasso(y ~ ., data = d$data, lambda = 8)
  expect_equal(coef(o1), coef(o2), tolerance = 1e-8)
  expect_equal(o1$se, o2$se, tolerance = 1e-8)
  expect_equal(o1$iter, o2$iter)
})
