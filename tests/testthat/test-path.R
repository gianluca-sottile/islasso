test_that("islasso.path: structure, ordering and agreement with islasso()", {
  d <- sim_gaus()
  o <- islasso.path(y ~ ., data = d$data, nlambda = 15)
  expect_s3_class(o, "islasso.path")
  lam <- o$Info[, "lambda"]
  expect_true(all(diff(lam) > 0))
  expect_equal(nrow(o$Coef), length(lam))
  expect_equal(nrow(o$Fitted.values), length(lam))
  expect_true(all(o$Info[, "conv"] == 0))
  # a single model of the path vs islasso() at the same lambda (warm start differs)
  k <- 8
  f <- islasso(y ~ ., data = d$data, lambda = lam[k])
  expect_equal(unname(o$Coef[k, ]), unname(coef(f)), tolerance = 2e-2)
  expect_equal(unname(o$Info[k, "df"]), f$rank, tolerance = 2e-2)
  # increasing lambda shrinks the effective number of parameters
  expect_true(o$Info[1, "df"] > o$Info[length(lam), "df"])
})

test_that("path methods, interpolation and boundaries", {
  d <- sim_gaus()
  o <- islasso.path(y ~ ., data = d$data, nlambda = 10)
  lam <- o$Info[, "lambda"]
  expect_equal(coef(o, lambda = lam[3]), o$Coef[3, ])
  expect_false(anyNA(coef(o, lambda = max(lam))))
  expect_false(anyNA(fitted(o, lambda = max(lam))))
  expect_error(coef(o, lambda = 2 * max(lam)), "out of bound")
  p <- predict(o, newdata = d$data[1:4, ], type = "response", lambda = c(lam[2], lam[5]))
  expect_equal(dim(p), c(4L, 2L))
  s <- summary(o, lambda = lam[4], pval = 0.1)
  expect_output(print(s), "Residual deviance")
  expect_equal(unname(s$dev), unname(o$Info[4, "deviance"]), tolerance = 1e-8)
  g <- GoF.islasso.path(o, plot = FALSE)
  expect_named(g$lambda.min, c("AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"))
  expect_length(deviance(o), length(lam))
  expect_output(print(o), "lambda")
})

test_that("information criteria are correctly defined", {
  v <- islasso:::.is_gof(ll = 100, dev = 50, df = 5, n = 60, p = 20, g = 0.5)
  expect_equal(unname(v["AICc"]), 100 + 10 + 2 * 5 * 6 / (60 - 5 - 1))
  expect_equal(unname(v["GCV"]), 60 * 50 / 55^2)
  expect_equal(unname(v["eBIC"]), 100 + (log(60) + log(20)) * 5)
})

test_that("GLM path and p > n path", {
  sim <- simulXy(n = 60, p = 10, beta = c(1, -1, rep(0, 8)), interc = 0.5, seed = 3,
                 family = poisson())
  o <- islasso.path(y ~ ., data = sim$data, family = poisson(), nlambda = 8)
  expect_true(all(is.finite(o$Coef)))
  sim2 <- simulXy(n = 25, p = 60, beta = c(2, -2, rep(0, 58)), seed = 4)
  o2 <- islasso.path(y ~ ., data = sim2$data, nlambda = 6)
  expect_true(all(is.finite(o2$SE)))
  expect_equal(ncol(o2$Coef), 61L)
})

test_that("each model of the path coincides with islasso() at the same lambda", {
  d <- sim_gaus()
  o <- islasso.path(y ~ ., data = d$data, nlambda = 8)
  lam <- o$Info[, "lambda"]
  for (k in c(2, 5, 7)) {
    f <- islasso(y ~ ., data = d$data, lambda = lam[k])
    expect_equal(unname(o$Info[k, "df"]), f$rank, tolerance = 5e-3)
    expect_equal(unname(o$Coef[k, ]), unname(coef(f)), tolerance = 5e-3)
  }
})

test_that("parallel path gives the same results", {
  skip_on_cran()
  skip_if_not(.Platform$OS.type == "unix")
  d <- sim_gaus()
  o1 <- islasso.path(y ~ ., data = d$data, nlambda = 6)
  o2 <- islasso.path(y ~ ., data = d$data, nlambda = 6, control = is.control(ncores = 2))
  expect_equal(o1$Coef, o2$Coef)
  expect_equal(o1$Info, o2$Info)
})
