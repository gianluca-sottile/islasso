test_that("tiny lambda reproduces glm() for several families", {
  set.seed(4)
  n <- 150
  X <- matrix(rnorm(n * 3), n); colnames(X) <- paste0("x", 1:3)
  eta <- 0.3 + drop(X %*% c(0.5, -0.4, 0.2))
  dats <- list(
    gaussian = data.frame(y = eta + rnorm(n), X),
    binomial = data.frame(y = rbinom(n, 1, plogis(eta)), X),
    poisson = data.frame(y = rpois(n, exp(eta)), X),
    Gamma = data.frame(y = rgamma(n, shape = 5, scale = exp(eta) / 5), X))
  fams <- list(gaussian = gaussian(), binomial = binomial(), poisson = poisson(),
               Gamma = Gamma(link = "log"))
  for (f in names(fams)) {
    o <- islasso(y ~ ., data = dats[[f]], family = fams[[f]], lambda = 1e-8)
    g <- glm(y ~ ., data = dats[[f]], family = fams[[f]])
    expect_equal(unname(coef(o)), unname(coef(g)), tolerance = 1e-5)
    # the damped fixed-point iterations stop when the change of the standard errors is
    # below 10 * tol: their accuracy is of a few percent at the default tolerance
    expect_equal(unname(o$se), unname(sqrt(diag(vcov(g)))), tolerance = 3e-2)
    expect_equal(o$deviance, g$deviance, tolerance = 1e-6)
  }
})

test_that("other links and families work through the family object", {
  set.seed(5)
  n <- 200
  X <- matrix(rnorm(n * 4), n); colnames(X) <- paste0("x", 1:4)
  eta <- -0.5 + drop(X %*% c(0.6, -0.5, 0, 0))
  d <- data.frame(y = rbinom(n, 1, pnorm(eta)), X)
  for (lk in c("probit", "cloglog", "cauchit")) {
    o <- islasso(y ~ ., data = d, family = binomial(link = lk), lambda = 1e-8)
    g <- glm(y ~ ., data = d, family = binomial(link = lk))
    expect_equal(unname(coef(o)), unname(coef(g)), tolerance = 1e-4)
  }
  d2 <- data.frame(y = rgamma(n, 3, scale = exp(1 + 0.3 * X[, 1]) / 3), X)
  o <- islasso(y ~ ., data = d2, family = inverse.gaussian(link = "log"), lambda = 1e-8)
  g <- glm(y ~ ., data = d2, family = inverse.gaussian(link = "log"))
  expect_equal(unname(coef(o)), unname(coef(g)), tolerance = 1e-4)
  o <- islasso(y ~ ., data = d2, family = quasipoisson(), lambda = 1e-8)
  g <- glm(y ~ ., data = d2, family = quasipoisson())
  expect_equal(unname(coef(o)), unname(coef(g)), tolerance = 1e-4)
  expect_equal(o$dispersion, summary(g)$dispersion, tolerance = 1e-3)
})

test_that("standardisation: penalty on standardised coefficients (glmnet semantics)", {
  d <- sim_gaus()
  o1 <- islasso(y ~ ., data = d$data, lambda = 8)
  d2 <- d$data; d2$X1 <- 10 * d2$X1
  o2 <- islasso(y ~ ., data = d2, lambda = 8)
  expect_equal(coef(o1)[["X1"]], 10 * coef(o2)[["X1"]], tolerance = 1e-5)
  expect_equal(unname(coef(o1)[-2]), unname(coef(o2)[-2]), tolerance = 1e-5)
  # without standardisation the estimates change
  o3 <- islasso(y ~ ., data = d2, lambda = 8, control = is.control(stand = FALSE))
  expect_false(isTRUE(all.equal(coef(o3)[["X1"]], coef(o2)[["X1"]], tolerance = 1e-3)))
})

test_that("models without intercept are not centred", {
  d <- sim_gaus()
  d$data$X1 <- d$data$X1 + 5
  o <- islasso(y ~ . - 1, data = d$data, lambda = 1e-8)
  g <- lm(y ~ . - 1, data = d$data)
  expect_equal(unname(coef(o)), unname(coef(g)), tolerance = 1e-5)
  expect_equal(unname(fitted(o)), unname(fitted(g)), tolerance = 1e-5)
})

test_that("weights, offset, unpenalized and subset", {
  d <- sim_gaus()$data
  w <- rep(1:2, length.out = nrow(d))
  o1 <- islasso(y ~ ., data = d, weights = w, lambda = 1e-8)
  g <- lm(y ~ ., data = d, weights = w)
  expect_equal(unname(coef(o1)), unname(coef(g)), tolerance = 1e-6)
  off <- rep(0.5, nrow(d))
  o3 <- islasso(y ~ ., data = d, offset = off, lambda = 6)
  o4 <- islasso(I(y - 0.5) ~ ., data = d, lambda = 6)
  expect_equal(unname(coef(o3)), unname(coef(o4)), tolerance = 1e-6)
  o5 <- islasso(y ~ ., data = d, lambda = 50, unpenalized = c("X4", "X5"))
  expect_true(all(o5$internal$unpenalized[c(1, 5, 6)]))
  expect_equal(sum(o5$internal$unpenalized), 3)
  expect_gt(abs(coef(o5)[["X4"]]) + abs(coef(o5)[["X5"]]), 0.01)
  o6 <- islasso(y ~ ., data = d, lambda = 6, subset = 1:60)
  expect_equal(nobs(o6), 60)
})

test_that("dense and Woodbury solvers give the same fit", {
  d <- simulXy(n = 30, p = 70, beta = c(2, -1, rep(0, 68)), seed = 2)
  o1 <- islasso(y ~ ., data = d$data, lambda = 4, control = is.control(solver = "dense"))
  o2 <- islasso(y ~ ., data = d$data, lambda = 4, control = is.control(solver = "woodbury"))
  expect_equal(coef(o1), coef(o2), tolerance = 1e-6)
  expect_equal(o1$se, o2$se, tolerance = 1e-5)
  expect_equal(vcov(o1), vcov(o2), tolerance = 1e-5)
})

test_that("MM and Newton inner solvers agree", {
  d <- sim_gaus(seed = 3)
  o1 <- islasso(y ~ ., data = d$data, lambda = 8, control = is.control(inner = "mm"))
  o2 <- islasso(y ~ ., data = d$data, lambda = 8, control = is.control(inner = "newton"))
  expect_equal(coef(o1), coef(o2), tolerance = 1e-3)
  z1 <- coef(o1) / o1$se; z2 <- coef(o2) / o2$se
  expect_lt(max(abs(z1 - z2)), 0.05)
})

test_that("methods work and are consistent", {
  d <- sim_gaus()
  o <- islasso(y ~ ., data = d$data, lambda = 8)
  s <- summary(o, pval = 0.05)
  expect_s3_class(s, "summary.islasso")
  expect_output(print(s), "Estimate")
  expect_output(print(o), "Coefficients")
  expect_equal(sqrt(diag(vcov(o))), o$se, tolerance = 1e-10, ignore_attr = TRUE)
  ci <- confint(o)
  expect_equal(unname(ci[, 2] - ci[, 1]), unname(2 * qnorm(0.975) * o$se), tolerance = 1e-8)
  expect_equal(sum(influence(o)$hat), o$rank, tolerance = 1e-3)
  expect_equal(unname(predict(o, type = "link")), unname(o$linear.predictors))
  nd <- d$data[1:5, ]
  expect_equal(unname(predict(o, newdata = nd)), unname(o$linear.predictors[1:5]), tolerance = 1e-10)
  pc <- predict(o, type = "coefficients", se.fit = TRUE)
  expect_true(is.matrix(pc))
  a <- anova(o, A = "X1 + X2 = 0")
  expect_s3_class(a, "anova.islasso")
  expect_equal(nobs(o), 80)
  expect_s3_class(logLik(o), "logLik")
  expect_length(residuals(o, type = "pearson"), 80)
  expect_length(cooks.distance(o), 80)
})

test_that("relax.islasso and aic.islasso", {
  d <- sim_gaus()
  o <- islasso(y ~ ., data = d$data, lambda = 8)
  r <- relax.islasso(o)
  expect_s3_class(r, "islasso")
  sel <- c("X1", "X2", "X3")
  expect_true(all(r$internal$unpenalized[sel]))
  g <- lm(y ~ X1 + X2 + X3, data = d$data)
  expect_equal(unname(coef(r)[c("(Intercept)", sel)]), unname(coef(g)), tolerance = 1e-3)
  o0 <- islasso(y ~ . - 1, data = d$data, lambda = 8)
  expect_s3_class(relax.islasso(o0, id = 1:2), "islasso")
  lb <- aic.islasso(o, method = "BIC", interval = c(1, 40), trace = FALSE)
  expect_true(lb >= 1 && lb <= 40)
})

test_that("errors are informative", {
  d <- sim_gaus()$data
  expect_error(islasso(y ~ ., data = d, lambda = -1), "non-negative")
  expect_error(islasso(y ~ ., data = d, lambda = 1, alpha = 2), "alpha")
  expect_error(islasso(y ~ ., data = d, lambda = 1, unpenalized = "nope"), "not columns")
  expect_error(islasso(y ~ ., data = d, family = poisson(), lambda = 1), "negative")
})
