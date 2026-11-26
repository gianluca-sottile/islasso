# Optimal Tuning Parameter by Information Criteria

Selects the tuning parameter `lambda` of an
[`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model by minimising an information criterion (AIC, BIC/eBIC, AICc, GCV
or GIC; see
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
for the definitions), as an alternative to cross-validation.

## Usage

``` r
aic.islasso(
  object,
  method = c("AIC", "BIC", "AICc", "GCV", "GIC"),
  interval,
  g = 0,
  y,
  X,
  intercept = FALSE,
  family = gaussian(),
  alpha = 1,
  offset,
  weights,
  unpenalized,
  control = is.control(),
  trace = TRUE
)
```

## Arguments

- object:

  a fitted model of class `"islasso"` (optional if `y` and `X` are
  given).

- method:

  the criterion to minimise.

- interval:

  the search interval for `lambda`; by default the range of the glmnet
  path used to obtain the starting values of `object`.

- g:

  numeric in \[0, 1\]: the parameter of the extended BIC (`g = 0` gives
  the BIC).

- y, X:

  response and model matrix (without the column of ones), used if
  `object` is missing.

- intercept, family, alpha, offset, weights, unpenalized, control:

  as in
  [`islasso.fit()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
  used if `object` is missing.

- trace:

  logical: print the progress of the optimisation?

## Value

The value of `lambda` minimising the criterion.

## Details

The criterion is minimised over \\\log\lambda\\ with
[`stats::optimize()`](https://rdrr.io/r/stats/optimize.html); the
degrees of freedom are the effective number of parameters (trace of the
pseudo hat matrix).

## See also

[`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`GoF.islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/GoF.islasso.path.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
# \donttest{
n <- 100; p <- 30
beta <- c(rep(1, 5), rep(0, p - 5))
sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
o <- islasso(y ~ ., data = sim1$data, family = gaussian(), lambda = 5)
lambda.bic <- aic.islasso(o, method = "BIC", interval = c(1, 50), trace = FALSE)
o2 <- update(o, lambda = lambda.bic)
# }
```
