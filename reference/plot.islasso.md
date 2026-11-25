# Diagnostic Plots for islasso Models

Produces standard diagnostic plots for a fitted
[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model to assess residuals, model fit, and variance structure.

## Usage

``` r
# S3 method for class 'islasso'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"islasso"`, typically created via
  [`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md).

- ...:

  Additional graphical parameters passed to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) functions.

## Details

Generates a 2x2 grid of diagnostic plots:

- Top-left: Deviance residuals vs fitted values.

- Top-right: Normal Q-Q plot of standardized deviance residuals (red
  line = reference).

- Bottom-left: Squared standardized Pearson residuals vs fitted values.

- Bottom-right: Working response vector vs linear predictor.

These plots help assess the assumptions of linearity, homoscedasticity,
and residual normality in penalized regression.

## See also

[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`summary.islasso`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[`residuals.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`logLik.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`predict.islasso`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.md),
[`deviance.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
if (FALSE) { # \dontrun{
  set.seed(1)
  n <- 100; p <- 100
  beta <- c(runif(20, -3, 3), rep(0, p - 20))
  sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
  fit <- islasso(y ~ ., data = sim$data, family = gaussian(), lambda = 2)
  plot(fit)
} # }
```
