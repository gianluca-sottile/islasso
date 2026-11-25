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

  Graphical parameters: `cex` (point size), `pch` (point shape), `col`
  (colour of the reference lines) and `lwd` (their width).

## Value

Invisibly, a list with the four `ggplot` objects.

## Details

Draws a 2x2 grid of ggplot2 diagnostic plots:

- Top-left: standardized deviance residuals (adjusted for leverage) vs
  linear predictor.

- Top-right: normal Q-Q plot of the standardized deviance residuals (red
  line = reference).

- Bottom-left: squared standardized Pearson residuals vs fitted values
  (variance function).

- Bottom-right: working response vs linear predictor (link function).

The leverages are the diagonal elements of \\W^{1/2} X \tilde U'^{-1}
X^\top W^{1/2}\\, computed from the derivative of the smoothed
estimating equation; their sum is the effective number of parameters.

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
# \donttest{
sim <- simulXy(n = 100, p = 20, beta = c(1, -1, 0.5, rep(0, 17)), seed = 1)
fit <- islasso(y ~ ., data = sim$data, family = gaussian(), lambda = 8)
plot(fit)

# }
```
