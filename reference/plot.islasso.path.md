# Coefficient Profile and Diagnostic Plots for islasso.path

Generates plots of coefficient profiles, standard errors, gradients,
weights, or goodness-of-fit criteria from a fitted
[`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
model.

## Usage

``` r
# S3 method for class 'islasso.path'
plot(
  x,
  yvar = c("coefficients", "se", "score", "gradient", "weight", "gof"),
  gof = c("none", "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"),
  label = FALSE,
  legend = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class `"islasso.path"`, typically created via
  [`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md).

- yvar:

  Character. Specifies what to display on the y-axis. Choices are:

  - `"coefficients"` - coefficient paths over `log(lambda)`,

  - `"se"` - standard errors over `log(lambda)`,

  - `"score"` - standardized scores \\x_j^\top W(z -
    X\hat\beta)/\lambda_j\\ over `log(lambda)`, the induced-smoothing
    analogue of the lasso optimality conditions (dashed lines at \\\pm
    1\\),

  - `"gradient"` - gradient of the smoothed estimating equation over
    `log(lambda)` (a convergence diagnostic),

  - `"weight"` - mixture weights used in smoothing,

  - `"gof"` - goodness-of-fit values.

- gof:

  Character. Criterion used for highlighting active variables. Choices:
  `"none"`, `"AIC"`, `"BIC"`, `"AICc"`, `"eBIC"`, `"GCV"`, `"GIC"`.

- label:

  Logical. Whether to annotate curves with variable names.

- legend:

  Logical. Whether to display a plot legend.

- ...:

  Additional graphical parameters, e.g. `main`, `xlab`, `ylab`, `xlim`,
  `ylim`, `lty`, `col`, `lwd`, `cex.axis`, `cex.lab`, `cex.main`,
  `gof_lty`, `gof_col`, `gof_lwd`.

## Value

A `ggplot` object (printed when not assigned).

## Details

This function visualizes the behavior of the solution path across a
sequence of lambda values, helping diagnose coefficient shrinkage,
influence of penalty, and variable selection stability. All the profiles
are drawn against \\\log\lambda\\.

The standardized score of a covariate (computed on the standardized
scale on which the penalty acts) equals, at the solution, the smoothed
sign \\P_j\\ (lasso) or \\\alpha P_j + (1-\alpha)\hat\beta_j\\ (elastic
net): it is close to \\\pm 1\\ for the covariates with sizeable effects
and strictly inside \\(-1, 1)\\ for those shrunk towards zero, as the
subgradient conditions of the lasso. The gradient, instead, is the
residual of the smoothed estimating equation: it is zero up to the
convergence tolerance, and its irregular values for large \\\lambda\\
reflect the very large curvature of the smoothed penalty when the
standard errors of null coefficients are tiny, not features of the
model.

## See also

[`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`GoF.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/GoF.islasso.path.md),
[`summary.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.path.md),
[`coef.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`fitted.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`predict.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.path.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
# \donttest{
set.seed(1)
n <- 100; p <- 30
beta <- c(runif(10, -2, 2), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso.path(y ~ ., data = sim$data, family = gaussian(), nlambda = 30)
plot(fit, yvar = "coefficients", gof = "AICc", label = TRUE)

plot(fit, yvar = "se", gof = "AICc")

plot(fit, yvar = "score", gof = "AICc")

plot(fit, yvar = "gradient", gof = "AICc")

plot(fit, yvar = "gof", gof = "AICc")

# }
```
