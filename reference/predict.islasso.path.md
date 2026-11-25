# Prediction Method for islasso.path Objects

Generates predictions from a fitted
[`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
model at one or more lambda values. Supports various output types
including linear predictors, response scale, class labels, and
coefficients.

## Usage

``` r
# S3 method for class 'islasso.path'
predict(
  object,
  newdata,
  type = c("link", "response", "coefficients", "class"),
  lambda,
  ...
)
```

## Arguments

- object:

  A fitted model object of class `"islasso.path"`.

- newdata:

  Optional data frame containing covariates for prediction. If omitted,
  returns fitted values from the original model.

- type:

  Character. Type of prediction:

  - `"link"` (default) - linear predictor scale,

  - `"response"` - original response scale,

  - `"coefficients"` - estimated coefficients,

  - `"class"` - predicted class labels (only for binomial models).

- lambda:

  Numeric value(s). Specific lambda value(s) at which predictions are
  required. If missing, predictions are computed for the full lambda
  sequence.

- ...:

  Additional arguments passed to lower-level methods.

## Value

A vector, matrix, or list depending on the `type` requested.

## See also

[`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`summary.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.path.md),
[`coef.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`GoF.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/GoF.islasso.path.md),
[`fitted.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`logLik.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`residuals.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`deviance.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
if (FALSE) { # \dontrun{
  set.seed(1)
  n <- 100; p <- 30
  beta <- c(runif(10, -3, 3), rep(0, p - 10))
  sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
  fit <- islasso.path(y ~ ., data = sim$data, family = gaussian())
  optimal <- GoF.islasso.path(fit)
  pred <- predict(fit, type = "response", lambda = optimal$lambda.min)
} # }
```
