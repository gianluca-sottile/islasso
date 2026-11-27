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

A vector or a matrix (one column per value of `lambda`) depending on
`type`.

## Details

For values of `lambda` not in the path, the coefficients are linearly
interpolated on the \\\log\lambda\\ scale.

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
set.seed(1)
n <- 100; p <- 30
beta <- c(runif(10, -3, 3), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso.path(y ~ ., data = sim$data, family = gaussian(), nlambda = 30)
optimal <- GoF.islasso.path(fit, plot = FALSE)
pred <- predict(fit, newdata = sim$data[1:5, ], type = "response",
                lambda = optimal$lambda.min[c("AIC", "BIC")])
pred
#>   lambda = 7.646 lambda = 7.646
#> 1       4.614027       4.614027
#> 2      -1.223396      -1.223396
#> 3       7.608375       7.608375
#> 4       1.551165       1.551165
#> 5       3.113812       3.113812
```
