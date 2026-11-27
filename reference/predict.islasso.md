# Prediction Method for islasso Objects

Computes predictions from a fitted
[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model object. Multiple output types supported, including response scale,
linear predictor, and coefficient values.

## Usage

``` r
# S3 method for class 'islasso'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "coefficients", "class", "terms"),
  se.fit = FALSE,
  ci = NULL,
  type.ci = "wald",
  level = 0.95,
  terms = NULL,
  na.action = na.pass,
  ...
)
```

## Arguments

- object:

  A fitted model of class `"islasso"`.

- newdata:

  Optional data frame containing predictors for prediction. If omitted,
  the fitted model matrix is used.

- type:

  Character. Specifies the prediction scale:

  - `"link"` (default): linear predictor scale;

  - `"response"`: original response scale;

  - `"coefficients"`: estimated coefficients;

  - `"class"`: predicted class (only for
    [`binomial()`](https://rdrr.io/r/stats/family.html) family);

  - `"terms"`: contribution of each term to the linear predictor.

- se.fit:

  Logical. If `TRUE`, Wald confidence intervals (at level `level`) are
  returned together with the predictions (on the scale of `type`).

- ci:

  Optional. Precomputed matrix of confidence intervals for the
  coefficients (2 columns), e.g. from
  [`confint.islasso`](https://gianluca-sottile.github.io/islasso/reference/confint.islasso.md).

- type.ci:

  Type of interval. Only `"wald"` is implemented.

- level:

  Confidence level for intervals. Default is `0.95`.

- terms:

  If `type = "terms"`, optionally specify which terms to extract.

- na.action:

  Function to handle missing values in `newdata`. Default: `na.pass`.

- ...:

  Additional arguments passed to downstream methods.

## Value

A numeric vector or, with `se.fit = TRUE` or `type = "terms"`, a matrix.

## See also

[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`summary.islasso`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[`logLik.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`residuals.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`deviance.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
set.seed(1)
n <- 100; p <- 100
beta <- c(runif(20, -3, 3), rep(0, p - 20))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso(y ~ ., data = sim$data, family = gaussian(), lambda = 2)
head(predict(fit, type = "response"))
#>         1         2         3         4         5         6 
#>  5.210184 -1.434233 -1.877905  3.387227 -1.921190  2.984015 
predict(fit, newdata = sim$data[1:3, ], se.fit = TRUE)
#>         Fit     2.5 %    97.5 %
#> 1  5.210184  3.391453 7.0289153
#> 2 -1.434233 -3.307983 0.4395162
#> 3 -1.877905 -4.064405 0.3085957
```
