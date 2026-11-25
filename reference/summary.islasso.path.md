# Summarize islasso.path Model at Specific Lambda

Extracts coefficient estimates, standard errors and p-values from an
[`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
fit at a given regularization level `lambda`.

## Usage

``` r
# S3 method for class 'islasso.path'
summary(object, pval = 1, use.t = FALSE, lambda, ...)
```

## Arguments

- object:

  A fitted object of class `"islasso.path"`.

- pval:

  Numeric threshold for displaying coefficients. Only variables with
  `p-value <= pval` are printed. Unpenalized coefficients (like the
  intercept) are always shown.

- use.t:

  Logical. If `TRUE`, p-values are computed using a t-distribution with
  residual degrees of freedom.

- lambda:

  Numeric. Value of the regularization parameter at which the summary
  should be extracted.

- ...:

  Currently unused.

## Value

An object of class `"summary.islasso.path"` containing filtered
estimates and significance metrics.

## See also

[`islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`GoF.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/GoF.islasso.path.md),
[`coef.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`fitted.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`predict.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.path.md),
[`residuals.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`logLik.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`deviance.islasso.path`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming object `o` is from islasso.path
summary(o, pval = 0.1, lambda = 5)
} # }
```
