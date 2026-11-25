# Summarize islasso Fitted Model

Provides a concise summary of a fitted
[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model, including p-values and optional filtering.

## Usage

``` r
# S3 method for class 'islasso'
summary(object, pval = 1, which, use.t = FALSE, type.pval = "wald", ...)
```

## Arguments

- object:

  A fitted model of class `"islasso"`.

- pval:

  Numeric threshold for displaying coefficients. Only those with \\p
  \le\\ `pval` are printed. Unpenalized coefficients (like intercepts)
  are always shown.

- which:

  Optional. Specifies a subset of coefficients to test. If missing, all
  parameters are evaluated.

- use.t:

  Logical. If `TRUE`, p-values are computed using the t-distribution and
  residual degrees of freedom.

- type.pval:

  Character. Type of p-value approximation. Only `"wald"` (default) is
  implemented.

- ...:

  Additional arguments (not currently used).

## Value

An object of class `"summary.islasso"` containing:

- coefficients:

  Coefficient estimates and related statistics

- pval:

  Threshold used to filter coefficients

- call:

  Original model call

## See also

[`islasso.fit`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`residuals.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`logLik.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`predict.islasso`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.md),
[`deviance.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming object `o` from an islasso fit
summary(o, pval = 0.1)  # Show coefficients with p <= 0.1
} # }
```
