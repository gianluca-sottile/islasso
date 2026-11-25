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

  Additional arguments; `digits` is used by the print method.

## Value

An object of class `"summary.islasso"` containing:

- coefficients:

  matrix with estimates, standard errors, pseudo hat values, Wald
  statistics and p-values;

- dispersion, df, res, aic, lambda, nulldev, dev, df.null, df.res,
  family, iter:

  model summaries;

- pval:

  threshold used to filter the printed coefficients;

- unpenalized:

  logical vector of the unpenalized coefficients (always printed);

- call:

  the original model call.

## Details

The column `Df` reports the pseudo hat values \\h_j\\, whose sum is the
effective number of parameters (see
[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)).
The Wald statistics are the ratios of the estimates to their standard
errors; with `use.t = TRUE` (Gaussian family only) the p-values use a t
distribution with the residual degrees of freedom.

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
sim <- simulXy(n = 100, p = 20, beta = c(1, -1, 0.5, rep(0, 17)), seed = 1)
o <- islasso(y ~ ., data = sim$data, lambda = 8)
summary(o, pval = 0.1)            # coefficients with p-value <= 0.1
#> 
#> Call:
#> islasso(formula = y ~ ., lambda = 8, data = sim$data)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.5941 -0.5941 -0.1104  0.5965  2.4541 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)  -0.1418     0.1057  1.000  -1.342     0.18    
#> X1            1.0200     0.1079  1.000   9.457  < 2e-16 ***
#> X2           -0.7371     0.1085  1.000  -6.794 1.09e-11 ***
#> X3            0.5672     0.1095  1.000   5.182 2.20e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 1.116669)
#> 
#>     Null deviance: 344.198  on 99.00  degrees of freedom
#> Residual deviance:  98.881  on 88.55  degrees of freedom
#> AIC: 307.56
#> Lambda: 8
#> 
#> Number of IS iterations: 77
#> 
summary(o, which = c("X1", "X4"))
#> 
#> Call:
#> islasso(formula = y ~ ., lambda = 8, data = sim$data)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.5941 -0.5941 -0.1104  0.5965  2.4541 
#> 
#>    Estimate Std. Error    Df z value Pr(>|z|)    
#> X1  1.02003    0.10786 1.000   9.457   <2e-16 ***
#> X4  0.01381    0.07504 0.679   0.184    0.854    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 1.116669)
#> 
#>     Null deviance: 344.198  on 99.00  degrees of freedom
#> Residual deviance:  98.881  on 88.55  degrees of freedom
#> AIC: 307.56
#> Lambda: 8
#> 
#> Number of IS iterations: 77
#> 
```
