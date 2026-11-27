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

## Details

Between two values of the path, estimates, standard errors and the other
quantities are linearly interpolated on the \\\log\lambda\\ scale.

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
sim <- simulXy(n = 100, p = 20, beta = c(1, -1, 0.5, rep(0, 17)), seed = 1)
o <- islasso.path(y ~ ., data = sim$data, nlambda = 20)
summary(o, pval = 0.1, lambda = 5)
#> 
#> Call:
#> islasso.path(formula = y ~ ., nlambda = 20, data = sim$data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.50784 -0.54899 -0.07702  0.58924  2.40280 
#> 
#>    Estimate Std. Error z value Pr(>|z|)    
#> X1   1.0502     0.1101   9.538  < 2e-16 ***
#> X2  -0.7539     0.1101  -6.850 7.41e-12 ***
#> X3   0.5986     0.1130   5.298 1.17e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 1.124881)
#> 
#>     Null deviance: 344.198  on 99.00  degrees of freedom
#> Residual deviance:  96.454  on 85.75  degrees of freedom
#> AIC: 310.67
#> Lambda: 5
#> 
#> Number of IS iterations: 59
#> 
```
