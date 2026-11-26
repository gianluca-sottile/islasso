# Relaxed islasso Fit

Refits an
[`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model leaving unpenalized a subset of covariates, chosen by their
p-values or given explicitly, while the other covariates are
(practically) removed by a large penalty.

## Usage

``` r
relax.islasso(fit, id = NULL, pval = 0.05, ...)
```

## Arguments

- fit:

  a fitted model of class `"islasso"`.

- id:

  optional names or indices (with respect to the covariates, intercept
  excluded) of the covariates to leave unpenalized; by default those
  with p-value not greater than `pval`.

- pval:

  the significance level used when `id` is missing (default 0.05).

- ...:

  further arguments replacing those of the original call.

## Value

An object of class `"islasso"`.

## Details

The refit uses a value of `lambda` equal to 100 times the smallest value
that shrinks all the penalized coefficients to zero, so that the
unselected covariates have (practically) null estimates.

## See also

[`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`summary.islasso()`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md)

## Examples

``` r
# \donttest{
n <- 100; p <- 30
beta <- c(rep(1, 5), rep(0, p - 5))
sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
o <- islasso(y ~ ., data = sim1$data, lambda = 10)
summary(relax.islasso(o), pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., lambda = 2167.66323428158, data = sim1$data, 
#>     unpenalized = c("X1", "X2", "X3", "X4", "X5"))
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.8817 -0.8952  0.0440  0.7494  2.6931 
#> 
#>             Estimate Std. Error    Df z value Pr(>|z|)    
#> (Intercept)   0.1419     0.1161 1.000   1.222    0.222    
#> X1            1.1201     0.1182 1.000   9.474   <2e-16 ***
#> X2            1.1426     0.1182 1.000   9.667   <2e-16 ***
#> X3            1.0686     0.1176 1.000   9.084   <2e-16 ***
#> X4            1.0519     0.1178 1.000   8.926   <2e-16 ***
#> X5            1.1069     0.1193 1.000   9.281   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 1.348562)
#> 
#>     Null deviance: 797.63  on 99  degrees of freedom
#> Residual deviance: 126.76  on 94  degrees of freedom
#> AIC: 321.5
#> Lambda: 2167.7
#> 
#> Number of IS iterations: 79
#> 
summary(relax.islasso(o, id = c("X1", "X2")), pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., lambda = 12113.1018937038, data = sim1$data, 
#>     unpenalized = c("X1", "X2"))
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#> -5.751 -1.524  0.285  1.625  5.078 
#> 
#>             Estimate Std. Error    Df z value Pr(>|z|)    
#> (Intercept)   0.1419     0.2284 1.000   0.621    0.534    
#> X1            1.2495     0.2295 1.000   5.443 5.23e-08 ***
#> X2            1.1779     0.2295 1.000   5.132 2.87e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 5.216484)
#> 
#>     Null deviance: 797.63  on 99  degrees of freedom
#> Residual deviance: 506.00  on 97  degrees of freedom
#> AIC: 453.92
#> Lambda: 12113
#> 
#> Number of IS iterations: 92
#> 
# }
```
