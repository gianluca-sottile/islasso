# confint method for `islasso` objects

Computes confidence intervals for `islasso` objects using a Wald-type
approach.

## Usage

``` r
# S3 method for class 'islasso'
confint(object, parm, level = 0.95, type.ci = "wald", trace = TRUE, ...)
```

## Arguments

- object:

  A fitted model object of class `"islasso"`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  The confidence level required.

- type.ci:

  Character. Only Wald-type confidence intervals are implemented yet!
  Set `type.ci = "wald"` to use estimates and standard errors to build
  the confidence interval.

- trace:

  Not used (kept for compatibility).

- ...:

  Additional arguments for methods.

## Details

confint method for islasso objects

## See also

[`islasso.fit`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`summary.islasso`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[`residuals.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`logLik.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`predict.islasso`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.md),
[`deviance.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)

## Author

Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
n <- 100; p <- 100; p1 <- 10
beta.veri <- sort(round(c(seq(0.5, 3, length.out = p1 / 2),
                          seq(-1, -2, length.out = p1 / 2)), 2))
beta <- c(beta.veri, rep(0, p - p1))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
o <- islasso(y ~ ., data = sim$data, family = gaussian())

ci <- confint(o, type.ci = "wald", parm = 1:11)
ci
#>              2.5 %    97.5 % 
#> (Intercept)  -0.1702   0.1876
#> X1           -2.1114  -1.7209
#> X2           -1.7961  -1.4116
#> X3           -1.6218  -1.2331
#> X4           -1.3961  -0.9959
#> X5           -0.9985  -0.6076
#> X6            0.2055   0.5953
#> X7            0.7388   1.1222
#> X8            1.4321   1.8268
#> X9            2.1133   2.5240
#> X10           2.6710   3.0764
plot(ci)

```
