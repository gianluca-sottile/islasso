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

  Logical. If `TRUE` (default), a bar shows the iterations status.

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

Maintainer: Gianluca Sottile \<gianluca.sottile@unipa.it\>

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
#> (Intercept)  -0.1669   0.1843
#> X1           -2.1210  -1.7335
#> X2           -1.8039  -1.4235
#> X3           -1.6206  -1.2367
#> X4           -1.4060  -1.0086
#> X5           -1.0058  -0.6194
#> X6            0.2127   0.5972
#> X7            0.7509   1.1291
#> X8            1.4426   1.8342
#> X9            2.1278   2.5356
#> X10           2.6815   3.0834
plot(ci)

```
