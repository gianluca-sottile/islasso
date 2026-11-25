# General Linear Hypotheses for islasso Models

Tests general linear hypotheses and computes confidence intervals for
linear combinations of coefficients from a fitted
[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model.

## Usage

``` r
# S3 method for class 'islasso'
anova(object, A, b = NULL, ci, ...)
```

## Arguments

- object:

  A fitted model object of class `"islasso"`.

- A:

  Hypothesis specification. Either:

  - A numeric matrix or vector with each row specifying a linear
    combination of coefficients,

  - Or a character vector with symbolic expressions (e.g.
    `"X1 + X2 = 3"`).

- b:

  Right-hand side vector for the null hypotheses `A %*% beta = b`. If
  omitted, defaults to zeros.

- ci:

  Optional 2-column matrix of confidence intervals for coefficients.

- ...:

  Currently unused.

## Value

An object of class `"anova.islasso"` containing:

- Estimate:

  Linear combination estimates

- SE:

  Standard errors

- Wald:

  Wald statistics

- p-value:

  Associated p-values

## Details

The method tests the null hypothesis \\H_0: A \beta = b\\, where \\A\\
and \\b\\ define a linear constraint on model coefficients.

Symbolic expressions support natural syntax: coefficients may be
added/subtracted, constants may be multiplied (e.g.
`"2 * X1 + 3 * X2 = 7"`). Equations with omitted `=` assume zero on the
right-hand side. See examples for syntax flexibility.

## See also

[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`summary.islasso`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[`confint.islasso`](https://gianluca-sottile.github.io/islasso/reference/confint.islasso.md),
[`predict.islasso`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.md),
[`logLik.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`residuals.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
n <- 100; p <- 100
beta <- c(runif(10, -2, 2), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso(y ~ . -1, data = sim$data, family = gaussian())

# Test if first 5 variables sum to -7.5
anova(fit, A = c("X1 + X2 + X3 + X4 + X5 = -7.5"))
#> 
#>   Simultaneous Tests for General Linear Combination
#> 
#> Fit: islasso(formula = y ~ . - 1, family = gaussian(), data = sim$data)
#> 
#> Linear Hypothesis:
#>                                         Estimate Std. Error Chisq Pr(>Chisq)
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5   4.9635     0.2072 573.9     <2e-16
#>                                            
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Test multiple hypotheses
anova(fit, A = c("X1 + X2 + X3 + X4 + X5", "X6 + X7 + X8 + X9 + X10"), b = c(-7.5, 8.75))
#> 
#>   Simultaneous Tests for General Linear Combination
#> 
#> Fit: islasso(formula = y ~ . - 1, family = gaussian(), data = sim$data)
#> 
#> Linear Hypothesis:
#>                                          Estimate Std. Error  Chisq Pr(>Chisq)
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5    4.9635     0.2072  573.9     <2e-16
#> 1*X6 + 1*X7 + 1*X8 + 1*X9 + 1*X10 = 8.75 -10.2064     0.2063 2447.6     <2e-16
#> Overall                                                      3172.8     <2e-16
#>                                             
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5  ***
#> 1*X6 + 1*X7 + 1*X8 + 1*X9 + 1*X10 = 8.75 ***
#> Overall                                  ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Full diagonal comparison to true coefficients
anova(fit, A = diag(p), b = beta)
#> Error in solve.default(vcov.hyp, tol = 1e-20): system is computationally singular: reciprocal condition number = 5.2108e-22
```
