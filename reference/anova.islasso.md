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
set.seed(1)
beta <- c(runif(10, -2, 2), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso(y ~ . -1, data = sim$data, family = gaussian(), lambda = 10)

# Test if first 5 variables sum to -7.5
anova(fit, A = c("X1 + X2 + X3 + X4 + X5 = -7.5"))
#> 
#>   Simultaneous Tests for General Linear Combination
#> 
#> Fit: islasso(formula = y ~ . - 1, family = gaussian(), lambda = 10, 
#>     data = sim$data)
#> 
#> Linear Hypothesis:
#>                                         Estimate Std. Error Chisq Pr(>Chisq)
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5    6.809      0.203  1125     <2e-16
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
#> Fit: islasso(formula = y ~ . - 1, family = gaussian(), lambda = 10, 
#>     data = sim$data)
#> 
#> Linear Hypothesis:
#>                                          Estimate Std. Error  Chisq Pr(>Chisq)
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5    6.8090     0.2030 1125.0     <2e-16
#> 1*X6 + 1*X7 + 1*X8 + 1*X9 + 1*X10 = 8.75  -6.3171     0.2058  942.4     <2e-16
#> Overall                                                      2162.0     <2e-16
#>                                             
#> 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 = -7.5  ***
#> 1*X6 + 1*X7 + 1*X8 + 1*X9 + 1*X10 = 8.75 ***
#> Overall                                  ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Full diagonal comparison to true coefficients
anova(fit, A = diag(p), b = beta)
#> Warning: the covariance matrix of the hypotheses is singular: the overall test uses 99 degrees of freedom instead of 100
#> 
#>   Simultaneous Tests for General Linear Combination
#> 
#> Fit: islasso(formula = y ~ . - 1, family = gaussian(), lambda = 10, 
#>     data = sim$data)
#> 
#> Linear Hypothesis:
#>                             Estimate Std. Error  Chisq Pr(>Chisq)  
#> 1*X1 = -0.9379653474316    4.205e-02  9.645e-02  0.190     0.6629  
#> 1*X2 = -0.511504401452839  1.225e-01  9.532e-02  1.652     0.1987  
#> 1*X3 = 0.291413453407586  -1.033e-01  8.990e-02  1.320     0.2506  
#> 1*X4 = 1.6328311599791    -2.495e-01  1.003e-01  6.190     0.0128 *
#> 1*X5 = -1.19327227585018   2.157e-01  9.647e-02  5.000     0.0253 *
#> 1*X6 = 1.59355873987079   -1.340e-01  9.662e-02  1.923     0.1656  
#> 1*X7 = 1.77870107442141   -1.883e-01  9.563e-02  3.879     0.0489 *
#> 1*X8 = 0.643191169947386  -1.336e-01  9.841e-02  1.843     0.1746  
#> 1*X9 = 0.516456175595522  -2.479e-02  9.926e-02  0.062     0.8028  
#> 1*X10 = -1.75285491812974  1.346e-01  9.991e-02  1.815     0.1780  
#> 1*X11 = 0                  4.254e-02  5.487e-02  0.601     0.4381  
#> 1*X12 = 0                  5.955e-07  1.204e-02  0.000     1.0000  
#> 1*X13 = 0                 -7.301e-03  3.562e-02  0.042     0.8376  
#> 1*X14 = 0                  2.460e-07  5.019e-03  0.000     1.0000  
#> 1*X15 = 0                  1.613e-03  3.928e-02  0.002     0.9672  
#> 1*X16 = 0                 -9.555e-03  4.032e-02  0.056     0.8127  
#> 1*X17 = 0                 -4.447e-07  8.972e-03  0.000     1.0000  
#> 1*X18 = 0                 -7.664e-07  1.530e-02  0.000     1.0000  
#> 1*X19 = 0                 -5.107e-07  1.026e-02  0.000     1.0000  
#> 1*X20 = 0                  8.008e-03  4.173e-02  0.037     0.8478  
#> 1*X21 = 0                  6.078e-07  1.219e-02  0.000     1.0000  
#> 1*X22 = 0                 -8.119e-07  1.620e-02  0.000     1.0000  
#> 1*X23 = 0                 -2.009e-03  3.822e-02  0.003     0.9581  
#> 1*X24 = 0                  6.325e-02  6.289e-02  1.012     0.3145  
#> 1*X25 = 0                  6.542e-07  1.316e-02  0.000     1.0000  
#> 1*X26 = 0                  1.922e-02  4.490e-02  0.183     0.6686  
#> 1*X27 = 0                 -4.477e-07  9.066e-03  0.000     1.0000  
#> 1*X28 = 0                  1.249e-07  2.543e-03  0.000     1.0000  
#> 1*X29 = 0                  1.949e-02  4.409e-02  0.195     0.6585  
#> 1*X30 = 0                  1.632e-02  4.431e-02  0.136     0.7127  
#> 1*X31 = 0                 -1.818e-02  4.552e-02  0.159     0.6897  
#> 1*X32 = 0                 -8.995e-03  3.962e-02  0.052     0.8204  
#> 1*X33 = 0                 -2.552e-08  9.566e-04  0.000     1.0000  
#> 1*X34 = 0                  6.742e-07  1.344e-02  0.000     1.0000  
#> 1*X35 = 0                  7.734e-07  1.557e-02  0.000     1.0000  
#> 1*X36 = 0                  9.233e-02  7.290e-02  1.604     0.2053  
#> 1*X37 = 0                  9.964e-07  1.977e-02  0.000     1.0000  
#> 1*X38 = 0                 -8.666e-08  1.760e-03  0.000     1.0000  
#> 1*X39 = 0                 -2.747e-07  5.574e-03  0.000     1.0000  
#> 1*X40 = 0                  7.342e-08  1.484e-03  0.000     1.0000  
#> 1*X41 = 0                  3.369e-02  4.984e-02  0.457     0.4990  
#> 1*X42 = 0                  1.605e-02  4.103e-02  0.153     0.6957  
#> 1*X43 = 0                  1.480e-03  4.036e-02  0.001     0.9707  
#> 1*X44 = 0                 -9.117e-07  1.822e-02  0.000     1.0000  
#> 1*X45 = 0                 -7.791e-03  4.252e-02  0.034     0.8546  
#> 1*X46 = 0                  4.041e-07  8.227e-03  0.000     1.0000  
#> 1*X47 = 0                  3.004e-02  4.875e-02  0.380     0.5378  
#> 1*X48 = 0                  2.235e-02  4.573e-02  0.239     0.6250  
#> 1*X49 = 0                 -2.050e-02  4.509e-02  0.207     0.6494  
#> 1*X50 = 0                  1.646e-02  4.196e-02  0.154     0.6949  
#> 1*X51 = 0                  9.772e-07  1.940e-02  0.000     1.0000  
#> 1*X52 = 0                  2.954e-03  3.702e-02  0.006     0.9364  
#> 1*X53 = 0                 -5.009e-03  4.103e-02  0.015     0.9028  
#> 1*X54 = 0                  2.837e-02  4.932e-02  0.331     0.5651  
#> 1*X55 = 0                 -8.911e-07  1.777e-02  0.000     1.0000  
#> 1*X56 = 0                 -4.016e-02  5.380e-02  0.557     0.4553  
#> 1*X57 = 0                 -4.042e-07  8.160e-03  0.000     1.0000  
#> 1*X58 = 0                  2.069e-02  4.519e-02  0.210     0.6471  
#> 1*X59 = 0                 -3.844e-07  7.864e-03  0.000     1.0000  
#> 1*X60 = 0                 -8.341e-07  1.666e-02  0.000     1.0000  
#> 1*X61 = 0                  6.386e-07  1.280e-02  0.000     1.0000  
#> 1*X62 = 0                  2.996e-03  3.997e-02  0.006     0.9402  
#> 1*X63 = 0                  3.103e-02  4.842e-02  0.411     0.5216  
#> 1*X64 = 0                  4.698e-07  9.475e-03  0.000     1.0000  
#> 1*X65 = 0                 -6.999e-07  1.391e-02  0.000     1.0000  
#> 1*X66 = 0                  1.169e-07  2.371e-03  0.000     1.0000  
#> 1*X67 = 0                  5.183e-02  5.816e-02  0.794     0.3728  
#> 1*X68 = 0                 -1.269e-07  2.568e-03  0.000     1.0000  
#> 1*X69 = 0                  1.736e-07  3.547e-03  0.000     1.0000  
#> 1*X70 = 0                 -1.491e-08  9.566e-04  0.000     1.0000  
#> 1*X71 = 0                  1.377e-07  2.807e-03  0.000     1.0000  
#> 1*X72 = 0                  1.320e-02  3.991e-02  0.109     0.7409  
#> 1*X73 = 0                 -5.443e-07  1.092e-02  0.000     1.0000  
#> 1*X74 = 0                  8.937e-07  1.772e-02  0.000     1.0000  
#> 1*X75 = 0                  9.857e-07  1.975e-02  0.000     1.0000  
#> 1*X76 = 0                  5.923e-07  1.180e-02  0.000     1.0000  
#> 1*X77 = 0                  3.030e-07  6.178e-03  0.000     1.0000  
#> 1*X78 = 0                  4.025e-03  3.955e-02  0.010     0.9189  
#> 1*X79 = 0                  1.944e-03  3.918e-02  0.002     0.9604  
#> 1*X80 = 0                 -2.316e-02  4.547e-02  0.259     0.6105  
#> 1*X81 = 0                  4.053e-08  9.568e-04  0.000     1.0000  
#> 1*X82 = 0                 -4.869e-08  9.981e-04  0.000     1.0000  
#> 1*X83 = 0                 -2.036e-02  4.572e-02  0.198     0.6560  
#> 1*X84 = 0                  3.114e-07  6.316e-03  0.000     1.0000  
#> 1*X85 = 0                  2.129e-02  4.728e-02  0.203     0.6525  
#> 1*X86 = 0                  3.462e-07  6.939e-03  0.000     1.0000  
#> 1*X87 = 0                 -7.605e-02  6.689e-02  1.293     0.2556  
#> 1*X88 = 0                 -9.728e-07  1.932e-02  0.000     1.0000  
#> 1*X89 = 0                 -8.369e-02  6.895e-02  1.473     0.2249  
#> 1*X90 = 0                 -3.837e-07  7.763e-03  0.000     1.0000  
#> 1*X91 = 0                 -4.337e-07  8.814e-03  0.000     1.0000  
#> 1*X92 = 0                 -6.018e-02  6.146e-02  0.959     0.3275  
#> 1*X93 = 0                  9.151e-07  1.797e-02  0.000     1.0000  
#> 1*X94 = 0                 -3.373e-02  5.126e-02  0.433     0.5106  
#> 1*X95 = 0                  2.608e-07  5.313e-03  0.000     1.0000  
#> 1*X96 = 0                  9.779e-07  1.939e-02  0.000     1.0000  
#> 1*X97 = 0                  7.952e-07  1.590e-02  0.000     1.0000  
#> 1*X98 = 0                 -2.740e-07  5.464e-03  0.000     1.0000  
#> 1*X99 = 0                 -7.658e-07  1.535e-02  0.000     1.0000  
#> 1*X100 = 0                 9.703e-07  1.938e-02  0.000     1.0000  
#> Overall                                         47.504     1.0000  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
