# High-dimensional data and computational aspects

The induced smoothed lasso works with full p \times p matrices (the
derivative of the smoothed estimating equation and the sandwich
covariance matrix), so its cost grows with the number of covariates.
This vignette summarises how **islasso** keeps the computations
manageable.

## Cost per iteration

| step | until 1.6.2 | from 1.7.0, p \le 2n | from 1.7.0, p \> 2n |
|:---|:--:|:--:|:--:|
| estimating equation at fixed s.e. | damped fixed point, one p\times p Cholesky per step | undamped MM (or Newton) steps, one Cholesky per step | n \times n factorisations (Woodbury), O(n^2p) |
| standard errors | full sandwich, \approx 4.3\\p^3 flops | diagonal only, \approx p^3 flops | diagonal only, O(n^2p) |
| X^\top W X (Gaussian path) | recomputed at every \lambda | computed once | not needed |
| Gaussian iterations | Fortran | entirely in C | entirely in C |
| full covariance matrix | every iteration | once, at convergence | once, O(np^2) |

The number of steps needed to solve the estimating equation is also
smaller, because the majorise–minimise iterations are not damped any
more (`inner = "mm"`, default) or use Newton’s method
(`inner = "newton"`).

## An example with p \gg n

``` r

library(islasso)
n <- 60; p <- 600
beta <- c(2, -2, 1.5, -1.5, 1, rep(0, p - 5))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, rho = 0.3)
system.time(fit <- islasso(y ~ ., data = sim$data, lambda = 10))
#>    user  system elapsed 
#>   6.843   0.103   2.521
fit$internal$mode   # 1 = Woodbury back end
#> [1] 1
summary(fit, pval = 0.01)
#> 
#> Call:
#> islasso(formula = y ~ ., lambda = 10, data = sim$data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.33508 -0.44950  0.03032  0.35238  1.51862 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)   0.1169     0.1122  1.000   1.042  0.29742    
#> X1            1.5168     0.1409  1.000  10.768  < 2e-16 ***
#> X2           -0.6826     0.1633  1.000  -4.179 2.93e-05 ***
#> X3            0.4499     0.1565  0.976   2.874  0.00405 ** 
#> X4           -0.7858     0.1540  1.000  -5.102 3.36e-07 ***
#> X5            0.5313     0.1359  0.999   3.911 9.21e-05 ***
#> X116         -0.3691     0.1301  0.982  -2.838  0.00454 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.755663)
#> 
#>     Null deviance: 695.983  on 59.00  degrees of freedom
#> Residual deviance:  19.513  on 25.82  degrees of freedom
#> AIC: 87.868
#> Lambda: 10
#> 
#> Number of IS iterations: 166
```

The same model with the dense back end (one 601 \times 601 factorisation
per step):

``` r

system.time(fit.d <- islasso(y ~ ., data = sim$data, lambda = 10,
                             control = is.control(solver = "dense")))
#>    user  system elapsed 
#>  35.696   0.178   9.170
max(abs(coef(fit) - coef(fit.d)))
#> [1] 2.706863e-13
```

## Practical advice

- Use
  [`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
  to explore the tuning parameter: with the Gaussian family X^\top W X
  is computed once, and the models for the different values of `lambda`
  are independent, so they can be fitted in parallel with
  `control = is.control(ncores = 4)`.
- `control = is.control(inner = "newton")` usually needs fewer
  factorisations when p is large.
- An optimised BLAS/LAPACK library (OpenBLAS, Apple Accelerate, Intel
  MKL) speeds up all the computations, since the numerical core relies
  on BLAS level 3 and LAPACK routines.
- The full covariance matrix is p \times p: with tens of thousands of
  covariates its storage (8 bytes per entry) becomes the limiting
  factor.
