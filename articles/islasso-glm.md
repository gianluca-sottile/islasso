# Generalised linear models and families in islasso

From version 1.7.0 the generalised linear models are fitted by
iteratively reweighted least squares (IRLS) driven by the `family`
object, exactly as in [`glm()`](https://rdrr.io/r/stats/glm.html): any
family with the components `linkinv`, `mu.eta`, `variance` and
`dev.resids` can be used, including all the links of the standard
families, the `quasi` families and, for instance,
[`MASS::negative.binomial()`](https://rdrr.io/pkg/MASS/man/negative.binomial.html)
with known \theta. At each iteration the working weights w_i =
\mu'(\eta_i)^2/V(\mu_i) and the working response z_i = \eta_i + (y_i -
\mu_i)/\mu'(\eta_i) define the weighted least squares problem solved by
the induced smoothed lasso; the step is halved whenever the linear
predictor leaves the valid region of the family.

The dispersion parameter is fixed to one for the binomial, Poisson and
negative binomial families and is otherwise estimated by the Pearson
statistic divided by the residual degrees of freedom, n - \sum_j h_j.

``` r

library(islasso)
set.seed(1)
n <- 200; p <- 30
X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
eta <- -0.3 + 0.8 * X[, 1] - 0.6 * X[, 2] + 0.4 * X[, 3]
```

## Binomial responses with different links

``` r

d <- data.frame(y = rbinom(n, 1, pnorm(eta)), X)
fit.probit <- islasso(y ~ ., data = d, family = binomial(link = "probit"), lambda = 4)
summary(fit.probit, pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., family = binomial(link = "probit"), 
#>     lambda = 4, data = d)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.4789 -0.7513 -0.2861  0.7826  2.3453 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)  -0.1715     0.1104  1.000  -1.553 0.120357    
#> X1            0.6401     0.1298  1.000   4.931 8.16e-07 ***
#> X2           -0.6386     0.1250  1.000  -5.109 3.23e-07 ***
#> X3            0.3796     0.1101  0.999   3.446 0.000569 ***
#> X4            0.2622     0.1044  0.992   2.512 0.012005 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 274.37  on 199.0  degrees of freedom
#> Residual deviance: 166.07  on 174.9  degrees of freedom
#> AIC: 216.26
#> Lambda: 4
#> 
#> Number of IS iterations: 193
fit.cloglog <- islasso(y ~ ., data = d, family = binomial(link = "cloglog"), lambda = 4)
coef(fit.cloglog)[1:5]
#> (Intercept)          X1          X2          X3          X4 
#>  -0.6236215   0.7818604  -0.6849379   0.4068574   0.2655697
```

With a (practically) null tuning parameter the estimates coincide with
those of [`glm()`](https://rdrr.io/r/stats/glm.html):

``` r

f0 <- islasso(y ~ X1 + X2 + X3, data = d, family = binomial(link = "probit"), lambda = 1e-8)
g0 <- glm(y ~ X1 + X2 + X3, data = d, family = binomial(link = "probit"))
cbind(islasso = coef(f0), glm = coef(g0))
#>                islasso        glm
#> (Intercept) -0.1906965 -0.1906967
#> X1           0.6653830  0.6653832
#> X2          -0.5901572 -0.5901570
#> X3           0.4283252  0.4283252
```

## Counts: Poisson with an offset, quasi-Poisson and negative binomial

``` r

expo <- runif(n, 1, 3)
mu <- expo * exp(eta)
d2 <- data.frame(y = rpois(n, mu), X, expo = expo)
fit.pois <- islasso(y ~ . - expo + offset(log(expo)), data = d2, family = poisson(), lambda = 5)
summary(fit.pois, pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ . - expo + offset(log(expo)), family = poisson(), 
#>     lambda = 5, data = d2)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.0393 -0.8853 -0.3327  0.5385  2.1021 
#> 
#>             Estimate Std. Error      Df z value Pr(>|z|)    
#> (Intercept) -0.35327    0.07236   1.000  -4.882 1.05e-06 ***
#> X1           0.88050    0.05745   1.000  15.327  < 2e-16 ***
#> X2          -0.62203    0.05015   1.000 -12.404  < 2e-16 ***
#> X3           0.34293    0.04445   1.000   7.715 1.21e-14 ***
#> X16          0.11240    0.05033   0.990   2.233   0.0255 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 843.26  on 199.0  degrees of freedom
#> Residual deviance: 182.22  on 173.7  degrees of freedom
#> AIC: 630.89
#> Lambda: 5
#> 
#> Number of IS iterations: 77

fit.qp <- islasso(y ~ . - expo + offset(log(expo)), data = d2, family = quasipoisson(), lambda = 5)
fit.qp$dispersion
#> [1] 1.034923
```

``` r

d3 <- data.frame(y = MASS::rnegbin(n, mu = exp(1 + eta), theta = 2), X)
fit.nb <- islasso(y ~ ., data = d3, family = MASS::negative.binomial(theta = 2), lambda = 5)
summary(fit.nb, pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., family = MASS::negative.binomial(theta = 2), 
#>     lambda = 5, data = d3)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.4373 -0.9239 -0.2181  0.4861  2.0752 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)  0.86000    0.07705  1.000  11.162  < 2e-16 ***
#> X1           0.63112    0.08059  1.000   7.831 4.83e-15 ***
#> X2          -0.56384    0.07471  1.000  -7.547 4.46e-14 ***
#> X3           0.30230    0.06703  1.000   4.510 6.49e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for Negative Binomial(2) family taken to be 1)
#> 
#>     Null deviance: 395.21  on 199.0  degrees of freedom
#> Residual deviance: 172.38  on 175.3  degrees of freedom
#> AIC: 838.09
#> Lambda: 5
#> 
#> Number of IS iterations: 74
```

## Positive continuous responses

``` r

d4 <- data.frame(y = rgamma(n, shape = 4, scale = exp(1 + 0.5 * eta) / 4), X)
fit.gamma <- islasso(y ~ ., data = d4, family = Gamma(link = "log"), lambda = 3)
fit.ig <- islasso(y ~ ., data = d4, family = inverse.gaussian(link = "log"), lambda = 3)
round(cbind(Gamma = coef(fit.gamma), inv.gaussian = coef(fit.ig))[1:6, ], 3)
#>              Gamma inv.gaussian
#> (Intercept)  0.819        0.794
#> X1           0.342        0.347
#> X2          -0.377       -0.378
#> X3           0.156        0.151
#> X4          -0.011       -0.029
#> X5          -0.049       -0.046
c(Gamma = fit.gamma$dispersion, inv.gaussian = fit.ig$dispersion)
#>        Gamma inv.gaussian 
#>    0.2596685    0.1258368
```

## Regularisation path

[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
accepts the same families; the information criteria use the
log-likelihood when the family provides one (`family$aic`), and the
deviance otherwise.

``` r

path <- islasso.path(y ~ ., data = d, family = binomial(link = "probit"), nlambda = 20)
GoF.islasso.path(path, plot = FALSE)$lambda.min
#>      AIC      BIC     AICc     eBIC      GCV      GIC 
#>  8.55741 22.56290  8.55741 36.63712  8.55741 22.56290
```
