# The Induced Smoothed Lasso

Fits regression models with a smoothed lasso (or elastic-net) penalty
under the induced smoothing paradigm, which yields point estimates,
standard errors and Wald tests for all the coefficients.

## Usage

``` r
islasso(
  formula,
  family = gaussian,
  lambda,
  alpha = 1,
  data,
  weights,
  subset,
  offset,
  unpenalized,
  contrasts = NULL,
  control = is.control()
)

islasso.fit(
  X,
  y,
  family = gaussian,
  lambda,
  alpha = 1,
  intercept = FALSE,
  weights = NULL,
  offset = NULL,
  unpenalized = NULL,
  control = is.control()
)
```

## Arguments

- formula:

  a symbolic description of the model (see
  [`stats::formula()`](https://rdrr.io/r/stats/formula.html)).

- family:

  the error distribution and link function: a family object, a family
  function or its name. Any family accepted by
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) can be used, for
  instance [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html) (links logit,
  probit, cauchit, log, cloglog),
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`Gamma()`](https://rdrr.io/r/stats/family.html),
  [`inverse.gaussian()`](https://rdrr.io/r/stats/family.html), the
  `quasi*` families and
  [`MASS::negative.binomial()`](https://rdrr.io/pkg/MASS/man/negative.binomial.html)
  with known `theta`.

- lambda:

  the tuning parameter (a non-negative number). If missing, it is chosen
  by cross-validation with
  [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  (see `nfolds` and `seed` in
  [`is.control()`](https://gianluca-sottile.github.io/islasso/reference/is.control.md)).
  With `lambda = 0` the model is not penalized (same estimates as
  [`glm()`](https://rdrr.io/r/stats/glm.html)) and the model matrix must
  have full column rank: aliased terms (`NA` coefficients in
  [`glm()`](https://rdrr.io/r/stats/glm.html)) give an error.

- alpha:

  the elastic-net mixing parameter, \\0 \le \alpha \le 1\\: `alpha = 1`
  is the lasso, `alpha = 0` the ridge penalty.

- data:

  an optional data frame containing the variables in the model.

- weights:

  optional vector of prior weights.

- subset:

  an optional vector specifying a subset of observations.

- offset:

  an optional offset in the linear predictor.

- unpenalized:

  names or indices (with respect to the columns of the model matrix
  without intercept) of the coefficients left unpenalized. The intercept
  is never penalized.

- contrasts:

  an optional list of contrasts for the factors (see
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)).

- control:

  a list of control parameters, see
  [`is.control()`](https://gianluca-sottile.github.io/islasso/reference/is.control.md).

- X, y:

  for `islasso.fit()`: the model matrix (without the column of ones) and
  the response.

- intercept:

  for `islasso.fit()`: logical, should an (unpenalized) intercept be
  included?

## Value

An object of class `"islasso"`, a list with components including

- coefficients, se:

  estimates and standard errors on the original scale;

- dispersion:

  the dispersion parameter;

- fitted.values, linear.predictors, residuals, weights, prior.weights:

  as in [`stats::glm()`](https://rdrr.io/r/stats/glm.html) (`residuals`
  are the working residuals and `weights` the working weights);

- deviance, null.deviance, aic, rank, df.residual, df.null:

  deviance-based summaries; `rank` is the effective number of
  parameters, the sum of the pseudo hat values;

- converged, iter:

  convergence flag and number of iterations;

- lambda, alpha, family, control, call, formula, terms, model, data,
  xlevels, contrasts, offset, y:

  the inputs;

- internal:

  other internal quantities (full covariance matrix, pseudo hat values
  `hi`, gradient, mixture weights, standardisation constants, ...).

## Details

For a fixed \\\lambda\\ the lasso estimating equation \\-X^\top W (z -
X\beta) + \lambda\\\mathrm{sign}(\beta) = 0\\ is not smooth. The induced
smoothing (IS) approach (Brown and Wang, 2005) replaces
\\\mathrm{sign}(\beta_j)\\ by its expectation under the approximate
sampling distribution of the estimator, modelled as the mixture
\\c\\N(0,1) + (1-c)\\N(0,10^{-10})\\ of a continuous component and a
component concentrated at zero: \$\$P_j = c\\\\2\Phi(\beta_j/s_j) -
1\\ + (1-c)\\\\2\Phi(10^5 \beta_j/s_j) - 1\\,\$\$ where \\s_j\\ is the
standard error of \\\hat\beta_j\\. The smoothed equation \\\tilde
U(\beta) = -X^\top W(z - X\beta) + \lambda\\\alpha P + (1-\alpha)\beta\\
= 0\\ admits the derivative \\\tilde U'(\beta) = X^\top W X +
\lambda\\\mathrm{diag}\\\alpha P' + 1 - \alpha\\\\, and the covariance
matrix is obtained by the sandwich formula \\V = \phi\\\tilde U'^{-1}
X^\top W X\\ \tilde U'^{-1}\\. Point estimates and standard errors are
computed alternately until the standard errors stabilise.

At fixed standard errors, \\\tilde U\\ is the gradient of a smooth
convex function, which is minimised by a monotone majorise-minimise
algorithm (default) or by Newton's method with line search
(`inner = "newton"` in
[`is.control()`](https://gianluca-sottile.github.io/islasso/reference/is.control.md)).
When the number of columns exceeds twice the number of observations, the
linear algebra uses the Woodbury identity, with cost \\O(n^2 p)\\
instead of \\O(p^3)\\. Generalised linear models are fitted by
iteratively reweighted least squares, so any family with the usual
components (`linkinv`, `mu.eta`, `variance`, `dev.resids`) is supported.

With `stand = TRUE` (default in
[`is.control()`](https://gianluca-sottile.github.io/islasso/reference/is.control.md))
the covariates are standardised (weighted mean and standard deviation
with divisor \\\sum w\\, as in glmnet) and the penalty acts on the
standardised coefficients; the results are reported on the original
scale. The scale of `lambda` is that of the deviance: \\\lambda = \sum_i
w_i\\\lambda\_{glmnet}\\.

Since the smoothed penalty is never exactly flat, no estimate is exactly
zero: the Wald statistics, rather than the point estimates, are used to
judge which covariates matter. For coefficients whose estimate is shrunk
to (practically) zero the induced-smoothing fixed point drives the
standard error towards zero as well, while the Wald ratio stays stable:
their standard errors depend on the convergence tolerance and should not
be used to build confidence intervals. See
[`vignette("islasso-methodology", package = "islasso")`](https://gianluca-sottile.github.io/islasso/articles/islasso-methodology.md).

## References

Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2020). The induced
smoothed lasso: a practical framework for hypothesis testing in high
dimensional regression. *Statistical Methods in Medical Research*,
**29**(3), 765–777.
[doi:10.1177/0962280219842890](https://doi.org/10.1177/0962280219842890)

Sottile G., Cilluffo G., Muggeo V.M.R. (2019). The R package islasso:
estimation and hypothesis testing in lasso regression. Technical report.
[doi:10.13140/RG.2.2.16360.11521](https://doi.org/10.13140/RG.2.2.16360.11521)

Brown B.M., Wang Y.G. (2005). Standard errors and covariance matrices
for smoothed rank estimators. *Biometrika*, **92**(1), 149–158.
[doi:10.1093/biomet/92.1.149](https://doi.org/10.1093/biomet/92.1.149)

## See also

[`summary.islasso()`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[`confint.islasso()`](https://gianluca-sottile.github.io/islasso/reference/confint.islasso.md),
[`predict.islasso()`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.md),
[`anova.islasso()`](https://gianluca-sottile.github.io/islasso/reference/anova.islasso.md),
[`plot.islasso()`](https://gianluca-sottile.github.io/islasso/reference/plot.islasso.md),
[`aic.islasso()`](https://gianluca-sottile.github.io/islasso/reference/aic.islasso.md),
[`relax.islasso()`](https://gianluca-sottile.github.io/islasso/reference/relax.islasso.md),
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`is.control()`](https://gianluca-sottile.github.io/islasso/reference/is.control.md),
[`simulXy()`](https://gianluca-sottile.github.io/islasso/reference/simulXy.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
n <- 100; p <- 50
beta <- c(rep(1, 5), rep(0, p - 5))
sim1 <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
o <- islasso(y ~ ., data = sim1$data, family = gaussian(), lambda = 10)
summary(o, pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., family = gaussian(), lambda = 10, data = sim1$data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.62155 -0.60848 -0.04476  0.51511  1.79945 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept) -0.11060    0.09886  1.000  -1.119    0.263    
#> X1           0.82242    0.10738  1.000   7.659 1.87e-14 ***
#> X2           0.95762    0.10829  1.000   8.843  < 2e-16 ***
#> X3           0.81484    0.10836  1.000   7.520 5.50e-14 ***
#> X4           0.99564    0.11054  1.000   9.007  < 2e-16 ***
#> X5           0.91294    0.10743  1.000   8.498  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.9772461)
#> 
#>     Null deviance: 657.721  on 99.00  degrees of freedom
#> Residual deviance:  70.928  on 72.58  degrees of freedom
#> AIC: 306.28
#> Lambda: 10
#> 
#> Number of IS iterations: 238
#> 
head(coef(o))
#> (Intercept)          X1          X2          X3          X4          X5 
#>  -0.1106006   0.8224221   0.9576181   0.8148399   0.9956382   0.9129428 
logLik(o)
#> 'log Lik.' -124.7188 (df=28.4205)

# \donttest{
# lambda chosen by 5-fold cross-validation (reproducible through 'seed')
o <- islasso(y ~ ., data = sim1$data, control = is.control(seed = 1))
summary(relax.islasso(o), pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., lambda = 1841.92420348731, data = sim1$data, 
#>     unpenalized = c("X1", "X2", "X3", "X4", "X5"), control = is.control(seed = 1))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.84475 -0.67994 -0.08082  0.65141  2.21448 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)  -0.1106     0.1011  1.000  -1.094    0.274    
#> X1            0.9211     0.1029  1.000   8.948   <2e-16 ***
#> X2            1.0819     0.1029  1.000  10.514   <2e-16 ***
#> X3            0.9715     0.1024  1.000   9.485   <2e-16 ***
#> X4            1.0924     0.1026  1.000  10.647   <2e-16 ***
#> X5            0.9702     0.1038  1.000   9.343   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 1.022312)
#> 
#>     Null deviance: 657.721  on 99  degrees of freedom
#> Residual deviance:  96.097  on 94  degrees of freedom
#> AIC: 293.81
#> Lambda: 1841.9
#> 
#> Number of IS iterations: 76
#> 

# binomial
beta <- c(1, 1, 1, rep(0, p - 3))
sim2 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1, size = 100,
                family = binomial())
o2 <- islasso(cbind(y.success, y.failure) ~ ., data = sim2$data, family = binomial(),
              lambda = 5)
summary(o2, pval = 0.05)
#> 
#> Call:
#> islasso(formula = cbind(y.success, y.failure) ~ ., family = binomial(), 
#>     lambda = 5, data = sim2$data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.82965 -0.29167  0.08265  0.44092  2.05016 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)  0.99044    0.03125  1.000  31.695  < 2e-16 ***
#> X1           1.02907    0.04334  1.000  23.747  < 2e-16 ***
#> X2           1.03592    0.04286  1.000  24.172  < 2e-16 ***
#> X3           1.00597    0.03859  1.000  26.071  < 2e-16 ***
#> X33          0.10786    0.03830  0.999   2.816  0.00486 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 3876.288  on 99.00  degrees of freedom
#> Residual deviance:   54.867  on 52.35  degrees of freedom
#> AIC: 574.12
#> Lambda: 5
#> 
#> Number of IS iterations: 50
#> 

# Poisson
sim3 <- simulXy(n = n, p = p, beta = beta, interc = 1, seed = 1, family = poisson())
o3 <- islasso(y ~ ., data = sim3$data, family = poisson(), lambda = 5)
summary(o3, pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., family = poisson(), lambda = 5, data = sim3$data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.69225 -0.74072  0.01484  0.28906  1.82577 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept)  1.08907    0.07652  1.000  14.232   <2e-16 ***
#> X1           0.81148    0.07744  1.000  10.479   <2e-16 ***
#> X2           1.05083    0.05831  1.000  18.022   <2e-16 ***
#> X3           0.87313    0.07075  1.000  12.341   <2e-16 ***
#> X45          0.15435    0.06264  0.993   2.464   0.0137 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 1734.954  on 99.00  degrees of freedom
#> Residual deviance:   52.619  on 59.94  degrees of freedom
#> AIC: 409.51
#> Lambda: 5
#> 
#> Number of IS iterations: 76
#> 

# Gamma with log link
sim4 <- simulXy(n = n, p = p, beta = beta, interc = -1, seed = 1, dispersion = 0.1,
                family = Gamma(link = "log"))
o4 <- islasso(y ~ ., data = sim4$data, family = Gamma(link = "log"), lambda = 5)
summary(o4, pval = 0.05)
#> 
#> Call:
#> islasso(formula = y ~ ., family = Gamma(link = "log"), lambda = 5, 
#>     data = sim4$data)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.72607 -0.23926 -0.01133  0.23229  0.67453 
#> 
#>             Estimate Std. Error     Df z value Pr(>|z|)    
#> (Intercept) -1.02700    0.03436   1.00  -29.89   <2e-16 ***
#> X1           0.98975    0.03547   1.00   27.90   <2e-16 ***
#> X2           0.94438    0.03542   1.00   26.66   <2e-16 ***
#> X3           1.01056    0.03567   1.00   28.33   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for Gamma family taken to be 0.1180596)
#> 
#>     Null deviance: 266.048  on 99.0  degrees of freedom
#> Residual deviance:  10.264  on 81.6  degrees of freedom
#> AIC: -119.02
#> Lambda: 5
#> 
#> Number of IS iterations: 56
#> 
# }
```
