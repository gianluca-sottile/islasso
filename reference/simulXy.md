# Simulate a Model Matrix and a Response Vector

Generates covariates (multivariate normal with AR(1) correlation) and a
response from a Gaussian, binomial, Poisson or Gamma regression model,
for simulation studies and examples.

## Usage

``` r
simulXy(
  n,
  p,
  interc = 0,
  beta,
  family = gaussian(),
  prop = 0.1,
  lim.b = c(-3, 3),
  sigma = 1,
  size = 1,
  rho = 0,
  scale.data = TRUE,
  seed = NULL,
  X = NULL,
  dispersion = 0.1
)
```

## Arguments

- n:

  number of observations.

- p:

  number of covariates.

- interc:

  the intercept of the linear predictor (default 0).

- beta:

  vector of `p` regression coefficients; if missing, a proportion `prop`
  of coefficients is drawn uniformly in `lim.b` and the others are zero.

- family:

  a family object, function or name: `gaussian`, `binomial`, `poisson`
  or `Gamma` (with any of their links).

- prop:

  proportion of non-null coefficients when `beta` is missing (default
  0.1).

- lim.b:

  range of the non-null coefficients when `beta` is missing.

- sigma:

  standard deviation of the Gaussian response.

- size:

  number of trials of the binomial response; with `size > 1` the
  response is a two-column matrix of failures and successes.

- rho:

  correlation parameter: \\cor(X_j, X_k) = \rho^{\|j-k\|}\\.

- scale.data:

  logical: should the columns of the model matrix be standardised?

- seed:

  optional seed; the state of the random number generator of the session
  is restored on exit.

- X:

  optional model matrix replacing the generated one.

- dispersion:

  dispersion parameter of the Gamma response (default 0.1).

## Value

A list with components

- data:

  a data frame with the response `y` (columns `y.failure`, `y.success`
  when `size > 1`) and the covariates `X1`, ..., `Xp`;

- beta0:

  the intercept;

- beta:

  the regression coefficients.

## See also

[`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
n <- 100; p <- 100
beta <- c(runif(10, -3, 3), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1234)
str(sim$data[, 1:5])
#> 'data.frame':    100 obs. of  5 variables:
#>  $ y : num  -3.411 7.985 -7.269 0.904 1.386 ...
#>  $ X1: num  -1.046 0.432 1.236 -2.179 0.583 ...
#>  $ X2: num  0.362 -0.5 0.024 -0.527 -0.84 ...
#>  $ X3: num  0.3443 0.5647 0.0322 0.5688 0.1636 ...
#>  $ X4: num  -0.5445 -0.8999 -0.1631 0.9692 0.0302 ...
```
