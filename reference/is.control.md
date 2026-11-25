# Control Parameters for islasso Fitting

Collects the parameters controlling the fitting process of \[islasso()\]
and \[islasso.path()\].

## Usage

``` r
is.control(
  sigma2 = -1,
  tol = 1e-05,
  itmax = 1000,
  stand = TRUE,
  trace = 0,
  nfolds = 5,
  seed = NULL,
  adaptive = FALSE,
  g = 0.5,
  b0 = NULL,
  V0 = NULL,
  c = 0.5,
  inner = c("mm", "newton"),
  solver = c("auto", "dense", "woodbury"),
  ncores = 1L,
  ...
)
```

## Arguments

- sigma2:

  fixed value of the dispersion parameter; if negative (default) it is
  estimated (Pearson statistic over the residual degrees of freedom).
  Ignored for the binomial, Poisson and negative binomial families,
  where it is 1.

- tol:

  tolerance of the convergence criteria (default \`1e-5\`).

- itmax:

  maximum number of iterations (default 1000).

- stand:

  logical: should the covariates be standardised (default \`TRUE\`)?
  With \`stand = TRUE\` the penalty acts on the standardised
  coefficients, as in glmnet; estimates are always returned on the
  original scale.

- trace:

  integer: 0 (default) silent, 1 compact and 2 detailed printing of the
  iterations.

- nfolds:

  number of folds for the cross-validation used when \`lambda\` is
  missing.

- seed:

  optional integer seed for the cross-validation folds (the random
  number generator state of the session is restored afterwards).

- adaptive:

  logical: adaptive lasso with weights based on the pseudo hat values
  (experimental).

- g:

  numeric in \\0, 1\\: the parameter of the extended BIC (\`g = 0\`
  gives the BIC).

- b0:

  optional starting values for the coefficients (original scale,
  intercept first).

- V0:

  optional starting covariance matrix (or a number for its diagonal) on
  the standardised scale; default \`0.01\`.

- c:

  the mixture weight of the induced smoothing (default 0.5), a number in
  \\0, 1\\ or a vector with one value per coefficient; a negative value
  means that the weights are re-estimated at each iteration.

- inner:

  the algorithm solving the smoothed estimating equation at fixed
  standard errors: \`"mm"\` (default, monotone majorise-minimise
  iterations, the undamped version of the fixed-point iterations of
  islasso 1.6) or \`"newton"\` (Newton's method with line search).

- solver:

  linear algebra: \`"auto"\` (default) uses \`"woodbury"\` when the
  number of columns exceeds twice the number of observations and
  \`"dense"\` otherwise.

- ncores:

  number of cores used by \[islasso.path()\] to fit the models for the
  different values of \`lambda\` in parallel (socket cluster of package
  'parallel'); default 1.

- ...:

  further arguments, ignored (allows \`do.call(is.control, control)\`).

## Value

A list with the control parameters.

## See also

\[islasso()\], \[islasso.path()\]

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
is.control(tol = 1e-6, seed = 1)
#> $sigma2
#> [1] -1
#> 
#> $tol
#> [1] 1e-06
#> 
#> $itmax
#> [1] 1000
#> 
#> $trace
#> [1] 0
#> 
#> $stand
#> [1] TRUE
#> 
#> $nfolds
#> [1] 5
#> 
#> $seed
#> [1] 1
#> 
#> $adaptive
#> [1] FALSE
#> 
#> $g
#> [1] 0.5
#> 
#> $b0
#> NULL
#> 
#> $V0
#> NULL
#> 
#> $c
#> [1] 0.5
#> 
#> $inner
#> [1] "mm"
#> 
#> $solver
#> [1] "auto"
#> 
#> $ncores
#> [1] 1
#> 
```
