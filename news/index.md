# Changelog

## islasso 1.7.0

### Numerical core

- The Fortran routines have been replaced by a new core written in C on
  top of the BLAS/LAPACK libraries of R (no Fortran compiler is needed
  any more to install from source).
- At fixed standard errors the smoothed estimating equation is the
  gradient of a smooth convex objective; it is now solved by undamped
  majorise-minimise iterations (default, the same fixed point as before
  with about half of the factorisations) or by Newton’s method with line
  search (`is.control(inner = "newton")`).
- Only the diagonal of the sandwich matrix and the pseudo hat values are
  computed during the iterations (about 4 times fewer flops per
  iteration); the full covariance matrix is computed once, at
  convergence.
- New Woodbury back end used when the number of columns exceeds twice
  the number of observations: cost O(n^2 p) instead of O(p^3) per step
  (`is.control(solver = )`).
- [`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
  fits every model as
  [`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
  would at the same `lambda` (glmnet starting values, standard errors
  started at their default values), so that path and single fits agree
  (in 1.6.x the standard errors were carried over along the path and the
  effective degrees of freedom, hence the information criteria, depended
  on the path); the models can be fitted in parallel with
  `is.control(ncores = )` (socket cluster); computes X’X once for the
  Gaussian family, keeps all the models (with a convergence flag, column
  `conv` of `Info`) and no longer drops covariates permanently from the
  path.
- For the Gaussian family the whole induced-smoothing iteration runs in
  compiled code.
- No more large arrays on the stack (possible crashes with large n \* p
  in 1.6.x).

### Models

- Any family object accepted by
  [`glm()`](https://rdrr.io/r/stats/glm.html) can be used (all the links
  of `binomial`, `poisson`, `Gamma`; `inverse.gaussian`; quasi families;
  [`MASS::negative.binomial()`](https://rdrr.io/pkg/MASS/man/negative.binomial.html),
  …): the IRLS loop is driven by the family object, with step halving
  when the linear predictor becomes invalid.
- With `stand = TRUE` the penalty now acts on the standardised
  coefficients, as in glmnet (whose cross-validation provides `lambda`
  when it is missing): before, the estimator penalised the coefficients
  on the original scale while glmnet worked on the standardised ones.
  Results change for covariates that are not standardised.
- The scale of `lambda` is `sum(weights)` times the glmnet scale (equal
  to `n` with unit weights).
- `is.control(seed = )` is now used (reproducible cross-validation; the
  RNG state is restored).

### Bug fixes

- Models without intercept and `stand = TRUE` were centred, changing the
  model.
- The probit inverse link had a wrong clamping threshold (8.1e6 instead
  of 8.13), producing infinite weights for large linear predictors.
- [`aic.islasso()`](https://gianluca-sottile.github.io/islasso/reference/aic.islasso.md)
  refitted the model with the IRLS working weights as prior weights; it
  now optimises over log(lambda).
- AICc had an extra factor n; GCV is now n \* deviance / (n - df)^2.
- [`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md):
  rows of fitted values/residuals could be misaligned with the
  coefficients when some fit failed.
- [`relax.islasso()`](https://gianluca-sottile.github.io/islasso/reference/relax.islasso.md)
  failed for models without intercept and used an infinite `lambda`.
- `predict(type = "coefficients", se.fit = TRUE)` failed;
  [`summary.islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.path.md)
  printed the log-likelihood as residual deviance; path methods returned
  `NA` at the largest `lambda`.
- An `offset` given as a covariate name removed that covariate from the
  model.
- [`logLik()`](https://rdrr.io/r/stats/logLik.html) is `NA` for families
  without likelihood (quasi).
- [`anova()`](https://rdrr.io/r/stats/anova.html): the overall Wald test
  uses a generalised inverse, so it no longer fails when the covariance
  matrix of the hypotheses is numerically singular (degrees of freedom =
  its rank).
- The dispersion is protected against non-positive residual degrees of
  freedom for all families.

### Other changes

- Influence measures and diagnostic plots use the leverages of the
  sandwich Hessian (no QR of the augmented design, no n x n hat matrix).
- glmnet moved from Depends to Imports; the startup banner and the
  dependencies on cli and gridExtra have been removed.
- New vignettes: methodology and algorithm, GLMs and families,
  high-dimensional data.
- Unit tests (testthat).

## islasso 1.6.2 (2025-11-15)

CRAN release: 2025-11-15

- Fixed the number of arguments in a call to the internal routine
  `inv_posdef`.

## islasso 1.6.1 (2025-08-17)

CRAN release: 2025-11-14

- Created GitHub page and website.
- Created introductory vignette.
- Fortran routines updated and speeded up (~2x).
- Manuals updated.
- Added new function
  [`relax.islasso()`](https://gianluca-sottile.github.io/islasso/reference/relax.islasso.md),
  which allows fitting a relaxed islasso model by selecting variables to
  remain unpenalized. Variables can be specified either by name or by
  index, or automatically selected according to a significance level
  (`alpha`). This extension provides additional flexibility in
  post-selection inference.
- Some bugs fixed.

## islasso 1.6.0 (2025-07-30)

CRAN release: 2025-07-31

### Performance & Refactoring

- Core computational routines have been cleaned up, and some bugs have
  been fixed.
- Legacy R routines have been revised, cleaned, and commented. Minor
  inconsistencies have been addressed.

### Documentation

- Help files and function manuals are now fully managed via `roxygen2`,
  with substantial updates to usage examples and descriptions.

### Visualization

- All plotting functions have been refactored to use the `ggplot2`
  framework for consistent and modern graphics.

### User Experience

- A custom ASCII startup banner has been added on package attach,
  providing a welcoming and informative message.

------------------------------------------------------------------------

## islasso 1.5.1

CRAN release: 2023-05-15

- Some bugs fixed.

## islasso 1.5.0

CRAN release: 2023-05-08

- Some bugs fixed.
- Other S3 methods implemented.

## islasso 1.4.3

CRAN release: 2022-03-01

- Some bugs fixed.

## islasso 1.4.2

CRAN release: 2021-12-06

- Some bugs for binomial family fixed.

## islasso 1.4.1

CRAN release: 2021-11-16

- Some bugs fixed.

## islasso 1.4.0

CRAN release: 2021-10-15

- New optimization algorithm for the ‘islasso’ method. The algorithm is
  now stable for all the implemented distributions.
- In
  [`aic.islasso()`](https://gianluca-sottile.github.io/islasso/reference/aic.islasso.md)
  function the available methods are “AIC”, “BIC”, “AICc”, “eBIC”,
  “GCV”, “GIC”.
- New class of functions named `islasso.path` created. The main function
  [`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
  builds the coefficient profile for a fixed sequence of lambda values.
- New function
  [`GoF.islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/GoF.islasso.path.md)
  extracts the optimal tuning parameter minimizing a fixed criterion.
  Available criteria are the same as in
  [`aic.islasso()`](https://gianluca-sottile.github.io/islasso/reference/aic.islasso.md).
- Some bugs fixed.

## islasso 1.3.1

CRAN release: 2021-06-16

- Some bugs fixed.

## islasso 1.3.0

CRAN release: 2021-05-24

- Vignette added to the package.
- Some bugs fixed.

## islasso 1.2.3

- Some bugs fixed.

## islasso 1.2.2

CRAN release: 2021-01-22

- Some bugs fixed.

## islasso 1.2.1

CRAN release: 2021-01-21

- Some bugs fixed.

## islasso 1.2.0

CRAN release: 2021-01-21

- New implementation of the estimating algorithm. Now islasso is much
  stabler and faster.
- New function: general linear hypotheses for linear combinations of the
  regression coefficients, including confidence intervals.
- Prediction function includes confidence intervals for the fitted
  values.
- Step halving with Armijo’s rule improved.
- Convergence criterion improved.
- Some bugs fixed.

## islasso 1.1.0

CRAN release: 2019-06-25

- New implementation of the estimating algorithm. Now islasso is much
  stabler and faster, reducing the number of iterations to reach
  convergence.
- Step halving with Armijo’s rule implemented.
- Elastic-net approach added via `alpha` parameter in the objective
  function (as in `glmnet`).
- Summary method now includes degrees of freedom for each covariate,
  with choice between t-test or z-test (only for Gaussian family).
- `optim.islasso` renamed to `aic.islasso`; interval specification no
  longer required.
- `islasso.control` renamed to `is.control`; control parameters
  modified.
- Two trace versions implemented in `is.control`: compact (`trace = 1`)
  and verbose (`trace = 2`).
- Some bugs fixed.
