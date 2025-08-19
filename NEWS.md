# islasso 1.6.1 (2025-08-17)

- Added a GitHub repository and project website.
- Introduced an introductory vignette.
- Optimized and modernized Fortran routines, improving computational speed (2x).
- Updated documentation and manuals.
- Fixed minor bugs.

# islasso 1.6.0 (2025-07-30)

## Performance & Refactoring
- Core computational routines have been cleaned up, and some bugs have been fixed.
- Legacy R routines have been revised, cleaned, and commented. Minor inconsistencies have been addressed.

## Documentation
- Help files and function manuals are now fully managed via `roxygen2`, with substantial updates to usage examples and descriptions.

## Visualization
- All plotting functions have been refactored to use the `ggplot2` framework for consistent and modern graphics.

## User Experience
- A custom ASCII startup banner has been added on package attach, providing a welcoming and informative message.

---

# islasso 1.5.1

- Some bugs fixed.

# islasso 1.5.0

- Some bugs fixed.
- Other S3 methods implemented.

# islasso 1.4.3

- Some bugs fixed.

# islasso 1.4.2

- Some bugs for binomial family fixed.

# islasso 1.4.1

- Some bugs fixed.

# islasso 1.4.0

- New optimization algorithm for the 'islasso' method. The algorithm is now stable for all the implemented distributions.
- In `aic.islasso()` function the available methods are "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC".
- New class of functions named `islasso.path` created. The main function `islasso.path()` builds the coefficient profile for a fixed sequence of lambda values.
- New function `GoF.islasso.path()` extracts the optimal tuning parameter minimizing a fixed criterion. Available criteria are the same as in `aic.islasso()`.
- Some bugs fixed.

# islasso 1.3.1

- Some bugs fixed.

# islasso 1.3.0

- Vignette added to the package.
- Some bugs fixed.

# islasso 1.2.3

- Some bugs fixed.

# islasso 1.2.2

- Some bugs fixed.

# islasso 1.2.1

- Some bugs fixed.

# islasso 1.2.0

- New implementation of the estimating algorithm. Now islasso is much stabler and faster.
- New function: general linear hypotheses for linear combinations of the regression coefficients, including confidence intervals.
- Prediction function includes confidence intervals for the fitted values.
- Step halving with Armijo's rule improved.
- Convergence criterion improved.
- Some bugs fixed.

# islasso 1.1.0

- New implementation of the estimating algorithm. Now islasso is much stabler and faster, reducing the number of iterations to reach convergence.
- Step halving with Armijo's rule implemented.
- Elastic-net approach added via `alpha` parameter in the objective function (as in `glmnet`).
- Summary method now includes degrees of freedom for each covariate, with choice between t-test or z-test (only for Gaussian family).
- `optim.islasso` renamed to `aic.islasso`; interval specification no longer required.
- `islasso.control` renamed to `is.control`; control parameters modified.
- Two trace versions implemented in `is.control`: compact (`trace = 1`) and verbose (`trace = 2`).
- Some bugs fixed.
