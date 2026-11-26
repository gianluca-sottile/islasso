# Select the Tuning Parameter by Information Criteria

Finds, along an
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
fit, the values of `lambda` minimising the information criteria AIC,
BIC, AICc, eBIC, GCV and GIC (see
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
for their definitions).

## Usage

``` r
GoF.islasso.path(object, plot = TRUE, ...)
```

## Arguments

- object:

  a fitted model of class `"islasso.path"`.

- plot:

  logical: should the criteria be plotted against \\\log\lambda\\?

- ...:

  not used.

## Value

A list with components

- gof:

  the matrix of criteria along the path;

- minimum:

  the positions of the minima;

- lambda.min:

  the values of `lambda` minimising each criterion;

- plot:

  the `ggplot` object (if `plot = TRUE`).

## Details

The degrees of freedom are the effective number of parameters (trace of
the pseudo hat matrix), which is not necessarily an integer.

## See also

[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md),
[`summary.islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.path.md),
[`predict.islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.path.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
n <- 100; p <- 30
beta <- c(runif(10, -2, 2), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso.path(y ~ ., data = sim$data, family = gaussian(), nlambda = 30)
GoF.islasso.path(fit, plot = FALSE)$lambda.min
#>       AIC       BIC      AICc      eBIC       GCV       GIC 
#>  3.217658 11.462096  8.343207 21.633467  6.072982 11.462096 
```
