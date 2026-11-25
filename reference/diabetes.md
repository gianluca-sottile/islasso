# Blood and other measurements in diabetics

The `diabetes` data frame contains 442 observations used in the Efron et
al. "Least Angle Regression" paper.

## Format

A data frame with 442 rows and 3 columns:

- x:

  Matrix with 10 numeric columns (standardized)

- y:

  Numeric response vector

- x2:

  Matrix with 64 columns including interactions

## Source

<https://web.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.ps>

## Details

The `x` matrix has been standardized to have unit L2 norm and zero mean
in each column. The `x2` matrix extends `x` by adding selected
interaction terms.

## References

Efron, Hastie, Johnstone and Tibshirani (2003). "Least Angle Regression"
(with discussion), *Annals of Statistics*.

## Examples

``` r
data(diabetes)
str(diabetes)
#> 'data.frame':    442 obs. of  3 variables:
#>  $ x : 'AsIs' num [1:442, 1:10] 0.03808 -0.00188 0.0853 -0.08906 0.00538 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:10] "age" "sex" "bmi" "map" ...
#>  $ y : num  151 75 141 206 135 97 138 63 110 310 ...
#>  $ x2: 'AsIs' num [1:442, 1:64] 0.03808 -0.00188 0.0853 -0.08906 0.00538 ...
#>   ..- attr(*, ".Names")= chr [1:28288] "age" "age" "age" "age" ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:442] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:64] "age" "sex" "bmi" "map" ...
summary(diabetes$y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    25.0    87.0   140.5   152.1   211.5   346.0 

if (FALSE) { # \dontrun{
  fit <- islasso(y ~ ., data = data.frame(y = diabetes$y, diabetes$x2),
                 family = gaussian())
  summary(fit, pval = 0.05)
  lambda.aic <- aic.islasso(fit, interval = c(1, 100))
  fit.aic <- update(fit, lambda = lambda.aic)
  summary(fit.aic, pval = 0.05)
} # }
```
