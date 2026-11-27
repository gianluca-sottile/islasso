# Prostate Cancer Data

This dataset originates from a study examining the correlation between
prostate-specific antigen levels and various clinical measures in men
scheduled for radical prostatectomy. It contains 97 rows and 9
variables.

## Format

A data frame with 97 observations and 9 variables:

- lcavol:

  Log of cancer volume

- lweight:

  Log of prostate weight

- age:

  Age of the patient

- lbph:

  Log of benign prostatic hyperplasia amount

- svi:

  Seminal vesicle invasion (binary)

- lcp:

  Log of capsular penetration

- gleason:

  Gleason score

- pgg45:

  Percentage of Gleason scores 4 or 5

- lpsa:

  Log of prostate-specific antigen

## Source

Stamey, T.A., et al. (1989). Prostate specific antigen in the diagnosis
and treatment of adenocarcinoma of the prostate: II. radical
prostatectomy treated patients. Journal of Urology, 141(5), 1076-1083.

## References

Stamey, T.A., Kabalin, J.N., McNeal, J.E., Johnstone, I.M., Freiha, F.,
Redwine, E.A., and Yang, N. (1989). Journal of Urology, 141(5),
1076-1083.

## Examples

``` r
data(Prostate)
summary(Prostate)
#>      lcavol           lweight           age             lbph        
#>  Min.   :-1.3471   Min.   :2.375   Min.   :41.00   Min.   :-1.3863  
#>  1st Qu.: 0.5128   1st Qu.:3.376   1st Qu.:60.00   1st Qu.:-1.3863  
#>  Median : 1.4469   Median :3.623   Median :65.00   Median : 0.3001  
#>  Mean   : 1.3500   Mean   :3.653   Mean   :63.87   Mean   : 0.1004  
#>  3rd Qu.: 2.1270   3rd Qu.:3.878   3rd Qu.:68.00   3rd Qu.: 1.5581  
#>  Max.   : 3.8210   Max.   :6.108   Max.   :79.00   Max.   : 2.3263  
#>       svi              lcp             gleason          pgg45       
#>  Min.   :0.0000   Min.   :-1.3863   Min.   :6.000   Min.   :  0.00  
#>  1st Qu.:0.0000   1st Qu.:-1.3863   1st Qu.:6.000   1st Qu.:  0.00  
#>  Median :0.0000   Median :-0.7985   Median :7.000   Median : 15.00  
#>  Mean   :0.2165   Mean   :-0.1794   Mean   :6.753   Mean   : 24.38  
#>  3rd Qu.:0.0000   3rd Qu.: 1.1787   3rd Qu.:7.000   3rd Qu.: 40.00  
#>  Max.   :1.0000   Max.   : 2.9042   Max.   :9.000   Max.   :100.00  
#>       lpsa        
#>  Min.   :-0.4308  
#>  1st Qu.: 1.7317  
#>  Median : 2.5915  
#>  Mean   : 2.4784  
#>  3rd Qu.: 3.0564  
#>  Max.   : 5.5829  
cor(Prostate$lpsa, Prostate$lcavol)
#> [1] 0.7344603
# \donttest{
fit <- islasso(lpsa ~ ., data = Prostate, family = gaussian(),
               control = is.control(seed = 1))
summary(fit, pval = 0.05)
#> 
#> Call:
#> islasso(formula = lpsa ~ ., family = gaussian(), data = Prostate, 
#>     control = is.control(seed = 1))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.65715 -0.39409  0.00099  0.39133  1.54191 
#> 
#>             Estimate Std. Error    Df z value Pr(>|z|)    
#> (Intercept)  0.47075    0.93597 1.000   0.503  0.61500    
#> lcavol       0.51990    0.07748 1.000   6.710 1.94e-11 ***
#> lweight      0.38162    0.16334 0.982   2.336  0.01947 *  
#> svi          0.58280    0.21558 0.992   2.703  0.00686 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.5085712)
#> 
#>     Null deviance: 127.918  on 96.00  degrees of freedom
#> Residual deviance:  45.746  on 89.95  degrees of freedom
#> AIC: 218.47
#> Lambda: 4.1676
#> 
#> Number of IS iterations: 70
#> 
lambda.aic <- aic.islasso(fit, method = "AIC", trace = FALSE)
fit.aic <- update(fit, lambda = lambda.aic)
summary(fit.aic, pval = 0.05)
#> 
#> Call:
#> islasso(formula = lpsa ~ ., family = gaussian(), lambda = lambda.aic, 
#>     data = Prostate, control = is.control(seed = 1))
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.72467 -0.38645  0.00495  0.39934  1.59690 
#> 
#>             Estimate Std. Error    Df z value Pr(>|z|)    
#> (Intercept)  0.62438    1.21792 1.000   0.513  0.60819    
#> lcavol       0.56756    0.08617 1.000   6.587  4.5e-11 ***
#> lweight      0.44002    0.16866 0.998   2.609  0.00908 ** 
#> svi          0.71444    0.23998 0.999   2.977  0.00291 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.5019419)
#> 
#>     Null deviance: 127.918  on 96.00  degrees of freedom
#> Residual deviance:  44.277  on 88.21  degrees of freedom
#> AIC: 218.78
#> Lambda: 0.80954
#> 
#> Number of IS iterations: 37
#> 
# }
```
