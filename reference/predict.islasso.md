# Prediction Method for islasso Objects

Computes predictions from a fitted
[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
model object. Multiple output types supported, including response scale,
linear predictor, and coefficient values.

## Usage

``` r
# S3 method for class 'islasso'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "coefficients", "class", "terms"),
  se.fit = FALSE,
  ci = NULL,
  type.ci = c("wald", "score"),
  level = 0.95,
  terms = NULL,
  na.action = na.pass,
  ...
)
```

## Arguments

- object:

  A fitted model of class `"islasso"`.

- newdata:

  Optional data frame containing predictors for prediction. If omitted,
  the fitted model matrix is used.

- type:

  Character. Specifies the prediction scale:

  - `"link"` (default): linear predictor scale;

  - `"response"`: original response scale;

  - `"coefficients"`: estimated coefficients;

  - `"class"`: predicted class (only for
    [`binomial()`](https://rdrr.io/r/stats/family.html) family);

  - `"terms"`: contribution of each term to the linear predictor.

- se.fit:

  Logical. Whether to compute standard errors/confidence intervals.

- ci:

  Optional. Precomputed matrix of confidence intervals (2 columns).

- type.ci:

  Type of interval. Only `"wald"` is implemented.

- level:

  Confidence level for intervals. Default is `0.95`.

- terms:

  If `type = "terms"`, optionally specify which terms to extract.

- na.action:

  Function to handle missing values in `newdata`. Default: `na.pass`.

- ...:

  Additional arguments passed to downstream methods.

## Value

A numeric vector, matrix, or list depending on `type`.

## See also

[`islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`summary.islasso`](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[`logLik.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`residuals.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md),
[`deviance.islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)

## Author

Gianluca Sottile <gianluca.sottile@unipa.it>

## Examples

``` r
set.seed(1)
n <- 100; p <- 100
beta <- c(runif(20, -3, 3), rep(0, p - 20))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso(y ~ ., data = sim$data, family = gaussian(), lambda = 2)
predict(fit, type = "response")
#>            1            2            3            4            5            6 
#>   5.21020664  -1.43421213  -1.87791088   3.38722076  -1.92121489   2.98399863 
#>            7            8            9           10           11           12 
#>  -9.43967416   4.15444608  -2.66866256  -4.14657527   3.46597081  -4.10865541 
#>           13           14           15           16           17           18 
#>   0.90571133  -2.76268169   1.89880109  -0.60790512  -6.55905785   0.06405997 
#>           19           20           21           22           23           24 
#>  -2.62465339   0.42162178  -4.82394986  -4.67200277  10.48924263  14.86563057 
#>           25           26           27           28           29           30 
#>  -4.48873490  -6.37536959   4.09580825  13.60024576   2.07323884   7.50932660 
#>           31           32           33           34           35           36 
#>  -6.73279096  -0.91592165   9.11033340  11.18301303  -4.52898177  -7.37401519 
#>           37           38           39           40           41           42 
#>  23.33784169   4.78118074  -1.57125946   5.41078578  -0.46321484  -8.41448565 
#>           43           44           45           46           47           48 
#>   6.61190530  -2.31717713  -7.40475331 -10.68565253  -0.01408795 -12.36278647 
#>           49           50           51           52           53           54 
#>  -0.89289804   9.39046946   3.35334858  -0.98578129  -5.27453464  -5.73864071 
#>           55           56           57           58           59           60 
#>   8.16540484  -5.46079767 -12.85745842  10.86630292 -11.87254466  11.34249952 
#>           61           62           63           64           65           66 
#>  -3.67390869   4.09016722  -8.00817873  -7.92823089   6.25605727  13.88266447 
#>           67           68           69           70           71           72 
#>  -6.17327077  14.24611341  -6.76382107   3.31284085 -15.90362475   8.00361199 
#>           73           74           75           76           77           78 
#>  -2.24689556  -4.95818129   5.03243898  -0.25807910   6.89749841   4.03794340 
#>           79           80           81           82           83           84 
#> -12.83216625  -0.93607050  -9.93800658  -0.36122027  -0.11142223  10.86635980 
#>           85           86           87           88           89           90 
#> -10.84552784  -3.02222837 -12.36530414  13.77558320  -2.55195678  -4.87463935 
#>           91           92           93           94           95           96 
#>   2.22065465  -6.34621525  -2.11717296  -3.08787669  -0.06717946   9.07336127 
#>           97           98           99          100 
#>  -9.20158960   6.35749136  13.40272158  -0.31340209 
```
