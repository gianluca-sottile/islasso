# **Induced Smoothed Lasso**

[![CRAN
version](https://www.r-pkg.org/badges/version/islasso)](https://cran.r-project.org/package=islasso)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/islasso)](https://cran.r-project.org/package=islasso)

## 📦 Installation

You can install the development version of `islasso` from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("gianluca-sottile/islasso")
```

Once installed, load the package:

``` r

library(islasso)
```

## 🔍 Description

`islasso` implements the **Induced Smoothed Lasso** (Cilluffo et al.,
2020): the non-smooth lasso (or elastic-net) penalty is replaced by its
expectation under the approximate sampling distribution of the
estimator, so that standard errors (sandwich formula), Wald tests,
confidence intervals and general linear hypotheses are available for
**all** the coefficients.

Key features:

- numerical core in C (BLAS/LAPACK): majorise–minimise or Newton
  iterations on a smooth convex objective, only the diagonal of the
  sandwich during the iterations, Woodbury identities when `p > 2n`
  (cost `O(n^2 p)` instead of `O(p^3)`);
- any GLM family and link accepted by
  [`glm()`](https://rdrr.io/r/stats/glm.html) (binomial with
  logit/probit/cloglog/cauchit/log, Poisson, Gamma, inverse Gaussian,
  quasi families, negative binomial, …);
- regularisation paths with warm starts and tuning-parameter selection
  by cross-validation or information criteria (AIC, BIC, AICc, eBIC,
  GCV, GIC);
- relaxed fits, general linear hypotheses, diagnostic plots (`ggplot2`).

## 🚀 Quick Example

``` r

sim <- simulXy(n = 100, p = 50, beta = c(rep(1, 5), rep(0, 45)), seed = 1)
mod <- islasso(y ~ ., data = sim$data, control = is.control(seed = 1))  # lambda by 5-fold CV
summary(mod, pval = 0.05)
confint(mod, parm = 2:6)

path <- islasso.path(y ~ ., data = sim$data, nlambda = 50)
GoF.islasso.path(path)$lambda.min
```

## 📚 Documentation

- 📘 Function reference:
  [`?islasso`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
- 📄 Vignettes:
  [`vignette("islasso-intro")`](https://gianluca-sottile.github.io/islasso/articles/islasso-intro.md),
  [`vignette("islasso-methodology")`](https://gianluca-sottile.github.io/islasso/articles/islasso-methodology.md),
  [`vignette("islasso-glm")`](https://gianluca-sottile.github.io/islasso/articles/islasso-glm.md),
  [`vignette("islasso-highdim")`](https://gianluca-sottile.github.io/islasso/articles/islasso-highdim.md)
- 🌐 Website: <https://gianluca-sottile.github.io/islasso/>

## 📖 References

> Cilluffo G, Sottile G, La Grutta S, Muggeo V (2020). *The Induced
> Smoothed lasso: A practical framework for hypothesis testing in high
> dimensional regression.* Statistical Methods in Medical Research\_,
> *29*(3), 765-777.
> [doi:10.1177/0962280219842890](https://doi.org/10.1177/0962280219842890)

## 🤝 Contributing

Feel free to open issues, suggest improvements, or submit pull
requests.  
Bug reports and feature requests are welcome!

## 📜 License

GPL (\>= 2)
