## **Induced Smoothed Lasso for R** <img src="man/figures/logo.png" align="right" height="120" alt="" />

## 📦 Installation

You can install the development version of `islasso` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("gianluca-sottile/islasso")
```

Once installed, load the package:

```r
library(islasso)
```

---

## 🔍 Description

`islasso` implements the **Induced Smoothed Lasso**, a robust and interpretable approach for hypothesis testing in high-dimensional linear and generalized linear models.

Key features include:

- Efficient Fortran backend for fast computation
- Support for Gaussian, Binomial, Poisson, and Gamma families
- Smoothed penalization for stable inference
- Automatic selection of active variables
- Visualization tools powered by `ggplot2`

---

## 🚀 Quick Example

```r
set.seed(123)
sim <- simulXy(n = 100, p = 20, family = "gaussian")
mod <- islasso(y ~ ., data = sim$data)
summary(mod)
plot(mod)
```

---

## 📚 Documentation

- 📘 Function reference: `?islasso`
- 📄 Vignette: `vignette("islasso-intro")`
- 🌐 Website: [https://gianluca-sottile.github.io/islasso](https://gianluca-sottile.github.io/islasso)

---

## 📖 References

> Cilluffo G, Sottile G, La Grutta S, Muggeo V (2020). *The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression.* 
> Statistical Methods in Medical Research_, *29*(3), 765-777. [doi:10.1177/0962280219842890](https://doi.org/10.1177/0962280219842890)

---

## 🤝 Contributing

Feel free to open issues, suggest improvements, or submit pull requests.  
Bug reports and feature requests are welcome!

---

## 📜 License

<a href="https://cran.r-project.org/web/packages/islasso/index.html">islasso</a> © 2019 by <a href="https://gianlucasottile.rbind.io">Gianluca Sottile</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">
