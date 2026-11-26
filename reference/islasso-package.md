# islasso: The Induced Smoothed Lasso

Estimation and inference for lasso and elastic-net regularised
regression models through the induced smoothing approach: the non-smooth
penalty is replaced by its expectation under the approximate sampling
distribution of the estimator, so that standard errors (sandwich
formula) and Wald tests are available for all the coefficients.

## Details

Main functions:
[`islasso()`](https://gianluca-sottile.github.io/islasso/reference/islasso.md)
fits a model at a given value of the tuning parameter (chosen by
cross-validation if missing),
[`islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/islasso.path.md)
fits a regularisation path, with
[`GoF.islasso.path()`](https://gianluca-sottile.github.io/islasso/reference/GoF.islasso.path.md)
and
[`aic.islasso()`](https://gianluca-sottile.github.io/islasso/reference/aic.islasso.md)
selecting the tuning parameter by information criteria;
[`relax.islasso()`](https://gianluca-sottile.github.io/islasso/reference/relax.islasso.md)
refits a relaxed model. Methods for
[summary()](https://gianluca-sottile.github.io/islasso/reference/summary.islasso.md),
[confint()](https://gianluca-sottile.github.io/islasso/reference/confint.islasso.md),
[anova()](https://gianluca-sottile.github.io/islasso/reference/anova.islasso.md)
(general linear hypotheses),
[predict()](https://gianluca-sottile.github.io/islasso/reference/predict.islasso.md)
and
[plot()](https://gianluca-sottile.github.io/islasso/reference/plot.islasso.md)
are available.

## References

Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2020). The induced
smoothed lasso: a practical framework for hypothesis testing in high
dimensional regression. *Statistical Methods in Medical Research*,
**29**(3), 765–777.
[doi:10.1177/0962280219842890](https://doi.org/10.1177/0962280219842890)

Sottile G., Cilluffo G., Muggeo V.M.R. (2019). The R package islasso:
estimation and hypothesis testing in lasso regression. Technical report.
[doi:10.13140/RG.2.2.16360.11521](https://doi.org/10.13140/RG.2.2.16360.11521)

## See also

Useful links:

- <https://gianluca-sottile.github.io/islasso/>

- <https://github.com/gianluca-sottile/islasso>

- Report bugs at <https://github.com/gianluca-sottile/islasso/issues>

## Author

**Maintainer**: Gianluca Sottile <gianluca.sottile@unipa.it>
([ORCID](https://orcid.org/0000-0001-9347-7251))

Authors:

- Gianluca Sottile <gianluca.sottile@unipa.it>
  ([ORCID](https://orcid.org/0000-0001-9347-7251))

- Giovanna Cilluffo

- Vito M. R. Muggeo
