\name{islasso-internal}
\title{Internal Functions}

\alias{qqNorm}
\alias{modelX}
\alias{lminfl}
\alias{is.influence}
\alias{islasso.diag}
\alias{islasso.diag.plots}
\alias{predislasso}

\alias{.checkinput}
\alias{.startpoint}
\alias{.islasso}

\alias{checkinput.islasso.path}
\alias{startpoint.islasso.path}
\alias{islasso.path.fit.glm}

\alias{interpolate}
\alias{create_coef_plot}
\alias{create_se_plot}
\alias{create_weight_plot}
\alias{calculate_gradient}
\alias{create_gradient_plot}
\alias{create_gof_plot}

\alias{makeHyp}
\alias{printHyp}
\alias{cislasso}
\alias{ci.fitted.islasso}

\description{Internal islasso functions.}
\usage{
qqNorm(x, probs = seq(0.005, 0.995, length.out = 200), centre = FALSE,
  scale = FALSE, leg = TRUE, mean = 0, sd = 1, dF = FALSE, ylab = NULL,
  color = "black", ...)
modelX(n, p, rho=.5, scale.data=TRUE)
lminfl(mod, tol = 1e-8)
is.influence(model, do.coef = TRUE)
islasso.diag(glmfit)
islasso.diag.plots(glmfit, glmdiag = islasso.diag(glmfit),
  subset = NULL, iden = FALSE, labels = NULL, ret = FALSE)
predislasso(object, newdata, type = c("response", "terms"),
  terms = NULL, na.action = na.pass, ...)

.checkinput(X, y, family, alpha, intercept, weights, offset,
  unpenalized, control)
.startpoint(X, y, lambda, alpha, weights, offset, mustart,
  family, intercept, setting)
.islasso(prep, start, Lambda, fam, link)

checkinput.islasso.path(X, y, family, lambda, nlambda, lambda.min.ratio,
  alpha, intercept, weights, offset, unpenalized, control)
startpoint.islasso.path(X, y, lambda, alpha, weights, offset, mustart,
  family, intercept, setting)
islasso.path.fit.glm(prep, start, lambda, fam, link)

interpolate(y1, y2, x1, x2, x.new)
create_coef_plot(coef1, loglambda, label, id.best, gof, dots, active,
  unactive, legend, nlambda)
create_se_plot(se1, coef1, loglambda, label, id.best, gof, dots,
  active, unactive, legend, nlambda)
create_weight_plot(weight1, coef1, loglambda, label, id.best, gof, dots,
  active, unactive, legend, nlambda)
calculate_gradient(object, lambda, nlambda, intercept)
create_gradient_plot(grad, coef1, lambda, label, id.best, gof, dots,
  active, unactive, legend, nlambda)
create_gof_plot(object, loglambda, id.best, gof, dots)

makeHyp(cnames, hypothesis, rhs = NULL)
printHyp(L, b, nms)
cislasso(object, a, ci)
ci.fitted.islasso(object, newx, ci = NULL, type.ci = "wald",
  conf.level=.95, only.ci = FALSE)
}
\author{Gianluca Sottile (\email{gianluca.sottile@unipa.it})}
\details{These functions are not intended for users.}
\keyword{internal}
