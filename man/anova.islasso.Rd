% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint_and_auxiliary.R
\name{anova.islasso}
\alias{anova.islasso}
\alias{print.anova.islasso}
\title{General Linear Hypotheses for islasso Models}
\usage{
\method{anova}{islasso}(object, A, b = NULL, ci, ...)
}
\arguments{
\item{object}{A fitted model object of class \code{"islasso"}.}

\item{A}{Hypothesis specification. Either:
\itemize{
  \item A numeric matrix or vector with each row specifying a linear combination of coefficients,
  \item Or a character vector with symbolic expressions (e.g. \code{"X1 + X2 = 3"}).
}}

\item{b}{Right-hand side vector for the null hypotheses \code{A \%*\% beta = b}. If omitted, defaults to zeros.}

\item{ci}{Optional 2-column matrix of confidence intervals for coefficients.}

\item{...}{Currently unused.}
}
\value{
An object of class \code{"anova.islasso"} containing:
\item{Estimate}{Linear combination estimates}
\item{SE}{Standard errors}
\item{Wald}{Wald statistics}
\item{p-value}{Associated p-values}
}
\description{
Tests general linear hypotheses and computes confidence intervals for linear combinations of coefficients
from a fitted \code{\link{islasso}} model.
}
\details{
The method tests the null hypothesis \eqn{H_0: A \beta = b}, where \eqn{A} and \eqn{b} define a linear constraint on model coefficients.

Symbolic expressions support natural syntax: coefficients may be added/subtracted, constants may be multiplied (e.g. \code{"2 * X1 + 3 * X2 = 7"}).
Equations with omitted \code{=} assume zero on the right-hand side. See examples for syntax flexibility.
}
\examples{
n <- 100; p <- 100
beta <- c(runif(10, -2, 2), rep(0, p - 10))
sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
fit <- islasso(y ~ . -1, data = sim$data, family = gaussian())

# Test if first 5 variables sum to -7.5
anova(fit, A = c("X1 + X2 + X3 + X4 + X5 = -7.5"))

# Test multiple hypotheses
anova(fit, A = c("X1 + X2 + X3 + X4 + X5", "X6 + X7 + X8 + X9 + X10"), b = c(-7.5, 8.75))

# Full diagonal comparison to true coefficients
anova(fit, A = diag(p), b = beta)

}
\seealso{
\code{\link{islasso}}, \code{\link{summary.islasso}}, \code{\link{confint.islasso}},
         \code{\link{predict.islasso}}, \code{\link{logLik.islasso}}, \code{\link{residuals.islasso}}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
