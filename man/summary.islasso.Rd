% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3methods_islasso.R
\name{summary.islasso}
\alias{summary.islasso}
\alias{print.summary.islasso}
\title{Summarize islasso Fitted Model}
\usage{
\method{summary}{islasso}(object, pval = 1, which, use.t = FALSE, type.pval = "wald", ...)
}
\arguments{
\item{object}{A fitted model of class \code{"islasso"}.}

\item{pval}{Numeric threshold for displaying coefficients. Only those with \eqn{p \le} \code{pval} are printed. Unpenalized coefficients (like intercepts) are always shown.}

\item{which}{Optional. Specifies a subset of coefficients to test. If missing, all parameters are evaluated.}

\item{use.t}{Logical. If \code{TRUE}, p-values are computed using the t-distribution and residual degrees of freedom.}

\item{type.pval}{Character. Type of p-value approximation. Only \code{"wald"} (default) is implemented.}

\item{...}{Additional arguments (not currently used).}
}
\value{
An object of class \code{"summary.islasso"} containing:
  \item{coefficients}{Coefficient estimates and related statistics}
  \item{pval}{Threshold used to filter coefficients}
  \item{call}{Original model call}
}
\description{
Provides a concise summary of a fitted \code{\link{islasso}} model, including p-values and optional filtering.
}
\examples{
\dontrun{
# Assuming object `o` from an islasso fit
summary(o, pval = 0.1)  # Show coefficients with p <= 0.1
}

}
\seealso{
\code{\link{islasso.fit}}, \code{\link{residuals.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}}, \code{\link{deviance.islasso}}
}
\author{
Gianluca Sottile \email{gianluca.sottile@unipa.it}
}
