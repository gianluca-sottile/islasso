#' @useDynLib islasso, .registration = TRUE, .fixes = "C_"
#' @import glmnet stats ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom cli col_cyan
#' @importFrom graphics abline legend lines par plot points segments axis text matplot identify
#' @importFrom utils setTxtProgressBar txtProgressBar tail head globalVariables
NULL

#' The Induced Smoothed Lasso: A practical framework for hypothesis testing in high dimensional regression
#'
#' This package implements an induced smoothed approach for hypothesis testing in Lasso regression.
#'
#' @details
#' \tabular{ll}{
#' Package: \tab islasso\cr
#' Type: \tab Package\cr
#' Version: \tab 1.6.0\cr
#' Date: \tab 2025-07-30\cr
#' License: \tab GPL-2\cr
#' }
#'
#' \code{\link{islasso}} fits generalized linear models with an L1 penalty on selected coefficients.
#' It returns both point estimates and full covariance matrices, enabling standard error-based inference.
#' Related methods include: \code{\link{summary.islasso}}, \code{\link{predict.islasso}}, \code{\link{logLik.islasso}}, \code{\link{deviance.islasso}}, and \code{\link{residuals.islasso}}.
#'
#' \code{\link{islasso.path}} fits regularization paths using the Induced Smoothed Lasso.
#' It computes coefficients and standard errors across a grid of \code{lambda} values.
#' Companion methods include: \code{\link{summary.islasso.path}}, \code{\link{predict.islasso.path}}, \code{\link{logLik.islasso.path}}, \code{\link{residuals.islasso.path}}, \code{\link{coef.islasso.path}}, and \code{\link{fitted.islasso.path}}.
#'
#' @author
#' Gianluca Sottile, based on preliminary work by Vito Muggeo.
#' Maintainer: \email{gianluca.sottile@unipa.it}
#'
#' @references
#' Cilluffo, G., Sottile, G., La Grutta, S., Muggeo, VMR (2019). *The Induced Smoothed lasso: A practical framework for hypothesis testing in high dimensional regression*, Statistical Methods in Medical Research.
#' DOI: \doi{10.1177/0962280219842890}
#'
#' Sottile, G., Cilluffo, G., Muggeo, VMR (2019). *The R package islasso: estimation and hypothesis testing in lasso regression*. Technical Report on ResearchGate.
#' DOI: \doi{10.13140/RG.2.2.16360.11521}
#'
#' @docType package
#' @name islasso
#' @keywords package
"_PACKAGE"
