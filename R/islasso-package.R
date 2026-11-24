#' islasso: The Induced Smoothed Lasso
#'
#' Estimation and inference for lasso and elastic-net regularised regression models through
#' the induced smoothing approach: the non-smooth penalty is replaced by its expectation under
#' the approximate sampling distribution of the estimator, so that standard errors (sandwich
#' formula) and Wald tests are available for all the coefficients.
#'
#' Main functions: [islasso()] fits a model at a given value of the tuning parameter (chosen by
#' cross-validation if missing), [islasso.path()] fits a regularisation path, with
#' [GoF.islasso.path()] and [aic.islasso()] selecting the tuning parameter by information
#' criteria; [relax.islasso()] refits a relaxed model. Methods for [summary()][summary.islasso],
#' [confint()][confint.islasso], [anova()][anova.islasso] (general linear hypotheses),
#' [predict()][predict.islasso] and [plot()][plot.islasso] are available.
#'
#' @references
#' Cilluffo G., Sottile G., La Grutta S., Muggeo V.M.R. (2020). The induced smoothed lasso: a
#' practical framework for hypothesis testing in high dimensional regression.
#' *Statistical Methods in Medical Research*, **29**(3), 765--777. \doi{10.1177/0962280219842890}
#'
#' Sottile G., Cilluffo G., Muggeo V.M.R. (2019). The R package islasso: estimation and
#' hypothesis testing in lasso regression. Technical report. \doi{10.13140/RG.2.2.16360.11521}
#'
#' @useDynLib islasso, .registration = TRUE, .fixes = "C_"
#' @import stats ggplot2
#' @importFrom graphics plot
#' @importFrom utils head
#' @keywords internal
"_PACKAGE"
