### S3 methods for islasso.path class

interpolate <- function(y1, y2, x1, x2, x.new) {
  if (x1 == x2) return(y1)
  m <- (y2 - y1) / (log(x2) - log(x1))
  y1 + m * (log(x.new) - log(x1))
}

#' @export
print.islasso.path <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (length(x$Coef)) {
    cat("Coefficients:\n")
    Info <- data.frame(x$Info[, 1:5, drop = FALSE])

    # Format all columns in one go with lapply
    formats <- list(
      lambda = list(format = "f", digits = 4),
      df = list(format = "f", digits = 4),
      phi = list(format = "f", digits = 4),
      deviance = list(format = "f", digits = 4),
      logLik = list(format = "f", digits = 4)
    )

    Info <- as.data.frame(lapply(names(formats), function(col) {
      fmt <- formats[[col]]
      width <- 4 + nchar(round(max(Info[[col]])))
      formatC(Info[[col]], format = fmt$format, digits = fmt$digits, width = width)
    }))

    names(Info) <- names(formats)
    print(Info)
  } else {
    cat("No coefficients\n\n")
  }

  cat("\n")
  invisible(x)
}

#' Summarize islasso.path Model at Specific Lambda
#'
#' Extracts coefficient estimates, standard errors and p-values from an \code{\link{islasso.path}} fit at a given regularization level \code{lambda}.
#'
#' @aliases summary.islasso.path print.summary.islasso.path
#'
#' @param object A fitted object of class \code{"islasso.path"}.
#' @param pval Numeric threshold for displaying coefficients. Only variables with \code{p-value <= pval} are printed.
#'        Unpenalized coefficients (like the intercept) are always shown.
#' @param use.t Logical. If \code{TRUE}, p-values are computed using a t-distribution with residual degrees of freedom.
#' @param lambda Numeric. Value of the regularization parameter at which the summary should be extracted.
#' @param ... Currently unused.
#'
#' @details Between two values of the path, estimates, standard errors and the other quantities
#' are linearly interpolated on the \eqn{\log\lambda} scale.
#'
#' @return An object of class \code{"summary.islasso.path"} containing filtered estimates and significance metrics.
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso.path}}, \code{\link{GoF.islasso.path}}, \code{\link{coef.islasso.path}},
#'          \code{\link{fitted.islasso.path}}, \code{\link{predict.islasso.path}},
#'          \code{\link{residuals.islasso.path}}, \code{\link{logLik.islasso.path}},
#'          \code{\link{deviance.islasso.path}}
#'
#' @examples
#' sim <- simulXy(n = 100, p = 20, beta = c(1, -1, 0.5, rep(0, 17)), seed = 1)
#' o <- islasso.path(y ~ ., data = sim$data, nlambda = 20)
#' summary(o, pval = 0.1, lambda = 5)
#'
#' @export
summary.islasso.path <- function(object, pval = 1, use.t = FALSE, lambda, ...) {
  temp <- list(...)

  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  # Handle single lambda case
  if (nlambda == 1) {
    lambda <- lambda.seq
    coef <- object$Coef[1, ]
    se <- object$SE[1, ]
    aic <- object$GoF[1, "AIC"]
    dispersion <- object$Info[1, "phi"]
    df0 <- object$Info[1, "df"]
    logLik <- object$Info[1, "logLik"]
    iter <- object$Info[1, "iter"]
  } else {
    # Handle missing lambda case
    if (missing(lambda)) return(print(object))

    # Validate lambda
    if (lambda < min(lambda.seq) || lambda > max(lambda.seq))
      stop("value of lambda out of bound")

    # Find indices for interpolation
    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda <= lambda.seq)[1]

    # Handle edge cases
    if (any(is.na(c(id1, id2)))) {
      id1 <- c(id1, id2)[!is.na(c(id1, id2))]
      coef <- object$Coef[id1, ]
      se <- object$SE[id1, ]
      aic <- object$GoF[id1, "AIC"]
      dispersion <- object$Info[id1, "phi"]
      df0 <- object$Info[id1, "df"]
      logLik <- object$Info[id1, "logLik"]
      iter <- object$Info[id1, "iter"]
    } else {
      # Use vectorized interpolation for efficiency
      params <- c("coef", "se", "aic", "dispersion", "df0", "logLik")
      values1 <- list(
        object$Coef[id1, ], object$SE[id1, ], object$GoF[id1, "AIC"],
        object$Info[id1, "phi"], object$Info[id1, "df"], object$Info[id1, "logLik"]
      )
      values2 <- list(
        object$Coef[id2, ], object$SE[id2, ], object$GoF[id2, "AIC"],
        object$Info[id2, "phi"], object$Info[id2, "df"], object$Info[id2, "logLik"]
      )

      results <- mapply(
        function(v1, v2) interpolate(v1, v2, lambda.seq[id1], lambda.seq[id2], lambda),
        values1, values2, SIMPLIFY = FALSE
      )

      coef <- results[[1]]
      se <- results[[2]]
      aic <- results[[3]]
      dispersion <- results[[4]]
      df0 <- results[[5]]
      logLik <- results[[6]]
      iter <- object$Info[id1, "iter"]
    }
  }

  # Prepare model statistics
  dev0 <- if (nlambda == 1) object$Info[1, "deviance"] else deviance.islasso.path(object, lambda = lambda)
  n <- object$Input$n
  p <- object$Input$p

  x <- model.matrix(object)
  y <- object$y
  offset <- object$offset
  family <- object$family

  eta <- drop(x %*% coef) + offset
  mu <- family$linkinv(eta)
  res <- drop(y - mu)
  rdf <- n - df0
  df <- c(p, rdf)

  # Calculate coefficient statistics
  chival <- (coef / se)
  type <- if (use.t) "t value" else "z value"
  pvalue <- if (use.t) 2 * pt(abs(chival), rdf, lower.tail = FALSE) else 2 * pnorm(-abs(chival))
  ptype <- if (use.t) "Pr(>|t|)" else "Pr(>|z|)"

  # Create coefficient matrix
  coefficients <- cbind(coef, se, chival, pvalue)
  dimnames(coefficients) <- list(colnames(object$Coef), c("Estimate", "Std. Error", type, ptype))

  # Calculate null model statistics
  intercept <- object$Input$intercept
  weights <- object$prior.weights
  wtdmu <- if (intercept) sum(weights * y) / sum(weights) else family$linkinv(offset)
  nulldev <- sum(family$dev.resids(y, wtdmu, weights))
  n.ok <- n - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)

  # Create output list
  out <- list(
    coefficients = coefficients,
    dispersion = dispersion,
    df = df,
    res = res,
    aic = aic,
    lambda = lambda,
    nulldev = nulldev,
    dev = dev0,
    df.null = nulldf,
    df.res = nulldf - (df0 - 1 * object$Input$intercept),
    family = object$family$family,
    iter = iter,
    pval = pval,
    temp = temp,
    call = object$call
  )

  class(out) <- "summary.islasso.path"
  return(out)
}

#' @export
print.summary.islasso.path <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  # Print call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  # Use custom digits if provided
  if (!is.null(x$temp$digits)) digits <- x$temp$digits

  # Print residuals summary
  resid <- x$res
  df <- x$df
  rdf <- df[2L]
  cat("Residuals:\n", sep = "")

  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    zz <- zapsmall(quantile(resid), digits + 1L)
    rq <- structure(zz, names = nam)
    print(rq, digits = digits, ...)
  } else {
    if (rdf > 0L) {
      print(resid, digits = digits, ...)
    } else {
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    }
  }
  cat("\n")

  # Filter and print coefficients
  coefs <- x$coefficients
  if (sum(coefs[, 4] <= x$pval) == 0) {
    warning("No coefficients lower than the selected p-value. The lowest p-value is printed.")
  }

  temp <- (coefs[, 4] <= x$pval)
  if (sum(temp) == 0) temp <- coefs[, 4] == min(coefs[, 4])
  coefs <- coefs[temp, , drop = FALSE]

  printCoefmat(coefs, digits = digits, signif.stars = TRUE,
               has.Pvalue = TRUE, na.print = "NA", cs.ind = 1L:2L, tst.ind = 3L, ...)

  # Format and print deviance info
  dev_info <- apply(cbind(
    paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
    format(c(x$nulldev, x$dev), digits = max(5L, digits + 1L)),
    " on",
    format(c(x$df.null, x$df.res), digits = digits),
    " degrees of freedom\n"
  ), 1L, paste, collapse = " ")

  cat(
    "\n(Dispersion parameter for ", x$family, " family taken to be ",
    format(x$dispersion), ")\n\n", dev_info, sep = ""
  )

  # Print AIC, lambda and iteration info
  cat(
    "AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
    "\nLambda: ", format(x$lambda, digits = max(4L, digits + 1L)),
    "\n\nNumber of IS iterations: ", x$iter, "\n\n",
    sep = ""
  )

  invisible(x)
}

#' Coefficient Profile and Diagnostic Plots for islasso.path
#'
#' Generates plots of coefficient profiles, standard errors, gradients, weights, or goodness-of-fit criteria
#' from a fitted \code{\link{islasso.path}} model.
#'
#' @param x An object of class \code{"islasso.path"}, typically created via \code{\link{islasso.path}}.
#' @param yvar Character. Specifies what to display on the y-axis. Choices are:
#'   \itemize{
#'     \item \code{"coefficients"} - coefficient paths over \code{log(lambda)},
#'     \item \code{"se"} - standard errors over \code{log(lambda)},
#'     \item \code{"score"} - standardized scores \eqn{x_j^\top W(z - X\hat\beta)/\lambda_j} over
#'       \code{log(lambda)}, the induced-smoothing analogue of the lasso optimality conditions
#'       (dashed lines at \eqn{\pm 1}),
#'     \item \code{"gradient"} - gradient of the smoothed estimating equation over \code{log(lambda)}
#'       (a convergence diagnostic),
#'     \item \code{"weight"} - mixture weights used in smoothing,
#'     \item \code{"gof"} - goodness-of-fit values.
#'   }
#' @param gof Character. Criterion used for highlighting active variables.
#'        Choices: \code{"none"}, \code{"AIC"}, \code{"BIC"}, \code{"AICc"}, \code{"eBIC"}, \code{"GCV"}, \code{"GIC"}.
#' @param label Logical. Whether to annotate curves with variable names.
#' @param legend Logical. Whether to display a plot legend.
#' @param ... Additional graphical parameters, e.g. \code{main}, \code{xlab}, \code{ylab}, \code{xlim}, \code{ylim},
#'        \code{lty}, \code{col}, \code{lwd}, \code{cex.axis}, \code{cex.lab}, \code{cex.main}, \code{gof_lty},
#'        \code{gof_col}, \code{gof_lwd}.
#'
#' @details
#' This function visualizes the behavior of the solution path across a sequence of lambda values,
#' helping diagnose coefficient shrinkage, influence of penalty, and variable selection stability.
#' All the profiles are drawn against \eqn{\log\lambda}.
#'
#' The standardized score of a covariate (computed on the standardized scale on which the penalty
#' acts) equals, at the solution, the smoothed sign \eqn{P_j} (lasso) or
#' \eqn{\alpha P_j + (1-\alpha)\hat\beta_j} (elastic net): it is close to \eqn{\pm 1} for the
#' covariates with sizeable effects and strictly inside \eqn{(-1, 1)} for those shrunk towards zero,
#' as the subgradient conditions of the lasso. The gradient, instead, is the residual of the smoothed
#' estimating equation: it is zero up to the convergence tolerance, and its irregular values for
#' large \eqn{\lambda} reflect the very large curvature of the smoothed penalty when the standard
#' errors of null coefficients are tiny, not features of the model.
#'
#' @return A \code{ggplot} object (printed when not assigned).
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso.path}}, \code{\link{GoF.islasso.path}}, \code{\link{summary.islasso.path}},
#'          \code{\link{coef.islasso.path}}, \code{\link{fitted.islasso.path}}, \code{\link{predict.islasso.path}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100; p <- 30
#' beta <- c(runif(10, -2, 2), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' fit <- islasso.path(y ~ ., data = sim$data, family = gaussian(), nlambda = 30)
#' plot(fit, yvar = "coefficients", gof = "AICc", label = TRUE)
#' plot(fit, yvar = "se", gof = "AICc")
#' plot(fit, yvar = "score", gof = "AICc")
#' plot(fit, yvar = "gradient", gof = "AICc")
#' plot(fit, yvar = "gof", gof = "AICc")
#' }
#'
#' @export
plot.islasso.path <- function(x,
                              yvar = c("coefficients", "se", "score", "gradient", "weight", "gof"),
                              gof = c("none", "AIC", "BIC", "AICc", "eBIC", "GCV", "GIC"),
                              label = FALSE,
                              legend = FALSE,
                              ...) {
  # Validate inputs
  if (!inherits(x, "islasso.path"))
    stop("x is not an object of islasso.path class")

  object <- x
  label <- abs(label)
  yvar <- match.arg(yvar)
  gof <- match.arg(gof)

  if (yvar == "gof" && gof == "none")
    stop("You have to select one criterion between AIC, BIC, AICc, eBIC, GCV, GIC")

  # Extract model information
  id.best <- apply(object$GoF, 2, which.min)
  lambda <- object$Info[, "lambda"]
  loglambda <- log(lambda)
  nlambda <- length(lambda)

  if (nlambda == 1)
    stop("no path in islasso.path fit!")

  # Process coefficients
  intercept <- object$Input$intercept
  coef1 <- object$Coef
  if (intercept) coef1 <- coef1[, -1]

  unactive <- abs(coef1[id.best[gof], ]) <= 1E-6
  active <- !unactive
  if (gof == "none") active <- rep(TRUE, length(active))

  # Process optional plotting parameters
  dots <- list(...)
  dots$lty <- if (is.null(dots$lty)) {
    if (gof != "none") ifelse(unactive, 2L, 1L) else rep(1L, length(unactive))
  } else dots$lty

  dots$col <- if (is.null(dots$col)) {
    if (gof != "none") ifelse(unactive, "gray80", "gray20") else rep("gray50", length(unactive))
  } else dots$col

  dots$lwd <- if (is.null(dots$lwd)) 0.75 else dots$lwd
  dots$gof_lty <- if (is.null(dots$gof_lty)) "dashed" else dots$gof_lty
  dots$gof_col <- if (is.null(dots$gof_col)) "red" else dots$gof_col
  dots$cex.axis <- if (is.null(dots$cex.axis)) 12L else dots$cex.axis
  dots$cex.lab <- if (is.null(dots$cex.lab)) 12L else dots$cex.lab

  # Create different plots based on the yvar parameter
  if (yvar == "coefficients") {
    p <- create_coef_plot(coef1, loglambda, label, id.best, gof, dots,
                          active, unactive, legend, nlambda)
  } else if (yvar == "se") {
    se1 <- object$SE
    if (intercept) se1 <- se1[, -1]
    p <- create_se_plot(se1, coef1, loglambda, label, id.best, gof, dots,
                        active, unactive, legend, nlambda)
  } else if (yvar == "weight") {
    weight1 <- object$Weight
    if (intercept) weight1 <- weight1[, -1]
    p <- create_weight_plot(weight1, coef1, loglambda, label, id.best, gof, dots,
                            active, unactive, legend, nlambda)
  } else if (yvar == "score") {
    if (is.null(object$Score))
      stop("the object has no standardized scores: refit it with islasso >= 1.7.0", call. = FALSE)
    score1 <- object$Score
    if (intercept) score1 <- score1[, -1, drop = FALSE]
    if (is.null(dots$ylab)) dots$ylab <- "Standardized score"
    p <- create_coef_plot(score1, loglambda, label, id.best, gof, dots,
                          active, unactive, legend, nlambda) +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", col = "gray50")
  } else if (yvar == "gradient") {
    grad <- object$Gradient
    if (intercept) grad <- grad[, -1, drop = FALSE]
    if (is.null(dots$ylab)) dots$ylab <- "Gradient"
    p <- create_coef_plot(grad, loglambda, label, id.best, gof, dots,
                          active, unactive, legend, nlambda)
  } else if (yvar == "gof") {
    p <- create_gof_plot(object, loglambda, id.best, gof, dots)
  }

  return(p)
}

create_coef_plot <- function(coef1, loglambda, label, id.best, gof, dots,
                             active, unactive, legend, nlambda) {
  # Set default axis labels if not provided
  if (is.null(dots$ylab)) dots$ylab <- "Coefficients"
  if (is.null(dots$xlab)) dots$xlab <- expression(log(lambda))
  if (is.null(dots$xlim)) dots$xlim <- range(loglambda) + c(-label, 0)

  # Create the plot
  p <- data.frame(
    lambda = loglambda,
    value = c(coef1),
    varname = factor(rep(names(unactive), each = nlambda), levels = names(unactive))
  ) |>
    ggplot() +
    geom_line(aes(x = lambda, y = value, linetype = varname),
              color = rep(dots$col, each = nlambda), linewidth = dots$lwd) +
    scale_linetype_manual(values = rep(dots$lty, each = nlambda)) +
    geom_hline(yintercept = 0.0, linetype = "dotted", col = "gray50") +
    theme_bw() +
    xlab(dots$xlab) +
    ylab(dots$ylab) +
    theme(legend.position = "none",
          axis.text = element_text(size = dots$cex.axis),
          axis.title = element_text(size = dots$cex.lab)) +
    xlim(dots$xlim)

  # Add labels if requested
  if (label) {
    p <- p + annotate("text",
                      x = rep(min(loglambda) - label, sum(active)),
                      y = coef1[1, active],
                      label = colnames(coef1)[active], size = 5L)
  }

  # Add vertical line for best model if gof criterion specified
  if (gof != "none") {
    p <- p + geom_vline(xintercept = loglambda[id.best[gof]],
                        col = dots$gof_col, linetype = dots$gof_lty)
  }

  # Add legend if requested
  if (legend && gof != "none") {
    p <- p + annotate("text",
                      x = quantile(loglambda, probs = .975),
                      y = max(coef1[1, ]),
                      label = paste0("min.", gof), size = 4L, color = dots$gof_col)
  }

  return(p)
}

create_se_plot <- function(se1, coef1, loglambda, label, id.best, gof, dots,
                           active, unactive, legend, nlambda) {
  # Set default axis labels if not provided
  if (is.null(dots$ylab)) dots$ylab <- "Std.Errs"
  if (is.null(dots$xlab)) dots$xlab <- expression(log(lambda))
  if (is.null(dots$xlim)) dots$xlim <- range(loglambda) + c(-label, 0)

  # Create the plot
  p <- data.frame(
    lambda = loglambda,
    value = c(se1),
    varname = factor(rep(names(unactive), each = nlambda), levels = names(unactive))
  ) |>
    ggplot() +
    geom_line(aes(x = lambda, y = value, linetype = varname),
              color = rep(dots$col, each = nlambda), linewidth = dots$lwd) +
    scale_linetype_manual(values = rep(dots$lty, each = nlambda)) +
    geom_hline(yintercept = 0.0, linetype = "dotted", col = "gray50") +
    theme_bw() +
    xlab(dots$xlab) +
    ylab(dots$ylab) +
    theme(legend.position = "none",
          axis.text = element_text(size = dots$cex.axis),
          axis.title = element_text(size = dots$cex.lab)) +
    xlim(dots$xlim)

  # Add labels if requested
  if (label) {
    p <- p + annotate("text",
                      x = rep(min(loglambda) - label, sum(active)),
                      y = se1[1, active],
                      label = colnames(coef1)[active], size = 5L)
  }

  # Add vertical line for best model if gof criterion specified
  if (gof != "none") {
    p <- p + geom_vline(xintercept = loglambda[id.best[gof]],
                        col = dots$gof_col, linetype = dots$gof_lty)
  }

  # Add legend if requested
  if (legend && gof != "none") {
    p <- p + annotate("text",
                      x = quantile(loglambda, probs = .975),
                      y = max(se1[1, ]),
                      label = paste0("min.", gof), size = 4L, color = dots$gof_col)
  }

  return(p)
}

create_weight_plot <- function(weight1, coef1, loglambda, label, id.best, gof, dots,
                               active, unactive, legend, nlambda) {
  # Set default axis labels if not provided
  if (is.null(dots$ylab)) dots$ylab <- "Mixture weight"
  if (is.null(dots$xlab)) dots$xlab <- expression(log(lambda))
  if (is.null(dots$xlim)) dots$xlim <- range(loglambda) + c(-label, 0)

  # Create the plot
  p <- data.frame(
    lambda = loglambda,
    value = c(weight1),
    varname = factor(rep(names(unactive), each = nlambda), levels = names(unactive))
  ) |>
    ggplot() +
    geom_line(aes(x = lambda, y = value, linetype = varname),
              color = rep(dots$col, each = nlambda), linewidth = dots$lwd) +
    scale_linetype_manual(values = rep(dots$lty, each = nlambda)) +
    geom_hline(yintercept = 0.0, linetype = "dotted", col = "gray50") +
    theme_bw() +
    xlab(dots$xlab) +
    ylab(dots$ylab) +
    theme(legend.position = "none",
          axis.text = element_text(size = dots$cex.axis),
          axis.title = element_text(size = dots$cex.lab)) +
    xlim(dots$xlim)

  # Add labels if requested
  if (label) {
    p <- p + annotate("text",
                      x = rep(min(loglambda) - label, sum(active)),
                      y = weight1[1, active],
                      label = colnames(coef1)[active], size = 5L)
  }

  # Add vertical line for best model if gof criterion specified
  if (gof != "none") {
    p <- p + geom_vline(xintercept = loglambda[id.best[gof]],
                        col = dots$gof_col, linetype = dots$gof_lty)
  }

  # Add legend if requested
  if (legend && gof != "none") {
    p <- p + annotate("text",
                      x = quantile(loglambda, probs = .975),
                      y = .95,
                      label = paste0("min.", gof), size = 4L, color = dots$gof_col)
  }

  return(p)
}


create_gof_plot <- function(object, loglambda, id.best, gof, dots) {
  # Set default axis labels if not provided
  if (is.null(dots$xlab)) dots$xlab <- expression(log(lambda))
  if (is.null(dots$ylab)) dots$ylab <- gof

  # Create the plot
  p <- data.frame(
    lambda = loglambda,
    value = object$GoF[, gof]
  ) |>
    ggplot() +
    geom_line(aes(x = lambda, y = value), linewidth = dots$lwd) +
    theme_bw() +
    xlab(dots$xlab) +
    ylab(dots$ylab) +
    theme(axis.text = element_text(size = dots$cex.axis),
          axis.title = element_text(size = dots$cex.lab)) +
    geom_vline(xintercept = loglambda[id.best[gof]],
               col = dots$gof_col, linetype = dots$gof_lty)

  return(p)
}

#' Prediction Method for islasso.path Objects
#'
#' Generates predictions from a fitted \code{\link{islasso.path}} model at one or more lambda values.
#' Supports various output types including linear predictors, response scale, class labels, and coefficients.
#'
#' @param object A fitted model object of class \code{"islasso.path"}.
#' @param newdata Optional data frame containing covariates for prediction.
#'                If omitted, returns fitted values from the original model.
#' @param type Character. Type of prediction:
#'   \itemize{
#'     \item \code{"link"} (default) - linear predictor scale,
#'     \item \code{"response"} - original response scale,
#'     \item \code{"coefficients"} - estimated coefficients,
#'     \item \code{"class"} - predicted class labels (only for binomial models).
#'   }
#' @param lambda Numeric value(s). Specific lambda value(s) at which predictions are required.
#'               If missing, predictions are computed for the full lambda sequence.
#' @param ... Additional arguments passed to lower-level methods.
#'
#' @details For values of \code{lambda} not in the path, the coefficients are linearly
#' interpolated on the \eqn{\log\lambda} scale.
#'
#' @return A vector or a matrix (one column per value of \code{lambda}) depending on \code{type}.
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso.path}}, \code{\link{summary.islasso.path}}, \code{\link{coef.islasso.path}},
#'          \code{\link{GoF.islasso.path}}, \code{\link{fitted.islasso.path}}, \code{\link{logLik.islasso.path}},
#'          \code{\link{residuals.islasso.path}}, \code{\link{deviance.islasso.path}}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 30
#' beta <- c(runif(10, -3, 3), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' fit <- islasso.path(y ~ ., data = sim$data, family = gaussian(), nlambda = 30)
#' optimal <- GoF.islasso.path(fit, plot = FALSE)
#' pred <- predict(fit, newdata = sim$data[1:5, ], type = "response",
#'                 lambda = optimal$lambda.min[c("AIC", "BIC")])
#' pred
#'
#' @export
predict.islasso.path <- function(object,
                                 newdata,
                                 type = c("link", "response", "coefficients", "class"),
                                 lambda,
                                 ...) {
  # Validate inputs
  type <- match.arg(type)
  family <- object$family

  # Check if 'class' type is being used with a non-binomial family
  if (type == "class" && family$family != "binomial") {
    stop(gettextf("Type 'class' is available only for the %s family",
                  sQuote(binomial()$family)),
         domain = NA)
  }

  # Get model terms
  tt <- terms(object)

  # Prepare model matrix and offset
  if (missing(newdata) || is.null(newdata)) {
    # Use original data if newdata not provided
    X <- model.matrix(object)
    offset <- object$offset
  } else {
    # Process newdata
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels)

    # Check model frame classes
    if (!is.null(cl <- attr(Terms, "dataClasses"))) {
      .checkMFClasses(cl, m)
    }

    # Create model matrix
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

    # Calculate offset
    offset <- rep(0, nrow(X))

    # Add offset from terms
    if (!is.null(off.num <- attr(tt, "offset"))) {
      for (i in off.num) {
        offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
      }
    }

    # Add offset from call
    if (!is.null(object$call$offset)) {
      offset <- offset + eval(object$call$offset, newdata)
    }
  }

  # Get lambda sequence from object
  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  # Handle lambda parameter
  if (missing(lambda) || nlambda == 1) {
    # Use all lambdas from the model or the single lambda if there's only one
    beta <- object$Coef
    predictor <- tcrossprod(X, beta)
    lambda <- lambda.seq
    nlmb <- nlambda
  } else {
    # Validate provided lambda values
    if (any(lambda < min(lambda.seq)) || any(lambda > max(lambda.seq))) {
      stop("value of lambda out of bound")
    }

    # Sort lambda values and interpolate coefficients
    lambda <- sort(lambda)
    nlmb <- length(lambda)
    beta <- matrix(0, nrow = ncol(X), ncol = nlmb)

    # Interpolate coefficients for each lambda value
    for (i in seq_len(nlmb)) {
      id1 <- rev(which(lambda[i] >= lambda.seq))[1]
      id2 <- which(lambda[i] <= lambda.seq)[1]
      beta[, i] <- interpolate(object$Coef[id1,],
                               object$Coef[id2,],
                               lambda.seq[id1],
                               lambda.seq[id2],
                               lambda[i])
    }

    predictor <- drop(X %*% beta)
  }

  # Add offset to predictor if available
  if (!is.null(offset)) {
    predictor <- predictor + offset
  }

  # Generate output based on requested type
  out <- switch(type,
                "link" = {
                  predictor
                },
                "response" = {
                  family$linkinv(predictor)
                },
                "coefficients" = {
                  t(beta)
                },
                "class" = {
                  mu <- family$linkinv(predictor)
                  1 * (mu >= 0.5)
                })

  # Add column names for matrix output
  if (is.matrix(out)) {
    colnames(out) <- paste0("lambda = ", formatC(lambda, format = "f", digits = 3))
  }

  return(out)
}

#' @export
coef.islasso.path <- function(object, lambda, ...) {
  if (missing(lambda)) {
    return(object$Coef)
  }

  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  if (nlambda == 1) {
    return(object$Coef[1, ])
  }

  if (lambda < min(lambda.seq) || lambda > max(lambda.seq)) {
    stop("value of lambda out of bound")
  }

  id1 <- rev(which(lambda >= lambda.seq))[1]
  id2 <- which(lambda <= lambda.seq)[1]
  interpolate(object$Coef[id1, ], object$Coef[id2, ],
              lambda.seq[id1], lambda.seq[id2], lambda)
}

#' @export
fitted.islasso.path <- function(object, lambda, ...) {
  if (missing(lambda)) {
    return(object$Fitted.values)
  }

  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  if (nlambda == 1) {
    return(object$Fitted.values[1, ])
  }

  if (lambda < min(lambda.seq) || lambda > max(lambda.seq)) {
    stop("value of lambda out of bound")
  }

  id1 <- rev(which(lambda >= lambda.seq))[1]
  id2 <- which(lambda <= lambda.seq)[1]
  interpolate(object$Fitted.values[id1, ], object$Fitted.values[id2, ],
              lambda.seq[id1], lambda.seq[id2], lambda)
}

#' @export
residuals.islasso.path <- function(object, type = c("working", "response", "deviance", "pearson"),
                                   lambda, ...) {
  type <- match.arg(type)

  if (missing(lambda)) {
    return(object$Residuals)
  }

  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  if (nlambda == 1) {
    res <- object$Residuals[1, ]
    mu <- object$Fitted.values[1, ]
  } else {
    if (lambda < min(lambda.seq) || lambda > max(lambda.seq)) {
      stop("value of lambda out of bound")
    }

    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda <= lambda.seq)[1]
    res <- interpolate(object$Residuals[id1, ], object$Residuals[id2, ],
                       lambda.seq[id1], lambda.seq[id2], lambda)
    mu <- interpolate(object$Fitted.values[id1, ], object$Fitted.values[id2, ],
                      lambda.seq[id1], lambda.seq[id2], lambda)
  }

  y <- object$y
  weights <- object$prior.weights
  family <- object$family

  switch(type,
         "working" = res,
         "response" = y - mu,
         "deviance" = sign(y - mu) * sqrt(family$dev.resids(y, mu, weights)),
         "pearson" = (y - mu) * sqrt(weights) / sqrt(family$variance(mu)))
}

#' @export
deviance.islasso.path <- function(object, lambda, ...) {
  if (missing(lambda)) {
    return(object$Info[, "deviance"])
  }

  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  if (nlambda == 1) {
    return(object$Info[1, "deviance"][[1]])
  }

  if (lambda < min(lambda.seq) || lambda > max(lambda.seq)) {
    stop("value of lambda out of bound")
  }

  id1 <- rev(which(lambda >= lambda.seq))[1]
  id2 <- which(lambda <= lambda.seq)[1]
  interpolate(object$Info[id1, "deviance"], object$Info[id2, "deviance"],
              lambda.seq[id1], lambda.seq[id2], lambda)[[1]]
}

#' @export
logLik.islasso.path <- function(object, lambda, ...) {
  if (missing(lambda)) {
    return(object$Info[, "logLik"])
  }

  lambda.seq <- object$Info[, "lambda"]
  nlambda <- length(lambda.seq)

  if (nlambda == 1) {
    ll <- object$Info[1, "logLik"]
    df <- object$Info[1, "df"]
    phi <- object$Info[1, "phi"]
  } else {
    if (lambda < min(lambda.seq) || lambda > max(lambda.seq)) {
      stop("value of lambda out of bound")
    }

    id1 <- rev(which(lambda >= lambda.seq))[1]
    id2 <- which(lambda <= lambda.seq)[1]
    ll <- interpolate(object$Info[id1, "logLik"], object$Info[id2, "logLik"],
                      lambda.seq[id1], lambda.seq[id2], lambda)
    df <- interpolate(object$Info[id1, "df"], object$Info[id2, "df"],
                      lambda.seq[id1], lambda.seq[id2], lambda)
    phi <- interpolate(object$Info[id1, "phi"], object$Info[id2, "phi"],
                       lambda.seq[id1], lambda.seq[id2], lambda)
  }

  out <- list(loglik = ll, df = df, object = object, phi = phi)
  class(out) <- "logLik.islasso.path"
  out
}


#' @export
print.logLik.islasso.path <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n'log Lik.' ", paste(format(c(x$loglik), digits = digits), collapse = ", "),
      " (df = ", format(x$df), ")\n\n", sep = "")
  invisible(x)
}
