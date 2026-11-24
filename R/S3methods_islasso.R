#' @export
print.islasso <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(round(x$coef, 6), digits = digits), print.gap = 2, quote = FALSE)
  } else {
    cat("No coefficients\n\n")
  }

  residual_df <- format(signif(x$df.null - (x$rank - 1 * x$internal$intercept * x$internal$hi[1]), digits))
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null);", residual_df, "Residual\n")

  cat("Null Deviance:", format(signif(x$null.deviance, digits)),
      "\nResidual Deviance:", format(signif(x$deviance, digits)),
      "\nAIC:", format(signif(x$aic, digits)),
      "\nLambda:", format(signif(x$lambda, digits)), "\n")

  invisible(x)
}

#' Summarize islasso Fitted Model
#'
#' Provides a concise summary of a fitted \code{\link{islasso}} model, including p-values and optional filtering.
#'
#' @aliases summary.islasso print.summary.islasso
#'
#' @param object A fitted model of class \code{"islasso"}.
#' @param pval Numeric threshold for displaying coefficients. Only those with \eqn{p \le} \code{pval} are printed. Unpenalized coefficients (like intercepts) are always shown.
#' @param which Optional. Specifies a subset of coefficients to test. If missing, all parameters are evaluated.
#' @param use.t Logical. If \code{TRUE}, p-values are computed using the t-distribution and residual degrees of freedom.
#' @param type.pval Character. Type of p-value approximation. Only \code{"wald"} (default) is implemented.
#' @param ... Additional arguments (not currently used).
#'
#' @return An object of class \code{"summary.islasso"} containing:
#'   \item{coefficients}{Coefficient estimates and related statistics}
#'   \item{pval}{Threshold used to filter coefficients}
#'   \item{call}{Original model call}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso.fit}}, \code{\link{residuals.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}}, \code{\link{deviance.islasso}}
#'
#' @examples
#' \dontrun{
#' # Assuming object `o` from an islasso fit
#' summary(o, pval = 0.1)  # Show coefficients with p <= 0.1
#' }
#'
#' @export
summary.islasso <- function(object, pval = 1, which, use.t = FALSE,
                            type.pval = "wald", ...) {
  # Process optional arguments
  dots <- list(...)
  type.pval <- match.arg(type.pval)

  # Extract model information
  n <- object$internal$n
  p <- object$internal$p
  coef <- object$coef

  # Handle coefficient selection
  if (missing(which)) {
    which <- seq_len(p)
  } else {
    if (is.character(which)) {
      which <- pmatch(which, names(coef))
    }
    if (length(which) == 0) {
      stop("Coefficient name not present!")
    }
  }

  # Extract specific coefficients and statistics
  coef <- coef[which]
  se <- object$se[which]
  h <- object$internal$hi[which]

  # Set family-specific options
  if (object$family$family != "gaussian") {
    use.t <- FALSE
  }

  # Calculate test statistics based on type
  if (type.pval == "wald") {
    chival <- (coef / se)
  }

  # Build coefficient table
  coefficients <- round(cbind(coef, se, h, chival), 6)

  # Set column names based on test type
  type <- if (use.t) "t value" else "z value"
  pvalue <- if (type == "t value") {
    2 * pt(abs(chival), object$df.residual, lower.tail = FALSE)
  } else {
    2 * pnorm(-abs(chival))
  }
  ptype <- if (type == "t value") "Pr(>|t|)" else "Pr(>|z|)"

  # Combine into final coefficient table
  coefficients <- cbind(coefficients, pvalue)
  dimnames(coefficients) <- list(
    object$internal$nms[which],
    c("Estimate", "Std. Error", "Df", type, ptype)
  )

  # Prepare residuals and other model information
  res <- residuals(object, type = "deviance")
  df_residual <- object$df.null - (object$rank - 1 * object$internal$intercept * object$internal$hi[1])

  # Create and return summary object
  out <- list(
    coefficients = coefficients,
    dispersion = object$dispersion,
    df = c(object$internal$p, object$df.residual),
    res = res,
    aic = object$aic,
    lambda = object$lambda,
    nulldev = object$null.deviance,
    dev = object$deviance,
    df.null = object$df.null,
    df.res = df_residual,
    family = object$family,
    iter = object$iter,
    pval = pval,
    unpenalized = object$internal$unpenalized[which],
    temp = dots,
    call = object$call
  )

  class(out) <- "summary.islasso"
  return(out)
}

#' @export
print.summary.islasso <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  # Print model call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  # Use custom digits if provided
  if (!is.null(x$temp$digits)) {
    digits <- x$temp$digits
  }

  # Print residual summary
  resid <- x$res
  df <- x$df
  rdf <- df[2L]

  cat("Residuals:\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  zz <- zapsmall(quantile(resid), digits + 1L)
  rq <- structure(zz, names = nam)
  print(rq, digits = digits, ...)
  cat("\n")

  # Filter and print coefficients
  coefs <- x$coefficients
  if (sum(coefs[, 5] <= x$pval) == 0) {
    warning("No coefficients lower than the selected p-value. The lowest p-value is printed.")
  }

  # Select coefficients to display
  temp <- (coefs[, 5] <= x$pval) | x$unpenalized
  if (sum(temp) == 0) {
    temp <- coefs[, 5] == min(coefs[, 5])
  }
  coefs <- coefs[temp, , drop = FALSE]

  # Print coefficient table
  printCoefmat(coefs,
               digits = digits,
               signif.stars = TRUE,
               has.Pvalue = TRUE,
               na.print = "NA",
               cs.ind = 1L:2L,
               tst.ind = 3L:4L,
               ...)

  # Print deviance information
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
      format(x$dispersion), ")\n\n",
      apply(cbind(
        paste(format(c("Null", "Residual"), justify = "right"), "deviance:"),
        format(c(x$nulldev, x$dev), digits = max(5L, digits + 1L)),
        " on",
        format(c(x$df.null, x$df.res), digits = digits),
        " degrees of freedom\n"),
        1L, paste, collapse = " "),
      sep = "")

  # Print model fit statistics
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
      "\nLambda: ", format(x$lambda, digits = max(4L, digits + 1L)),
      "\n\nNumber of Newton-Raphson iterations: ", x$iter, "\n\n",
      sep = "")

  invisible(x)
}

#' Diagnostic Plots for islasso Models
#'
#' Produces standard diagnostic plots for a fitted \code{\link{islasso}} model to assess residuals,
#' model fit, and variance structure.
#'
#' @param x An object of class \code{"islasso"}, typically created via \code{\link{islasso}}.
#' @param ... Additional graphical parameters passed to the underlying \code{plot()} functions.
#'
#' @details
#' Generates a 2x2 grid of diagnostic plots:
#' \itemize{
#'   \item Top-left: Deviance residuals vs fitted values.
#'   \item Top-right: Normal Q-Q plot of standardized deviance residuals (red line = reference).
#'   \item Bottom-left: Squared standardized Pearson residuals vs fitted values.
#'   \item Bottom-right: Working response vector vs linear predictor.
#' }
#'
#' These plots help assess the assumptions of linearity, homoscedasticity, and residual normality in penalized regression.
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}},
#' \code{\link{logLik.islasso}}, \code{\link{predict.islasso}}, \code{\link{deviance.islasso}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 100; p <- 100
#'   beta <- c(runif(20, -3, 3), rep(0, p - 20))
#'   sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#'   fit <- islasso(y ~ ., data = sim$data, family = gaussian(), lambda = 2)
#'   plot(fit)
#' }
#'
#' @export
plot.islasso <- function(x, ...) {
  # Process optional parameters with defaults
  opts <- list(...)
  # ggplot equivalents for base R parameters
  point_size <- if (is.null(opts$cex)) 1 else opts$cex * 3  # Scaling factor for ggplot
  point_shape <- if (is.null(opts$pch)) 19 else opts$pch
  line_color <- if (is.null(opts$col)) "red" else opts$col
  line_width <- if (is.null(opts$lwd)) 1 else opts$lwd * 0.5  # Scaling factor for ggplot

  # Extract model components
  family <- x$family
  nX <- model.matrix(x)
  y <- x$y
  eta <- x$linear.predictors
  mu <- x$fitted.values
  dev <- residuals(x, type = "deviance")
  sdev <- dev / sqrt(x$dispersion)

  pea <- residuals(x, type = "pearson")
  spea <- pea / sqrt(x$dispersion)

  # Calculate hat matrix and related statistics
  w2 <- x$weights
  invH <- x$internal$invH

  w <- sqrt(w2)
  WX <- nX * w

  H <- WX %*% invH %*% t(WX)
  h <- diag(H)
  p <- x$rank
  rp <- spea / sqrt(1 - h)
  rd <- sdev / sqrt(1 - h)
  res <- sign(sdev) * sqrt(sdev^2 + h * rp^2)

  # For QQ plot
  rd_clean <- na.omit(rd)
  qq_data <- data.frame(
    theoretical = qnorm(ppoints(length(rd_clean))),
    sample = sort(rd_clean)
  )
  qq_line <- data.frame(
    theoretical = c(min(qq_data$theoretical), max(qq_data$theoretical))
  )
  qq_line$sample <- quantile(rd_clean, c(0.25, 0.75))[1] +
    diff(quantile(rd_clean, c(0.25, 0.75))) /
    diff(qnorm(c(0.25, 0.75))) *
    (qq_line$theoretical - qnorm(0.25))

  # For link function
  zmod <- eta + (y - mu) / w2
  lm_result <- lm(zmod ~ eta)
  link_line <- data.frame(
    eta = c(min(eta), max(eta))
  )
  link_line$predicted <- predict(lm_result, newdata = data.frame(eta = link_line$eta))

  # Create data frames for plots
  plot_data <- data.frame(
    eta = eta,
    res = res,
    mu = mu,
    spea = spea,
    zmod = zmod,
    rd = rd
  )

  # Plot 1: Residuals vs Linear predictor
  p1 <- ggplot(plot_data, aes(x = eta, y = res)) +
    geom_point(size = point_size, shape = point_shape) +
    labs(x = "Linear predictor", y = "Residuals") +
    theme_bw()

  # Plot 2: QQ-plot of deviance residuals
  p2 <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_point(size = point_size, shape = point_shape) +
    geom_line(data = qq_line, aes(x = theoretical, y = sample),
              color = line_color, linewidth = line_width) +
    labs(x = "Theoretical Quantiles", y = "Ordered deviance residuals") +
    theme_bw()

  # Plot 3: Variance function
  p3 <- ggplot(plot_data, aes(x = mu, y = spea^2)) +
    geom_point(size = point_size, shape = point_shape) +
    labs(x = "Fitted values", y = "Squared stand Pearson residuals",
         title = "Variance function") +
    theme_bw()

  # Plot 4: Link function
  p4 <- ggplot(plot_data, aes(x = eta, y = zmod)) +
    geom_point(size = point_size, shape = point_shape) +
    geom_line(data = link_line, aes(x = eta, y = predicted),
              color = line_color, linewidth = line_width) +
    labs(x = "Linear predictor", y = "Working vector",
         title = "Link function") +
    theme_bw()

  # Arrange the four plots in a 2x2 grid using gridExtra
  grid_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

  # Return the arranged plot
  return(grid_plot)
}

#' Prediction Method for islasso Objects
#'
#' Computes predictions from a fitted \code{\link{islasso}} model object. Multiple output types supported, including response scale, linear predictor, and coefficient values.
#'
#' @param object A fitted model of class \code{"islasso"}.
#' @param newdata Optional data frame containing predictors for prediction. If omitted, the fitted model matrix is used.
#' @param type Character. Specifies the prediction scale:
#'   \itemize{
#'     \item \code{"link"} (default): linear predictor scale;
#'     \item \code{"response"}: original response scale;
#'     \item \code{"coefficients"}: estimated coefficients;
#'     \item \code{"class"}: predicted class (only for \code{binomial()} family);
#'     \item \code{"terms"}: contribution of each term to the linear predictor.
#'   }
#' @param se.fit Logical. Whether to compute standard errors/confidence intervals.
#' @param ci Optional. Precomputed matrix of confidence intervals (2 columns).
#' @param type.ci Type of interval. Only \code{"wald"} is implemented.
#' @param level Confidence level for intervals. Default is \code{0.95}.
#' @param terms If \code{type = "terms"}, optionally specify which terms to extract.
#' @param na.action Function to handle missing values in \code{newdata}. Default: \code{na.pass}.
#' @param ... Additional arguments passed to downstream methods.
#'
#' @return A numeric vector, matrix, or list depending on \code{type}.
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso}}, \code{\link{summary.islasso}}, \code{\link{logLik.islasso}}, \code{\link{residuals.islasso}}, \code{\link{deviance.islasso}}
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 100
#' beta <- c(runif(20, -3, 3), rep(0, p - 20))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' fit <- islasso(y ~ ., data = sim$data, family = gaussian(), lambda = 2)
#' predict(fit, type = "response")
#'
#' @export
predict.islasso <- function(object, newdata = NULL,
                            type = c("link", "response", "coefficients", "class", "terms"),
                            se.fit = FALSE, ci = NULL,
                            type.ci = c("wald", "score"),
                            level = 0.95,
                            terms = NULL,
                            na.action = na.pass, ...) {
  # Match arguments
  type <- match.arg(type)
  type.ci <- match.arg(type.ci)

  # Initialize variables
  ci.fit <- NULL
  na.act <- object$na.action
  object$na.action <- NULL

  # Validate type for binomial family
  if (type == "class" && family(object)$family != "binomial") {
    stop(gettextf("Type 'class' is available only for the %s family",
                  sQuote(binomial()$family)),
         domain = NA)
  }

  # Handle coefficients request
  if (type == "coefficients") {
    fit <- coef(object)
    if (se.fit) {
      ci.fit <- ci.fitted.islasso(object, X, ci,
                                  type.ci = type.ci,
                                  conf.level = level,
                                  only.ci = TRUE)
    }
  } else {
    # Handle prediction with or without new data
    if (missing(newdata)) {
      fit <- switch(type,
                    "response" = ,
                    "class" = ,
                    "link" = object$linear.predictors,
                    "terms" = predislasso(object, type = "terms", terms = terms))
      X <- model.matrix(object)
    } else {
      fit <- predislasso(object, newdata,
                         type = if (type %in% c("link", "class")) "response" else type,
                         terms = terms,
                         na.action = na.action)
      X <- attr(fit, "X")
    }

    # Clean up attributes
    attr(fit, "X") <- NULL

    # Handle NA values
    if (missing(newdata) && !is.null(na.act)) {
      fit <- napredict(na.act, fit)
    }

    # Calculate confidence intervals if requested
    if (se.fit && !(type %in% c("class", "terms"))) {
      ci.fit <- ci.fitted.islasso(object, X, ci,
                                  type.ci = type.ci,
                                  conf.level = level)

      # Handle NA values in confidence intervals
      if (missing(newdata) && !is.null(na.act)) {
        ci.fit[, 1L] <- napredict(na.act, ci.fit[, 1L])
        ci.fit[, 2L] <- napredict(na.act, ci.fit[, 2L])
      }
    }
  }

  # Combine fit and confidence intervals
  out <- if (type != "terms") cbind("Fit" = fit, ci.fit) else fit

  # Transform outputs based on prediction type
  family <- family(object)
  out <- switch(type,
                "response" = family$linkinv(out),
                "class" = {
                  mu <- family$linkinv(out)
                  1 * (mu >= .5)
                },
                "terms" = out,
                "coefficients" = ,
                "link" = drop(out))

  drop(out)
}

#' @export
model.matrix.islasso <- function(object, ...) {
  # Extract model frame
  data <- model.frame(object, xlev = object$xlevels, ...)

  # Use appropriate method based on context
  if (exists(".GenericCallEnv", inherits = FALSE)) {
    NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
  } else {
    # Build arguments for model.matrix.default
    dots <- list(...)
    dots$data <- dots$contrasts.arg <- NULL
    do.call("model.matrix.default",
            c(list(object = object,
                   data = data,
                   contrasts.arg = object$contrasts),
              dots))
  }
}

#' @export
vcov.islasso <- function(object, ...) {
  object$internal$vcov
}

#' @export
residuals.islasso <- function(object,
                              type = c("deviance", "pearson", "working", "response", "partial"),
                              ...) {
  # Match type argument
  type <- match.arg(type)

  # Extract model components
  y <- object$y
  r <- object$residuals
  mu <- object$fitted.values
  wts <- object$prior.weights

  # Reconstruct response if needed
  if (is.null(y) && type %in% c("deviance", "pearson", "response")) {
    mu.eta <- object$family$mu.eta
    eta <- object$linear.predictors
    y <- mu + r * mu.eta(eta)
  }

  # Calculate requested residual type
  res <- switch(type,
                deviance = if (object$df.residual > 0) {
                  d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts), 0))
                  ifelse(y > mu, d.res, -d.res)
                } else {
                  rep.int(0, length(mu))
                },
                pearson = (y - mu) * sqrt(wts) / sqrt(object$family$variance(mu)),
                working = r,
                response = y - mu,
                partial = r)

  # Handle NA values
  if (!is.null(object$na.action)) {
    res <- naresid(object$na.action, res)
  }

  # Add terms to partial residuals
  if (type == "partial") {
    res <- res + predict(object, type = "terms")
  }

  res
}

#' @export
deviance.islasso <- function(object, ...) {
  object$deviance
}

#' @export
logLik.islasso <- function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }

  # Determine degrees of freedom based on family
  fam <- family(object)$family
  p <- object$rank
  if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) {
    p <- p + 1
  }

  # Calculate log-likelihood from AIC
  val <- p - object$aic/2

  # Add attributes
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}

#' @export
extractAIC.islasso <- function(fit, scale = 0, k = 2, ...) {
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  aic <- fit$aic
  c(edf, aic + (k - 2) * edf)
}

#' @export
family.islasso <- function(object, ...) {
  object$family
}

#' @export
formula.islasso <- function(x, ...) {
  form <- x$formula
  if (!is.null(form)) {
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
  } else {
    formula(x$terms)
  }
}

#' @export
model.frame.islasso <- function(formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]

  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(stats::glm)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms) %||% parent.frame()
    eval(fcall, env)
  } else {
    formula$model
  }
}

#' @export
nobs.islasso <- function(object, ...) {
  if (!is.null(w <- object$prior.weights)) {
    sum(w != 0)
  } else {
    length(object$residuals)
  }
}

#' @export
weights.islasso <- function(object, type = c("prior", "working"), ...) {
  type <- match.arg(type)

  # Get weights based on type
  res <- if (type == "prior") {
    object$prior.weights
  } else {
    object$weights
  }

  # Handle NA values
  if (is.null(object$na.action)) {
    res
  } else {
    naresid(object$na.action, res)
  }
}

#' @export
variable.names.islasso <- function(object, ...) {
  object$internal$nms
}

#' @export
influence.islasso <- function(model, do.coef = TRUE, ...) {
  # Calculate influence measures
  res <- is.influence(model, do.coef = do.coef, ...)

  # Extract Pearson residuals
  pRes <- na.omit(residuals(model, type = "pearson"))[model$prior.weights != 0]
  pRes <- naresid(model$na.action, pRes)

  # Rename fields for consistency
  names(res)[names(res) == "wt.res"] <- "dev.res"

  # Return combined results
  c(res, list(pear.res = pRes))
}

#' @export
cooks.distance.islasso <- function(model,
                                   infl = influence(model, do.coef = FALSE),
                                   res = infl$pear.res,
                                   dispersion = summary(model)$dispersion,
                                   hat = infl$hat, ...) {
  # Calculate Cook's distance
  p <- model$rank
  res <- (res/(1 - hat))^2 * hat/(dispersion * p)

  # Replace infinite values with NaN
  res[is.infinite(res)] <- NaN
  res
}

#' @export
rstandard.islasso <- function(model,
                              infl = influence(model, do.coef = FALSE),
                              type = c("deviance", "pearson"), ...) {
  # Match type argument
  type <- match.arg(type)

  # Extract residuals based on type
  res <- switch(type,
                pearson = infl$pear.res,
                infl$dev.res)

  # Standardize residuals
  res <- res/sqrt(summary(model)$dispersion * (1 - infl$hat))

  # Replace infinite values with NaN
  res[is.infinite(res)] <- NaN
  res
}

#' @export
rstudent.islasso <- function(model, infl = influence(model, do.coef = FALSE), ...) {
  # Calculate studentized residuals
  r <- infl$dev.res
  r <- sign(r) * sqrt(r^2 + (infl$hat * infl$pear.res^2)/(1 - infl$hat))

  # Replace infinite values with NaN
  r[is.infinite(r)] <- NaN

  # Scale by sigma for non-exponential family models
  if (any(family(model)$family == c("binomial", "poisson"))) {
    r
  } else {
    r/infl$sigma
  }
}
