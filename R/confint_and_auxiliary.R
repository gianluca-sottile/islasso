#' confint method for islasso objects
#'
#' @name confint.islasso
#' @aliases confint.islasso print.confint.islasso plot.confint.islasso
#' @title confint method for \code{islasso} objects
#'
#' @param object A fitted model object of class \code{"islasso"}.
#' @param parm A specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level The confidence level required.
#' @param type.ci Character. Only Wald-type confidence intervals are implemented yet! Set \code{type.ci = "wald"} to use estimates and standard errors to build the confidence interval.
#' @param trace Logical. If \code{TRUE} (default), a bar shows the iterations status.
#' @param ... Additional arguments for methods.
#'
#' @description
#' Computes confidence intervals for \code{islasso} objects using a Wald-type approach.
#'
#' @author
#' Maintainer: Gianluca Sottile <gianluca.sottile@unipa.it>
#'
#' @seealso
#' \code{\link{islasso.fit}}, \code{\link{summary.islasso}}, \code{\link{residuals.islasso}}, \code{\link{logLik.islasso}}, \code{\link{predict.islasso}}, \code{\link{deviance.islasso}}
#'
#' @examples
#' n <- 100; p <- 100; p1 <- 10
#' beta.veri <- sort(round(c(seq(0.5, 3, length.out = p1 / 2),
#'                           seq(-1, -2, length.out = p1 / 2)), 2))
#' beta <- c(beta.veri, rep(0, p - p1))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' o <- islasso(y ~ ., data = sim$data, family = gaussian())
#'
#' ci <- confint(o, type.ci = "wald", parm = 1:11)
#' ci
#' plot(ci)
#'
#' @export
confint.islasso <- function(object, parm, level = 0.95, type.ci = "wald", trace = TRUE, ...) {
  type.ci <- match.arg(type.ci)

  # Input validation
  if (!inherits(object, "islasso"))
    stop("'object' must be of class 'islasso'")

  if (level <= 0 || level >= 1)
    stop("'level' must be between 0 and 1")

  # Parameter extraction
  pnames <- names(B0 <- coef(object))
  if (missing(parm))
    parm <- seq_along(pnames)
  else if (is.character(parm))
    parm <- match(parm, pnames, nomatch = 0L)

  # Wald-type confidence intervals
  if(type.ci == "wald") {
    alpha <- (1 - level)/2
    a <- c(alpha, 1 - alpha)
    fac <- qnorm(a)
    pct <- paste(round(100 * a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(pnames[parm], pct))
    ses <- object$se
    ci[] <- B0[parm] + ses[parm] %o% fac
  }

  attr(ci, "coefficients") <- B0[parm]
  attr(ci, "level") <- level
  class(ci) <- "confint.islasso"
  ci
}

#' @export
print.confint.islasso <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.default(format(round(x, 6), digits = digits), print.gap = 2, quote = FALSE)
  invisible(x)
}

#' @export
plot.confint.islasso <- function(x, parm, ...) {
  # Extract coefficients
  beta <- attr(x, "coefficients")

  # Handle parm parameter
  if(missing(parm)) parm <- seq_along(beta)
  nms <- names(beta)

  # Default plot parameters with sensible defaults
  defaults <- list(
    ylab = "Confidence intervals",
    xlab = "",
    srt = 45,
    lwd = 1.1,
    col = "grey",
    cex = 2.5,
    text.cex = 12,
    pch = 16
  )

  # Merge user-provided parameters with defaults
  plot_params <- list(...)
  for (param in names(defaults)) {
    if (is.null(plot_params[[param]])) {
      plot_params[[param]] <- defaults[[param]]
    }
  }

  # Calculate confidence level for y-axis label
  conf_level <- attr(x, "level")
  if (!is.null(conf_level)) {
    plot_params$ylab <- sprintf("Estimates and %d%% Confidence Intervals", round(100 * conf_level))
  }

  # Create data frame for plotting in one step
  plot_data <- data.frame(
    "Estimate" = beta,
    "Low" = unname(x[, 1]),
    "Up" = x[, 2],
    "varName" = factor(nms, levels = nms)
  )[parm, , drop = FALSE]

  # Create plot with ggplot2
  p <- ggplot(plot_data) +
    geom_segment(
      aes(x = varName, y = Low, xend = varName, yend = Up),
      colour = plot_params$col,
      linewidth = plot_params$lwd
    ) +
    geom_point(
      aes(x = varName, y = Estimate),
      size = plot_params$cex,
      shape = plot_params$pch
    ) +
    geom_hline(
      yintercept = 0.0,
      colour = plot_params$col,
      linewidth = .75 * plot_params$lwd,
      linetype = "dashed"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = plot_params$text.cex * 1.2),
      axis.text.y = element_text(size = plot_params$text.cex),
      axis.text.x = element_text(
        size = plot_params$text.cex * .8,
        angle = plot_params$srt,
        hjust = 1
      )
    ) +
    xlab(plot_params$xlab) +
    ylab(plot_params$ylab)

  return(p)
}

makeHyp <- function (cnames, hypothesis, rhs = NULL) {
  parseTerms <- function(terms) {
    component <- gsub("^[-\\ 0-9\\.]+", "", terms)
    component <- gsub(" ", "", component, fixed = TRUE)
    component
  }
  stripchars <- function(x) {
    x <- gsub("\\n", " ", x)
    x <- gsub("\\t", " ", x)
    x <- gsub(" ", "", x, fixed = TRUE)
    x <- gsub("*", "", x, fixed = TRUE)
    x <- gsub("-", "+-", x, fixed = TRUE)
    x <- strsplit(x, "+", fixed = TRUE)[[1]]
    x <- x[x != ""]
    x
  }
  char2num <- function(x) {
    x[x == ""] <- "1"
    x[x == "-"] <- "-1"
    as.numeric(x)
  }
  constants <- function(x, y) {
    with.coef <- unique(unlist(sapply(y, function(z) which(z == parseTerms(x)))))
    if (length(with.coef) > 0) x <- x[-with.coef]
    x <- if (is.null(x)) 0 else sum(as.numeric(x))
    if (any(is.na(x))) stop("The hypothesis \"", hypothesis, "\" is not well formed: contains bad coefficient/variable names.")
    x
  }
  coefvector <- function(x, y) {
    rv <- gsub(" ", "", x, fixed = TRUE) == parseTerms(y)
    if (!any(rv)) return(0)
    if (sum(rv) > 1) stop("The hypothesis \"", hypothesis, "\" is not well formed.")
    rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed = TRUE))))
    if (is.na(rv)) stop("The hypothesis \"", hypothesis, "\" is not well formed: contains non-numeric coefficients.")
    rv
  }
  if (!is.null(rhs)) rhs <- rep(rhs, length.out = length(hypothesis))
  if (length(hypothesis) > 1)
    return(rbind(Recall(cnames, hypothesis[1], rhs[1]), Recall(cnames, hypothesis[-1], rhs[-1])))
  cnames_symb <- sapply(c("@", "#", "~"), function(x) length(grep(x, cnames)) < 1)
  if (any(cnames_symb)) {
    cnames_symb <- head(c("@", "#", "~")[cnames_symb], 1)
    cnames_symb <- paste(cnames_symb, seq_along(cnames), cnames_symb, sep = "")
    hypothesis_symb <- hypothesis
    for (i in order(nchar(cnames), decreasing = TRUE))
      hypothesis_symb <- gsub(cnames[i], cnames_symb[i], hypothesis_symb, fixed = TRUE)
  }
  else {
    stop("The hypothesis \"", hypothesis, "\" is not well formed: contains non-standard coefficient names.")
  }
  lhs <- strsplit(hypothesis_symb, "=", fixed = TRUE)[[1]]
  if (is.null(rhs)) {
    if (length(lhs) < 2)
      rhs <- "0"
    else if (length(lhs) == 2) {
      rhs <- lhs[2]
      lhs <- lhs[1]
    }
    else stop("The hypothesis \"", hypothesis, "\" is not well formed: contains more than one = sign.")
  }
  else {
    if (length(lhs) < 2)
      as.character(rhs)
    else stop("The hypothesis \"", hypothesis, "\" is not well formed: contains a = sign although rhs was specified.")
  }
  lhs <- stripchars(lhs)
  rhs <- stripchars(rhs)
  rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb, coefvector, y = rhs)
  rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs, cnames_symb))
  names(rval) <- c(cnames, "*rhs*")
  rval
}

printHyp <- function (L, b, nms) {
  nomi <- apply(L, 1, function(l){
    id <- which(l != 0)
    temp <- paste(l[id], nms[id], sep = "*", collapse = " + ")
    if(nchar(temp) > min(getOption("width"), 50)) temp <- paste(strtrim(temp, min(getOption("width"), 50)), "...")
    temp
  })
  paste0(nomi, " = ", b)
}

#' General Linear Hypotheses for islasso Models
#'
#' Tests general linear hypotheses and computes confidence intervals for linear combinations of coefficients
#' from a fitted \code{\link{islasso}} model.
#'
#' @aliases anova.islasso print.anova.islasso
#'
#' @param object A fitted model object of class \code{"islasso"}.
#' @param A Hypothesis specification. Either:
#'   \itemize{
#'     \item A numeric matrix or vector with each row specifying a linear combination of coefficients,
#'     \item Or a character vector with symbolic expressions (e.g. \code{"X1 + X2 = 3"}).
#'   }
#' @param b Right-hand side vector for the null hypotheses \code{A \%*\% beta = b}. If omitted, defaults to zeros.
#' @param ci Optional 2-column matrix of confidence intervals for coefficients.
#' @param ... Currently unused.
#'
#' @details
#' The method tests the null hypothesis \eqn{H_0: A \beta = b}, where \eqn{A} and \eqn{b} define a linear constraint on model coefficients.
#'
#' Symbolic expressions support natural syntax: coefficients may be added/subtracted, constants may be multiplied (e.g. \code{"2 * X1 + 3 * X2 = 7"}).
#' Equations with omitted \code{=} assume zero on the right-hand side. See examples for syntax flexibility.
#'
#' @return An object of class \code{"anova.islasso"} containing:
#' \item{Estimate}{Linear combination estimates}
#' \item{SE}{Standard errors}
#' \item{Wald}{Wald statistics}
#' \item{p-value}{Associated p-values}
#'
#' @author Gianluca Sottile \email{gianluca.sottile@unipa.it}
#'
#' @seealso \code{\link{islasso}}, \code{\link{summary.islasso}}, \code{\link{confint.islasso}},
#'          \code{\link{predict.islasso}}, \code{\link{logLik.islasso}}, \code{\link{residuals.islasso}}
#'
#' @examples
#' n <- 100; p <- 100
#' beta <- c(runif(10, -2, 2), rep(0, p - 10))
#' sim <- simulXy(n = n, p = p, beta = beta, seed = 1, family = gaussian())
#' fit <- islasso(y ~ . -1, data = sim$data, family = gaussian())
#'
#' # Test if first 5 variables sum to -7.5
#' anova(fit, A = c("X1 + X2 + X3 + X4 + X5 = -7.5"))
#'
#' # Test multiple hypotheses
#' anova(fit, A = c("X1 + X2 + X3 + X4 + X5", "X6 + X7 + X8 + X9 + X10"), b = c(-7.5, 8.75))
#'
#' # Full diagonal comparison to true coefficients
#' anova(fit, A = diag(p), b = beta)
#'
#' @export
anova.islasso <- function(object, A, b = NULL, ci, ...) {
  # Extract coefficients and variance-covariance matrix
  beta <- coef(object)
  V <- vcov(object)
  nms <- names(beta)

  # Check for aliased coefficients
  if (any(aliased <- is.na(beta)))
    stop("There are aliased coefficients in the model")

  beta <- beta[!aliased]
  if (is.null(beta))
    stop(paste("There is no coef() method for models of class",
               paste(class(object), collapse = ", ")))

  # Handle hypotheses specification
  if (is.character(A)) {
    L <- makeHyp(nms, A, b)
    if (is.null(dim(L))) L <- t(L)
    b <- L[, NCOL(L)]
    L <- L[, -NCOL(L), drop = FALSE]
    rownames(L) <- A
  } else {
    L <- if (is.null(dim(A))) t(A) else A
    if (is.null(b)) b <- rep(0, nrow(L))
  }

  q <- NROW(L)

  # Calculate hypothesis values and variance-covariance matrix
  value.hyp <- L %*% beta - b
  vcov.hyp <- L %*% V %*% t(L)

  # Create result matrix
  rval <- matrix(NA, nrow = q + 1L, ncol = 4L)
  colnames(rval) <- c("Estimate", "Std. Error", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- c(printHyp(L, b, nms), "Overall")

  # Fill in individual hypothesis results
  rval[1:q, 1L] <- value.hyp
  rval[1:q, 2L] <- sqrt(diag(vcov.hyp))
  rval[1:q, 3L] <- value.hyp^2 / diag(vcov.hyp)
  rval[1:q, 4L] <- pchisq(rval[1:q, 3L], 1L, lower.tail = FALSE)

  # Calculate overall test statistic if multiple hypotheses
  if (q > 1) {
    # Use solve with tolerance parameter for better numerical stability
    statistic2 <- as.vector(t(value.hyp) %*% solve(vcov.hyp, tol = 1e-20) %*% value.hyp)
    pval2 <- pchisq(statistic2, q, lower.tail = FALSE)
    rval[q + 1L, 3:4] <- c(statistic2, pval2)
  } else {
    rval <- rval[-c(q + 1L), , drop = FALSE]
  }

  # Handle confidence intervals if provided
  ic <- NULL
  if (!missing(ci)) {
    # Use lapply instead of apply for better performance
    ic_list <- lapply(1:nrow(L), function(i) cislasso(object, a = L[i,], ci = ci))
    ic <- do.call(rbind, ic_list)
    colnames(ic) <- colnames(ci)
    rownames(ic) <- rownames(rval)[1:q]
  }

  # Set results and class
  object$anova <- list(result = as.data.frame(rval), ci = ic)
  class(object) <- c("anova.islasso", class(object))

  return(object)
}

#' @export
print.anova.islasso <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  # Print header
  cat("\n\t", "Simultaneous Tests for General Linear Combination\n\n")

  # Print model call if available
  call <- x$call
  if (!is.null(call)) {
    cat("Fit: ")
    print(call)
    cat("\n")
  }

  # Extract anova results
  pq <- x$anova

  # Print linear hypothesis test results
  cat("Linear Hypothesis:\n")
  printCoefmat(
    pq$result,
    digits = digits,
    has.Pvalue = TRUE,
    P.values = TRUE,
    eps.Pvalue = .Machine$double.eps,
    na.print = ""
  )
  cat("\n")

  # Print confidence intervals if available
  if (!is.null(pq$ci)) {
    cat("Confidence Interval Estimation for Linear Hypothesis:\n")
    # Format confidence intervals with consistent precision
    formatted_ci <- format(round(pq$ci, 6), digits = digits)
    print.default(formatted_ci, print.gap = 2, quote = FALSE)
    cat("\n")
  }

  # Return invisibly
  invisible(x)
}

cislasso <- function(object, a, ci) {
  # Check for required arguments
  if(missing(a)) stop("Vector of linear combination is missing")
  if(missing(ci)) stop("Confidence intervals are missing")

  # Extract necessary components
  V <- vcov(object)
  est <- coef(object)

  # Calculate the linear combination estimate
  th <- sum(a * est)

  # Extract confidence interval bounds
  low <- ci[, 1]
  up <- ci[, 2]

  # Calculate differences for lower and upper bounds
  difmin <- a * est - pmin(a * low, a * up)
  difmax <- a * est - pmax(a * low, a * up)

  # Initialize correlation terms
  Cl <- Cu <- 0

  # Calculate correlation matrix and extract unique correlations
  R <- cov2cor(V)
  r <- R[col(R) > row(R)]

  # Calculate lower bound correlation contribution
  DMIN <- outer(difmin, difmin)
  dmin <- DMIN[col(R) > row(R)]
  Cl <- 2 * sum(r * dmin)

  # Calculate upper bound correlation contribution
  DMAX <- outer(difmax, difmax)
  dmax <- DMAX[col(R) > row(R)]
  Cu <- 2 * sum(r * dmax)

  # Calculate confidence interval bounds
  L <- th - sqrt(sum(difmin^2) + Cl)
  U <- th + sqrt(sum(difmax^2) + Cu)

  # Return confidence interval as a vector
  return(c(L, U))
}

ci.fitted.islasso <- function(object, newx, ci = NULL, type.ci = "wald",
                              conf.level = .95, only.ci = FALSE) {
  # Match argument for CI type
  type.ci <- match.arg(type.ci)

  # Use model matrix if newx is missing
  if(missing(newx)) newx <- model.matrix(object)

  # Calculate confidence intervals if not provided
  if(is.null(ci)) {
    ci <- confint.islasso(object, level = conf.level, type.ci = type.ci, trace = FALSE)
  }

  # Return CI object if only.ci is TRUE
  if(only.ci) return(ci)

  # Get number of observations
  n <- nrow(newx)

  # Pre-allocate result matrix for better performance
  ris <- matrix(NA, nrow = n, ncol = 2)
  colnames(ris) <- colnames(ci)
  rownames(ris) <- rownames(newx)

  # Calculate confidence intervals for each observation
  for(i in 1:n) {
    ris[i, ] <- cislasso(object, newx[i, ], ci)
  }

  return(ris)
}

