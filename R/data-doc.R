#' Prostate Cancer Data
#'
#' This dataset originates from a study examining the correlation between prostate-specific antigen levels and various clinical measures in men scheduled for radical prostatectomy.
#' It contains 97 rows and 9 variables.
#'
#' @format A data frame with 97 observations and 9 variables:
#' \describe{
#'   \item{lcavol}{Log of cancer volume}
#'   \item{lweight}{Log of prostate weight}
#'   \item{age}{Age of the patient}
#'   \item{lbph}{Log of benign prostatic hyperplasia amount}
#'   \item{svi}{Seminal vesicle invasion (binary)}
#'   \item{lcp}{Log of capsular penetration}
#'   \item{gleason}{Gleason score}
#'   \item{pgg45}{Percentage of Gleason scores 4 or 5}
#'   \item{lpsa}{Log of prostate-specific antigen}
#' }
#'
#' @source Stamey, T.A., et al. (1989).
#' Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate: II. radical prostatectomy treated patients.
#' Journal of Urology, 141(5), 1076-1083.
#'
#' @references
#' Stamey, T.A., Kabalin, J.N., McNeal, J.E., Johnstone, I.M., Freiha, F., Redwine, E.A., and Yang, N. (1989).
#' Journal of Urology, 141(5), 1076-1083.
#'
#' @examples
#' data(Prostate)
#' summary(Prostate)
#' cor(Prostate$lpsa, Prostate$lcavol)
#' \dontrun{
#'   fit <- islasso(lpsa ~ ., data = Prostate, family = gaussian())
#'   summary(fit, pval = 0.05)
#'   lambda.aic <- aic.islasso(fit, method = "AIC")
#'   fit.aic <- update(fit, lambda = lambda.aic)
#'   summary(fit.aic, pval = 0.05)
#' }
#'
#' @keywords datasets
#' @name Prostate
#' @docType data
NULL


#' Breast Cancer microarray experiment
#'
#' This data set details a microarray experiment for 52 breast cancer patients. The binary variable \code{status} indicates whether or not the patient died of breast cancer (\code{status = 0}: did not die, \code{status = 1}: died). The other variables represent amplification or deletion of specific genes.
#'
#' Unlike gene expression studies, this experiment focuses on measuring gene amplification or deletion-the number of DNA copies for a given genomic sequence. The goal is to identify key genomic markers distinguishing aggressive from non-aggressive breast cancer.
#'
#' The experiment was conducted by Dr. John Bartlett and Dr. Caroline Witton in the Division of Cancer Sciences and Molecular Pathology at the University of Glasgow's Royal Infirmary.
#'
#' @format A data frame with 52 rows and multiple variables, including a binary \code{status} and gene-level measurements.
#'
#' @source Dr. John Bartlett and Dr. Caroline Witton, Division of Cancer Sciences and Molecular Pathology, University of Glasgow, Glasgow Royal Infirmary.
#'
#' @references
#' Augugliaro L., Mineo A.M. and Wit E.C. (2013). \emph{dgLARS: a differential geometric approach to sparse generalized linear models}, Journal of the Royal Statistical Society. Series B, Vol 75(3), 471-498.
#' Wit E.C. and McClure J. (2004). \emph{Statistics for Microarrays: Design, Analysis and Inference}, Chichester: Wiley.
#'
#' @examples
#' data(breast)
#' str(breast)
#' table(breast$status)
#'
#' \dontrun{
#'   fit <- islasso.path(status ~ ., data = breast, family = binomial(),
#'                       alpha = 0, control = is.control(trace = 2L))
#'   temp <- GoF.islasso.path(fit)
#'   lambda.aic <- temp$lambda.min["AIC"]
#'   fit.aic <- islasso(status ~ ., data = breast, family = binomial(),
#'                      alpha = 0, lambda = lambda.aic)
#'   summary(fit.aic, pval = 0.05)
#' }
#'
#' @keywords datasets
#' @name breast
#' @docType data
NULL

#' Blood and other measurements in diabetics
#'
#' The \code{diabetes} data frame contains 442 observations used in the Efron et al. "Least Angle Regression" paper.
#'
#' @details
#' The \code{x} matrix has been standardized to have unit L2 norm and zero mean in each column.
#' The \code{x2} matrix extends \code{x} by adding selected interaction terms.
#'
#' @format A data frame with 442 rows and 3 columns:
#' \describe{
#'   \item{x}{Matrix with 10 numeric columns (standardized)}
#'   \item{y}{Numeric response vector}
#'   \item{x2}{Matrix with 64 columns including interactions}
#' }
#'
#' @source \url{https://web.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.ps}
#'
#' @references
#' Efron, Hastie, Johnstone and Tibshirani (2003). "Least Angle Regression" (with discussion), \emph{Annals of Statistics}.
#'
#' @examples
#' data(diabetes)
#' str(diabetes)
#' summary(diabetes$y)
#'
#' \dontrun{
#'   fit <- islasso(y ~ ., data = data.frame(y = diabetes$y, diabetes$x2),
#'                  family = gaussian())
#'   summary(fit, pval = 0.05)
#'   lambda.aic <- aic.islasso(fit, interval = c(1, 100))
#'   fit.aic <- update(fit, lambda = lambda.aic)
#'   summary(fit.aic, pval = 0.05)
#' }
#'
#' @keywords datasets
#' @name diabetes
#' @docType data
NULL
