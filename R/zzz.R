utils::globalVariables(c(
  "lambda", "value", "opt", "varname",
  "theoretical", "predicted", "varName",
  "Low", "Up", "Estimate", "Empirical",
  "Theoretical"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(cli::col_cyan("
+----------------------------------------------+
|          Welcome to *islasso*                |
|  The Induced Smoothed Lasso for R            |
|  Hypothesis testing in high-dimensional data |
+----------------------------------------------+
"))
}
