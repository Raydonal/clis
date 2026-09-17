## ===========================================================================
##  clis-package.R -- package-level documentation and namespace imports.
##  The imports live here so that roxygen2 regenerates NAMESPACE completely;
##  before this file they were maintained by hand and were silently lost on
##  the first `devtools::document()` run.
## ===========================================================================

#' @keywords internal
#' @importFrom grDevices adjustcolor
#' @importFrom utils read.csv
#' @importFrom graphics abline barplot legend par points text lines mtext plot
#' @importFrom stats runif sd qnorm pnorm dnorm quantile median p.adjust
#' @importFrom stats residuals vcov AIC BIC as.formula complete.cases
#' @importFrom stats rbinom rnorm rpois qqnorm qqline model.frame model.matrix
#' @importFrom gamlss.dist pBEZI pBEOI
"_PACKAGE"
