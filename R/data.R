## ===========================================================================
##  data.R -- documentation for bundled datasets
## ===========================================================================

#' National DTP3 vaccination coverage (2022)
#'
#' Reads country-level proportions of diphtheria-tetanus-pertussis (DTP3)
#' vaccine coverage for 2022, with development covariates, from the CSV
#' shipped in the package's `extdata` directory. The response places a point
#' mass at one (countries achieving complete coverage), making it a natural
#' application of one-inflated beta (BEOI) regression.
#'
#' @return A data frame with one row per country and the columns
#'   \describe{
#'     \item{iso3c}{ISO 3166-1 alpha-3 country code.}
#'     \item{dtp3}{DTP3 coverage proportion in \eqn{[0,1]} (response).}
#'     \item{ln_gdp}{Natural log of GDP per capita (PPP, constant USD).}
#'     \item{urb}{Urbanisation rate (proportion of urban population).}
#'     \item{ln_pop}{Natural log of total population.}
#'     \item{hdi}{Human Development Index in \eqn{[0,1]}.}
#'   }
#'
#' @source WHO/UNICEF Joint Monitoring Programme (coverage); World Bank
#'   (GDP, urbanisation, population); UNDP Human Development Report (HDI).
#'
#' @examples
#' vaccination <- load_vaccination()
#' summary(vaccination$dtp3)
#' mean(vaccination$dtp3 == 1)   # fraction at the boundary
#'
#' @export
load_vaccination <- function() {
  f <- system.file("extdata", "vaccination.csv", package = "clis")
  if (!nzchar(f))
    stop("vaccination.csv not found in the installed package.", call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}
