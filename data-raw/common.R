## ===========================================================================
##  common.R -- shared generator, contamination and helpers for the paper's
##  simulation experiments.  Sourced by the sim-*.R scripts.
## ===========================================================================
suppressPackageStartupMessages({
  library(gamlss); library(gamlss.dist); library(clis); library(parallel)
})
CORES <- max(1L, min(4L, parallel::detectCores() - 1L))

## Data-generating process of equation (eq:simmodel):
##   logit(alpha) = -1.4 + 0.8 z,  logit(mu) = 0.3 + x,  log(phi) = 1.4 + 0.5 v
gen <- function(n) {
  x <- rnorm(n); v <- rnorm(n); z <- rnorm(n)
  y <- rBEZI(n, mu = plogis(0.3 + x), sigma = exp(1.4 + 0.5 * v),
             nu = plogis(-1.4 + 0.8 * z))
  data.frame(y = pmin(y, 1 - 1e-8), x = x, v = v, z = z)   # BEZI support is [0,1)
}

## The single planting mechanism used in every experiment: the mean covariate
## is displaced by Delta standard deviations and the response is redrawn from
## the model with mu replaced by mu_out, so it contradicts the displacement.
contaminate <- function(df, bad, Delta = 5, mu_out = 0.05) {
  df$x[bad] <- df$x[bad] + Delta * sd(df$x)
  df$y[bad] <- pmin(rBEZI(length(bad), mu = mu_out,
                          sigma = exp(1.4 + 0.5 * df$v[bad]),
                          nu = plogis(-1.4 + 0.8 * df$z[bad])), 1 - 1e-8)
  df
}

fitit <- function(df) suppressWarnings(tryCatch({
  f <- gamlss(y ~ x, sigma.formula = ~ v, nu.formula = ~ z, family = BEZI,
              data = df, control = gamlss.control(trace = FALSE, n.cyc = 100))
  if (!isTRUE(f$converged)) NULL else f
}, error = function(e) NULL))

scores_of <- function(f) {
  io <- bic_info(f)
  cnc_scores_linear(delta_caseweights(f)$Delta, io$info_inv)
}
