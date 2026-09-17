# clis ![](reference/figures/logo.png)

**clis** provides scalable, statistically calibrated influence
diagnostics for zero-or-one inflated beta (BIc) regression models with
variable dispersion. It turns the conformal normal curvature of an
observation into a non-conformity score, then wraps it in a
split-conformal testing procedure so that the declared influential set
has its **false discovery rate controlled** at a level you choose.

## Why clis?

Classical local influence analysis (Cook, 1986; Poon & Poon, 1999) ranks
observations by curvature and asks you to eyeball an index plot. That
does not scale to large data and gives no error guarantee. **clis**
fixes both:

- **Finite-sample FDR control** via conformal $`p`$-values and
  Benjamini-Hochberg (Bates et al., 2023).
- **Linear time per observation** after a single model fit, so it scales
  to large $`n`$.
- **Block decomposition** that attributes each observation’s influence
  to the inflation-probability submodel or the
  conditional-mean/precision submodel.

## Installation

``` r

install.packages("clis_0.3.6.tar.gz", repos = NULL, type = "source")
```

## Quick start

``` r

library(clis)
library(gamlss)

vaccination <- load_vaccination()

fit <- gamlss(
  dtp3 ~ ln_gdp + urb,
  sigma.formula = ~ ln_gdp + ln_pop,
  nu.formula    = ~ hdi,
  family  = gamlss.dist::BEOI,
  data    = vaccination,
  control = gamlss.control(trace = FALSE)
)

res <- clis_screen(fit, alpha = 0.1, seed = 1)
res
plot_clis(res)
```

## Learn more

See [`vignette("clis-intro")`](articles/clis-intro.md) for a full
walkthrough. The scripts under `data-raw/` reproduce every table and
figure of the accompanying paper; `data-raw/README.md` is the runbook.

## References

- Bates, S., Candès, E., Lei, L., Romano, Y., & Sesia, M. (2023).
  Testing for outliers with conformal p-values. *The Annals of
  Statistics*, 51(1), 149–178.
- Ospina, R., & Ferrari, S. L. P. (2012). A general class of zero-or-one
  inflated beta regression models. *Computational Statistics & Data
  Analysis*, 56(6), 1609–1623.
- Poon, W.-Y., & Poon, Y. S. (1999). Conformal normal curvature and
  assessment of local influence. *JRSS-B*, 61(1), 51–61.
