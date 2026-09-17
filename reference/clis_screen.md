# Conformal Local Influence Screening

Performs scalable, FDR-controlled influence screening for a fitted
zero-or-one inflated beta regression model. The conformal normal
curvature score of each observation is used as a non-conformity score
within a split-conformal procedure: the data are partitioned into a
calibration set (presumed clean) and a screening set; conformal p-values
are computed for the screening set; and Benjamini-Hochberg adjustment
declares an influential subset with false discovery rate controlled at
level `alpha`.

## Usage

``` r
clis_screen(
  object,
  scheme = "caseweights",
  p = 1L,
  alpha = 0.1,
  calib_frac = 0.5,
  calib_idx = NULL,
  use_fisher = TRUE,
  r_max = 4L,
  score = c("B_Et", "m_r"),
  penalised = FALSE,
  penalty = NULL,
  seed = NULL
)
```

## Arguments

- object:

  A fitted `gamlss` model of family `BEZI` or `BEOI`.

- scheme:

  Character; the perturbation scheme. One of `"caseweights"` (default),
  `"disccovar"`, `"meancovar"`, or `"preccovar"`.

- p:

  Integer; covariate index for the covariate-perturbation schemes.

- alpha:

  Numeric in (0, 1); target false discovery rate. Default 0.1.

- calib_frac:

  Numeric in (0, 1); fraction of observations used for conformal
  calibration when `calib_idx` is not supplied. Default 0.5.

- calib_idx:

  Optional integer vector of observation indices to use as the
  calibration set. The conformal false discovery rate guarantee requires
  the calibration set to be (nearly) free of influential points; when a
  trusted clean subset is known, supply it here. If `NULL` (default) a
  random subset of size `floor(calib_frac * n)` is drawn, which is
  appropriate only when influential points are rare.

- use_fisher:

  Logical; use the Fisher information (default `TRUE`).

- r_max:

  Integer; order for aggregate contributions used as the score.

- score:

  Character; which CNC quantity to use as the non-conformity score:
  `"B_Et"` (basic perturbation, default) or `"m_r"` (aggregate
  contribution at order `r_max`).

- penalised:

  Logical; if `TRUE`, use the penalised information `J + S` obtained
  from the smooth terms of an additive `gamlss` fit, implementing the
  semiparametric extension. The penalty is extracted with
  [`bic_penalty()`](bic_penalty.md). Default `FALSE`.

- penalty:

  Optional penalty matrix, of the dimension of the full coefficient
  vector, added to the information before inversion. Use this when the
  smooth basis is supplied explicitly as design columns.

- seed:

  Optional integer for reproducible calibration splitting.

## Value

An object of class `clis` with components:

- influential:

  Integer indices (into the screening set) declared influential at FDR
  level `alpha`.

- influential_global:

  Integer indices into the original data.

- pvalues:

  Conformal p-values for the screening set.

- padj:

  Benjamini-Hochberg adjusted p-values.

- scores:

  Non-conformity scores for all observations.

- calib_idx, screen_idx:

  Calibration and screening indices.

- alpha, scheme:

  Inputs echoed back.

- decomp:

  Block decomposition (case-weights scheme only).

## Details

Classical local influence diagnostics rank observations by a curvature
measure and rely on visual inspection of an index plot against an ad-hoc
threshold. This neither scales to large \\n\\ nor provides any control
of the type-I error rate. CLIS addresses both limitations: the conformal
wrapper supplies a finite-sample FDR guarantee, and the per-observation
score is computed in linear time after a single model fit.

The calibration set is assumed to be predominantly free of influential
observations. Because influential points are typically rare, a random
split satisfies this approximately; for adversarial settings, a robust
pre-filter can be applied before calibration.

## Examples

``` r
# \donttest{
if (requireNamespace("gamlss", quietly = TRUE) &&
    requireNamespace("betareg", quietly = TRUE)) {
  data("ReadingSkills", package = "betareg")
  ReadingSkills$dys <- as.numeric(ReadingSkills$dyslexia == "yes")
  fit <- gamlss::gamlss(accuracy1 ~ dys * iq,
                        sigma.formula = ~ dys, nu.formula = ~ iq,
                        family = gamlss.dist::BEOI, data = ReadingSkills,
                        control = gamlss::gamlss.control(trace = FALSE))
  res <- clis_screen(fit, alpha = 0.1, seed = 1)
  print(res)
}
#> Conformal Local Influence Screening
#>   Scheme:        caseweights
#>   Score:         B_Et
#>   Target FDR:    0.1
#>   Observations:  44 (22 calibration, 22 screening)
#>   Declared influential: 0 (FDR <= 0.1)
# }
```
