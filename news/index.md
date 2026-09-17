# Changelog

## clis 0.3.6

- `tests/testthat/test-screening.R`: the helper that plants outliers
  drew responses from
  [`rBEZI()`](https://rdrr.io/pkg/gamlss.dist/man/BEZI.html) without
  clamping them strictly below one, so `gamlss` rejected the response
  with “response variable out of range” and two screening tests failed.
  The clamp applied to the simulation scripts in 0.2.3 is now applied
  here too.
- New `data-raw/` directory with the scripts that regenerate every table
  and figure of the accompanying paper, and a runbook in
  `data-raw/README.md`.

## clis 0.3.5

- [`bic_penalty()`](../reference/bic_penalty.md) now stops with an
  explanatory message when the fit contains a
  [`gamlss::pb()`](https://rdrr.io/pkg/gamlss/man/ps.html)-style smooth.
  Such terms keep their basis coefficients outside the model coefficient
  vector, so no penalty block can be recovered and
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) returns
  only the linear part of the term. The previous behaviour was to return
  a zero penalty, which silently turned penalised screening into
  unpenalised screening on a design that omitted the basis. Penalised
  screening is supported when the basis is supplied explicitly as design
  columns and the penalty matrix is passed to
  [`bic_info()`](../reference/bic_info.md).

## clis 0.3.4

- Bug fix: the observed information for the inflation block omitted the
  term in the second derivative of the inverse link. That term vanishes
  in expectation, which justifies dropping it from the Fisher weight but
  not from the observed one. `bic_info(use_fisher = FALSE)` now agrees
  with a numerical Hessian of the log-likelihood; `use_fisher = TRUE`,
  the default, was already correct.
- Bug fix: the three covariate perturbation schemes omitted the direct
  term by which the perturbed design column enters the score for its own
  coefficient. [`delta_disccovar()`](../reference/delta_disccovar.md),
  [`delta_meancovar()`](../reference/delta_meancovar.md) and
  [`delta_preccovar()`](../reference/delta_preccovar.md) now agree with
  numerical mixed partial derivatives.
- Bug fix: the delta block of
  [`delta_preccovar()`](../reference/delta_preccovar.md) had the wrong
  sign.
- The covariate schemes now stop with an informative message when the
  selected column is constant, which happens with the default `p = 1`
  for a design whose first column is the intercept and previously
  returned a matrix of zeros without warning.

## clis 0.3.3

- Bug fix: the benchmark for the aggregate contribution `m[r]` returned
  by [`cnc_scores()`](../reference/cnc_scores.md) was computed as the
  square root of twice the mean of the selected eigenvalues, which is on
  a different scale from `m[r]` itself. The threshold was therefore far
  too large and the rule never flagged an observation. It is now
  `sqrt(2)` times the mean of `m[r]` across observations, as intended.

## clis 0.3.2

- Documentation: regenerated `man/` from the roxygen sources. All 20
  exported functions now have help pages;
  [`clis_screen()`](../reference/clis_screen.md) documentation includes
  the `calib_idx` and `penalised` arguments, which had been missing.
- Packaging: removed the empty `data/` directory and the unused
  `LazyData` field; the bundled vaccination data is loaded from
  `inst/extdata` via
  [`load_vaccination()`](../reference/load_vaccination.md). The README
  example was corrected accordingly.
- Removed the unused `zoib` dependency from Suggests and updated the
  vignette to refer to the reading-accuracy and lung-function
  applications.

## clis 0.3.1

- Fix: added
  [`bic_quantile_residuals()`](../reference/bic_quantile_residuals.md),
  a correct randomized quantile residual for BEZI/BEOI. The generic
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) mishandled the
  inflation atom for one-inflated (BEOI) fits, making every residual
  share the same sign;
  [`plot_residuals()`](../reference/plot_residuals.md) and
  [`envelope_bic()`](../reference/envelope_bic.md) now use the corrected
  version.

## clis 0.3.0

- Scalability: [`clis_screen()`](../reference/clis_screen.md) with the
  default `B_Et` score now uses a linear-time algorithm that never forms
  the n x n curvature matrix. On a 67,000-observation fit the influence
  scores compute in under a second, where the dense n x n construction
  needs ~33 GB and fails. New exported
  [`cnc_scores_linear()`](../reference/cnc_scores_linear.md) exposes
  this path directly, and
  [`cnc_block_decomp()`](../reference/cnc_block_decomp.md) was rewritten
  to be linear-time as well.
- The `m_r` aggregate score still uses the dense eigendecomposition and
  is now guarded by `options(clis.max_dense_n=)` (default 5000) to avoid
  an accidental out-of-memory build on large samples.

## clis 0.2.3

- Fix: `print.clis` and `summary.clis` are now registered as S3 methods,
  so `print(res)` and `summary(res)` dispatch correctly.
- [`clis_screen()`](../reference/clis_screen.md) gains a `calib_idx`
  argument to supply an explicit, trusted calibration set. The conformal
  guarantee requires the calibration set to be (nearly) free of
  influential points; a random split can let outliers leak into
  calibration and suppress power, so when a clean subset is known it
  should be supplied via `calib_idx`.
- Simulation scripts revised: the response is clamped strictly below one
  (the zero-inflated beta admits \[0, 1) only, and draws at exactly one
  made gamlss reject the response), influential points are planted so
  that the response contradicts the covariate (genuine influence rather
  than accommodated leverage), and screening uses a large clean
  calibration set.
- All of the above validated by running the package under R 4.3.3 with
  gamlss on the real AlcoholUse data and on simulated data.

## clis 0.2.2

- Fix: influence computations no longer require the model to be fitted
  with `x = TRUE`. The design matrices are now recovered via the gamlss
  `model.matrix` method (with a manual fallback), so
  [`bic_info()`](../reference/bic_info.md),
  [`clis_screen()`](../reference/clis_screen.md), and the plots work on
  a standard `gamlss` fit. This was the cause of silent per-replication
  failures in the simulation scripts.
- The simulation scripts now surface the first real error when an entire
  cell fails, instead of returning a non-numeric result downstream.

## clis 0.2.1

- New `data-raw/sim-fdr-power.R` reproduces the linear-model false
  discovery rate table and the detection-power curves.
- New `data-raw/sim-classical.R` reproduces the classical sensitivity
  analysis and its figure.
- New `data-raw/README.md` gives a runbook for reproducing every table
  and figure in the paper.
- Worked examples now use the real `AlcoholUse` data from the `zoib`
  package, matching the paper’s application.

## clis 0.2.0

- Semiparametric extension: penalised additive submodels.
- [`bic_info()`](../reference/bic_info.md) gains a `penalty` argument
  for the penalised information `J + S`, and reports effective degrees
  of freedom and their block split.
- New [`bic_penalty()`](../reference/bic_penalty.md) assembles the
  block-diagonal penalty matrix from the smooth terms of a fitted
  additive `gamlss` model.
- [`clis_screen()`](../reference/clis_screen.md) gains a `penalised`
  argument; the false discovery rate guarantee carries over to penalised
  fits.
- New reproducibility script `data-raw/sim-semiparametric.R` for the
  semiparametric simulation (FDR control, curve recovery, and effective
  degrees of freedom under REML and GCV smoothing).
- New [`plot_influence()`](../reference/plot_influence.md) reproduces
  the classical local-influence index plot (curvature vs. observation
  index with cutoff and labels), in the style of the beta-regression
  diagnostics literature.
- New `data-raw/application.R` reproduces the paper’s application on the
  real `AlcoholUse` data (from the `zoib` package): fit, classical index
  plot, conformal screening, and the semiparametric refit.

## clis 0.1.0

- Initial release.
- [`clis_screen()`](../reference/clis_screen.md): conformal local
  influence screening with finite-sample false discovery rate control
  for zero-or-one inflated beta regression.
- Conformal normal curvature scores via
  [`cnc_matrix()`](../reference/cnc_matrix.md) and
  [`cnc_scores()`](../reference/cnc_scores.md).
- Block decomposition of influence into inflation vs. mean/precision
  components via
  [`cnc_block_decomp()`](../reference/cnc_block_decomp.md).
- Four perturbation schemes: case-weights, discrete-covariate,
  mean-covariate, and precision-covariate.
- Diagnostic plots: [`plot_clis()`](../reference/plot_clis.md),
  [`plot_cnc_panels()`](../reference/plot_cnc_panels.md),
  [`plot_residuals()`](../reference/plot_residuals.md),
  [`envelope_bic()`](../reference/envelope_bic.md).
- Bundled `vaccination` dataset (national DTP3 coverage, 2022).
