# Reproduction scripts

Every table and figure in *Conformal Local Influence Screening for
Bounded-Response Regression* is produced by the scripts in this directory.
Run them from inside `data-raw/`, with the `clis` package installed:

```r
install.packages(c("gamlss", "gamlss.dist", "gamlss.data", "betareg"))
install.packages("clis", repos = NULL, type = "source")
```

| Script | Produces |
|---|---|
| `application-reading.R` | Tables 6 and 7; Figures 6, 7, 8 (reading accuracy, `betareg::ReadingSkills`) |
| `application-lung.R`    | Table 8; Figures 9, 10 (lung function, `gamlss.data::lungFunction`) |
| `sim-classical.R`       | Table 3 (classical diagnostics) |
| `sim-exchangeability.R` | the numbers quoted in Section 6.3 |
| `sim-fdr-power.R`       | Table 4; Figures 4 and 5 |
| `sim-semiparametric.R`  | Table 5 (penalised screening) |

`common.R` holds the data-generating process, the planting mechanism and the
fitting helpers shared by the simulation scripts. `fakefit.R` builds the
minimal fitted-model object that lets `clis` read a penalised fit whose
coefficients were obtained outside `gamlss`; it is validated against a real
`gamlss` fit (the information matrices agree to `1e-13`).

Figures are written to `../../artigo1_clis/figuras` by default; set the
environment variable `CLIS_FIGDIR` to change that.

## Notes

* The applications are deterministic apart from the randomisation in the
  quantile residuals and the conformal calibration split, both of which are
  seeded. Screening results are reported over a range of splits rather than
  for one seed, because a single split is not a stable summary at these
  sample sizes.
* The simulation scripts are seeded per replication and parallelised over
  `parallel::detectCores() - 1` cores, capped at four. Monte Carlo error is
  roughly `0.004` for the realised false discovery rates at `R = 400`.
* The penalised lung-function and semiparametric fits maximise the penalised
  log-likelihood directly, with an explicit B-spline basis in the design and
  the penalty in the information, because `gamlss::pb()` keeps its basis
  coefficients outside the coefficient vector and `bic_penalty()` therefore
  refuses such fits (see `NEWS.md` for 0.3.5).
