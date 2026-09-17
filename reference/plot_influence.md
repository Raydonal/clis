# Classical local-influence index plot

Reproduces the standard influence display of the beta-regression
diagnostics literature: the per-observation conformal normal curvature
\\B\_{E_t}\\ (or the aggregate contribution \\m\[r\]\_t\\) plotted
against the observation index, with the reference cutoff drawn and the
most influential observations labelled. This is the plot practitioners
of local influence expect to see; the conformal screening in
[`clis_screen()`](clis_screen.md) and [`plot_clis()`](plot_clis.md)
complements it with an error-controlled declaration.

## Usage

``` r
plot_influence(
  object,
  scheme = "caseweights",
  p = 1L,
  measure = c("B_Et", "m_r"),
  r_plot = 3L,
  use_fisher = TRUE,
  penalised = FALSE,
  label_top = 5L,
  labels = NULL
)
```

## Arguments

- object:

  A fitted BEZI/BEOI `gamlss` model.

- scheme:

  Perturbation scheme: "caseweights" (default), "disccovar",
  "meancovar", or "preccovar".

- p:

  Covariate index for the covariate-perturbation schemes.

- measure:

  Which quantity to plot: "B_Et" (default) or "m_r".

- r_plot:

  Aggregate-contribution order when `measure = "m_r"`.

- use_fisher:

  Logical; use the Fisher information (default `TRUE`).

- penalised:

  Logical; use the penalised information (default `FALSE`).

- label_top:

  Integer; number of most-influential points to label.

- labels:

  Optional character vector of observation labels (for example country
  codes); defaults to the observation index.

## Value

Invisibly, a list with the plotted scores and the cutoff.
