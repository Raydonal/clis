# Plot the classical CNC diagnostic panels

Plot the classical CNC diagnostic panels

## Usage

``` r
plot_cnc_panels(
  cnc,
  scores,
  decomp = NULL,
  r_plot = 3L,
  label_top = 5L,
  suffix = ""
)
```

## Arguments

- cnc:

  A list from [`cnc_matrix()`](cnc_matrix.md).

- scores:

  A list from [`cnc_scores()`](cnc_scores.md).

- decomp:

  Optional list from [`cnc_block_decomp()`](cnc_block_decomp.md) for a
  third panel.

- r_plot:

  Integer; which aggregate-contribution order to plot.

- label_top:

  Integer; number of points to label.

- suffix:

  Character; appended to panel titles.

## Value

Invisibly returns `scores`.
