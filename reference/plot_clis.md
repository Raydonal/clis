# Plot a CLIS screening result

Produces a two-panel diagnostic plot for a
[`clis_screen()`](clis_screen.md) result: a plot of conformal p-values
(with the Benjamini-Hochberg rejection boundary) and a plot of the
non-conformity scores with the declared influential observations
highlighted.

## Usage

``` r
plot_clis(x, label_top = 8L, ...)
```

## Arguments

- x:

  A `clis` object from [`clis_screen()`](clis_screen.md).

- label_top:

  Integer; how many top influential points to label.

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.
