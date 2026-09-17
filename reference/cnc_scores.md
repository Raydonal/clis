# Per-observation conformal normal curvature scores

Computes the basic-perturbation CNC scores \\B\_{E_t}\\ and the
aggregate contribution measures \\m\[r\]\_t\\ from a fitted CNC matrix.

## Usage

``` r
cnc_scores(cnc, r_max = 4L)
```

## Arguments

- cnc:

  A list returned by [`cnc_matrix()`](cnc_matrix.md).

- r_max:

  Integer; the maximum order of `r`-influential eigenvectors to consider
  for the aggregate contributions.

## Value

A list with per-observation scores `B_Et`, the matrix of aggregate
contributions `m_r` (one column per `r`), the eigenvalue thresholds
`threshold_r`, the aggregate-contribution thresholds `threshold_mt`, the
number of `r`-influential eigenvectors `k_r`, and the cutoff `b2` for
`B_Et`.
