# Conformal normal curvature matrix

Computes the Frobenius-normalised conformal normal curvature (CNC)
matrix of Poon and Poon (1999) for a given perturbation matrix and
information matrix inverse.

## Usage

``` r
cnc_matrix(Delta, info_inv)
```

## Arguments

- Delta:

  A perturbation matrix of dimension `(M+m+p) x n`, e.g. from
  [`delta_caseweights()`](delta_caseweights.md).

- info_inv:

  The inverse information matrix from [`bic_info()`](bic_info.md).

## Value

A list with the `n x n` CNC matrix `B`, its eigenvalues `lambda`
(normalised so that the sum of squares is one), eigenvectors `vectors`,
and the Frobenius norm `normF` of the unnormalised curvature.

## References

Poon, W.-Y. and Poon, Y. S. (1999). Conformal normal curvature and
assessment of local influence. *Journal of the Royal Statistical
Society: Series B*, 61(1), 51-61.
