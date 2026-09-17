# Linear-time per-observation CNC scores

Computes the basic-perturbation CNC scores \\B\_{E_t}\\ without ever
forming the \\n \times n\\ curvature matrix. The score is the scaled
diagonal of \\F_0 = \Delta^\top \mathcal{I}^{-1}\Delta\\, and both the
diagonal and the Frobenius norm \\\\F_0\\\_F\\ are obtained from
quantities of dimension \\(M+m+p)\\, so the cost is linear in `n` and
the memory footprint is negligible. This is the scalable path used for
large samples, where the dense \\n\times n\\ eigenproblem of
[`cnc_matrix()`](cnc_matrix.md) is infeasible.

## Usage

``` r
cnc_scores_linear(Delta, info_inv)
```

## Arguments

- Delta:

  A perturbation matrix of dimension `(M+m+p) x n`.

- info_inv:

  The inverse information matrix from [`bic_info()`](bic_info.md).

## Value

A list with per-observation scores `B_Et`, the cutoff `b2`, the
Frobenius norm `normF`, and `n`.
