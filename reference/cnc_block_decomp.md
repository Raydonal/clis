# Block decomposition of conformal normal curvature

Decomposes the per-observation CNC scores into a contribution from the
inflation-probability submodel and a contribution from the
conditional-mean/precision submodel, exploiting the information
orthogonality of the BIc model. The two contributions sum exactly to the
total, with no cross terms.

## Usage

``` r
cnc_block_decomp(delta_out, info_out)
```

## Arguments

- delta_out:

  A perturbation-matrix list (e.g. from
  [`delta_caseweights()`](delta_caseweights.md)) with a full `Delta`
  component.

- info_out:

  An information list from [`bic_info()`](bic_info.md) with `info_inv`
  and `idx`.

## Value

A list with the discrete contribution `B_gamma`, the continuous
contribution `B_betadelta`, their sum `B_total`, and the discrete
fraction `ratio_gamma`.
