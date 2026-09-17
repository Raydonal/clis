# Extract the penalty matrix from a fitted additive BIc model

Assembles the block-diagonal penalty matrix `S(lambda)` from the smooth
terms of a fitted `gamlss` model whose submodels use penalised additive
terms (for example [`pb()`](https://rdrr.io/pkg/gamlss/man/ps.html)
P-spline terms). The result is suitable as the `penalty` argument of
[`bic_info()`](bic_info.md).

## Usage

``` r
bic_penalty(object)
```

## Arguments

- object:

  A fitted `gamlss` model of family `BEZI` or `BEOI` with penalised
  additive terms in one or more submodels.

## Value

A block-diagonal penalty matrix of dimension `(M+m+p) x (M+m+p)`, block
diagonal across the inflation, mean, and precision coefficient groups.

## Details

The function reads the smoothing structure stored by `gamlss` for each
parameter (`mu`, `sigma`, `nu`) and places the corresponding penalty
contributions on the diagonal blocks. Terms without a penalty contribute
a zero block, so a model with a mix of linear and smooth terms is
handled transparently. By construction the returned matrix is block
diagonal across the three coefficient groups, so it preserves the
separability required by the semiparametric theory.

## See also

[`bic_info()`](bic_info.md)
