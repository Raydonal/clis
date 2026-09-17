# Observed or Fisher information matrix for a BIc model

Computes the block-structured information matrix of a fitted zero-or-one
inflated beta regression with variable dispersion. By construction the
matrix is block diagonal between the inflation parameters and the
mean/precision parameters (information orthogonality).

## Usage

``` r
bic_info(object, use_fisher = TRUE, penalty = NULL)
```

## Arguments

- object:

  A fitted `gamlss` model of family `BEZI` or `BEOI`.

- use_fisher:

  Logical; if `TRUE` (default) returns the expected (Fisher)
  information, which is guaranteed positive definite; if `FALSE` returns
  the observed information.

- penalty:

  Optional penalty matrix `S` of dimension `(M+m+p) x (M+m+p)` for
  penalised additive submodels. When supplied, the returned information
  is the penalised information `J + S` (or `I + S`), following the
  semiparametric extension. The penalty must be block diagonal across
  the `gamma`, `beta`, `delta` groups so that separability is preserved;
  see [`bic_penalty()`](bic_penalty.md).

## Value

A list with the information matrix (`info`), its inverse (`info_inv`),
the parameter-block indices (`idx`), the weight sequences (`weights`),
and (if a penalty was supplied) the effective degrees of freedom (`edf`)
and their block split (`edf_blocks`).

## References

Ospina, R. and Ferrari, S. L. P. (2012). A general class of zero-or-one
inflated beta regression models. *Computational Statistics & Data
Analysis*, 56(6), 1609-1623.
