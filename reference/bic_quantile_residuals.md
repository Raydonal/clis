# Randomized quantile residuals for BEZI/BEOI models

Computes randomized quantile residuals for a fitted zero- or
one-inflated beta model. For the continuous observations the beta CDF is
scaled by the non-inflation probability; for the boundary observations
the residual is randomized uniformly over the probability atom, which is
the standard construction for a mixed discrete-continuous distribution.
This is provided because the generic
[`residuals()`](https://rdrr.io/r/stats/residuals.html) method does not
always handle the inflation atom correctly, which can make every
residual share the same sign.

## Usage

``` r
bic_quantile_residuals(object, seed = NULL)
```

## Arguments

- object:

  A fitted `gamlss` model of family `BEZI` or `BEOI`.

- seed:

  Optional integer seed for the randomization at the boundary.

## Value

A numeric vector of randomized quantile residuals.
