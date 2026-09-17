# Simulated envelope for quantile residuals

Simulated envelope for quantile residuals

## Usage

``` r
envelope_bic(object, B = 200L, level = 0.95, seed = NULL)
```

## Arguments

- object:

  A fitted BEZI/BEOI `gamlss` model.

- B:

  Integer; number of simulated samples.

- level:

  Numeric; envelope coverage (default 0.95).

- seed:

  Optional integer seed.

## Value

Invisibly returns the observed residuals and envelope bands.
