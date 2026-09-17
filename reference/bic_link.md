# Link function and its derivatives for BIc submodels

Constructs a list with the link function, its inverse, and the first and
second derivatives of the inverse link, used throughout the influence
computations.

## Usage

``` r
bic_link(link)
```

## Arguments

- link:

  Character string naming the link. One of `"logit"`, `"probit"`,
  `"cloglog"`, or `"log"`.

## Value

A list with components `name`, `linkfun`, `linkinv`, `mu.eta` (first
derivative of the inverse link), and `mu.eta2` (second derivative).

## Examples

``` r
lk <- bic_link("logit")
lk$mu.eta(0.3)   # d mu / d eta at mu = 0.3
#> [1] 0.21
```
