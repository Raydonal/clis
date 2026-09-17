# National DTP3 vaccination coverage (2022)

Reads country-level proportions of diphtheria-tetanus-pertussis (DTP3)
vaccine coverage for 2022, with development covariates, from the CSV
shipped in the package's `extdata` directory. The response places a
point mass at one (countries achieving complete coverage), making it a
natural application of one-inflated beta (BEOI) regression.

## Usage

``` r
load_vaccination()
```

## Source

WHO/UNICEF Joint Monitoring Programme (coverage); World Bank (GDP,
urbanisation, population); UNDP Human Development Report (HDI).

## Value

A data frame with one row per country and the columns

- iso3c:

  ISO 3166-1 alpha-3 country code.

- dtp3:

  DTP3 coverage proportion in \\\[0,1\]\\ (response).

- ln_gdp:

  Natural log of GDP per capita (PPP, constant USD).

- urb:

  Urbanisation rate (proportion of urban population).

- ln_pop:

  Natural log of total population.

- hdi:

  Human Development Index in \\\[0,1\]\\.

## Examples

``` r
vaccination <- load_vaccination()
summary(vaccination$dtp3)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.3300  0.8200  0.9100  0.8638  0.9700  1.0000 
mean(vaccination$dtp3 == 1)   # fraction at the boundary
#> [1] 0.09803922
```
