# clis: Conformal Local Influence Screening for Bounded-Response Regression

Provides scalable, statistically calibrated influence diagnostics for
zero-or-one inflated beta (BIc) regression models with variable
dispersion. The core idea is to use the conformal normal curvature of
Poon and Poon (1999)
[doi:10.1111/1467-9868.00162](https://doi.org/10.1111/1467-9868.00162)
as a non-conformity score within a split-conformal testing procedure,
yielding per-observation conformal p-values whose Benjamini-Hochberg
adjustment controls the false discovery rate at a user-specified level
(Bates and others, 2023)
[doi:10.1214/22-AOS2244](https://doi.org/10.1214/22-AOS2244) . Unlike
classical local influence diagnostics, which rely on visual inspection
of index plots and do not scale beyond a few hundred observations,
'clis' provides a finite-sample error guarantee and runs in linear time
per observation after a single model fit. Methods for four perturbation
schemes, block decomposition of influence into the inflation-probability
and conditional-mean/precision components, penalised additive
(semiparametric) submodels, and a full suite of diagnostic plots are
included.

## See also

Useful links:

- <https://github.com/Raydonal/clis>

- <https://raydonal.github.io/clis/>

- Report bugs at <https://github.com/Raydonal/clis/issues>

## Author

**Maintainer**: Raydonal Ospina <raydonal@de.ufpe.br>
([ORCID](https://orcid.org/0000-0002-9884-9090))

Authors:

- Raydonal Ospina <raydonal@de.ufpe.br>
  ([ORCID](https://orcid.org/0000-0002-9884-9090))
