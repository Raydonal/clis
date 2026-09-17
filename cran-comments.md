## R CMD check results

0 errors | 0 warnings | 2 notes

Checked with `R CMD check --as-cran` under R 4.6.1 on Ubuntu.

* The first NOTE is the usual one for a new submission.
* The second reports `data-raw` as a non-standard top-level directory. It is
  intentional: it holds the scripts that reproduce every table and figure of
  the accompanying manuscript, and they are referenced from the paper. They
  are never run at build or check time.

## Test environments

* local: R 4.6.1 on Ubuntu 24.04

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Additional comments

The package implements the methods described in an accompanying manuscript
submitted for peer review. Examples that depend on suggested packages
('betareg', 'gamlss.data') are wrapped in \donttest{} and guarded by
requireNamespace(), so they are skipped when the suggested package is
unavailable. The vignette executes only chunks that use synthetic data or the
'gamlss' and 'gamlss.dist' packages, guarded by requireNamespace().

The scripts under data-raw/ reproduce every table and figure of the
accompanying manuscript; they are not run at build or check time.
