NOTE (2026-09-18): this file is kept for the author's own reference only.
It is excluded from the built tarball via `.Rbuildignore` and is not sent
to CRAN — per CRAN's own reply, submission comments belong in the comment
box of the CRAN web submission form, not as a file in the package
(https://cran.r-project.org/submit.html). `data-raw/` was also moved to
`inst/data-raw/` (installed at the top level of the installed package) so
it no longer triggers the "non-standard top-level directory" NOTE.

## R CMD check results

0 errors | 0 warnings | 1 note

Checked with `R CMD check --as-cran` under R 4.6.1 on Ubuntu. The one NOTE
is the usual "New submission" one. CRAN's own win-builder pretest has
additionally and consistently flagged "possibly misspelled words in
DESCRIPTION" for BIc (the model class defined in the Description text
itself) and for the surnames Poon, Benjamini and Hochberg (cited alongside
their publications); these are not misspellings and are not fixable
without removing legitimate content.

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
