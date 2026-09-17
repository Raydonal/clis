# Conformal p-values from non-conformity scores

Given calibration scores (assumed to come from non-influential, "clean"
observations) and test scores, computes marginal conformal p-values in
the sense of Bates and others (2023). Larger scores indicate stronger
evidence of influence, so the p-value for a test point is the calibrated
rank of its score among the calibration scores.

## Usage

``` r
conformal_pvalues(cal_scores, test_scores)
```

## Arguments

- cal_scores:

  Numeric vector of calibration non-conformity scores.

- test_scores:

  Numeric vector of test non-conformity scores.

## Value

A numeric vector of conformal p-values, one per test score.

## Details

For a test score \\s\\ and calibration scores \\c_1, \ldots, c_n\\, the
conformal p-value is \$\$p = \frac{1 + \\\\i : c_i \ge s\\}{n + 1}.\$\$
These p-values are marginally valid (super-uniform under the null that
the test point is exchangeable with the calibration set) and, by the
positive-dependence result of Bates and others (2023), permit
Benjamini-Hochberg FDR control.

## References

Bates, S., Candes, E., Lei, L., Romano, Y. and Sesia, M. (2023). Testing
for outliers with conformal p-values. *The Annals of Statistics*, 51(1),
149-178.

## Examples

``` r
set.seed(1)
cal  <- rnorm(100)
test <- c(rnorm(8), 4, 5)   # last two are outliers
conformal_pvalues(cal, test)
#>  [1] 0.80198020 0.52475248 0.88118812 0.49504950 0.82178218 0.03960396
#>  [7] 0.23762376 0.17821782 0.00990099 0.00990099
```
