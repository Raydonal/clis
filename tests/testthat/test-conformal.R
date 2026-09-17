test_that("conformal p-values are in (0, 1] and super-uniform under null", {
  set.seed(1)
  cal  <- rnorm(500)
  test <- rnorm(500)                 # exchangeable with calibration
  p <- conformal_pvalues(cal, test)
  expect_true(all(p > 0 & p <= 1))
  # super-uniformity: mean should be near 0.5, not systematically small
  expect_gt(mean(p), 0.4)
  expect_lt(mean(p), 0.6)
})

test_that("conformal p-values flag clear outliers", {
  set.seed(2)
  cal  <- rnorm(200)
  test <- c(rnorm(8), 6, 7)          # last two far in the tail
  p <- conformal_pvalues(cal, test)
  expect_true(all(p[9:10] < 0.05))   # large scores -> small p-values
})

test_that("smallest possible p-value is 1/(n+1)", {
  cal  <- rnorm(99)
  test <- 1e6                        # larger than all calibration scores
  expect_equal(conformal_pvalues(cal, test), 1 / (99 + 1))
})
