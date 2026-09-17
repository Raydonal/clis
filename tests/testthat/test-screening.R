make_fit_with_outliers <- function(n = 300, n_out = 15, seed = 1) {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  set.seed(seed)
  x1 <- rnorm(n); v1 <- rnorm(n)
  mu  <- 1 / (1 + exp(-(-0.4 + 0.8 * x1)))
  phi <- exp(1.4 + 0.5 * v1)
  ## BEZI has support [0, 1); a draw at exactly 1 makes gamlss reject the
  ## response with "response variable out of range".
  y <- pmin(gamlss.dist::rBEZI(n, mu = mu, sigma = phi, nu = rep(0.18, n)),
            1 - 1e-8)
  out_idx <- sample.int(n, n_out)
  x1[out_idx] <- x1[out_idx] + 4         # shift leverage
  mu2 <- 1 / (1 + exp(-(-0.4 + 0.8 * x1[out_idx])))
  y[out_idx] <- pmin(gamlss.dist::rBEZI(n_out, mu = mu2, sigma = phi[out_idx],
                                        nu = rep(0.18, n_out)), 1 - 1e-8)
  df <- data.frame(y = y, x1 = x1, v1 = v1)
  fit <- suppressWarnings(gamlss::gamlss(
    y ~ x1, sigma.formula = ~ v1, nu.formula = ~ 1,
    family = gamlss.dist::BEZI, data = df,
    control = gamlss::gamlss.control(trace = FALSE, n.cyc = 80)))
  list(fit = fit, out_idx = out_idx)
}

test_that("clis_screen returns a valid clis object", {
  obj <- make_fit_with_outliers()
  res <- clis_screen(obj$fit, alpha = 0.1, seed = 1)
  expect_s3_class(res, "clis")
  expect_true(all(res$pvalues > 0 & res$pvalues <= 1))
  expect_true(all(res$influential_global %in% seq_len(res$n)))
})

test_that("empirical FDR is approximately controlled across replications", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  fdrs <- replicate(20, {
    obj <- make_fit_with_outliers(seed = sample.int(1e6, 1))
    res <- clis_screen(obj$fit, alpha = 0.2, seed = 1)
    decl <- res$influential_global
    if (length(decl) == 0) return(0)
    # false discoveries: declared but not actually outliers
    mean(!(decl %in% obj$out_idx))
  })
  # mean empirical FDR should not greatly exceed the nominal 0.2
  expect_lt(mean(fdrs), 0.35)
})

test_that("higher alpha declares at least as many influential points", {
  obj <- make_fit_with_outliers()
  r1 <- clis_screen(obj$fit, alpha = 0.05, seed = 1)
  r2 <- clis_screen(obj$fit, alpha = 0.30, seed = 1)
  expect_gte(length(r2$influential), length(r1$influential))
})
