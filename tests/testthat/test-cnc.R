make_test_fit <- function(n = 120, seed = 1) {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  set.seed(seed)
  x1 <- rnorm(n); v1 <- rnorm(n)
  eta <- -0.4 + 0.8 * x1
  xi  <- 1.4 + 0.5 * v1
  mu  <- 1 / (1 + exp(-eta)); phi <- exp(xi)
  al  <- rep(0.18, n)
  y <- gamlss.dist::rBEZI(n, mu = mu, sigma = phi, nu = al)
  df <- data.frame(y = y, x1 = x1, v1 = v1)
  suppressWarnings(gamlss::gamlss(
    y ~ x1, sigma.formula = ~ v1, nu.formula = ~ 1,
    family = gamlss.dist::BEZI, data = df,
    control = gamlss::gamlss.control(trace = FALSE, n.cyc = 80)))
}

test_that("bic_info returns a block-diagonal, symmetric matrix", {
  fit <- make_test_fit()
  io <- bic_info(fit, use_fisher = TRUE)
  J <- io$info
  expect_true(isSymmetric(unname(J), tol = 1e-8))
  ig <- io$idx$gamma; ibd <- c(io$idx$beta, io$idx$delta)
  # gamma block must be orthogonal to (beta, delta) block
  expect_true(all(abs(J[ig, ibd]) < 1e-10))
})

test_that("Fisher information is positive definite", {
  fit <- make_test_fit()
  io <- bic_info(fit, use_fisher = TRUE)
  ev <- eigen(io$info, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("CNC eigenvalues are non-negative and normalised", {
  fit <- make_test_fit()
  io <- bic_info(fit)
  dd <- delta_caseweights(fit)
  cnc <- cnc_matrix(dd$Delta, io$info_inv)
  expect_true(all(cnc$lambda >= -1e-12))
  expect_equal(sum(cnc$lambda^2), 1, tolerance = 1e-8)
})

test_that("block decomposition sums exactly to the total", {
  fit <- make_test_fit()
  io <- bic_info(fit)
  dd <- delta_caseweights(fit)
  dec <- cnc_block_decomp(dd, io)
  err <- max(abs(dec$B_gamma + dec$B_betadelta - dec$B_total))
  expect_lt(err, 1e-10)
})

test_that("cnc_scores returns finite, non-negative scores", {
  fit <- make_test_fit()
  io <- bic_info(fit)
  dd <- delta_caseweights(fit)
  cnc <- cnc_matrix(dd$Delta, io$info_inv)
  sc  <- cnc_scores(cnc, r_max = 4L)
  expect_true(all(is.finite(sc$B_Et)))
  expect_true(all(sc$B_Et >= 0))
  expect_true(all(is.finite(sc$m_r)))
})
