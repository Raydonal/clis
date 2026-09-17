test_that("bic_info accepts a zero penalty and matches unpenalised", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  set.seed(1)
  n <- 120
  x1 <- rnorm(n); v1 <- rnorm(n)
  mu <- 1 / (1 + exp(-(-0.4 + 0.8 * x1))); phi <- exp(1.4 + 0.5 * v1)
  y <- gamlss.dist::rBEZI(n, mu = mu, sigma = phi, nu = rep(0.18, n))
  df <- data.frame(y = y, x1 = x1, v1 = v1)
  fit <- suppressWarnings(gamlss::gamlss(
    y ~ x1, sigma.formula = ~ v1, nu.formula = ~ 1,
    family = gamlss.dist::BEZI, data = df,
    control = gamlss::gamlss.control(trace = FALSE)))

  io0 <- bic_info(fit)
  np  <- nrow(io0$info)
  io1 <- bic_info(fit, penalty = matrix(0, np, np))
  expect_equal(io0$info, io1$info, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("a positive penalty inflates the information and lowers edf", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  set.seed(2)
  n <- 120
  x1 <- rnorm(n); v1 <- rnorm(n)
  mu <- 1 / (1 + exp(-(-0.4 + 0.8 * x1))); phi <- exp(1.4 + 0.5 * v1)
  y <- gamlss.dist::rBEZI(n, mu = mu, sigma = phi, nu = rep(0.18, n))
  df <- data.frame(y = y, x1 = x1, v1 = v1)
  fit <- suppressWarnings(gamlss::gamlss(
    y ~ x1, sigma.formula = ~ v1, nu.formula = ~ 1,
    family = gamlss.dist::BEZI, data = df,
    control = gamlss::gamlss.control(trace = FALSE)))

  np <- nrow(bic_info(fit)$info)
  S  <- diag(2, np)                       # ridge penalty on all coefficients
  io <- bic_info(fit, penalty = S)
  # effective df must be positive and below the nominal parameter count
  expect_true(io$penalised)
  expect_gt(io$edf, 0)
  expect_lt(io$edf, np)
  # edf blocks sum to total edf
  expect_equal(sum(io$edf_blocks), io$edf, tolerance = 1e-8)
})

test_that("penalty cross-block between gamma and (beta,delta) is zeroed", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  set.seed(3)
  n <- 120
  x1 <- rnorm(n); v1 <- rnorm(n)
  mu <- 1 / (1 + exp(-(-0.4 + 0.8 * x1))); phi <- exp(1.4 + 0.5 * v1)
  y <- gamlss.dist::rBEZI(n, mu = mu, sigma = phi, nu = rep(0.18, n))
  df <- data.frame(y = y, x1 = x1, v1 = v1)
  fit <- suppressWarnings(gamlss::gamlss(
    y ~ x1, sigma.formula = ~ v1, nu.formula = ~ 1,
    family = gamlss.dist::BEZI, data = df,
    control = gamlss::gamlss.control(trace = FALSE)))

  io0 <- bic_info(fit)
  np  <- nrow(io0$info)
  ig  <- io0$idx$gamma; ibd <- c(io0$idx$beta, io0$idx$delta)
  Sbad <- matrix(0, np, np)
  Sbad[ig[1], ibd[1]] <- Sbad[ibd[1], ig[1]] <- 5   # illegal cross term
  expect_warning(io <- bic_info(fit, penalty = Sbad), "cross-block")
  # separability preserved: gamma-(beta,delta) cross block still zero
  expect_true(all(abs(io$info[ig, ibd]) < 1e-10))
})

test_that("penalised block decomposition still sums to the total", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  set.seed(4)
  n <- 150
  x1 <- rnorm(n); v1 <- rnorm(n)
  mu <- 1 / (1 + exp(-(-0.4 + 0.8 * x1))); phi <- exp(1.4 + 0.5 * v1)
  y <- gamlss.dist::rBEZI(n, mu = mu, sigma = phi, nu = rep(0.18, n))
  df <- data.frame(y = y, x1 = x1, v1 = v1)
  fit <- suppressWarnings(gamlss::gamlss(
    y ~ x1, sigma.formula = ~ v1, nu.formula = ~ 1,
    family = gamlss.dist::BEZI, data = df,
    control = gamlss::gamlss.control(trace = FALSE)))

  np <- nrow(bic_info(fit)$info)
  io <- bic_info(fit, penalty = diag(1, np))
  dd <- delta_caseweights(fit)
  dec <- cnc_block_decomp(dd, io)
  err <- max(abs(dec$B_gamma + dec$B_betadelta - dec$B_total))
  expect_lt(err, 1e-10)
})
