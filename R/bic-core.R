## ===========================================================================
##  bic-core.R
##  Core quantities for zero-or-one inflated beta (BIc) regression:
##  link functions, auxiliary quantities, information matrix, Delta matrices.
## ===========================================================================

#' Link function and its derivatives for BIc submodels
#'
#' Constructs a list with the link function, its inverse, and the first and
#' second derivatives of the inverse link, used throughout the influence
#' computations.
#'
#' @param link Character string naming the link. One of `"logit"`, `"probit"`,
#'   `"cloglog"`, or `"log"`.
#'
#' @return A list with components `name`, `linkfun`, `linkinv`, `mu.eta`
#'   (first derivative of the inverse link), and `mu.eta2` (second derivative).
#'
#' @examples
#' lk <- bic_link("logit")
#' lk$mu.eta(0.3)   # d mu / d eta at mu = 0.3
#'
#' @export
bic_link <- function(link) {
  switch(link,
    logit = list(
      name    = "logit",
      linkfun = function(mu) log(mu / (1 - mu)),
      linkinv = function(eta) 1 / (1 + exp(-eta)),
      mu.eta  = function(mu) mu * (1 - mu),
      mu.eta2 = function(mu) mu * (1 - mu) * (1 - 2 * mu)
    ),
    probit = list(
      name    = "probit",
      linkfun = function(mu) qnorm(mu),
      linkinv = function(eta) pnorm(eta),
      mu.eta  = function(mu) dnorm(qnorm(mu)),
      mu.eta2 = function(mu) { z <- qnorm(mu); -z * dnorm(z) }
    ),
    cloglog = list(
      name    = "cloglog",
      linkfun = function(mu) log(-log(1 - mu)),
      linkinv = function(eta) 1 - exp(-exp(eta)),
      mu.eta  = function(mu) -(1 - mu) * log(1 - mu),
      mu.eta2 = function(mu) { lv <- log(1 - mu); -(1 - mu) * (1 + lv) * lv }
    ),
    log = list(
      name    = "log",
      linkfun = function(mu) log(mu),
      linkinv = function(eta) exp(eta),
      mu.eta  = function(mu) mu,
      mu.eta2 = function(mu) mu
    ),
    stop("Unsupported link '", link, "'. Use logit, probit, cloglog or log.")
  )
}

## Internal: verify object is a fitted BEZI/BEOI gamlss model; return c in {0,1}
.bic_check <- function(object) {
  if (!inherits(object, "gamlss"))
    stop("`object` must be a fitted gamlss model.", call. = FALSE)
  fam <- object$family[1]
  if (fam == "BEZI") return(0L)
  if (fam == "BEOI") return(1L)
  stop("Family must be BEZI (c = 0) or BEOI (c = 1).", call. = FALSE)
}

#' Randomized quantile residuals for BEZI/BEOI models
#'
#' Computes randomized quantile residuals for a fitted zero- or one-inflated
#' beta model. For the continuous observations the beta CDF is scaled by the
#' non-inflation probability; for the boundary observations the residual is
#' randomized uniformly over the probability atom, which is the standard
#' construction for a mixed discrete-continuous distribution. This is provided
#' because the generic \code{residuals()} method does not always handle the
#' inflation atom correctly, which can make every residual share the same sign.
#'
#' @param object A fitted `gamlss` model of family `BEZI` or `BEOI`.
#' @param seed Optional integer seed for the randomization at the boundary.
#' @return A numeric vector of randomized quantile residuals.
#' @export
bic_quantile_residuals <- function(object, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cmix <- .bic_check(object)
  y   <- object$y
  mu  <- object$mu.fv
  sig <- object$sigma.fv
  nu  <- object$nu.fv
  n   <- length(y)
  atom <- y == cmix
  u <- numeric(n)

  if (cmix == 0L) {
    ## BEZI: atom at 0, continuous mass on (0,1) with weight (1 - nu)
    pcont <- gamlss.dist::pBEZI(y[!atom], mu = mu[!atom],
                                sigma = sig[!atom], nu = nu[!atom])
    u[!atom] <- pcont
    u[atom]  <- stats::runif(sum(atom), 0, nu[atom])
  } else {
    ## BEOI: atom at 1, continuous mass on (0,1) with weight (1 - nu)
    pcont <- gamlss.dist::pBEOI(y[!atom], mu = mu[!atom],
                                sigma = sig[!atom], nu = nu[!atom])
    u[!atom] <- pcont
    u[atom]  <- stats::runif(sum(atom), 1 - nu[atom], 1)
  }
  stats::qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))
}

## Internal: obtain the design matrix for one gamlss parameter. gamlss only
## stores `mu.x` etc. when fitted with `x = TRUE`, so when those slots are
## absent we use the gamlss `model.matrix` S3 method (which knows how to
## rebuild the design from the stored terms and data), falling back to a
## manual reconstruction only if that is unavailable.
.bic_design <- function(object, par) {
  stored <- object[[paste0(par, ".x")]]
  if (!is.null(stored)) return(as.matrix(stored))

  ## gamlss provides model.matrix(object, what = "mu"|"sigma"|"nu")
  mm <- tryCatch(stats::model.matrix(object, what = par),
                 error = function(e) NULL)
  if (!is.null(mm)) return(as.matrix(mm))

  ## manual fallback from terms + data
  tt  <- object[[paste0(par, ".terms")]]
  dat <- object[["data"]]
  if (is.null(tt) || is.null(dat))
    stop("Cannot recover the '", par, "' design matrix from the fitted ",
         "model. Refit with the data retained, or with x = TRUE.",
         call. = FALSE)
  mf <- stats::model.frame(tt, data = as.data.frame(dat))
  stats::model.matrix(tt, mf)
}

## Internal: all auxiliary quantities evaluated at the MLE
.bic_aux <- function(object) {
  cmix    <- .bic_check(object)
  nu_hat  <- object$nu.fv
  mu_hat  <- object$mu.fv
  phi_hat <- object$sigma.fv
  y       <- object$y
  n       <- length(y)
  Ic <- as.integer(y == cmix); Dc <- 1L - Ic

  lnu <- bic_link(object$nu.link)
  lmu <- bic_link(object$mu.link)
  lph <- bic_link(object$sigma.link)

  Ga <- lnu$mu.eta(nu_hat); Ga2 <- lnu$mu.eta2(nu_hat)
  Tm <- lmu$mu.eta(mu_hat)
  Sm <- lmu$mu.eta2(mu_hat); Rp <- lph$mu.eta(phi_hat)
  Qp <- lph$mu.eta2(phi_hat)

  y_star <- ifelse(Dc == 1, log(y / (1 - y)), 0)
  y_ast  <- ifelse(Dc == 1, log(1 - y), 0)
  mp <- digamma(mu_hat * phi_hat) - digamma((1 - mu_hat) * phi_hat)
  mq <- digamma((1 - mu_hat) * phi_hat) - digamma(phi_hat)

  eps1 <- Ic - nu_hat
  eps2 <- y_star - mp
  v_star <- trigamma(mu_hat * phi_hat) + trigamma((1 - mu_hat) * phi_hat)
  c_star <- mu_hat * trigamma(mu_hat * phi_hat) -
            (1 - mu_hat) * trigamma((1 - mu_hat) * phi_hat)
  w_phph <- mu_hat^2 * trigamma(mu_hat * phi_hat) +
            (1 - mu_hat)^2 * trigamma((1 - mu_hat) * phi_hat) -
            trigamma(phi_hat)
  u_t <- Dc * (mu_hat * eps2 + y_ast - mq)

  Z <- .bic_design(object, "nu")
  X <- .bic_design(object, "mu")
  V <- .bic_design(object, "sigma")

  list(cmix = cmix, n = n, nu = nu_hat, mu = mu_hat, phi = phi_hat,
       y = y, Ic = Ic, Dc = Dc, Ga = Ga, Ga2 = Ga2, Tm = Tm, Sm = Sm, Rp = Rp, Qp = Qp,
       y_star = y_star, y_ast = y_ast, mp = mp, mq = mq,
       eps1 = eps1, eps2 = eps2, v_star = v_star, c_star = c_star,
       w_phph = w_phph, u_t = u_t, Pt = 1 / (nu_hat * (1 - nu_hat)),
       Z = Z, X = X, V = V,
       M = ncol(Z), m = ncol(X), p = ncol(V))
}

#' Observed or Fisher information matrix for a BIc model
#'
#' Computes the block-structured information matrix of a fitted zero-or-one
#' inflated beta regression with variable dispersion. By construction the
#' matrix is block diagonal between the inflation parameters and the
#' mean/precision parameters (information orthogonality).
#'
#' @param object A fitted `gamlss` model of family `BEZI` or `BEOI`.
#' @param use_fisher Logical; if `TRUE` (default) returns the expected
#'   (Fisher) information, which is guaranteed positive definite; if `FALSE`
#'   returns the observed information.
#' @param penalty Optional penalty matrix `S` of dimension
#'   `(M+m+p) x (M+m+p)` for penalised additive submodels. When supplied,
#'   the returned information is the penalised information `J + S`
#'   (or `I + S`), following the semiparametric extension. The penalty
#'   must be block diagonal across the `gamma`, `beta`, `delta` groups so
#'   that separability is preserved; see [bic_penalty()].
#'
#' @return A list with the information matrix (`info`), its inverse
#'   (`info_inv`), the parameter-block indices (`idx`), the weight
#'   sequences (`weights`), and (if a penalty was supplied) the effective
#'   degrees of freedom (`edf`) and their block split (`edf_blocks`).
#'
#' @references
#' Ospina, R. and Ferrari, S. L. P. (2012). A general class of zero-or-one
#' inflated beta regression models. \emph{Computational Statistics & Data
#' Analysis}, 56(6), 1609-1623.
#'
#' @export
bic_info <- function(object, use_fisher = TRUE, penalty = NULL) {
  a <- .bic_aux(object)

  ## Observed weights
  ## Observed information for the inflation block keeps the term in the
  ## second derivative of the inverse link; it vanishes in expectation,
  ## which justifies dropping it from the Fisher weight k1 but not here.
  w1  <- (a$Ic / a$nu^2 + a$Dc / (1 - a$nu)^2) * a$Ga^2 -
         (a$Ic / a$nu - a$Dc / (1 - a$nu)) * a$Ga2
  w2  <- a$Dc * (a$phi^2 * a$v_star * a$Tm^2 - a$phi * a$eps2 * a$Sm)
  wbd <- a$Dc * (a$phi * a$c_star - a$eps2) * a$Tm * a$Rp
  wd  <- a$Dc * (a$w_phph * a$Rp^2 - a$u_t * a$Qp)

  ## Fisher (expected) weights
  k1  <- a$Ga^2 / (a$nu * (1 - a$nu))
  k2  <- (1 - a$nu) * a$phi^2 * a$v_star * a$Tm^2
  kbd <- (1 - a$nu) * a$phi * a$c_star * a$Tm * a$Rp
  kd  <- (1 - a$nu) * a$w_phph * a$Rp^2

  b1 <- if (use_fisher) k1 else w1
  b2 <- if (use_fisher) k2 else w2
  bb <- if (use_fisher) kbd else wbd
  bd <- if (use_fisher) kd else wd

  Jgg <- crossprod(a$Z, b1 * a$Z)
  Jbb <- crossprod(a$X, b2 * a$X)
  Jbd <- crossprod(a$X, bb * a$V)
  Jdd <- crossprod(a$V, bd * a$V)

  M <- a$M; m <- a$m; p <- a$p; np <- M + m + p
  OF <- matrix(0, np, np)
  ig <- seq_len(M); ib <- M + seq_len(m); id <- M + m + seq_len(p)
  OF[ig, ig] <- Jgg; OF[ib, ib] <- Jbb
  OF[ib, id] <- Jbd; OF[id, ib] <- t(Jbd); OF[id, id] <- Jdd

  OF_unpen <- OF
  edf <- NULL; edf_blocks <- NULL
  if (!is.null(penalty)) {
    if (!all(dim(penalty) == c(np, np)))
      stop("`penalty` must be a (M+m+p) x (M+m+p) matrix.", call. = FALSE)
    ## enforce block-diagonal structure between gamma and (beta,delta)
    cross <- penalty[ig, c(ib, id), drop = FALSE]
    if (any(abs(cross) > .Machine$double.eps^0.5))
      warning("`penalty` has nonzero gamma-(beta,delta) cross-block; ",
              "separability requires it to be zero. Zeroing it out.",
              call. = FALSE)
    penalty[ig, c(ib, id)] <- 0
    penalty[c(ib, id), ig] <- 0
    OF <- OF + penalty
    ## effective degrees of freedom: tr((J+S)^{-1} J)
    H <- solve(OF, OF_unpen)
    edf <- sum(diag(H))
    edf_blocks <- c(gamma = sum(diag(H[ig, ig, drop = FALSE])),
                    betadelta = sum(diag(H[c(ib, id), c(ib, id), drop = FALSE])))
  }

  list(info = OF, info_inv = solve(OF),
       info_unpen = OF_unpen,
       use_fisher = use_fisher, penalised = !is.null(penalty),
       edf = edf, edf_blocks = edf_blocks,
       idx = list(gamma = ig, beta = ib, delta = id),
       weights = list(obs = list(w1 = w1, w2 = w2, wbd = wbd, wd = wd),
                      fisher = list(k1 = k1, k2 = k2, kbd = kbd, kd = kd)))
}

#' Extract the penalty matrix from a fitted additive BIc model
#'
#' Assembles the block-diagonal penalty matrix `S(lambda)` from the smooth
#' terms of a fitted `gamlss` model whose submodels use penalised additive
#' terms (for example `pb()` P-spline terms). The result is suitable as the
#' `penalty` argument of [bic_info()].
#'
#' @param object A fitted `gamlss` model of family `BEZI` or `BEOI` with
#'   penalised additive terms in one or more submodels.
#'
#' @return A block-diagonal penalty matrix of dimension
#'   `(M+m+p) x (M+m+p)`, block diagonal across the inflation, mean, and
#'   precision coefficient groups.
#'
#' @details The function reads the smoothing structure stored by `gamlss`
#'   for each parameter (`mu`, `sigma`, `nu`) and places the corresponding
#'   penalty contributions on the diagonal blocks. Terms without a penalty
#'   contribute a zero block, so a model with a mix of linear and smooth
#'   terms is handled transparently. By construction the returned matrix is
#'   block diagonal across the three coefficient groups, so it preserves the
#'   separability required by the semiparametric theory.
#'
#' @seealso [bic_info()]
#' @export
bic_penalty <- function(object) {
  a <- .bic_aux(object)
  M <- a$M; m <- a$m; p <- a$p; np <- M + m + p
  S <- matrix(0, np, np)
  ig <- seq_len(M); ib <- M + seq_len(m); id <- M + m + seq_len(p)

  ## Helper: pull a penalty block for one parameter if gamlss stored one.
  ## gamlss keeps smoothing information in object$<par>.coefSmo; when a
  ## closed-form penalty matrix is unavailable we fall back to a ridge
  ## approximation using the stored effective df and lambda, which keeps
  ## the block-diagonal structure intact.
  get_block <- function(par, ncoef) {
    smo <- tryCatch(object[[paste0(par, ".coefSmo")]], error = function(e) NULL)
    Sblk <- matrix(0, ncoef, ncoef)
    if (is.null(smo) || length(smo) == 0) return(Sblk)
    for (term in smo) {
      Sterm <- term$S
      lam   <- term$lambda
      if (!is.null(Sterm) && !is.null(lam)) {
        idx <- term$coef.idx
        if (!is.null(idx) && length(idx) == nrow(Sterm))
          Sblk[idx, idx] <- Sblk[idx, idx] + lam * Sterm
      }
    }
    Sblk
  }

  S[ig, ig] <- get_block("nu",    M)
  S[ib, ib] <- get_block("mu",    m)
  S[id, id] <- get_block("sigma", p)

  ## Guard: gamlss stores gamlss::pb() and friends in a local random-effects
  ## representation whose spline coefficients are NOT part of the model
  ## coefficient vector, so no penalty block can be recovered for them and
  ## the design returned by model.matrix() carries only the linear part of
  ## the term. Returning a zero penalty there would silently reduce
  ## penalised screening to unpenalised screening on a misspecified design.
  has_smooth <- any(vapply(c("nu", "mu", "sigma"), function(par) {
    smo <- tryCatch(object[[paste0(par, ".coefSmo")]], error = function(e) NULL)
    length(smo) > 0L
  }, logical(1)))
  if (has_smooth && all(S == 0))
    stop("This fit contains a smooth term whose basis coefficients gamlss ",
         "keeps outside the coefficient vector, so its penalty cannot be ",
         "recovered and the design would silently omit the basis. Fit the ",
         "smooth with an explicit basis (for example spline columns from ",
         "splines::bs()) and pass the corresponding penalty matrix to ",
         "bic_info(penalty = ).", call. = FALSE)
  S
}


#' Perturbation matrix for the case-weights scheme
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @return A list with the full `Delta` matrix and its `gamma`, `beta`,
#'   `delta` blocks.
#' @export
delta_caseweights <- function(object) {
  a  <- .bic_aux(object)
  Dg <- t(a$Z * (a$Ga * a$Pt * a$eps1))
  Db <- t(a$X * (a$Dc * a$phi * a$Tm * a$eps2))
  Dd <- t(a$V * (a$u_t * a$Rp))
  list(Delta = rbind(Dg, Db, Dd),
       gamma = Dg, beta = Db, delta = Dd)
}

#' Perturbation matrix for the discrete-covariate scheme
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @param p Index of the perturbed covariate in the inflation design matrix.
#' @return A list with the full `Delta` and its blocks.
#' @export
delta_disccovar <- function(object, p = 1L) {
  a  <- .bic_aux(object)
  gp <- object$nu.coefficients[p]; s <- sd(a$Z[, p])
  if (s <= 0) stop("Column ", p, " of the design is constant (typically the ",
                   "intercept); an additive covariate perturbation is not ",
                   "defined for it. Choose a non-constant column.", call. = FALSE)
  ## chain factor: d/d(zeta) of the inflation score, which is minus the
  ## observed weight w1 of bic_info()
  w1s <- -((a$Ic / a$nu^2 + a$Dc / (1 - a$nu)^2) * a$Ga^2 -
           (a$Ic / a$nu - a$Dc / (1 - a$nu)) * a$Ga2)
  Dg  <- gp * s * t(a$Z * w1s)
  ## direct term: the perturbed covariate enters the score for its own
  ## coefficient, contributing s * dl/d(eta) in row p
  Dg[p, ] <- Dg[p, ] + s * a$eps1 * a$Ga * a$Pt
  list(Delta = rbind(Dg, matrix(0, a$m + a$p, a$n)),
       gamma = Dg, beta = NULL, delta = NULL)
}

#' Perturbation matrix for the mean-covariate scheme
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @param p Index of the perturbed covariate in the mean design matrix.
#' @return A list with the full `Delta` and its blocks.
#' @export
delta_meancovar <- function(object, p = 1L) {
  a  <- .bic_aux(object)
  bp <- object$mu.coefficients[p]; s <- sd(a$X[, p])
  if (s <= 0) stop("Column ", p, " of the design is constant (typically the ",
                   "intercept); an additive covariate perturbation is not ",
                   "defined for it. Choose a non-constant column.", call. = FALSE)
  w2p  <- a$Dc * (-a$phi^2 * a$v_star * a$Tm^2 + a$phi * a$eps2 * a$Sm)
  wbdp <- a$Dc * (a$eps2 - a$phi * a$c_star) * a$Tm * a$Rp
  Db   <- bp * s * t(a$X * w2p)
  Dd   <- bp * s * t(a$V * wbdp)
  Db[p, ] <- Db[p, ] + s * a$Dc * a$phi * a$eps2 * a$Tm
  list(Delta = rbind(matrix(0, a$M, a$n), Db, Dd),
       gamma = NULL, beta = Db, delta = Dd)
}

#' Perturbation matrix for the precision-covariate scheme
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @param p Index of the perturbed covariate in the precision design matrix.
#' @return A list with the full `Delta` and its blocks.
#' @export
delta_preccovar <- function(object, p = 1L) {
  a  <- .bic_aux(object)
  dp <- object$sigma.coefficients[p]; s <- sd(a$V[, p])
  if (s <= 0) stop("Column ", p, " of the design is constant (typically the ",
                   "intercept); an additive covariate perturbation is not ",
                   "defined for it. Choose a non-constant column.", call. = FALSE)
  wbp <- a$Dc * (a$eps2 - a$phi * a$c_star) * a$Tm * a$Rp
  wdp <- -a$Dc * (a$w_phph * a$Rp^2 - a$u_t * a$Qp)
  Db  <- dp * s * t(a$X * wbp)
  Dd  <- dp * s * t(a$V * wdp)
  Dd[p, ] <- Dd[p, ] + s * a$u_t * a$Rp
  list(Delta = rbind(matrix(0, a$M, a$n), Db, Dd),
       gamma = NULL, beta = Db, delta = Dd)
}
