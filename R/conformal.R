## ===========================================================================
##  conformal.R
##  Conformal Local Influence Screening (CLIS).
##
##  The central contribution of the package: wrap the conformal normal
##  curvature score in a split-conformal testing procedure that yields
##  per-observation conformal p-values, then apply Benjamini-Hochberg to
##  control the false discovery rate (FDR) of the declared influential set.
## ===========================================================================

#' Conformal p-values from non-conformity scores
#'
#' Given calibration scores (assumed to come from non-influential, "clean"
#' observations) and test scores, computes marginal conformal p-values in the
#' sense of Bates and others (2023). Larger scores indicate stronger evidence
#' of influence, so the p-value for a test point is the calibrated rank of its
#' score among the calibration scores.
#'
#' @param cal_scores Numeric vector of calibration non-conformity scores.
#' @param test_scores Numeric vector of test non-conformity scores.
#'
#' @return A numeric vector of conformal p-values, one per test score.
#'
#' @details For a test score \eqn{s} and calibration scores
#'   \eqn{c_1, \ldots, c_n}, the conformal p-value is
#'   \deqn{p = \frac{1 + \#\{i : c_i \ge s\}}{n + 1}.}
#'   These p-values are marginally valid (super-uniform under the null that
#'   the test point is exchangeable with the calibration set) and, by the
#'   positive-dependence result of Bates and others (2023), permit
#'   Benjamini-Hochberg FDR control.
#'
#' @references
#' Bates, S., Candes, E., Lei, L., Romano, Y. and Sesia, M. (2023). Testing
#' for outliers with conformal p-values. \emph{The Annals of Statistics},
#' 51(1), 149-178.
#'
#' @examples
#' set.seed(1)
#' cal  <- rnorm(100)
#' test <- c(rnorm(8), 4, 5)   # last two are outliers
#' conformal_pvalues(cal, test)
#'
#' @export
conformal_pvalues <- function(cal_scores, test_scores) {
  n <- length(cal_scores)
  vapply(test_scores, function(s) (1 + sum(cal_scores >= s)) / (n + 1),
         numeric(1))
}

#' Conformal Local Influence Screening
#'
#' Performs scalable, FDR-controlled influence screening for a fitted
#' zero-or-one inflated beta regression model. The conformal normal curvature
#' score of each observation is used as a non-conformity score within a
#' split-conformal procedure: the data are partitioned into a calibration set
#' (presumed clean) and a screening set; conformal p-values are computed for
#' the screening set; and Benjamini-Hochberg adjustment declares an
#' influential subset with false discovery rate controlled at level `alpha`.
#'
#' @param object A fitted `gamlss` model of family `BEZI` or `BEOI`.
#' @param scheme Character; the perturbation scheme. One of `"caseweights"`
#'   (default), `"disccovar"`, `"meancovar"`, or `"preccovar"`.
#' @param p Integer; covariate index for the covariate-perturbation schemes.
#' @param alpha Numeric in (0, 1); target false discovery rate. Default 0.1.
#' @param calib_frac Numeric in (0, 1); fraction of observations used for
#'   conformal calibration when `calib_idx` is not supplied. Default 0.5.
#' @param calib_idx Optional integer vector of observation indices to use as
#'   the calibration set. The conformal false discovery rate guarantee
#'   requires the calibration set to be (nearly) free of influential points;
#'   when a trusted clean subset is known, supply it here. If `NULL` (default)
#'   a random subset of size `floor(calib_frac * n)` is drawn, which is
#'   appropriate only when influential points are rare.
#' @param use_fisher Logical; use the Fisher information (default `TRUE`).
#' @param r_max Integer; order for aggregate contributions used as the score.
#' @param score Character; which CNC quantity to use as the non-conformity
#'   score: `"B_Et"` (basic perturbation, default) or `"m_r"` (aggregate
#'   contribution at order `r_max`).
#' @param penalty Optional penalty matrix, of the dimension of the full
#'   coefficient vector, added to the information before inversion. Use this
#'   when the smooth basis is supplied explicitly as design columns.
#' @param penalised Logical; if `TRUE`, use the penalised information
#'   `J + S` obtained from the smooth terms of an additive `gamlss` fit,
#'   implementing the semiparametric extension. The penalty is extracted
#'   with [bic_penalty()]. Default `FALSE`.
#' @param seed Optional integer for reproducible calibration splitting.
#'
#' @return An object of class `clis` with components:
#'   \describe{
#'     \item{influential}{Integer indices (into the screening set) declared
#'       influential at FDR level `alpha`.}
#'     \item{influential_global}{Integer indices into the original data.}
#'     \item{pvalues}{Conformal p-values for the screening set.}
#'     \item{padj}{Benjamini-Hochberg adjusted p-values.}
#'     \item{scores}{Non-conformity scores for all observations.}
#'     \item{calib_idx, screen_idx}{Calibration and screening indices.}
#'     \item{alpha, scheme}{Inputs echoed back.}
#'     \item{decomp}{Block decomposition (case-weights scheme only).}
#'   }
#'
#' @details
#' Classical local influence diagnostics rank observations by a curvature
#' measure and rely on visual inspection of an index plot against an ad-hoc
#' threshold. This neither scales to large \eqn{n} nor provides any control of
#' the type-I error rate. CLIS addresses both limitations: the conformal
#' wrapper supplies a finite-sample FDR guarantee, and the per-observation
#' score is computed in linear time after a single model fit.
#'
#' The calibration set is assumed to be predominantly free of influential
#' observations. Because influential points are typically rare, a random
#' split satisfies this approximately; for adversarial settings, a robust
#' pre-filter can be applied before calibration.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("gamlss", quietly = TRUE) &&
#'     requireNamespace("betareg", quietly = TRUE)) {
#'   data("ReadingSkills", package = "betareg")
#'   ReadingSkills$dys <- as.numeric(ReadingSkills$dyslexia == "yes")
#'   fit <- gamlss::gamlss(accuracy1 ~ dys * iq,
#'                         sigma.formula = ~ dys, nu.formula = ~ iq,
#'                         family = gamlss.dist::BEOI, data = ReadingSkills,
#'                         control = gamlss::gamlss.control(trace = FALSE))
#'   res <- clis_screen(fit, alpha = 0.1, seed = 1)
#'   print(res)
#' }
#' }
#'
#' @export
clis_screen <- function(object, scheme = "caseweights", p = 1L,
                         alpha = 0.1, calib_frac = 0.5, calib_idx = NULL,
                         use_fisher = TRUE,
                         r_max = 4L, score = c("B_Et", "m_r"),
                         penalised = FALSE, penalty = NULL, seed = NULL) {
  score <- match.arg(score)
  if (!is.null(seed)) set.seed(seed)

  ## An explicit penalty always wins: it is the route to penalised screening
  ## for fits whose basis is supplied as design columns.
  S <- if (!is.null(penalty)) penalty
       else if (isTRUE(penalised)) bic_penalty(object) else NULL
  info_out <- bic_info(object, use_fisher = use_fisher, penalty = S)
  delta_out <- switch(scheme,
    caseweights = delta_caseweights(object),
    disccovar   = delta_disccovar(object, p = p),
    meancovar   = delta_meancovar(object, p = p),
    preccovar   = delta_preccovar(object, p = p),
    stop("Unknown scheme: ", scheme, call. = FALSE))

  ## Score computation. The B_Et score has a linear-time form that never
  ## builds the n x n curvature matrix, which is essential at scale (the dense
  ## matrix is n x n and its eigendecomposition is infeasible for large n).
  ## The aggregate m_r score genuinely needs the eigenvectors, so it falls
  ## back to the dense path; a threshold guards against an out-of-memory dense
  ## build for very large n.
  n <- ncol(delta_out$Delta)
  if (score == "B_Et") {
    sc_lin <- cnc_scores_linear(delta_out$Delta, info_out$info_inv)
    scores <- sc_lin$B_Et
    cnc <- NULL; sc <- sc_lin
  } else {
    max_dense <- getOption("clis.max_dense_n", 5000L)
    if (n > max_dense)
      stop("The 'm_r' score requires the dense ", n, " x ", n,
           " curvature matrix, which is too large. Use score = 'B_Et' ",
           "(linear time) for large samples, or raise ",
           "options(clis.max_dense_n=).", call. = FALSE)
    cnc <- cnc_matrix(delta_out$Delta, info_out$info_inv)
    sc  <- cnc_scores(cnc, r_max = r_max)
    scores <- sc$m_r[, r_max]
  }
  ## Calibration set: either supplied explicitly (recommended when a clean
  ## subset is known) or drawn at random. The conformal guarantee requires
  ## the calibration set to be (nearly) free of influential points; when
  ## outliers may be present, supply a trusted clean subset via `calib_idx`.
  if (!is.null(calib_idx)) {
    calib_idx <- as.integer(calib_idx)
    if (any(calib_idx < 1L | calib_idx > n))
      stop("`calib_idx` out of range.", call. = FALSE)
    cal_idx <- calib_idx
  } else {
    n_cal   <- floor(calib_frac * n)
    cal_idx <- sample.int(n, n_cal)
  }
  screen_idx <- setdiff(seq_len(n), cal_idx)

  cal_scores    <- scores[cal_idx]
  screen_scores <- scores[screen_idx]

  pvals <- conformal_pvalues(cal_scores, screen_scores)
  padj  <- p.adjust(pvals, method = "BH")
  reject <- which(padj <= alpha)

  decomp <- if (scheme == "caseweights")
    cnc_block_decomp(delta_out, info_out) else NULL

  structure(
    list(influential = reject,
         influential_global = screen_idx[reject],
         pvalues = pvals, padj = padj,
         scores = scores, cnc = cnc, cnc_scores = sc,
         calib_idx = cal_idx, screen_idx = screen_idx,
         alpha = alpha, scheme = scheme, score = score,
         penalised = penalised, edf = info_out$edf,
         edf_blocks = info_out$edf_blocks,
         decomp = decomp, n = n),
    class = "clis")
}

#' @export
print.clis <- function(x, ...) {
  cat("Conformal Local Influence Screening\n")
  cat("  Scheme:        ", x$scheme, "\n", sep = "")
  cat("  Score:         ", x$score, "\n", sep = "")
  cat("  Target FDR:    ", x$alpha, "\n", sep = "")
  cat("  Observations:  ", x$n,
      " (", length(x$calib_idx), " calibration, ",
      length(x$screen_idx), " screening)\n", sep = "")
  cat("  Declared influential: ", length(x$influential),
      " (FDR <= ", x$alpha, ")\n", sep = "")
  if (isTRUE(x$penalised) && !is.null(x$edf))
    cat("  Penalised fit, effective df: ", round(x$edf, 2),
        " (gamma ", round(x$edf_blocks[["gamma"]], 2),
        ", beta+delta ", round(x$edf_blocks[["betadelta"]], 2), ")\n",
        sep = "")
  if (length(x$influential_global) > 0) {
    idx <- x$influential_global
    show <- if (length(idx) > 15) c(idx[1:15], NA) else idx
    cat("  Global indices: ",
        paste(ifelse(is.na(show), "...", show), collapse = " "), "\n",
        sep = "")
  }
  invisible(x)
}

#' @export
summary.clis <- function(object, ...) {
  cat("Conformal Local Influence Screening -- summary\n\n")
  print(object)
  if (!is.null(object$decomp) && length(object$influential_global) > 0) {
    idx <- object$influential_global
    rg  <- object$decomp$ratio_gamma[idx]
    cat("\n  Block attribution of declared influential observations:\n")
    cat("    Inflation-driven (ratio_gamma > 0.5): ", sum(rg > 0.5), "\n",
        sep = "")
    cat("    Mean/precision-driven:                ", sum(rg <= 0.5), "\n",
        sep = "")
  }
  invisible(object)
}
