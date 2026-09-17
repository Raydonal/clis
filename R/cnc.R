## ===========================================================================
##  cnc.R
##  Conformal normal curvature (CNC) matrix, per-observation scores, and the
##  block decomposition into inflation vs. mean/precision contributions.
## ===========================================================================

#' Conformal normal curvature matrix
#'
#' Computes the Frobenius-normalised conformal normal curvature (CNC) matrix
#' of Poon and Poon (1999) for a given perturbation matrix and information
#' matrix inverse.
#'
#' @param Delta A perturbation matrix of dimension `(M+m+p) x n`, e.g. from
#'   [delta_caseweights()].
#' @param info_inv The inverse information matrix from [bic_info()].
#'
#' @return A list with the `n x n` CNC matrix `B`, its eigenvalues
#'   `lambda` (normalised so that the sum of squares is one), eigenvectors
#'   `vectors`, and the Frobenius norm `normF` of the unnormalised curvature.
#'
#' @references
#' Poon, W.-Y. and Poon, Y. S. (1999). Conformal normal curvature and
#' assessment of local influence. \emph{Journal of the Royal Statistical
#' Society: Series B}, 61(1), 51-61.
#'
#' @export
cnc_matrix <- function(Delta, info_inv) {
  F0 <- crossprod(Delta, info_inv %*% Delta)   # t(Delta) %*% info_inv %*% Delta
  normF <- sqrt(sum(F0 * F0))
  if (normF < .Machine$double.eps)
    stop("Frobenius norm of the curvature is zero.", call. = FALSE)
  B   <- F0 / normF
  eig <- eigen(B, symmetric = TRUE)
  lam <- pmax(eig$values, 0)
  ss  <- sqrt(sum(lam^2))
  if (ss > 0) lam <- lam / ss
  ord <- order(lam, decreasing = TRUE)
  list(B = B, lambda = lam[ord], vectors = eig$vectors[, ord, drop = FALSE],
       normF = normF, F0 = F0)
}

#' Per-observation conformal normal curvature scores
#'
#' Computes the basic-perturbation CNC scores \eqn{B_{E_t}} and the aggregate
#' contribution measures \eqn{m[r]_t} from a fitted CNC matrix.
#'
#' @param cnc A list returned by [cnc_matrix()].
#' @param r_max Integer; the maximum order of `r`-influential eigenvectors to
#'   consider for the aggregate contributions.
#'
#' @return A list with per-observation scores `B_Et`, the matrix of aggregate
#'   contributions `m_r` (one column per `r`), the eigenvalue thresholds
#'   `threshold_r`, the aggregate-contribution thresholds `threshold_mt`, the
#'   number of `r`-influential eigenvectors `k_r`, and the cutoff `b2` for
#'   `B_Et`.
#'
#' @export
cnc_scores <- function(cnc, r_max = 4L) {
  lam <- cnc$lambda
  V   <- cnc$vectors
  n   <- length(lam)
  aij2 <- V^2

  B_Et <- pmax(as.vector(aij2 %*% lam), 0)
  b2   <- 2 * sum(lam) / n
  thr_r <- seq_len(r_max) / sqrt(n)

  m_r <- matrix(0, n, r_max)
  k_r <- integer(r_max)
  thr_mt <- numeric(r_max)
  for (ri in seq_len(r_max)) {
    idx <- which(lam > thr_r[ri])
    k_r[ri] <- length(idx)
    if (length(idx) > 0L) {
      lsel <- lam[idx]
      m_r[, ri] <- sqrt(pmax(as.vector(aij2[, idx, drop = FALSE] %*% lsel), 0))
      ## Benchmark for the aggregate contribution: sqrt(2) times its mean
      ## across observations. Using the mean of the selected eigenvalues
      ## instead puts the threshold on a different scale and makes the rule
      ## effectively never fire.
      thr_mt[ri] <- sqrt(2) * mean(m_r[, ri])
    }
  }
  list(B_Et = B_Et, m_r = m_r, threshold_r = thr_r,
       threshold_mt = thr_mt, k_r = k_r, b2 = b2, n = n)
}

#' Linear-time per-observation CNC scores
#'
#' Computes the basic-perturbation CNC scores \eqn{B_{E_t}} without ever
#' forming the \eqn{n \times n} curvature matrix. The score is the scaled
#' diagonal of \eqn{F_0 = \Delta^\top \mathcal{I}^{-1}\Delta}, and both the
#' diagonal and the Frobenius norm \eqn{\|F_0\|_F} are obtained from
#' quantities of dimension \eqn{(M+m+p)}, so the cost is linear in `n` and
#' the memory footprint is negligible. This is the scalable path used for
#' large samples, where the dense \eqn{n\times n} eigenproblem of
#' [cnc_matrix()] is infeasible.
#'
#' @param Delta A perturbation matrix of dimension `(M+m+p) x n`.
#' @param info_inv The inverse information matrix from [bic_info()].
#'
#' @return A list with per-observation scores `B_Et`, the cutoff `b2`, the
#'   Frobenius norm `normF`, and `n`.
#'
#' @export
cnc_scores_linear <- function(Delta, info_inv) {
  LiD    <- info_inv %*% Delta                 # (M+m+p) x n
  diagF0 <- colSums(Delta * LiD)               # diag(F0) = per-column q-form
  A      <- info_inv %*% tcrossprod(Delta)     # (M+m+p) x (M+m+p), small
  normF  <- sqrt(sum(A * t(A)))                # ||F0||_F = sqrt(tr(A^2))
  if (normF < .Machine$double.eps)
    stop("Frobenius norm of the curvature is zero.", call. = FALSE)
  B_Et <- pmax(diagF0 / normF, 0)
  n    <- length(B_Et)
  ## b2 cutoff = 2 * sum(eigenvalues) / n = 2 * tr(B) / n = 2 * (tr F0/normF)/n
  b2 <- 2 * (sum(diagF0) / normF) / n
  list(B_Et = B_Et, b2 = b2, normF = normF, n = n)
}

#' Block decomposition of conformal normal curvature
#'
#' Decomposes the per-observation CNC scores into a contribution from the
#' inflation-probability submodel and a contribution from the
#' conditional-mean/precision submodel, exploiting the information
#' orthogonality of the BIc model. The two contributions sum exactly to the
#' total, with no cross terms.
#'
#' @param delta_out A perturbation-matrix list (e.g. from
#'   [delta_caseweights()]) with a full `Delta` component.
#' @param info_out An information list from [bic_info()] with `info_inv` and
#'   `idx`.
#'
#' @return A list with the discrete contribution `B_gamma`, the continuous
#'   contribution `B_betadelta`, their sum `B_total`, and the discrete
#'   fraction `ratio_gamma`.
#'
#' @export
cnc_block_decomp <- function(delta_out, info_out) {
  ig  <- info_out$idx$gamma
  ibd <- c(info_out$idx$beta, info_out$idx$delta)
  Li  <- info_out$info_inv
  D   <- delta_out$Delta

  Dg   <- D[ig, , drop = FALSE]
  Dbd  <- D[ibd, , drop = FALSE]
  Lig  <- Li[ig, ig, drop = FALSE]
  Libd <- Li[ibd, ibd, drop = FALSE]

  ## Linear-time block scores: each block score is the scaled diagonal of its
  ## own F0 block, computed without forming any n x n matrix. The shared
  ## normalising constant is the Frobenius norm of the total F0, obtained from
  ## small (M+m+p)-dimensional traces.
  diag_g  <- colSums(Dg  * (Lig  %*% Dg))
  diag_bd <- colSums(Dbd * (Libd %*% Dbd))

  ## ||F0||_F^2 = tr(F0^2), F0 = F0g + F0bd, via small traces only:
  ##   tr(F0g^2)  = tr(Ag^2),  Ag  = Lig  (Dg Dg^T)
  ##   tr(F0bd^2) = tr(Abd^2), Abd = Libd (Dbd Dbd^T)
  ##   tr(F0g F0bd) = tr( Lig (Dg Dbd^T) Libd (Dbd Dg^T) )
  Ag  <- Lig  %*% tcrossprod(Dg)
  Abd <- Libd %*% tcrossprod(Dbd)
  G   <- tcrossprod(Dg, Dbd)                     # (|g|) x (|bd|)
  M1  <- (Lig %*% G) %*% Libd                    # (|g|) x (|bd|)
  cross <- sum(M1 * G)                           # tr(F0g F0bd)
  normF <- sqrt(sum(Ag * t(Ag)) + sum(Abd * t(Abd)) + 2 * cross)
  if (normF < .Machine$double.eps) normF <- 1

  Bg  <- pmax(diag_g  / normF, 0)
  Bbd <- pmax(diag_bd / normF, 0)
  Bt  <- Bg + Bbd
  list(B_gamma = Bg, B_betadelta = Bbd, B_total = Bt,
       ratio_gamma = Bg / pmax(Bt, .Machine$double.eps))
}
