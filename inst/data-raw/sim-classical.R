## Table 3 of the paper: classical diagnostics under 5% planted influence.
## Usage:  Rscript sim-classical.R
source("common.R")
R <- 250L
one <- function(seed, n) {
  set.seed(seed)
  df <- gen(n); idx <- sample.int(n, round(0.05 * n)); df <- contaminate(df, idx)
  f <- fitit(df); if (is.null(f)) return(NULL)
  io <- bic_info(f); dd <- delta_caseweights(f)
  cn <- tryCatch(cnc_matrix(dd$Delta, io$info_inv), error = function(e) NULL)
  if (is.null(cn)) return(NULL)
  sc <- cnc_scores(cn, r_max = 4L)
  truth <- logical(n); truth[idx] <- TRUE
  d1 <- sc$B_Et > sc$b2
  m3 <- sc$m_r[, 3]; d2 <- m3 > sqrt(2) * mean(m3)
  c(share1 = cn$lambda[1] / sum(cn$lambda),
    sens1 = mean(d1[truth]), fpr1 = mean(d1[!truth]),
    sens2 = mean(d2[truth]), fpr2 = mean(d2[!truth]))
}
cat(sprintf("%5s %10s %10s %8s %10s %8s %9s\n",
            "n", "lambda1", "sens(2b)", "FPR", "sens(m3)", "FPR", "nonconv"))
for (n in c(50L, 150L, 300L)) {
  rr <- mclapply(seq_len(R), function(s) one(s + n * 1000L, n), mc.cores = CORES)
  ok <- !vapply(rr, is.null, logical(1)); M <- do.call(rbind, rr[ok])
  cat(sprintf("%5d %10.3f %10.2f %8.3f %10.2f %8.3f %9.3f\n", n,
              mean(M[, "share1"]), mean(M[, "sens1"]), mean(M[, "fpr1"]),
              mean(M[, "sens2"]), mean(M[, "fpr2"]), mean(!ok)))
}
