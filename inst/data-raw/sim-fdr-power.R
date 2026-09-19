## Table 4, Figure 5 (power curves) and Figure 4 (sensitivity to inflation).
## Usage:  Rscript sim-fdr-power.R
source("common.R")
FIGDIR <- Sys.getenv("CLIS_FIGDIR", "../../artigo1_clis/figuras")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
nb <- 200L

screen_once <- function(seed, n, Delta = 5, classical = FALSE, g0 = -1.4) {
  set.seed(seed)
  df <- if (g0 == -1.4) gen(n) else {
    x <- rnorm(n); v <- rnorm(n); z <- rnorm(n)
    y <- rBEZI(n, mu = plogis(0.3 + x), sigma = exp(1.4 + 0.5 * v),
               nu = plogis(g0 + 0.8 * z))
    data.frame(y = pmin(y, 1 - 1e-8), x = x, v = v, z = z)
  }
  scr <- (n - nb + 1):n; bad <- sample(scr, round(0.15 * nb))
  df$x[bad] <- df$x[bad] + Delta * sd(df$x)
  df$y[bad] <- pmin(rBEZI(length(bad), mu = 0.05, sigma = exp(1.4 + 0.5 * df$v[bad]),
                          nu = plogis(g0 + 0.8 * df$z[bad])), 1 - 1e-8)
  f <- fitit(df); if (is.null(f)) return(NULL)
  s <- scores_of(f); cal <- setdiff(seq_len(n), scr)
  if (classical) { rej <- scr[s$B_Et[scr] > s$b2]
    return(c(fdp = if (length(rej)) mean(!(rej %in% bad)) else 0, pow = mean(bad %in% rej))) }
  p <- conformal_pvalues(s$B_Et[cal], s$B_Et[scr]); padj <- p.adjust(p, "BH")
  out <- c()
  for (q in c(0.05, 0.10)) { rej <- scr[padj <= q]
    out <- c(out, if (length(rej)) mean(!(rej %in% bad)) else 0, mean(bad %in% rej)) }
  c(setNames(out, c("fdr05", "pow05", "fdr10", "pow10")), abar = mean(df$y == 0))
}
avg <- function(R, ...) {
  M <- do.call(rbind, Filter(Negate(is.null), mclapply(seq_len(R), ..., mc.cores = CORES)))
  colMeans(M)
}
cat("== Table 4 ==\n")
cat(sprintf("%5s %8s %7s %8s %7s %8s %7s\n","n","FDR.05","pow.05","FDR.10","pow.10","FDP.2b","pow.2b"))
for (n in c(600L, 1000L, 2000L)) {
  A <- avg(300L, function(s) screen_once(s + n, n))
  B <- avg(120L, function(s) screen_once(s + n + 5e5, n, classical = TRUE))
  cat(sprintf("%5d %8.3f %7.2f %8.3f %7.2f %8.3f %7.2f\n", n,
              A[1], A[2], A[3], A[4], B[1], B[2]))
}
cat("\n== Figure 5: power curves (R = 200) ==\n")
M <- matrix(NA_real_, 3, 6, dimnames = list(c("600","1000","2000"), paste0("D", 1:6)))
for (n in c(600L, 1000L, 2000L)) for (D in 1:6)
  M[as.character(n), D] <- avg(200L, function(s) screen_once(s + n*13L + D*101L, n, Delta = D))["pow10"]
print(round(M, 3))
cols <- c("gray45","steelblue4","firebrick"); pchs <- c(15,17,19); ltys <- c(3,2,1)
pdf(file.path(FIGDIR, "fig_sim_power.pdf"), width = 7, height = 5)
par(mar = c(4.4,4.4,1.6,1.1), mgp = c(2.6,0.8,0))
matplot(1:6, t(M), type="n", ylim=c(0,1), xaxt="n",
        xlab=expression("Outlier magnitude  "*Delta*"  (standard deviations)"),
        ylab=expression("Detection power at  "*q*" = 0.10"))
axis(1, at=1:6); abline(h=seq(0,1,.2), col="gray88")
for (i in 1:3) { lines(1:6, M[i,], col=cols[i], lwd=2.1, lty=ltys[i])
                 points(1:6, M[i,], col=cols[i], pch=pchs[i], cex=1.05) }
legend("bottomright", bty="n", legend=paste("n =", rownames(M)),
       col=cols, lty=ltys, pch=pchs, lwd=2.1, cex=.95)
dev.off()

cat("\n== Figure 4: sensitivity to the inflation level (n = 1000, R = 200) ==\n")
targets <- c(.05,.10,.15,.20,.25,.30,.35,.40)
g0s <- vapply(targets, function(a) uniroot(function(g) mean(plogis(g + 0.8*rnorm(20000))) - a,
                                           c(-8,4))$root, numeric(1))
S <- t(vapply(seq_along(targets), function(i) {
  r <- avg(200L, function(s) screen_once(s + i*331L, 1000L, g0 = g0s[i]))
  c(abar = r["abar"], fdr = r["fdr10"], pow = r["pow10"]) }, numeric(3)))
colnames(S) <- c("abar","fdr","pow"); print(round(S, 3))
pdf(file.path(FIGDIR, "fig_sim_sensitivity_alpha.pdf"), width = 6, height = 4.5)
par(mar = c(4.4,4.4,1.6,1.1), mgp = c(2.6,0.8,0))
plot(S[,"abar"], S[,"pow"], type="n", ylim=c(0,1),
     xlab=expression("Average inflation level  "*bar(alpha)),
     ylab="Power (solid) and false discovery rate (dashed)")
abline(h=seq(0,1,.2), col="gray88"); abline(h=0.10, col="gray50", lty=3)
lines(S[,"abar"], S[,"pow"], col="firebrick", lwd=2.2)
points(S[,"abar"], S[,"pow"], col="firebrick", pch=19)
lines(S[,"abar"], S[,"fdr"], col="steelblue4", lwd=2.2, lty=2)
points(S[,"abar"], S[,"fdr"], col="steelblue4", pch=17)
legend("right", bty="n", cex=.88, lty=c(1,2,3), lwd=c(2.2,2.2,1),
       col=c("firebrick","steelblue4","gray50"),
       legend=c("power","false discovery rate","nominal q = 0.10"))
dev.off()
cat("figures written to", FIGDIR, "\n")
