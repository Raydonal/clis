## Table 8 and Figures 9-10 for the lung-function application.
## Usage:  Rscript application-lung.R      (needs 'gamlss.data')
suppressPackageStartupMessages({library(gamlss); library(gamlss.dist); library(clis); library(splines)})
source("fakefit.R")
FIGDIR <- Sys.getenv("CLIS_FIGDIR", "../../artigo1_clis/figuras")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
data("lungFunction", package = "gamlss.data"); L <- lungFunction
n <- nrow(L); y <- L$slf
hc <- L$height - mean(L$height); ac <- L$age - mean(L$age)
cat(sprintf("n = %d; at the maximum = %d (%.1f%%); median age = %.2f; max age = %g (%d subjects)\n",
    n, sum(y == 1), 100*mean(y == 1), median(L$age), max(L$age), sum(L$age == max(L$age))))
cat(sprintf("upper quartile of age = %.1f; share over 18 = %.1f%%\n",
    quantile(L$age, .75), 100*mean(L$age > 18)))

## mean design: intercept + 8-column cubic B-spline in age + height;
## penalty: lambda times the second-difference operator on the spline block.
B <- bs(L$age, df = 8); Xmu <- cbind(1, B, hc); Xsi <- cbind(1, hc); Xnu <- cbind(1, ac)
pm <- ncol(Xmu); np <- pm + 4
Pmu <- matrix(0, pm, pm); Pmu[2:9, 2:9] <- crossprod(diff(diag(8), differences = 2))
loglik <- function(p) {
  mu <- plogis(Xmu %*% p[1:pm]); ph <- exp(Xsi %*% p[pm + 1:2]); al <- plogis(Xnu %*% p[pm + 3:4])
  at <- y == 1
  sum(log(al[at])) + sum(log(1 - al[!at]) +
      dbeta(y[!at], mu[!at]*ph[!at], (1 - mu[!at])*ph[!at], log = TRUE))
}
dat <- data.frame(y = y, as.data.frame(B), hc = hc, ac = ac)
names(dat)[2:9] <- paste0("b", 1:8)
f0 <- suppressWarnings(gamlss(as.formula(paste("y ~", paste(c(paste0("b",1:8), "hc"), collapse = "+"))),
      sigma.formula = ~hc, nu.formula = ~ac, family = BEOI, data = dat,
      control = gamlss.control(trace = FALSE, n.cyc = 200)))
start <- c(coef(f0, "mu"), coef(f0, "sigma"), coef(f0, "nu"))
Sof <- function(lam) { S <- matrix(0, np, np); S[1:pm, 1:pm] <- lam * Pmu; S }
fit_lambda <- function(lam) {
  op <- nlminb(start, function(p) -(loglik(p) - 0.5*lam*drop(t(p[1:pm]) %*% Pmu %*% p[1:pm])),
               control = list(iter.max = 3000, eval.max = 6000, rel.tol = 1e-12))
  f <- mkfit("BEOI", y, Xmu, Xsi, Xnu, op$par[1:pm], op$par[pm + 1:2], op$par[pm + 3:4])
  list(par = op$par, fit = f, io = suppressWarnings(bic_info(f, penalty = Sof(lam))))
}
cat("\n== Table 8 ==\n")
lin <- suppressWarnings(gamlss(y ~ ac + hc, sigma.formula = ~hc, nu.formula = ~ac,
       family = BEOI, data = data.frame(y = y, ac = ac, hc = hc),
       control = gamlss.control(trace = FALSE, n.cyc = 200)))
cat(sprintf("%-16s %8s %8s %10s %10s\n", "mean in age", "lambda", "edf", "AIC", "BIC"))
cat(sprintf("%-16s %8s %8.2f %10.1f %10.1f\n", "linear", "---", lin$df.fit,
            AIC(lin), AIC(lin, k = log(n))))
res <- list()
for (lam in c(1, 10, 100, 1e3, 1e6)) {
  r <- fit_lambda(lam); res[[as.character(lam)]] <- r
  m2 <- -2*loglik(r$par)
  cat(sprintf("%-16s %8g %8.2f %10.1f %10.1f\n", "penalised spline", lam,
              r$io$edf, m2 + 2*r$io$edf, m2 + log(n)*r$io$edf))
}
sel <- res[["1"]]
cat(sprintf("\nselected fit (lambda = 1): k(phi) = %.3f + %.4f height ; h(alpha) = %.3f %+.3f age\n",
    sel$par[pm+1], sel$par[pm+2], sel$par[pm+3], sel$par[pm+4]))

sc <- cnc_scores_linear(delta_caseweights(sel$fit)$Delta, sel$io$info_inv)
dec <- cnc_block_decomp(delta_caseweights(sel$fit), sel$io)
top <- which.max(sc$B_Et); o <- order(sc$B_Et, decreasing = TRUE)
cat(sprintf("above the cutoff: %d of %d | case %d: B_Et = %.4f, %.1f times the next, inflation share %.4f%%\n",
    sum(sc$B_Et > sc$b2), n, top, sc$B_Et[top], sc$B_Et[top]/sc$B_Et[o[2]],
    100*dec$ratio_gamma[top]))
r <- bic_quantile_residuals(sel$fit, seed = 1)
cat(sprintf("residual of case %d = %.3f (the largest in absolute value)\n", top, r[top]))
mm <- sapply(1:20, function(s) min(clis_screen(sel$fit, alpha = .10, seed = s, penalty = Sof(1))$padj))
nn <- sapply(1:20, function(s) length(clis_screen(sel$fit, alpha = .10, seed = s, penalty = Sof(1))$influential))
cat(sprintf("conformal screen over 20 splits: declared %d..%d; min adjusted p %.2f..%.2f\n",
            min(nn), max(nn), min(mm), max(mm)))
keep <- setdiff(seq_len(n), top)
Xmu2 <- Xmu[keep,]; Xsi2 <- Xsi[keep,]; Xnu2 <- Xnu[keep,]; y2 <- y[keep]
ll2 <- function(p) { mu <- plogis(Xmu2 %*% p[1:pm]); ph <- exp(Xsi2 %*% p[pm+1:2])
  al <- plogis(Xnu2 %*% p[pm+3:4]); at <- y2 == 1
  sum(log(al[at])) + sum(log(1-al[!at]) + dbeta(y2[!at], mu[!at]*ph[!at], (1-mu[!at])*ph[!at], log=TRUE)) }
op2 <- nlminb(start, function(p) -(ll2(p) - 0.5*drop(t(p[1:pm]) %*% Pmu %*% p[1:pm])),
              control = list(iter.max = 6000, eval.max = 12000, rel.tol = 1e-14))
cat(sprintf("refit without case %d: precision-on-height %.6f -> %.6f ; inflation slope %.4f -> %.4f\n",
    top, sel$par[pm+2], op2$par[pm+2], sel$par[pm+4], op2$par[pm+4]))

GREY <- "gray55"; RED <- "firebrick"
pdf(file.path(FIGDIR, "fig_lung_influence.pdf"), width = 8, height = 4.2)
par(mar = c(4.4,4.5,2.6,1.1), mgp = c(2.6,0.8,0))
plot(sc$B_Et, type = "h", col = ifelse(seq_len(n) == top, RED, GREY),
     lwd = ifelse(seq_len(n) == top, 2, 1), xlab = "Observation index",
     ylab = expression(B[E[t]]),
     main = "Conformal normal curvature, lung-function data (n = 3,164)")
abline(h = sc$b2, lty = 2)
text(top, sc$B_Et[top], labels = top, pos = 1, offset = .7, cex = .85, col = RED, font = 2)
legend("topleft", bty = "n", cex = .82, lty = c(1,2), lwd = c(2,1), col = c(RED,"black"),
       legend = c(paste("case", top), "reference cutoff 2b"))
dev.off()
pdf(file.path(FIGDIR, "fig_lung_resid.pdf"), width = 8, height = 4.2)
par(mfrow = c(1,2), mar = c(4.4,4.3,2.8,1.1), mgp = c(2.5,0.8,0))
qq <- qqnorm(r, plot.it = FALSE)
plot(qq$x, qq$y, pch = 19, cex = .35, col = GREY, xlab = "Theoretical quantiles",
     ylab = "Randomised quantile residual", main = "(a) Normal Q-Q plot")
abline(0, 1, col = RED, lwd = 1.8)
points(qq$x[order(r)[n]], max(r), pch = 1, cex = 2, col = RED, lwd = 2)
text(qq$x[order(r)[n]], max(r), labels = top, pos = 2, cex = .85, col = RED, font = 2)
plot(L$age, r, pch = 19, cex = .35, col = GREY, xlab = "Age (years)",
     ylab = "Randomised quantile residual", main = "(b) Residuals against age")
abline(h = c(-3,0,3), lty = c(2,1,2), col = c(GREY,"black",GREY))
points(L$age[top], r[top], pch = 1, cex = 2, col = RED, lwd = 2)
text(L$age[top], r[top], labels = top, pos = 2, cex = .85, col = RED, font = 2)
dev.off()
cat("figures written to", FIGDIR, "\n")
