## Tables 6 and 7, and Figures 6-8, for the reading-accuracy application.
## Usage:  Rscript application-reading.R      (needs the 'betareg' package)
suppressPackageStartupMessages({library(gamlss); library(gamlss.dist); library(clis)})
FIGDIR <- Sys.getenv("CLIS_FIGDIR", "../../artigo1_clis/figuras")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
data("ReadingSkills", package = "betareg")
d <- ReadingSkills; d$dys <- as.numeric(d$dyslexia == "yes"); d$y <- d$accuracy1
ctl <- gamlss.control(trace = FALSE, n.cyc = 200)
cat(sprintf("n = %d, dyslexic = %d, at the boundary = %d (%.1f%%), all non-dyslexic = %s\n",
    nrow(d), sum(d$dys), sum(d$y == 1), 100*mean(d$y == 1),
    all(d$dys[d$y == 1] == 0)))

fitit <- function(mu, si, nu, data = d) suppressWarnings(tryCatch(
  gamlss(as.formula(paste("y ~", mu)), sigma.formula = as.formula(paste("~", si)),
         nu.formula = as.formula(paste("~", nu)), family = BEOI, data = data,
         control = ctl), error = function(e) NULL))

cat("\n== Table 6: model selection ==\n")
cat(sprintf("%-14s %-6s %-7s %9s %9s\n", "mean", "prec", "infl", "AIC", "BIC"))
for (r in list(c("dys*iq","dys","iq"), c("dys","dys","iq"), c("dys+iq","dys","iq"),
               c("dys","dys","1"), c("dys*iq","dys","1"), c("dys*iq","1","1"),
               c("dys*iq","dys","dys"), c("dys+pb(iq)","dys","iq"))) {
  f <- fitit(r[1], r[2], r[3])
  cat(sprintf("%-14s %-6s %-7s %9.2f %9.2f\n", r[1], r[2], r[3],
              AIC(f), AIC(f, k = log(nrow(d)))))
}
cat("separated fit: SE of the dyslexia coefficient in the inflation submodel =",
    round(sqrt(diag(vcov(fitit("dys*iq","dys","dys"))))[8], 1), "\n")

fit <- fitit("dys*iq", "dys", "iq")
cat("\n== Table 7: estimates (p-values from the t distribution on 36 df) ==\n")
co <- summary(fit); print(round(co, 4))
g <- coef(fit, what = "nu")
cat(sprintf("P(boundary) at iq = mean -/+ 1 SD: %.3f  %.3f\n",
    plogis(g[1] - g[2]*sd(d$iq)), plogis(g[1] + g[2]*sd(d$iq))))

io <- bic_info(fit); dd <- delta_caseweights(fit)
sc <- cnc_scores_linear(dd$Delta, io$info_inv); dec <- cnc_block_decomp(dd, io)
above <- which(sc$B_Et > sc$b2); above <- above[order(sc$B_Et[above], decreasing = TRUE)]
cat("\n== Influence ==\ncutoff 2b =", round(sc$b2, 4), "| above cutoff:", above, "\n")
for (i in above) cat(sprintf("  child %2d  B_Et = %.4f  inflation share = %6.2f%%\n",
                             i, sc$B_Et[i], 100*dec$ratio_gamma[i]))
m <- sapply(1:50, function(s) min(clis_screen(fit, alpha = 0.10, seed = s)$padj))
n0 <- sapply(1:50, function(s) length(clis_screen(fit, alpha = 0.10, seed = s)$influential))
cat(sprintf("conformal screen over 50 splits: declared %d..%d; min adjusted p %.2f..%.2f (median %.2f)\n",
            min(n0), max(n0), min(m), max(m), median(m)))
f2 <- fitit("dys*iq", "dys", "iq", data = d[-above, ])
cat(sprintf("refit without the flagged children: beta_dys %.4f -> %.4f ; gamma_iq %.4f -> %.4f\n",
    coef(fit)[2], coef(f2)[2], coef(fit, "nu")[2], coef(f2, "nu")[2]))

## ---- figures ----
GREY <- "gray55"; RED <- "firebrick"
r <- bic_quantile_residuals(fit, seed = 1); atom <- d$y == 1
pdf(file.path(FIGDIR, "fig_app_resid.pdf"), width = 8, height = 4.2)
par(mfrow = c(1,2), mar = c(4.4,4.3,2.8,1.1), mgp = c(2.5,0.8,0))
qq <- qqnorm(r, plot.it = FALSE)
plot(qq$x, qq$y, type="n", xlab="Theoretical quantiles",
     ylab="Randomised quantile residual", main="(a) Normal Q-Q plot")
abline(0, 1, col = RED, lwd = 1.8)
points(qq$x[!atom], qq$y[!atom], pch=19, cex=.8, col=GREY)
points(qq$x[atom], qq$y[atom], pch=1, cex=.95, lwd=1.4)
legend("topleft", bty="n", cex=.8, pch=c(19,1), col=c(GREY,"black"),
       legend=c("interior observation","boundary observation"))
plot(fitted(fit,"mu"), r, type="n", xlab=expression("Fitted mean  "*hat(mu)[t]),
     ylab="Randomised quantile residual", main="(b) Residuals against fitted mean")
abline(h=c(-2,0,2), lty=c(2,1,2), col=c(GREY,"black",GREY))
points(fitted(fit,"mu")[!atom], r[!atom], pch=19, cex=.8, col=GREY)
points(fitted(fit,"mu")[atom], r[atom], pch=1, cex=.95, lwd=1.4)
dev.off()

pdf(file.path(FIGDIR, "fig_app_influence.pdf"), width = 8, height = 4.2)
par(mar = c(4.4,4.5,2.6,1.1), mgp = c(2.6,0.8,0))
plot(sc$B_Et, type="h", lwd=2, col=ifelse(seq_len(nrow(d)) %in% above, RED, GREY),
     xlab="Child index", ylab=expression(B[E[t]]),
     main="Local-influence index plot, reading-accuracy data")
abline(h = sc$b2, lty = 2)
text(above, sc$B_Et[above], labels=above, pos=3, cex=.8, col=RED, font=2)
legend("topright", bty="n", cex=.85, lty=2, legend="reference cutoff 2b")
dev.off()

pdf(file.path(FIGDIR, "fig_app_block.pdf"), width = 8, height = 4.2)
par(mar = c(4.4,4.5,2.6,1.1), mgp = c(2.6,0.8,0))
mat <- rbind(dec$B_betadelta[above], dec$B_gamma[above])
bp <- barplot(mat, names.arg = above, col = c("gray35","gray80"), border = NA,
              xlab = "Child", ylab = expression(B[E[t]]),
              main = "Block decomposition of the curvature")
text(bp, colSums(mat), labels = sprintf("%.0f%%", 100*dec$ratio_gamma[above]),
     pos = 3, cex = .8, col = RED, font = 2)
legend("topright", bty="n", cex=.85, fill=c("gray35","gray80"), border=NA,
       legend=c("mean and precision","inflation"))
dev.off()
cat("figures written to", FIGDIR, "\n")
