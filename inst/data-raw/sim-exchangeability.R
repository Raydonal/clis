source("common.R")
R <- 250L
## (a) uniformity of conformal p-values under no contamination
one <- function(seed, n) {
  set.seed(seed); df <- gen(n); f <- fitit(df); if (is.null(f)) return(NULL)
  s <- scores_of(f)$B_Et
  nc <- n %/% 2L; cal <- sample.int(n, nc); scr <- setdiff(seq_len(n), cal)
  p <- conformal_pvalues(s[cal], s[scr])
  list(p = p, rank1 = sum(s[cal] < s[scr[1]]) + 1L, nc = nc)
}
for (n in c(300L, 150L)) {
  rr <- mclapply(seq_len(R), function(s) one(s + n*7919L, n), mc.cores=CORES)
  rr <- rr[!vapply(rr, is.null, logical(1))]
  P <- unlist(lapply(rr, `[[`, "p")); rk <- vapply(rr, `[[`, integer(1), "rank1"); nc <- rr[[1]]$nc
  cat(sprintf("\n== n = %d, %d reps, %d pooled p-values ==\n", n, length(rr), length(P)))
  cat(sprintf("mean p = %.3f\n", mean(P)))
  for (a in c(.05,.10,.20,.50)) cat(sprintf("  P(p <= %.2f) = %.3f\n", a, mean(P <= a)))
  brk <- floor(seq(0, nc+1, length.out = 11))
  obs <- table(cut(rk, breaks = brk, include.lowest=TRUE))
  ct <- chisq.test(obs)
  cat(sprintf("rank test: chi2 = %.1f on %d df, p = %.2f | mean rank = %.1f (expected %.1f)\n",
      ct$statistic, ct$parameter, ct$p.value, mean(rk), (nc+2)/2))
}
