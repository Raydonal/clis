suppressPackageStartupMessages({library(gamlss.dist); library(clis); library(splines); library(parallel)})
source("fakefit.R")
CORES <- 4L; R <- 400L; n <- 800L; nb <- 150L; q <- 0.10
D2 <- diff(diag(8), differences = 2); P8 <- crossprod(D2)

one <- function(seed, lambda) {
  set.seed(seed)
  x1 <- runif(n,-1,1); v1 <- runif(n,-1,1); z1 <- runif(n,-1,1)
  mu <- plogis(0.2 + sin(2*pi*x1)); phi <- exp(1.5 + 0.75*v1); al <- plogis(-1.4 + 0.5*z1)
  y <- pmin(rBEZI(n, mu=mu, sigma=phi, nu=al), 1-1e-8)
  scr <- (n-nb+1):n; bad <- sample(scr, round(0.15*nb))
  y[bad] <- pmin(rBEZI(length(bad), mu=0.05, sigma=phi[bad], nu=al[bad]), 1-1e-8)

  B <- bs(x1, df=8); Xmu <- cbind(1,B); Xsi <- cbind(1,v1); Xnu <- cbind(1,z1)
  pm <- ncol(Xmu); ps <- 2; pn <- 2; np <- pm+ps+pn
  Pmu <- matrix(0,pm,pm); Pmu[2:9,2:9] <- P8
  S <- matrix(0,np,np); S[1:pm,1:pm] <- lambda*Pmu
  npll <- function(p) {
    m <- plogis(Xmu%*%p[1:pm]); ph <- exp(Xsi%*%p[pm+1:ps]); a <- plogis(Xnu%*%p[pm+ps+1:pn])
    m <- pmin(pmax(m,1e-10),1-1e-10); a <- pmin(pmax(a,1e-12),1-1e-12); ph <- pmin(pmax(ph,1e-8),1e8)
    at <- y==0
    ll <- sum(log(a[at])) + sum(log(1-a[!at]) + dbeta(y[!at], m[!at]*ph[!at], (1-m[!at])*ph[!at], log=TRUE))
    -(ll - 0.5*lambda*drop(t(p[1:pm])%*%Pmu%*%p[1:pm]))
  }
  st <- c(qlogis(mean(pmin(pmax(y,.01),.99))), rep(0,8), 1.5, 0, -1.4, 0)
  op <- try(nlminb(st, npll, control=list(iter.max=2000, eval.max=4000)), silent=TRUE)
  if (inherits(op,"try-error") || op$convergence != 0) return(NULL)
  p <- op$par
  f <- mkfit("BEZI", y, Xmu, Xsi, Xnu, p[1:pm], p[pm+1:ps], p[pm+ps+1:pn])
  io <- try(suppressWarnings(if (lambda==0) bic_info(f) else bic_info(f, penalty=S)), silent=TRUE)
  if (inherits(io,"try-error")) return(NULL)
  edf <- if (lambda==0) np else io$edf
  s <- try(cnc_scores_linear(delta_caseweights(f)$Delta, io$info_inv)$B_Et, silent=TRUE)
  if (inherits(s,"try-error")) return(NULL)
  cal <- setdiff(seq_len(n), scr)
  pv <- conformal_pvalues(s[cal], s[scr])
  nulls <- pv[!(scr %in% bad)]
  r_bh <- scr[p.adjust(pv,"BH") <= q]; r_by <- scr[p.adjust(pv,"BY") <= q]
  c(edf = edf,
    fdr = if (length(r_bh)) mean(!(r_bh %in% bad)) else 0, pow = mean(bad %in% r_bh),
    fdrBY = if (length(r_by)) mean(!(r_by %in% bad)) else 0, powBY = mean(bad %in% r_by),
    e05 = mean(nulls <= .05), e10 = mean(nulls <= .10))
}
cat(sprintf("%-7s %7s %7s %7s %8s %8s %7s %7s\n","lambda","edf","FDR","power","FDR_BY","pow_BY","P<=.05","P<=.10"))
for (lam in c(0,1,10,100)) {
  rr <- mclapply(seq_len(R), function(s) one(s + 977L, lam), mc.cores=CORES)
  M <- do.call(rbind, Filter(Negate(is.null), rr))
  cat(sprintf("%-7g %7.2f %7.3f %7.2f %8.3f %8.2f %7.3f %7.3f   (%d reps)\n", lam,
      mean(M[,"edf"]), mean(M[,"fdr"]), mean(M[,"pow"]), mean(M[,"fdrBY"]),
      mean(M[,"powBY"]), mean(M[,"e05"]), mean(M[,"e10"]), nrow(M)))
}
