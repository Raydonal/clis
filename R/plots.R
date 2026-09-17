## ===========================================================================
##  plots.R
##  Diagnostic plots: CLIS p-value plots, CNC panels, residuals, envelopes.
## ===========================================================================

#' Plot a CLIS screening result
#'
#' Produces a two-panel diagnostic plot for a [clis_screen()] result: a plot
#' of conformal p-values (with the Benjamini-Hochberg rejection boundary) and
#' a plot of the non-conformity scores with the declared influential
#' observations highlighted.
#'
#' @param x A `clis` object from [clis_screen()].
#' @param label_top Integer; how many top influential points to label.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns `x`.
#' @export
plot_clis <- function(x, label_top = 8L, ...) {
  stopifnot(inherits(x, "clis"))
  op <- par(mfrow = c(1, 2), mar = c(4.4, 4.2, 2.6, 1))
  on.exit(par(op), add = TRUE)

  ## (a) sorted p-values with BH line
  m  <- length(x$pvalues)
  ord <- order(x$pvalues)
  ps  <- x$pvalues[ord]
  bh  <- x$alpha * seq_len(m) / m
  plot(seq_len(m), ps, pch = 20, cex = .7, col = "gray50",
       xlab = "Ordered screening index", ylab = "Conformal p-value",
       main = "(a) Conformal p-values and BH boundary")
  lines(seq_len(m), bh, col = "firebrick", lwd = 2)
  rej_sorted <- which(ps <= bh)
  if (length(rej_sorted) > 0)
    points(rej_sorted, ps[rej_sorted], pch = 20, cex = .9, col = "firebrick")
  legend("topleft", bty = "n", cex = .8,
         legend = c("p-value", "BH boundary", "rejected"),
         pch = c(20, NA, 20), lty = c(NA, 1, NA),
         col = c("gray50", "firebrick", "firebrick"))

  ## (b) scores with influential highlighted
  sc  <- x$scores
  col <- rep("gray60", length(sc))
  col[x$influential_global] <- "firebrick"
  plot(sc, type = "h", col = col, lwd = 1.2,
       xlab = "Observation index", ylab = "CNC score",
       main = "(b) Influence scores (FDR-controlled)")
  ig <- x$influential_global
  if (length(ig) > 0) {
    top <- ig[order(sc[ig], decreasing = TRUE)][seq_len(min(label_top, length(ig)))]
    text(top, sc[top], labels = top, pos = 3, cex = .72,
         col = "firebrick", font = 2)
  }
  invisible(x)
}

#' Plot the classical CNC diagnostic panels
#'
#' @param cnc A list from [cnc_matrix()].
#' @param scores A list from [cnc_scores()].
#' @param decomp Optional list from [cnc_block_decomp()] for a third panel.
#' @param r_plot Integer; which aggregate-contribution order to plot.
#' @param label_top Integer; number of points to label.
#' @param suffix Character; appended to panel titles.
#'
#' @return Invisibly returns `scores`.
#' @export
plot_cnc_panels <- function(cnc, scores, decomp = NULL, r_plot = 3L,
                            label_top = 5L, suffix = "") {
  ncols <- if (!is.null(decomp)) 3L else 2L
  op <- par(mfrow = c(1, ncols), mar = c(4.5, 4.2, 2.5, 1.2))
  on.exit(par(op), add = TRUE)

  lam <- cnc$lambda
  plot(seq_along(lam), lam, pch = 20, cex = .8, col = "steelblue",
       ylim = c(0, max(lam) * 1.1), xlab = "Index",
       ylab = expression(lambda[i]^"*"),
       main = paste0("(a) Normalized eigenvalues", suffix))
  for (ri in seq_along(scores$threshold_r))
    abline(h = scores$threshold_r[ri], lty = 2, col = "gray50")

  mr  <- scores$m_r[, r_plot]; thr <- scores$threshold_mt[r_plot]
  col <- ifelse(mr > thr, "firebrick", "gray60")
  plot(mr, type = "h", lwd = 1.4, col = col, xlab = "Observation index",
       ylab = bquote(m*"["*.(r_plot)*"]"[t]),
       main = paste0("(b) Aggregate contribution", suffix))
  abline(h = thr, lty = 2, col = "firebrick")
  top <- order(mr, decreasing = TRUE)[seq_len(min(label_top, length(mr)))]
  top <- top[mr[top] > thr]
  if (length(top) > 0)
    text(top, mr[top], labels = top, pos = 3, cex = .72,
         col = "firebrick", font = 2)

  if (!is.null(decomp)) {
    Bg <- decomp$B_gamma; Bbd <- decomp$B_betadelta; Bt <- Bg + Bbd
    bp <- barplot(rbind(Bg, Bbd), col = c("steelblue", "coral"),
                  border = NA, xlab = "Observation index",
                  ylab = expression(B[E[t]]),
                  main = paste0("(c) Block decomposition", suffix),
                  ylim = c(0, max(Bt) * 1.15))
    legend("topright", bty = "n", cex = .85, fill = c("steelblue", "coral"),
           legend = c(expression(B^(gamma)), expression(B^(beta*","*delta))))
  }
  invisible(scores)
}

#' Diagnostic residual plots for a BIc model
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @param label_top Integer; number of extreme points to label.
#' @return Invisibly returns a list of residual vectors.
#' @export
plot_residuals <- function(object, label_top = 3L) {
  a <- .bic_aux(object)
  op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
  on.exit(par(op), add = TRUE)

  rp_d <- a$eps1 / sqrt(a$nu * (1 - a$nu))
  rp_c <- ifelse(a$Dc == 1, a$eps2 / sqrt(a$v_star), NA)
  rq   <- bic_quantile_residuals(object)
  lab  <- function(r) {
    top <- order(abs(r), decreasing = TRUE, na.last = NA)[seq_len(label_top)]
    text(top, r[top], labels = top, pos = 3, cex = .7, col = "firebrick")
  }
  plot(rp_d, pch = 20, cex = .8, ylim = c(-4, 4), xlab = "t",
       ylab = expression(r[Pt]^D), main = "(a) Pearson D vs index")
  abline(h = c(-2, 2), lty = 2); lab(rp_d)
  plot(a$nu, rp_d, pch = 20, cex = .8, ylim = c(-4, 4),
       xlab = expression(hat(alpha)[t]), ylab = expression(r[Pt]^D),
       main = "(b) Pearson D vs alpha"); abline(h = c(-2, 2), lty = 2); lab(rp_d)
  plot(rp_c, pch = 20, cex = .8, ylim = c(-4, 4), xlab = "t",
       ylab = expression(r[Pt]^C), main = "(c) Pearson C vs index")
  abline(h = c(-2, 2), lty = 2); lab(rp_c)
  plot(a$mu, rp_c, pch = 20, cex = .8, ylim = c(-4, 4),
       xlab = expression(hat(mu)[t]), ylab = expression(r[Pt]^C),
       main = "(d) Pearson C vs mu"); abline(h = c(-2, 2), lty = 2); lab(rp_c)
  plot(rq, pch = 20, cex = .8, xlab = "t", ylab = expression(r[t]^q),
       main = "(e) Quantile residual"); abline(h = c(-2, 2), lty = 2)
  tq <- order(abs(rq), decreasing = TRUE)[seq_len(label_top)]
  text(tq, rq[tq], labels = tq, pos = 3, cex = .7, col = "firebrick")
  stats::qqnorm(rq, pch = 20, cex = .8, main = "(f) Normal Q-Q")
  stats::qqline(rq, col = "steelblue", lwd = 1.5)
  invisible(list(rp_d = rp_d, rp_c = rp_c, rq = rq))
}

#' Simulated envelope for quantile residuals
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @param B Integer; number of simulated samples.
#' @param level Numeric; envelope coverage (default 0.95).
#' @param seed Optional integer seed.
#' @return Invisibly returns the observed residuals and envelope bands.
#' @export
envelope_bic <- function(object, B = 200L, level = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fam <- object$family[1]; n <- length(object$y)
  rq_obs <- sort(bic_quantile_residuals(object))
  q_norm <- qnorm((seq_len(n) - .5) / n)
  rfun <- match.fun(paste0("r", fam))
  a <- (1 - level) / 2
  sim <- matrix(NA_real_, n, B)
  for (b in seq_len(B)) {
    ys <- rfun(n, mu = object$mu.fv, sigma = object$sigma.fv, nu = object$nu.fv)
    ft <- tryCatch(suppressWarnings(gamlss::gamlss(
      ys ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1, family = fam,
      data = data.frame(ys = ys),
      control = gamlss::gamlss.control(trace = FALSE, n.cyc = 30))),
      error = function(e) NULL)
    if (!is.null(ft)) sim[, b] <- sort(bic_quantile_residuals(ft))
  }
  lo  <- apply(sim, 1, quantile, probs = a, na.rm = TRUE)
  hi  <- apply(sim, 1, quantile, probs = 1 - a, na.rm = TRUE)
  med <- apply(sim, 1, median, na.rm = TRUE)
  op <- par(mar = c(4.4, 4.2, 2.5, 1)); on.exit(par(op), add = TRUE)
  plot(q_norm, rq_obs, pch = 16, cex = .8,
       ylim = range(c(rq_obs, lo, hi), na.rm = TRUE),
       xlab = "Normal quantiles", ylab = expression(r[t]^q),
       main = "Simulated envelope")
  lines(q_norm, sort(lo), lty = 2, col = "gray40")
  lines(q_norm, sort(hi), lty = 2, col = "gray40")
  lines(q_norm, sort(med), lty = 1, col = "gray60")
  abline(0, 1, lty = 3, col = "steelblue")
  invisible(list(obs = rq_obs, lower = lo, upper = hi, median = med))
}

#' Classical local-influence index plot
#'
#' Reproduces the standard influence display of the beta-regression
#' diagnostics literature: the per-observation conformal normal curvature
#' \eqn{B_{E_t}} (or the aggregate contribution \eqn{m[r]_t}) plotted against
#' the observation index, with the reference cutoff drawn and the most
#' influential observations labelled. This is the plot practitioners of local
#' influence expect to see; the conformal screening in [clis_screen()] and
#' [plot_clis()] complements it with an error-controlled declaration.
#'
#' @param object A fitted BEZI/BEOI `gamlss` model.
#' @param scheme Perturbation scheme: "caseweights" (default), "disccovar",
#'   "meancovar", or "preccovar".
#' @param p Covariate index for the covariate-perturbation schemes.
#' @param measure Which quantity to plot: "B_Et" (default) or "m_r".
#' @param r_plot Aggregate-contribution order when `measure = "m_r"`.
#' @param use_fisher Logical; use the Fisher information (default `TRUE`).
#' @param penalised Logical; use the penalised information (default `FALSE`).
#' @param label_top Integer; number of most-influential points to label.
#' @param labels Optional character vector of observation labels (for example
#'   country codes); defaults to the observation index.
#'
#' @return Invisibly, a list with the plotted scores and the cutoff.
#' @export
plot_influence <- function(object, scheme = "caseweights", p = 1L,
                           measure = c("B_Et", "m_r"), r_plot = 3L,
                           use_fisher = TRUE, penalised = FALSE,
                           label_top = 5L, labels = NULL) {
  measure <- match.arg(measure)
  S <- if (isTRUE(penalised)) bic_penalty(object) else NULL
  info  <- bic_info(object, use_fisher = use_fisher, penalty = S)
  delta <- switch(scheme,
    caseweights = delta_caseweights(object),
    disccovar   = delta_disccovar(object, p = p),
    meancovar   = delta_meancovar(object, p = p),
    preccovar   = delta_preccovar(object, p = p),
    stop("Unknown scheme: ", scheme, call. = FALSE))
  cnc <- cnc_matrix(delta$Delta, info$info_inv)
  sc  <- cnc_scores(cnc, r_max = max(4L, r_plot))

  if (measure == "B_Et") {
    val <- sc$B_Et; cut <- sc$b2
    ylab <- expression(B[E[t]])
  } else {
    val <- sc$m_r[, r_plot]; cut <- sc$threshold_mt[r_plot]
    ylab <- bquote(m*"["*.(r_plot)*"]"[t])
  }
  n <- length(val)
  if (is.null(labels)) labels <- as.character(seq_len(n))

  op <- par(mar = c(4.4, 4.4, 2.6, 1)); on.exit(par(op), add = TRUE)
  cols <- ifelse(val > cut, "firebrick", "gray55")
  plot(val, type = "h", lwd = 1.3, col = cols,
       xlab = "Observation index", ylab = ylab,
       main = paste0("Local influence (", scheme, ")"))
  abline(h = cut, lty = 2, col = "firebrick")
  top <- order(val, decreasing = TRUE)[seq_len(min(label_top, n))]
  top <- top[val[top] > cut]
  if (length(top) > 0)
    text(top, val[top], labels = labels[top], pos = 3, cex = .72,
         col = "firebrick", font = 2)
  invisible(list(scores = val, cutoff = cut, flagged = which(val > cut)))
}
