# PROTOTYPE study: leave-block-out cross-validation (B1) for (lambda1, lambda2)
#
# Compares, on the simulation design of Section 4 (m = 50, FAR(1) with
# Hilbert-Schmidt norm 0.5, Brownian-motion innovations) plus trend T2b
# (deviation from the s-mean nonlinear in t), the ISE_T of
#   oracle    lambda minimising the true ISE (grid, then Nelder-Mead),
#   CV        leave-block-out CV, gap l from the residual autocorrelation of
#             the "bspline" pilot fit,
#   bspline   current default (marginal REML on the B-spline bases),
#   supp      marginal REML as in the supplementary code,
#   Lin       linear trend of Kokoszka and Young (2017).
#
# Usage: Rscript blockcv_study.R [reps] [cores] [output.rds] [mode]
#   mode = "main" (default) or "gap" (sensitivity to the gap l)

library(ftrend)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 40L
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
out_file <- if (length(args) >= 3) args[3] else "blockcv_study.rds"
mode <- if (length(args) >= 4) args[4] else "main"

bs <- bspline_basis(10); bt <- bspline_basis(15)
m <- 50

fit_ise <- function(prep, lam, sim) {
  Th <- ftrend:::penalised_fit(prep, bs, bt, lam)$Theta
  ise(sim$trend, basis_eval(bs, sim$s) %*% Th %*% t(basis_eval(bt, sim$t)), sim$s, sim$t)
}

oracle <- function(prep, sim) {
  f <- function(lr) fit_ise(prep, 10^lr, sim)
  gr <- expand.grid(l1 = seq(-10, 4, by = 1), l2 = seq(-8, 4, by = 1))
  v <- apply(gr, 1, f)
  o <- optim(unlist(gr[which.min(v), ]), f, method = "Nelder-Mead",
             control = list(reltol = 1e-6, maxit = 150))
  c(ise = min(o$value, min(v)), lambda1 = unname(10^o$par[1]), lambda2 = unname(10^o$par[2]))
}

one_rep <- function(trend, N, r) {
  sim <- simulate_fts(N = N, m = m, trend = trend, norm = 0.5,
                      seed = 7e5 + 1e3 * match(trend, c("T1", "T2", "T3", "T4", "T5", "T2b")) + 10 * N + r)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  fb <- ftrend(sim$Y, s = sim$s, select = "reml-bspline")    # bspline (pilot)
  gap_hat <- ftrend:::cv_gap_from_residuals(prep, fb$residuals_internal)
  if (mode == "gap") {
    out <- sapply(c(0, 2, 5, 10, 20), function(l)
      fit_ise(prep, ftrend:::select_lambda_blockcv(prep, bs, bt, as.integer(l))$lambda, sim))
    names(out) <- paste0("CV_l", c(0, 2, 5, 10, 20))
    return(c(out, gap_hat = gap_hat, oracle = oracle(prep, sim)[["ise"]]))
  }
  fs <- ftrend(sim$Y, s = sim$s, select = "reml-supplement",
                control = ftrend_control(supp_k = "te-default"))
  cvs <- ftrend:::select_lambda_blockcv(prep, bs, bt, gap_hat)
  orc <- oracle(prep, sim)
  c(oracle = orc[["ise"]], CV = fit_ise(prep, cvs$lambda, sim),
    bspline = ise(sim$trend, predict(fb, s = sim$s), sim$s, sim$t),
    supp = ise(sim$trend, predict(fs, s = sim$s), sim$s, sim$t),
    Lin = ise(sim$trend, linear_trend(sim$Y)$fitted, sim$s, sim$t),
    gap = gap_hat,
    l1_oracle = orc[["lambda1"]], l2_oracle = orc[["lambda2"]],
    l1_CV = cvs$lambda[1], l2_CV = cvs$lambda[2],
    l1_bspline = unname(fb$lambda[1]), l2_bspline = unname(fb$lambda[2]),
    l1_supp = unname(fs$lambda[1]), l2_supp = unname(fs$lambda[2]))
}

settings <- if (mode == "gap") {
  expand.grid(trend = c("T1", "T5", "T2b"), N = 300, stringsAsFactors = FALSE)
} else {
  expand.grid(trend = c("T1", "T2", "T3", "T4", "T5", "T2b"), N = c(100, 300, 500),
              stringsAsFactors = FALSE)
}

results <- list()
for (i in seq_len(nrow(settings))) {
  tr <- settings$trend[i]; N <- settings$N[i]
  t0 <- Sys.time()
  res <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) one_rep(tr, N, r),
                                           mc.cores = cores))
  results[[i]] <- data.frame(trend = tr, N = N, rep = seq_len(reps), res, check.names = FALSE)
  cm <- colMeans(res)
  cat(sprintf("%-3s N=%4d | %s | %.0fs\n", tr, N,
              paste(sprintf("%s %.4f", names(cm)[seq_len(if (mode == "gap") 7 else 6)],
                            cm[seq_len(if (mode == "gap") 7 else 6)]), collapse = " | "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(do.call(rbind, results), out_file)
}
