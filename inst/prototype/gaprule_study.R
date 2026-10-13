# PROTOTYPE study: rules for the gap l in leave-block-out CV (B1)
#
# Rules compared (ISE_T relative to the oracle lambda):
#   R1  fixed l0 = max(5, ceiling(N^(1/3)))
#   R2  plug-in: tail fraction tau(l) <= 0.05 from R1 residuals (flat-top lag),
#       bounded by 2 <= l, 2l + 1 <= N / (2 df_t)
#   R3  R2 + one-standard-error rule along lambda2 (dependence-robust SE)
#   R4  fixed l = 10
# All lambdas are capped where the fit stops changing (lambda_caps), which
# avoids the ill-conditioning of the penalised system at very large lambda.
# Design: Section 4 (m = 50, FAR(1), Hilbert-Schmidt norm 0.5) plus T2b.
#
#   R1se (mode "r1se") fixed l0 + one-standard-error rule, without the plug-in step
#
# Usage: Rscript gaprule_study.R [reps] [cores] [output.rds] [mode] [Ns]
#   mode = "main" (R1-R4, default) or "r1se"; Ns = comma-separated sample sizes

library(ftrend)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 30L
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
out_file <- if (length(args) >= 3) args[3] else "gaprule_study.rds"
mode <- if (length(args) >= 4) args[4] else "main"
Ns <- if (length(args) >= 5) as.numeric(strsplit(args[5], ",")[[1]]) else c(100, 300, 500, 1000)

bs <- bspline_basis(10); bt <- bspline_basis(15)
m <- 50
trends <- c("T1", "T2", "T3", "T4", "T5", "T2b")

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
  min(o$value, min(v))
}

one_rep <- function(trend, N, r) {
  sim <- simulate_fts(N = N, m = m, trend = trend, norm = 0.5,
                      seed = 9e5 + 1e3 * match(trend, trends) + 10 * N + r)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  # R1
  l0 <- as.integer(max(5, ceiling(N^(1/3))))
  s0 <- ftrend:::select_lambda_blockcv(prep, bs, bt, l0)
  if (mode == "r1se") {
    M <- max(2L, 2L * l0)                                   # lag for the SE of the paired differences
    lam <- ftrend:::one_se_lambda2(s0$pre, s0$lambda, M, cap2 = s0$caps[2])
    return(c(oracle = oracle(prep, sim), R1 = fit_ise(prep, s0$lambda, sim),
             R1se = fit_ise(prep, lam, sim), l0 = l0))
  }
  # R2: plug-in gap from R1 residuals
  Th <- ftrend:::penalised_fit(prep, bs, bt, s0$lambda)$Theta
  V <- basis_eval(bs, sim$s); Z <- basis_eval(bt, prep$u_t)
  res <- lapply(seq_len(N), function(n) sim$Y[, n] - as.vector(V %*% Th %*% Z[n, ]))
  gp <- ftrend:::cv_gap_plugin(prep, res, ftrend:::time_df(prep, bt, s0$lambda[2]))
  s2 <- if (gp$l == l0) s0 else
    ftrend:::select_lambda_blockcv(prep, bs, bt, gp$l, start = s0$lambda, maxit = 100L)
  # R3: one-SE rule along lambda2
  lam3 <- ftrend:::one_se_lambda2(s2$pre, s2$lambda, gp$M, cap2 = s0$caps[2])
  # R4: fixed l = 10
  s4 <- ftrend:::select_lambda_blockcv(prep, bs, bt, 10L, start = s0$lambda, maxit = 100L)
  c(oracle = oracle(prep, sim),
    R1 = fit_ise(prep, s0$lambda, sim), R2 = fit_ise(prep, s2$lambda, sim),
    R3 = fit_ise(prep, lam3, sim), R4 = fit_ise(prep, s4$lambda, sim),
    l0 = l0, l_plugin = gp$l, l_dep = gp$l_dep, l_max = gp$l_max,
    bound_binds = gp$bound_binds, M = gp$M,
    l2_R2 = s2$lambda[2], l2_R3 = lam3[2])
}

settings <- expand.grid(trend = trends, N = Ns, stringsAsFactors = FALSE)
rules <- if (mode == "r1se") c("R1", "R1se") else c("R1", "R2", "R3", "R4")
results <- list()
for (i in seq_len(nrow(settings))) {
  tr <- settings$trend[i]; N <- settings$N[i]
  t0 <- Sys.time()
  res <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) one_rep(tr, N, r),
                                           mc.cores = cores))
  results[[i]] <- data.frame(trend = tr, N = N, rep = seq_len(reps), res, check.names = FALSE)
  rat <- res[, rules, drop = FALSE] / res[, "oracle"]
  cat(sprintf("%-3s N=%4d | median ratio %s | q90 %s | %.0fs\n",
              tr, N, paste(sprintf("%.2f", apply(rat, 2, median)), collapse = " "),
              paste(sprintf("%.2f", apply(rat, 2, quantile, 0.9)), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(do.call(rbind, results), out_file)
}
