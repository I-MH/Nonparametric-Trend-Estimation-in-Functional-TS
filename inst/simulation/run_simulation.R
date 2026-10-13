# Simulation study of Section 4 in Martínez-Hernández and Genton (2021)
#
# Y_n(s) = T(s, n/N) + X_n(s), s on an equispaced grid of m = 50 points,
# X_n FAR(1) with kernel beta(u, v) = C exp{-(u^2 + v^2)/2}, Hilbert-Schmidt
# norm 0.5, Brownian-motion innovations. For each replicate the trend is
# estimated by
#   TPS      ftrend() with defaults (leave-block-out cross-validation),
#   TPS-supp ftrend(select = "reml-supplement", supp_k = "te-default")
#            (lambda selection of the supplementary code),
#   Lin      the linear least-squares trend of Kokoszka and Young (2017),
# and ISE_T^2 = int int {T - T_hat}^2 ds dt is computed on the grid
# s x {n/N}.
#
# Usage:  Rscript run_simulation.R [reps] [cores] [output.rds]
#   e.g.  Rscript run_simulation.R 1000 4 sim_results.rds

library(ftrend)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 1000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
out_file <- if (length(args) >= 3) args[3] else "sim_results.rds"

settings <- expand.grid(trend = c("T1", "T2", "T3", "T4", "T5"),
                        N = c(100, 300, 500, 1000), stringsAsFactors = FALSE)
m <- 50
norm <- 0.5

one_rep <- function(trend, N, r) {
  sim <- simulate_fts(N = N, m = m, trend = trend, norm = norm,
                      seed = 1e5 * match(trend, paste0("T", 1:5)) + 10 * N + r)
  f_tps <- ftrend(sim$Y, s = sim$s)
  f_sup <- ftrend(sim$Y, s = sim$s, select = "reml-supplement",
                  control = ftrend_control(supp_k = "te-default"))
  lin <- linear_trend(sim$Y)
  c(TPS = ise(sim$trend, predict(f_tps, s = sim$s), sim$s, sim$t),
    `TPS-supp` = ise(sim$trend, predict(f_sup, s = sim$s), sim$s, sim$t),
    Lin = ise(sim$trend, lin$fitted, sim$s, sim$t),
    lambda1 = unname(f_tps$lambda[1]), lambda2 = unname(f_tps$lambda[2]),
    edf = f_tps$edf, edf_supp = f_sup$edf)
}

results <- vector("list", nrow(settings))
for (i in seq_len(nrow(settings))) {
  tr <- settings$trend[i]; N <- settings$N[i]
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(reps), function(r) one_rep(tr, N, r),
                            mc.cores = cores)
  res <- do.call(rbind, res)
  results[[i]] <- data.frame(trend = tr, N = N, rep = seq_len(reps), res,
                             check.names = FALSE)
  cat(sprintf("%s N = %4d: mean ISE TPS %.4f | TPS-supp %.4f | Lin %.4f  (%.1f s)\n",
              tr, N, mean(res[, "TPS"]), mean(res[, "TPS-supp"]), mean(res[, "Lin"]),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
results <- do.call(rbind, results)
saveRDS(results, out_file)

summ <- aggregate(cbind(TPS, `TPS-supp`, Lin, edf, edf_supp) ~ trend + N,
                  data = results, FUN = function(v) c(mean = mean(v), sd = sd(v)))
print(summ)
