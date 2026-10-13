# PROTOTYPE study: variants of select = "reml-supplement"
#
# Basis dimension of the marginal te(., bs = "cr") fits: "te-default" (5, as in
# the supplementary code) or "matched" (k1 = 10, k2 = 15); error model: none,
# AR(1) in time (series int Y_n), AR(1) in both directions (also along s for
# the mean curve). Also "reml-bspline" for reference. Same simulated data as
# gaprule_study.R (same seeds), so the oracle and the CV default (R3) are taken
# from that study's results and merged.
#
# Usage: Rscript supplement_ar_study.R [reps] [cores] [output.rds] [gaprule.rds]

library(ftrend)
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 30L
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
out_file <- if (length(args) >= 3) args[3] else "supplement_ar_study.rds"
gap_file <- if (length(args) >= 4) args[4] else "gaprule_study.rds"

m <- 50
trends <- c("T1", "T2", "T3", "T4", "T5", "T2b")
variants <- expand.grid(k = c("te-default", "matched"), ar = c("none", "time", "both"),
                        stringsAsFactors = FALSE)

one_rep <- function(trend, N, r) {
  sim <- simulate_fts(N = N, m = m, trend = trend, norm = 0.5,
                      seed = 9e5 + 1e3 * match(trend, trends) + 10 * N + r)
  isef <- function(f) ise(sim$trend, predict(f, s = sim$s), sim$s, sim$t)
  out <- c(bspline = isef(ftrend(sim$Y, s = sim$s, select = "reml-bspline")))
  for (i in seq_len(nrow(variants))) {
    f <- suppressWarnings(ftrend(sim$Y, s = sim$s, select = "reml-supplement",
                                 control = ftrend_control(supp_k = variants$k[i],
                                                          supp_ar = variants$ar[i])))
    nm <- paste0("supp_", ifelse(variants$k[i] == "matched", "kM", "k5"), "_", variants$ar[i])
    out[nm] <- isef(f)
    out[paste0(nm, "_phi_t")] <- f$selection$ar_phi[["t"]]
  }
  out
}

settings <- expand.grid(trend = trends, N = c(100, 300, 500, 1000), stringsAsFactors = FALSE)
gp <- readRDS(gap_file)
results <- list()
for (i in seq_len(nrow(settings))) {
  tr <- settings$trend[i]; N <- settings$N[i]
  t0 <- Sys.time()
  res <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) one_rep(tr, N, r),
                                           mc.cores = cores))
  d <- data.frame(trend = tr, N = N, rep = seq_len(reps), res, check.names = FALSE)
  d <- merge(d, gp[gp$trend == tr & gp$N == N, c("trend", "N", "rep", "oracle", "R3")],
             by = c("trend", "N", "rep"))
  names(d)[names(d) == "R3"] <- "cv"
  results[[i]] <- d
  cols <- c("cv", "bspline", grep("^supp_k[5M]_[a-z]+$", names(d), value = TRUE))
  rat <- d[, cols] / d$oracle
  cat(sprintf("%-3s N=%4d | %s | %.0fs\n", tr, N,
              paste(sprintf("%s %.2f/%.2f", cols, apply(rat, 2, median),
                            apply(rat, 2, quantile, 0.9)), collapse = " | "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(do.call(rbind, results), out_file)
}
