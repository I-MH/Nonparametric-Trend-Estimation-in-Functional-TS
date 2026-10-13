test_that("fast leave-block-out CV equals brute-force refits", {
  sim <- simulate_fts(N = 40, m = 15, trend = "T2b", seed = 21)
  bs <- bspline_basis(6); bt <- bspline_basis(8)
  lam <- c(1e-3, 1e-2)
  for (common in c(TRUE, FALSE)) {
    prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
    prep$common <- common
    for (gap in c(0L, 3L)) {
      pre <- ftrend:::cv_precompute(prep, bs, bt, gap)
      fast <- ftrend:::cv_score(lam, pre)
      brute <- 0
      for (n in seq_len(prep$N)) {
        keep <- setdiff(seq_len(prep$N), max(1, n - gap):min(prep$N, n + gap))
        pk <- prep
        pk$curves <- prep$curves[keep]; pk$u_t <- prep$u_t[keep]; pk$N <- length(keep)
        Th <- ftrend:::penalised_fit(pk, bs, bt, lam)$Theta
        cv <- prep$curves[[n]]
        fit_n <- basis_eval(bs, cv$u) %*% Th %*% t(basis_eval(bt, prep$u_t[n]))
        brute <- brute + sum(cv$w * (cv$y - fit_n)^2)
      }
      expect_equal(fast, brute, tolerance = 1e-8)
    }
  }
})

test_that("block CV selection returns finite lambdas", {
  sim <- simulate_fts(N = 60, m = 20, trend = "T5", seed = 22)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  sel <- ftrend:::select_lambda_blockcv(prep, bspline_basis(10), bspline_basis(15), gap = 3L,
                                        refine = FALSE)
  expect_true(all(is.finite(sel$lambda) & sel$lambda > 0))
})

test_that("per-curve CV contributions sum to the criterion", {
  sim <- simulate_fts(N = 50, m = 20, trend = "T1", seed = 23)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  pre <- ftrend:::cv_precompute(prep, bspline_basis(8), bspline_basis(10), 3L)
  e <- ftrend:::cv_score(c(1e-2, 1e-1), pre, contributions = TRUE)
  expect_length(e, 50)
  expect_equal(sum(e), ftrend:::cv_score(c(1e-2, 1e-1), pre))
})

test_that("time degrees of freedom have the right limits", {
  sim <- simulate_fts(N = 80, m = 20, trend = "T1", seed = 24)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  bt <- bspline_basis(15)
  expect_equal(ftrend:::time_df(prep, bt, 1e12), 2, tolerance = 1e-3)
  expect_equal(ftrend:::time_df(prep, bt, 1e-14), 15, tolerance = 1e-3)
})

test_that("plug-in gap recovers the dependence range of the FAR(1) design", {
  sim <- simulate_fts(N = 2000, m = 30, trend = "T1", seed = 25)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$X, s = sim$s), "trapezoid")
  gp <- ftrend:::cv_gap_plugin(prep, lapply(seq_len(2000), function(n) sim$X[, n]),
                               df_t = 5)
  # tau(l) for this design is about 0.08 at l = 3 and 0.05 at l = 4
  expect_gte(gp$l, 3L)
  expect_lte(gp$l, 8L)
  expect_false(gp$bound_binds)
  # white noise gives the minimum gap
  wn <- lapply(1:2000, function(n) rnorm(30))
  expect_equal(ftrend:::cv_gap_plugin(prep, wn, df_t = 5)$l, 2L)
})

test_that("one-SE rule never decreases lambda2", {
  sim <- simulate_fts(N = 60, m = 20, trend = "T1", seed = 26)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  sel <- ftrend:::select_lambda_blockcv(prep, bspline_basis(10), bspline_basis(15), 5L,
                                        refine = FALSE)
  lam <- ftrend:::one_se_lambda2(sel$pre, sel$lambda, M = 4L)
  expect_equal(lam[1], sel$lambda[1])
  expect_gte(lam[2], sel$lambda[2])
})

test_that("lambda caps sit where the fit has stopped changing", {
  sim <- simulate_fts(N = 300, m = 50, trend = "T1", seed = 27)
  prep <- ftrend:::prepare_fts(as_fts_irreg(sim$Y, s = sim$s), "trapezoid")
  bs <- bspline_basis(10); bt <- bspline_basis(15)
  caps <- ftrend:::lambda_caps(prep, bs, bt)
  expect_lt(caps[2], 1e9)
  expect_equal(ftrend:::time_df(prep, bt, caps[2]), 2 + 1e-3, tolerance = 1e-6)
  V <- basis_eval(bs, sim$s); Z <- basis_eval(bt, prep$u_t)
  f <- function(l2) V %*% ftrend:::penalised_fit(prep, bs, bt, c(1, l2))$Theta %*% t(Z)
  expect_lt(max(abs(f(caps[2]) - f(10 * caps[2]))), 0.05)
})
