kron_setup <- function(m = 12, N = 9, k1 = 6, k2 = 7, seed = 1) {
  set.seed(seed)
  list(V = matrix(rnorm(m * k1), m, k1), Z = matrix(rnorm(N * k2), N, k2),
       Y = matrix(rnorm(m * N), m, N), w = runif(m))
}

test_that("Kronecker identities used for the common grid hold", {
  k <- kron_setup()
  ZV <- kronecker(k$Z, k$V)
  expect_equal(crossprod(ZV), kronecker(crossprod(k$Z), crossprod(k$V)))
  expect_equal(as.vector(crossprod(ZV, as.vector(k$Y))),
               as.vector(crossprod(k$V, k$Y) %*% k$Z))
  # weighted versions, W = I_N (x) diag(w)
  Wd <- rep(k$w, ncol(k$Y))
  expect_equal(crossprod(ZV, Wd * ZV),
               kronecker(crossprod(k$Z), crossprod(k$V, k$w * k$V)))
  expect_equal(as.vector(crossprod(ZV, Wd * as.vector(k$Y))),
               as.vector(crossprod(k$V, k$w * k$Y) %*% k$Z))
})

dense_solution <- function(prep, bs, bt, lambda) {
  V <- basis_eval(bs, prep$curves[[1]]$u)
  Z <- basis_eval(bt, prep$u_t)
  X <- kronecker(Z, V)
  w <- rep(prep$curves[[1]]$w, prep$N)
  y <- unlist(lapply(prep$curves, `[[`, "y"))
  S <- lambda[1] * kronecker(bt$J, bs$P) + lambda[2] * kronecker(bt$P, bs$J)
  matrix(solve(crossprod(X, w * X) + S, crossprod(X, w * y)), ncol(V), ncol(Z))
}

test_that("common-grid solver equals the dense solution of equation (8)", {
  sim <- simulate_fts(N = 40, m = 25, trend = "T2", seed = 3)
  x <- as_fts_irreg(sim$Y, s = sim$s)
  bs <- bspline_basis(8); bt <- bspline_basis(9)
  for (wt in c("trapezoid", "equal", "none")) {
    prep <- ftrend:::prepare_fts(x, wt)
    fit <- ftrend:::penalised_fit(prep, bs, bt, c(0.3, 2))
    expect_equal(fit$Theta, dense_solution(prep, bs, bt, c(0.3, 2)), tolerance = 1e-9)
  }
})

test_that("common-grid and curve-specific code paths agree when grids coincide", {
  sim <- simulate_fts(N = 30, m = 20, trend = "T5", seed = 4)
  x <- as_fts_irreg(sim$Y, s = sim$s)
  bs <- bspline_basis(10); bt <- bspline_basis(12)
  prep <- ftrend:::prepare_fts(x, "trapezoid")
  prep_irr <- prep
  prep_irr$common <- FALSE
  a <- ftrend:::penalised_fit(prep, bs, bt, c(0.05, 0.7))
  b <- ftrend:::penalised_fit(prep_irr, bs, bt, c(0.05, 0.7))
  expect_equal(a$Theta, b$Theta, tolerance = 1e-10)
  expect_equal(a$edf, b$edf, tolerance = 1e-10)
  # and the same through the user interface (data frame input)
  df <- as.data.frame(x)
  df <- df[nrow(df):1, ]
  f1 <- ftrend(sim$Y, s = sim$s, lambda = c(0.05, 0.7), k = c(10, 12))
  f2 <- ftrend(df, lambda = c(0.05, 0.7), k = c(10, 12))
  expect_true(f1$common_grid)
  expect_equal(coef(f1), coef(f2), tolerance = 1e-10)
})

test_that("the supplementary smooth.bibasis2 algebra is reproduced", {
  # smooth.bibasis2: Phi = kronecker(tbasismat, sbasismat), unit weights,
  # Bmat = Phi'Phi + lambdas * kronecker(inprod_t, penmat_s)
  #                + lambdat * kronecker(penmat_t, inprod_s), t = seq(0, 1, len = N)
  set.seed(5)
  m <- 21; N <- 30
  s <- seq(0, 1, length.out = m); tt <- seq(0, 1, length.out = N)
  Y <- outer(s, tt, function(s, t) 25 * t * sin(2 * pi * s)) + matrix(rnorm(m * N), m)
  bs <- bspline_basis(10); bt <- bspline_basis(15)
  Phi <- kronecker(basis_eval(bt, tt), basis_eval(bs, s))
  lam <- c(0.8, 3)
  Bmat <- crossprod(Phi) + lam[1] * kronecker(bt$J, bs$P) + lam[2] * kronecker(bt$P, bs$J)
  ref <- matrix(solve(Bmat, crossprod(Phi, as.vector(Y))), 10, 15)
  prep <- ftrend:::prepare_fts(as_fts_irreg(Y, s = s), "none")
  prep$u_t <- tt
  expect_equal(ftrend:::penalised_fit(prep, bs, bt, lam)$Theta, ref, tolerance = 1e-9)
})

test_that("mgcv paraPen uses the penalty unscaled", {
  set.seed(6)
  b <- bspline_basis(10)
  x <- seq(0, 1, length.out = 60)
  X <- basis_eval(b, x)
  y <- sin(3 * x) + rnorm(60, sd = 0.2)
  w <- ftrend:::quad_weights(x, "trapezoid")
  dat <- list(y = y, X = X, w = w)
  g <- mgcv::gam(y ~ X - 1, data = dat, weights = w,
                 paraPen = list(X = list(b$P)), sp = 0.01)
  direct <- solve(crossprod(X, w * X) + 0.01 * b$P, crossprod(X, w * y))
  expect_equal(unname(coef(g)), as.vector(direct), tolerance = 1e-7)
})

test_that("REML smoothing parameter matches an independent implementation", {
  reml_crit <- function(loglam, y, X, P, w) {
    lam <- exp(loglam)
    ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    r <- sum(ev > 1e-10 * max(ev)); M <- ncol(X) - r
    A <- crossprod(X, w * X) + lam * P
    beta <- solve(A, crossprod(X, w * y))
    D <- sum(w * (y - X %*% beta)^2) + lam * drop(t(beta) %*% P %*% beta)
    (length(y) - M) * log(D) + as.numeric(determinant(A)$modulus) - r * log(lam)
  }
  set.seed(7)
  b <- bspline_basis(10)
  x <- seq(0, 1, length.out = 80)
  X <- basis_eval(b, x)
  y <- cos(4 * x) + rnorm(80, sd = 0.3)
  for (w in list(rep(1, 80), ftrend:::quad_weights(x, "trapezoid"))) {
    lam_mgcv <- ftrend:::reml_lambda(y, X, b$P, w)
    opt <- optimize(reml_crit, c(-20, 20), y = y, X = X, P = b$P, w = w, tol = 1e-10)
    expect_equal(log(lam_mgcv), opt$minimum, tolerance = 1e-3)
  }
})

test_that("a trend in the spline space is reproduced exactly without noise", {
  s <- seq(0, 2, length.out = 30)
  tt <- 1:40
  Tst <- function(s, t) 1 + s - 0.5 * s^2 + (t / 40)^2 * s
  Y <- outer(s, tt, Tst)
  fit <- ftrend(Y, s = s, t = tt, lambda = c(1e-10, 1e-10))
  expect_equal(unname(fitted(fit)), Y, tolerance = 1e-6)
  expect_equal(unname(residuals(fit)), Y * 0, tolerance = 1e-6)
  expect_equal(unname(predict(fit, s = c(0.3, 1.7), t = c(5, 33))),
               outer(c(0.3, 1.7), c(5, 33), Tst), tolerance = 1e-6)
})

test_that("large lambdas give a trend that is linear in s and in t", {
  sim <- simulate_fts(N = 50, m = 30, trend = "T5", seed = 8)
  fit <- ftrend(sim$Y, s = sim$s, lambda = c(1e8, 1e8))
  M <- predict(fit, s = seq(0, 1, length.out = 21), t = fit$data$t[seq(1, 50, 5)])
  expect_lt(max(abs(diff(M, differences = 2))), 1e-4)
  expect_lt(max(abs(diff(t(M), differences = 2))), 1e-4)
})

test_that("ftrend selects lambda by cross-validation by default and records how", {
  sim <- simulate_fts(N = 100, m = 40, trend = "T2", seed = 9)
  fit <- ftrend(sim$Y, s = sim$s)
  expect_s3_class(fit, "ftrend")
  expect_true(all(fit$lambda > 0))
  expect_equal(unname(fit$lambda_method), rep("CV (leave-block-out)", 2))
  expect_gte(fit$selection$gap, 2L)
  # one-SE rule never decreases lambda2; lambda1 is the CV minimiser
  expect_gte(unname(fit$lambda[2]), fit$selection$lambda_cv_min[2])
  expect_equal(unname(fit$lambda[1]), fit$selection$lambda_cv_min[1])
  expect_output(print(summary(fit)), "Cross-validation: gap l =")
  ise_fit <- ise(sim$trend, predict(fit, s = sim$s, t = sim$t), sim$s, sim$t)
  expect_lt(ise_fit, 0.5)
  # one lambda supplied: only the other is selected
  fit2 <- ftrend(sim$Y, s = sim$s, lambda = c(NA, 5))
  expect_equal(unname(fit2$lambda_method), c("CV (leave-block-out)", "supplied"))
  expect_equal(unname(fit2$lambda[2]), 5)
  fit3 <- ftrend(sim$Y, s = sim$s, lambda = c(1e-3, NA))
  expect_equal(unname(fit3$lambda[1]), 1e-3)
  expect_error(ftrend(sim$Y, s = sim$s, reml = "bspline"), "replaced by 'select'")
})

test_that("marginal REML options record their scale and settings", {
  sim <- simulate_fts(N = 100, m = 40, trend = "T2", seed = 9)
  fb <- ftrend(sim$Y, s = sim$s, select = "reml-bspline")
  # lambda1 = N * lambda_s, lambda2 = lambda_t
  expect_equal(unname(fb$lambda), unname(c(100, 1) * fb$marginal_lambda))
  expect_equal(unname(fb$lambda_method), rep("REML (bspline)", 2))
  fs <- ftrend(sim$Y, s = sim$s, select = "reml-supplement")
  expect_match(fs$lambda_method[1], "REML \\(supplement, k te-default")
  fm <- ftrend(sim$Y, s = sim$s, select = "reml-supplement",
               control = ftrend_control(supp_k = "matched"))
  expect_match(fm$lambda_method[1], "k matched")
  # supplementary-code lambdas, unweighted scale -> trapezoid weights (sum 1, m = 40)
  expect_equal(unname(fs$lambda), unname(fs$marginal_lambda) / 40, tolerance = 1e-12)
  fn <- ftrend(sim$Y, s = sim$s, select = "reml-supplement", weights = "none",
               control = ftrend_control(supp_k = "te-default"))
  expect_equal(unname(fn$lambda), unname(fn$marginal_lambda), tolerance = 1e-12)
  # te-default reproduces the supplementary code exactly
  ybar <- rowMeans(sim$Y); s0 <- sim$s
  expect_equal(unname(fn$marginal_lambda[1]),
               unname(mgcv::gam(ybar ~ te(s0, bs = "cr"), method = "REML")$sp))
  fa <- ftrend(sim$Y, s = sim$s, select = "reml-supplement",
               control = ftrend_control(supp_ar = "time"))
  expect_true(is.na(fa$selection$ar_phi["s"]))
  expect_true(abs(fa$selection$ar_phi["t"]) < 1)
  expect_error(ftrend(sim$Y, s = sim$s, control = list()), "ftrend_control")
})

test_that("default selection works for curve-specific grids", {
  set.seed(10)
  sim <- simulate_fts(N = 60, m = 60, trend = "T2", seed = 10)
  df <- do.call(rbind, lapply(1:60, function(n) {
    keep <- sort(sample(60, sample(25:45, 1)))
    keep <- unique(c(1, keep, 60))
    data.frame(t = n, s = sim$s[keep], y = sim$Y[keep, n])
  }))
  fit <- ftrend(df)
  expect_false(fit$common_grid)
  expect_true(all(is.finite(fit$lambda)))
  r <- residuals(fit)
  expect_s3_class(r, "data.frame")
  expect_equal(r$t, df$t)
  expect_equal(r$s, df$s)
  expect_equal(fitted(fit)$y + r$y, df$y, tolerance = 1e-12)
  That <- predict(fit, s = sim$s, t = 1:60)
  expect_lt(ise(sim$trend, That, sim$s, sim$t), 1)
})
