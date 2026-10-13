# Leave-block-out cross-validation for (lambda1, lambda2) --------------------------
#
# Used by ftrend(select = "cv"), the default. Studies in inst/prototype/.
#
# Criterion (modified cross-validation for dependent data; Chu and Marron 1991,
# Burman, Chow and Nolan 1994):
#   CV_l(lambda) = sum_n sum_j w_nj {Y_n(s_nj) - T_hat^{(-n,l)}(s_nj, t_n)}^2,
# where T_hat^{(-n,l)} is the estimator (the same penalised criterion, eq. 4
# and eq. 6) fitted without the curves at times n' with |n' - n| <= l. When l
# exceeds the dependence range of {X_n}, E CV_l = const + E ISE(T_hat^{(-n,l)})
# up to negligible cross terms, so CV_l targets the integrated squared error of
# the whole surface. Only the selection of lambda changes; the estimator is
# the one of the paper.
#
# Fast computation. Write A = A_data + S_lambda for the penalised normal
# matrix and b for the right-hand side. Removing the block B = {n - l, ..., n + l}
# subtracts
#   D_B = sum_{n' in B} (z_{n'} z_{n'}^T) (x) G_{n'},  d_B = sum_{n' in B} z_{n'} (x) g_{n'},
# with z_n = eta(t_n), G_n = V_n^T W_n V_n, g_n = V_n^T W_n y_n. Because the
# B-splines eta are local, D_B is supported on the coefficient indices idx of
# the few time-basis functions that are non-zero on B. With U = I[, idx] and
# Q = U^T A^{-1} U, the Woodbury identity gives
#   (A - U D U^T)^{-1} = A^{-1} + A^{-1} U (I - D Q)^{-1} D U^T A^{-1},
# so, with w = theta[idx] - Q d_B,
#   theta^{(-B)}[idx] = w + Q (I - D Q)^{-1} D w.
# Only theta^{(-B)}[idx] is needed to predict curve n. The cost per held-out
# block is that of a (|idx| x |idx|) solve, with |idx| about 5 * k1.

cv_precompute <- function(prep, bs, bt, gap) {
  N <- prep$N
  k1 <- ncol(basis_gram(bs))
  k2 <- ncol(basis_gram(bt))
  Z <- basis_eval(bt, prep$u_t)
  G <- vector("list", N); g <- vector("list", N); cc <- numeric(N)
  Vcommon <- if (prep$common) basis_eval(bs, prep$curves[[1L]]$u) else NULL
  for (n in seq_len(N)) {
    cv <- prep$curves[[n]]
    Vn <- if (prep$common) Vcommon else basis_eval(bs, cv$u)
    G[[n]] <- crossprod(Vn, cv$w * Vn)
    g[[n]] <- as.vector(crossprod(Vn, cv$w * cv$y))
    cc[n] <- sum(cv$w * cv$y^2)
  }
  A_data <- matrix(0, k1 * k2, k1 * k2)
  Bm <- matrix(0, k1, k2)
  for (n in seq_len(N)) {
    A_data <- A_data + kronecker(tcrossprod(Z[n, ]), G[[n]])
    Bm <- Bm + tcrossprod(g[[n]], Z[n, ])
  }
  supp <- lapply(seq_len(N), function(n) which(Z[n, ] > 0))
  blocks <- lapply(seq_len(N), function(n) {
    B <- max(1L, n - gap):min(N, n + gap)
    S <- sort(unique(unlist(supp[B])))
    idx <- as.vector(outer(seq_len(k1), (S - 1L) * k1, `+`))
    D <- matrix(0, length(idx), length(idx))
    d <- numeric(length(idx))
    for (np in B) {
      D <- D + kronecker(tcrossprod(Z[np, S]), G[[np]])
      d <- d + kronecker(Z[np, S], g[[np]])
    }
    list(S = S, idx = idx, D = D, d = d, zS = Z[n, S])
  })
  list(N = N, k1 = k1, k2 = k2, A_data = A_data, b = as.vector(Bm),
       G = G, g = g, cc = cc, blocks = blocks, gap = gap,
       S1 = kronecker(basis_gram(bt), basis_penalty(bs)),
       S2 = kronecker(basis_penalty(bt), basis_gram(bs)))
}

# CV_l(lambda) for one value of lambda = c(lambda1, lambda2). With
# contributions = TRUE, returns the vector of per-curve terms
#   e_n(lambda) = sum_j w_nj {Y_n(s_nj) - T_hat^{(-n,l)}(s_nj, t_n)}^2.
cv_score <- function(lambda, pre, centres = seq_len(pre$N), contributions = FALSE) {
  bad <- if (contributions) rep(Inf, length(centres)) else Inf
  A <- pre$A_data + lambda[1L] * pre$S1 + lambda[2L] * pre$S2
  R <- tryCatch(chol(symmetrise(A)), error = function(e) NULL)
  if (is.null(R)) return(bad)
  Ainv <- chol2inv(R)
  theta <- as.vector(Ainv %*% pre$b)
  k1 <- pre$k1
  e <- numeric(length(centres))
  for (i in seq_along(centres)) {
    n <- centres[i]
    bl <- pre$blocks[[n]]
    Q <- Ainv[bl$idx, bl$idx, drop = FALSE]
    w <- theta[bl$idx] - as.vector(Q %*% bl$d)
    DQ <- bl$D %*% Q
    corr <- tryCatch(solve(diag(nrow(DQ)) - DQ, bl$D %*% w), error = function(e) NULL)
    if (is.null(corr)) return(bad)
    th <- w + as.vector(Q %*% corr)
    beta <- as.vector(matrix(th, k1) %*% bl$zS)          # Theta^{(-B)} eta(t_n)
    e[i] <- pre$cc[n] - 2 * sum(beta * pre$g[[n]]) + sum(beta * (pre$G[[n]] %*% beta))
  }
  if (contributions) e else sum(e)
}

# Choosing the gap l ------------------------------------------------------------------
#
# E CV_l = const + E ISE_l(lambda) + b_l(lambda), with dependence bias
#   b_l(lambda) = -2 sum_n sum_{|h| > l} tr(L_{n,n+h}(lambda) C_h),
# whose size relative to the variance term is governed by the tail fraction of
# the long-run trace variance,
#   tau(l) = sum_{|h| > l} gamma(h) / sum_h gamma(h),  gamma(h) = tr C_h.
# Plug-in rule: l = min{l : |tau_hat(l')| <= eps for all l' >= l}, with gamma
# estimated from pilot residuals and truncated at the flat-top lag of Politis
# (2003); bounded below by l_min and above so that 2l + 1 <= N / (2 df_t),
# where df_t is the time degrees of freedom of the pilot fit.

# Residual curves as a matrix on a common grid with quadrature weights
# (curve-specific grids: linear interpolation onto a reference grid).
resid_on_grid <- function(prep, resid, m_ref = NULL) {
  if (prep$common) {
    return(list(R = matrix(unlist(resid), ncol = prep$N), w = prep$curves[[1L]]$w))
  }
  lo <- max(vapply(prep$curves, function(cv) min(cv$u), numeric(1L)))
  hi <- min(vapply(prep$curves, function(cv) max(cv$u), numeric(1L)))
  if (is.null(m_ref))
    m_ref <- stats::median(vapply(prep$curves, function(cv) length(cv$u), numeric(1L)))
  u <- seq(lo, hi, length.out = max(5L, as.integer(m_ref)))
  R <- vapply(seq_len(prep$N), function(n)
    stats::approx(prep$curves[[n]]$u, resid[[n]], xout = u)$y, numeric(length(u)))
  list(R = R, w = quad_weights(u, "trapezoid"))
}

# Trace autocovariances gamma(h), h = 0, ..., H.
trace_acov <- function(R, w, H) {
  N <- ncol(R)
  Rc <- R - rowMeans(R)
  vapply(0:H, function(h) sum(w * Rc[, 1:(N - h), drop = FALSE] *
                                Rc[, (1 + h):N, drop = FALSE]) / N, numeric(1L))
}

# Flat-top lag (Politis 2003): q = smallest q with |rho(q + k)| < c sqrt(log10 N / N)
# for k = 1, ..., K; returns M = 2q.
flat_top_lag <- function(rho, N, c = 2, K = 5L) {
  thr <- c * sqrt(log10(N) / N)
  H <- length(rho) - 1L
  for (q in 0:(H - K)) if (all(abs(rho[q + 1L + seq_len(K)]) < thr)) return(2L * max(q, 1L))
  H
}

# Plug-in gap. Returns l and diagnostics.
cv_gap_plugin <- function(prep, resid, df_t, eps = 0.05, l_min = 2L, H = NULL) {
  N <- prep$N
  if (is.null(H)) H <- min(N - 2L, max(30L, as.integer(4 * sqrt(N))))
  rg <- resid_on_grid(prep, resid)
  g <- trace_acov(rg$R, rg$w, H)
  M <- min(flat_top_lag(g / g[1L], N), H)
  x <- (0:M) / M
  lam <- ifelse(x <= 0.5, 1, 2 * (1 - x))                   # flat-top taper
  gt <- g[1:(M + 1L)] * lam
  lr <- gt[1L] + 2 * sum(gt[-1L])
  tau <- vapply(0:M, function(l) if (l >= M) 0 else 2 * sum(gt[(l + 2L):(M + 1L)]) / lr,
                numeric(1L))
  ok <- rev(cumprod(rev(abs(tau) <= eps)))                   # all larger l also satisfy
  l_dep <- (0:M)[which(ok == 1)[1L]]
  l_max <- max(l_min, floor((N / (2 * df_t) - 1) / 2))
  l <- min(max(l_dep, l_min), l_max)
  list(l = as.integer(l), l_dep = l_dep, l_max = as.integer(l_max), M = M,
       bound_binds = l_dep > l_max, rho = g / g[1L], tau = tau, long_run = lr)
}

# Degrees of freedom in time of a fit with lambda2: by the exact decoupling,
# the s-mean of the trend is the penalised spline of z_n with lambda2, so
# df_t = tr{(Z^T W Z + lambda2 P_2)^{-1} Z^T W Z}, W = diag(W_n).
time_df <- function(prep, bt, lambda2) {
  Z <- basis_eval(bt, prep$u_t)
  Wn <- vapply(prep$curves, function(cv) sum(cv$w), numeric(1L))
  ZtZ <- crossprod(Z, Wn * Z)
  # df = sum_i 1 / (1 + lambda2 d_i), d_i eigenvalues of R^{-T} P_2 R^{-1}, ZtZ = R^T R
  Ri <- backsolve(chol(symmetrise(ZtZ)), diag(ncol(Z)))
  d <- eigen(symmetrise(crossprod(Ri, basis_penalty(bt) %*% Ri)), symmetric = TRUE,
             only.values = TRUE)$values
  d[d < 1e-9 * max(d)] <- 0                                 # null space of P_2
  sum(1 / (1 + lambda2 * d))
}

# Bartlett long-run variance of a series, lag M.
lrv_bartlett <- function(x, M) {
  x <- x - mean(x)
  N <- length(x)
  M <- min(M, N - 1L)
  g <- vapply(0:M, function(h) sum(x[1:(N - h)] * x[(1 + h):N]) / N, numeric(1L))
  g[1L] + 2 * sum((1 - seq_len(M) / (M + 1)) * g[-1L])
}

# One-standard-error rule along lambda2, lambda1 fixed: the largest lambda2
# whose CV exceeds the minimum by at most one standard error of the paired
# difference sum_n {e_n(lambda) - e_n(lambda_hat)}, the standard error
# allowing for serial dependence (Bartlett, lag M).
one_se_lambda2 <- function(pre, lambda_hat, M, cap2 = Inf, step = 0.25, max_up = 8) {
  e0 <- cv_score(lambda_hat, pre, contributions = TRUE)
  best <- lambda_hat[2L]
  for (k in seq_len(max_up / step)) {
    if (best >= cap2) break
    l2 <- min(lambda_hat[2L] * 10^(k * step), cap2)
    d <- cv_score(c(lambda_hat[1L], l2), pre, contributions = TRUE) - e0
    if (!all(is.finite(d))) break
    se <- sqrt(length(d) * max(lrv_bartlett(d, M), 0))
    if (sum(d) <= se) best <- l2 else break
  }
  c(lambda_hat[1L], best)
}

# Upper limits for (lambda1, lambda2). Beyond these the fit no longer changes
# (the penalised direction is already in the null space of the penalty to
# within 'tol' degrees of freedom), while the penalised normal matrix becomes
# ill-conditioned (for k2 = 15, lambda2 >= 1e10 already corrupts the solution).
#   df_t(lambda2) = sum_i 1 / (1 + lambda2 d_i),       d_i: eigenvalues for (Z^T W Z, P_2)
#   df_s(lambda1) = sum_i 1 / (1 + (lambda1 / N) e_i), e_i: eigenvalues for (Gbar, P_1)
# The cap is the smallest lambda with df - (null-space dimension) <= tol.
lambda_caps <- function(prep, bs, bt, tol = 1e-3) {
  gen_eig <- function(G, P) {
    Ri <- backsolve(chol(symmetrise(G)), diag(ncol(G)))
    d <- eigen(symmetrise(crossprod(Ri, P %*% Ri)), symmetric = TRUE, only.values = TRUE)$values
    d[d > 1e-9 * max(d)]
  }
  cap <- function(d, scale) {
    f <- function(lr) sum(1 / (1 + 10^lr * d / scale)) - tol
    10^stats::uniroot(f, c(-20, 30))$root
  }
  Z <- basis_eval(bt, prep$u_t)
  Wn <- vapply(prep$curves, function(cv) sum(cv$w), numeric(1L))
  Gbar <- Reduce(`+`, lapply(prep$curves, function(cv) {
    V <- basis_eval(bs, cv$u); crossprod(V, cv$w * V) })) / prep$N
  c(cap(gen_eig(Gbar, basis_penalty(bs)), prep$N),
    cap(gen_eig(crossprod(Z, Wn * Z), basis_penalty(bt)), 1))
}

# Select (lambda1, lambda2) by minimising CV_l over log-lambda: coarse grid
# (skipped if 'start' is given), then Nelder-Mead.
select_lambda_blockcv <- function(prep, bs, bt, gap, grid1 = seq(-10, 4, by = 2),
                                  grid2 = seq(-8, 4, by = 2), refine = TRUE,
                                  start = NULL, maxit = 200L, pre = NULL, caps = NULL,
                                  fixed = c(NA_real_, NA_real_)) {
  if (is.null(pre)) pre <- cv_precompute(prep, bs, bt, gap)
  if (is.null(caps)) caps <- lambda_caps(prep, bs, bt)
  lcap <- log10(caps)
  free <- is.na(fixed)
  if (!any(free)) stop("Nothing to select: both lambdas are fixed.", call. = FALSE)
  if (!all(free)) {
    # one smoothing parameter fixed: grid, then a one-dimensional search
    full <- function(lr) { v <- log10(fixed); v[free] <- lr; v }
    f1 <- function(lr) cv_score(10^pmin(full(lr), lcap), pre)
    g <- if (free[1L]) grid1 else grid2
    vals <- vapply(g, f1, numeric(1L))
    par <- g[which.min(vals)]; value <- min(vals)
    if (refine) {
      o <- stats::optimize(f1, c(par - 2, par + 2), tol = 1e-4)
      if (o$objective < value) { par <- o$minimum; value <- o$objective }
    }
    lam <- 10^full(par)
    lam[free] <- pmin(lam[free], caps[free])
    return(list(lambda = unname(lam), cv = value, gap = gap, pre = pre, caps = caps))
  }
  f <- function(lr) cv_score(10^pmin(lr, lcap), pre)
  if (is.null(start)) {
    gr <- expand.grid(l1 = grid1, l2 = grid2)
    vals <- apply(gr, 1L, f)
    start <- unlist(gr[which.min(vals), ])
    best <- list(par = start, value = min(vals))
  } else {
    start <- log10(start)
    best <- list(par = start, value = f(start))
  }
  if (refine) {
    o <- stats::optim(start, f, method = "Nelder-Mead",
                      control = list(reltol = 1e-6, maxit = maxit))
    if (o$value < best$value) best <- o
  }
  list(lambda = unname(10^pmin(best$par, lcap)), cv = best$value, gap = gap, pre = pre,
       caps = caps)
}

# Full rule: pilot with l0 = max(5, ceiling(N^(1/3))), plug-in gap, refit, and
# optionally the one-standard-error rule along lambda2.
select_lambda_cv_rule <- function(prep, bs, bt, eps = 0.05, one_se = TRUE,
                                  fixed = c(NA_real_, NA_real_)) {
  N <- prep$N
  l0 <- as.integer(max(5, ceiling(N^(1/3))))
  s0 <- select_lambda_blockcv(prep, bs, bt, l0, fixed = fixed)
  Th <- penalised_fit(prep, bs, bt, s0$lambda)$Theta
  Z <- basis_eval(bt, prep$u_t)
  resid <- lapply(seq_len(N), function(n) {
    cv <- prep$curves[[n]]
    cv$y - as.vector(basis_eval(bs, cv$u) %*% Th %*% Z[n, ])
  })
  gp <- cv_gap_plugin(prep, resid, time_df(prep, bt, s0$lambda[2L]), eps = eps)
  s1 <- if (gp$l == l0) s0 else
    select_lambda_blockcv(prep, bs, bt, gp$l, start = s0$lambda, maxit = 100L,
                          fixed = fixed, caps = s0$caps)
  lam <- s1$lambda
  if (one_se && is.na(fixed[2L]))
    lam <- one_se_lambda2(s1$pre, s1$lambda, gp$M, cap2 = s1$caps[2L])
  list(lambda = lam, lambda_min = s1$lambda, lambda_pilot = s0$lambda, l0 = l0,
       gap = gp, caps = s0$caps)
}

# Earlier rule used in the first prototype study (kept for reproducibility):
# l = smallest h >= 1 with |rho(h)| < 2 / sqrt(N), from pilot residuals.
cv_gap_from_residuals <- function(prep, resid, l_max = 20L) {
  N <- prep$N
  rg <- resid_on_grid(prep, resid)
  H <- min(l_max, N - 2L)
  rho <- trace_acov(rg$R, rg$w, H)
  rho <- rho / rho[1L]
  h <- which(abs(rho[-1L]) < 2 / sqrt(N))[1L]
  if (is.na(h)) H else h
}
