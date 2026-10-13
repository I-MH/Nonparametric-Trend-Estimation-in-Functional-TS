# Comparison with the supplementary code of Martínez-Hernández and Genton (2021)
#
# Reproduces the toy example of the Supporting Information (Section 1):
#   set.seed(1985); sim.FTS(N = 200, m = 50, norm.op = 0.5)
# and compares
#   (a) the original estimator: lambdas from gam(. ~ te(., bs = "cr")) on the
#       marginal means, then smooth.bibasis2() with t = seq(0, 1, len = N);
#   (b) ftrend with select = "reml-supplement" (supp_k = "te-default"),
#       weights = "none" and the same time
#       grid, which must reproduce (a) to rounding error;
#   (c) ftrend with its defaults (marginal REML on the B-spline bases with the
#       lambda conversion lambda1 = N * lambda_s, lambda2 = lambda_t,
#       trapezoidal weights, t_n = n/N).
#
# If the 'fda' package is installed and 'supp_dir' points to a copy of
# https://github.com/Martinez-Hernandez/Nonparametric-Trend-Estimation-in-Functional-TS
# at tag biometrics-2021 (git checkout biometrics-2021),
# the original smooth.bibasis2() is also run. Otherwise its algebra is
# emulated below with base R (the fda pieces it uses are eval.basis, inprod
# and eval.penalty, which equal basis_eval(), J and P of ftrend; see the unit
# tests).

library(ftrend)
library(mgcv)
supp_dir <- Sys.getenv("FTREND_SUPP_DIR", unset = NA)

# --- Data generator of the supplementary code (simFTS.R), verbatim ----------
TTrend <- function(s, t) 25 * t * sin(2 * pi * s)
brow <- function(TT, h, N, sigma = 1) {
  m <- length(seq(0, TT, h))
  Aux <- matrix(rnorm(n = N * (m - 1), sd = sqrt(h * sigma)), m - 1, N)
  rbind(rep(0, N), apply(Aux, 2, cumsum))
}
data.ind <- function(N = 100, m = 50, TT = 1) brow(TT, TT / m, N)
Tkernel <- function(t, u, norm.op = 0.5) norm.op / 0.746824 * exp(-(t^2 + u^2) / 2)
sim.FTS <- function(N = 100, m = 50, norm.op = 0.5, rangeval = c(0, 1), burn = 50) {
  h <- 1 / m
  s <- seq(rangeval[1], rangeval[2], by = h)
  KK <- outer(s, s, function(t, s) Tkernel(t, s, norm.op = norm.op))
  fw <- data.ind(N = N + burn + 1, m = m, TT = rangeval[2] - rangeval[1])
  X0 <- fw[, 1]; fw <- fw[, -1]
  D <- matrix(0, length(s), N + burn); D[, 1] <- X0
  for (i in 2:(N + burn)) D[, i] <- (t(KK) %*% D[, i - 1]) * h + fw[, i]
  D <- D[, (burn + 1):(N + burn)]
  t <- seq(0, 1, length.out = N)
  trendM <- outer(s, t, TTrend)
  list(s = s, DataM = D + trendM, TrendM = trendM, DataX = D)
}

set.seed(1985)
Data <- sim.FTS(N = 200, m = 50, norm.op = 0.5)
YY <- Data$DataM
m <- nrow(YY); N <- ncol(YY)                 # m = 51 (grid by = 1/50), N = 200
s <- seq(0, 1, length.out = m)
t <- seq(0, 1, length.out = N)               # (n - 1)/(N - 1), as in TrendEstimator.R

# --- (a) Original estimator --------------------------------------------------
auxs <- rowMeans(YY); auxt <- colMeans(YY)
lambdas <- gam(auxs ~ te(s, bs = "cr"), method = "REML")$sp
lambdat <- gam(auxt ~ te(t, bs = "cr"), method = "REML")$sp

bs <- bspline_basis(10); bt <- bspline_basis(15)
Phi <- kronecker(basis_eval(bt, t), basis_eval(bs, s))
Bmat <- crossprod(Phi) + lambdas * kronecker(bt$J, bs$P) + lambdat * kronecker(bt$P, bs$J)
edf_a <- sum(diag(solve(Bmat, crossprod(Phi))))
if (!is.na(supp_dir) && requireNamespace("fda", quietly = TRUE)) {
  library(fda)
  source(file.path(supp_dir, "smooth.bifda2.R"))
  fdPs <- fdPar(create.bspline.basis(c(0, 1), 10, 4), 2, lambda = lambdas)
  fdPt <- fdPar(create.bspline.basis(c(0, 1), 15, 4), 2, lambda = lambdat)
  Theta_a <- smooth.bibasis2(s, t, YY, fdPs, fdPt)$coef
  how_a <- "original smooth.bibasis2 (fda)"
} else {
  Theta_a <- matrix(solve(Bmat, crossprod(Phi, as.vector(YY))), 10, 15)
  how_a <- "emulated smooth.bibasis2 algebra"
}
That_a <- basis_eval(bs, s) %*% Theta_a %*% t(basis_eval(bt, t))

# --- (b) ftrend reproducing the supplementary code -----------------------------
x <- as_fts_irreg(YY, s = s)
prep <- ftrend:::prepare_fts(x, "none")
prep$u_t <- t                                  # time grid of the supplementary code
md <- ftrend:::marginal_data(prep, "none")
md$u_ref <- s
sel <- ftrend:::select_lambda_supplement(md, prep)
pf <- ftrend:::penalised_fit(prep, bs, bt, sel$lambda)
That_b <- basis_eval(bs, s) %*% pf$Theta %*% t(basis_eval(bt, t))

# --- (c) ftrend defaults ---------------------------------------------------------
fit_c <- ftrend(YY, s = s)                     # t_n = n/N
That_c <- predict(fit_c, s = s, t = fit_c$data$t)

# --- Report -------------------------------------------------------------------------
cat("Original estimator computed with:", how_a, "\n\n")
cat(sprintf("(a) original:            lambda1 = %.4g, lambda2 = %.4g\n", lambdas, lambdat))
cat(sprintf("(b) ftrend, supplement:  lambda1 = %.4g, lambda2 = %.4g\n", sel$lambda[1], sel$lambda[2]))
cat(sprintf("    max |Theta_a - Theta_b| = %.2e\n", max(abs(Theta_a - pf$Theta))))
cat(sprintf("(c) ftrend defaults:     lambda1 = %.4g, lambda2 = %.4g  (functional scale)\n",
            fit_c$lambda[1], fit_c$lambda[2]))
cat(sprintf("    on the unweighted scale of eq. (8): lambda1 = %.4g, lambda2 = %.4g\n",
            fit_c$lambda[1] * (m - 1), fit_c$lambda[2] * (m - 1)))
cat(sprintf("    effective degrees of freedom: original %.1f, default %.1f\n", edf_a, fit_c$edf))
cat("\nISE_T^2 against the true trend (trapezoidal rule on the data grid):\n")
cat(sprintf("  (a) original  %.4f\n", ise(Data$TrendM, That_a, s, t)))
cat(sprintf("  (b) ftrend supplement mode  %.4f\n", ise(Data$TrendM, That_b, s, t)))
# curve n was generated with T(s, (n - 1)/(N - 1)); ftrend places it at n/N
cat(sprintf("  (c) ftrend defaults  %.4f\n", ise(Data$TrendM, That_c, s, (1:N) / N)))
