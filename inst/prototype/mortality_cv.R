# PROTOTYPE: leave-block-out CV (B1) on the French mortality data.
# Needs demography::fr.mort (no download). Run: Rscript mortality_cv.R

library(ftrend)
mort <- demography::fr.mort
ages <- 0:100; years <- 1816:2006
Ylog <- log(mort$rate$total[match(ages, mort$age), match(years, mort$year)])
train <- years <= 2002; test <- Ylog[, !train]
fpca_forecast <- function(Z, h, r = 4, d = 0) {
  mu <- rowMeans(Z); pc <- prcomp(t(Z - mu), center = FALSE)
  phi <- pc$rotation[, 1:r]; scores <- pc$x[, 1:r]
  orders <- expand.grid(p = 0:3, q = 0:2)
  fc <- sapply(1:r, function(j) {
    fits <- lapply(seq_len(nrow(orders)), function(i) tryCatch(suppressWarnings(arima(scores[, j], order = c(orders$p[i], d, orders$q[i]), include.mean = (d == 0))), error = function(e) NULL))
    best <- fits[[which.min(sapply(fits, function(f) if (is.null(f)) Inf else AIC(f)))]]
    predict(best, n.ahead = h)$pred })
  mu + phi %*% t(matrix(fc, h))
}
l1 <- function(A, B) sum(colSums(abs(A - B) * ftrend:::quad_weights(ages / 100)))
fb <- ftrend(Ylog[, train], s = ages, t = years[train], select = "reml-bspline")
fs <- ftrend(Ylog[, train], s = ages, t = years[train], select = "reml-supplement",
             control = ftrend_control(supp_k = "te-default"))
prep <- ftrend:::prepare_fts(fb$data, "trapezoid")
bs <- fb$basis_s; bt <- fb$basis_t
rule <- ftrend:::select_lambda_cv_rule(prep, bs, bt, one_se = TRUE)
gp <- rule$gap
cat(sprintf("N = %d, pilot gap l0 = %d, flat-top lag M = %d\n", prep$N, rule$l0, gp$M))
cat("rho(h), h = 0..10:", sprintf("%.2f", gp$rho[1:11]), "\n")
cat("tau(l), l = 0..10:", sprintf("%.3f", gp$tau[1:min(11, length(gp$tau))]), "\n")
cat(sprintf("l from dependence = %d, upper bound = %d, bound binds: %s, l used = %d\n",
            gp$l_dep, gp$l_max, gp$bound_binds, gp$l))
fits <- list(
  bspline = fb, supplement = fs,
  R1 = ftrend(Ylog[, train], s = ages, t = years[train], lambda = rule$lambda_pilot),
  R2 = ftrend(Ylog[, train], s = ages, t = years[train], lambda = rule$lambda_min),
  R3 = ftrend(Ylog[, train], s = ages, t = years[train], lambda = rule$lambda))
for (nm in names(fits)) {
  f <- fits[[nm]]
  yh <- forecast_trend(f, 4, s = ages) + fpca_forecast(residuals(f), 4)
  cat(sprintf("%-10s l1=%.2e l2=%.2e edf=%6.1f  mean resid in [%.2f, %.2f]  L1 forecast %.3f\n", nm,
              f$lambda[1], f$lambda[2], f$edf, min(rowMeans(residuals(f))),
              max(rowMeans(residuals(f))), l1(test, yh)))
}
