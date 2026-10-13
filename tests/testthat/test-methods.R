sim <- simulate_fts(N = 80, m = 30, trend = "T3", seed = 11)
years <- 1901:1980
fit <- ftrend(sim$Y, s = sim$s * 100, t = years, domain = c(0, 100))

test_that("output keeps the input format", {
  expect_true(is.matrix(fitted(fit)))
  expect_equal(dim(residuals(fit)), dim(sim$Y))
  expect_equal(unname(fitted(fit) + residuals(fit)), sim$Y, tolerance = 1e-12)
  expect_equal(dim(coef(fit)), c(10L, 15L))
  x <- as_fts_irreg(sim$Y, s = sim$s)
  fx <- ftrend(x, lambda = c(1, 1))
  expect_s3_class(residuals(fx), "fts_irreg")
  expect_equal(residuals(fx)$curves[[5]]$s, sim$s)
})

test_that("predict reproduces fitted values and handles extrapolation", {
  P <- predict(fit, s = sim$s * 100, t = years)
  expect_equal(unname(P), unname(fitted(fit)), tolerance = 1e-10)
  v <- predict(fit, s = c(10, 50), t = c(1920, 1960), grid = FALSE)
  expect_equal(unname(v), unname(diag(predict(fit, c(10, 50), c(1920, 1960)))))
  expect_warning(predict(fit, s = 50, t = 1985), "outside the observed time range")
  expect_warning(predict(fit, s = 50, t = 1900), "outside the observed time range")
  expect_error(predict(fit, s = 101, t = 1950), "domain")
})

test_that("forecast_trend is the linear Taylor forecast of equation (5)", {
  s0 <- c(0, 25, 50, 75, 100)
  fc <- forecast_trend(fit, h = 4, s = s0)
  expect_equal(dim(fc), c(5L, 4L))
  # derivative at u = 1 from the basis, step 1/N on the rescaled time axis
  bt <- fit$basis_t
  dT <- basis_eval(fit$basis_s, s0 / 100) %*% coef(fit) %*% t(basis_eval(bt, 1, deriv = 1))
  T1 <- predict(fit, s = s0, t = 1980)
  for (h in 1:4) expect_equal(unname(fc[, h]), as.vector(T1 + h / 80 * dT), tolerance = 1e-10)
  expect_warning(p <- predict(fit, s = s0, t = 1982))
  expect_equal(unname(p[, 1]), unname(fc[, 2]), tolerance = 1e-10)
})

test_that("print, summary and plot methods run", {
  expect_output(print(fit), "lambda1 = ")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.ftrend")
  expect_gt(sm$r_squared, 0.5)
  expect_output(print(sm), "Effective degrees of freedom")
  pdf(NULL)
  on.exit(dev.off())
  for (ty in c("surface", "curves", "data", "residuals"))
    expect_invisible(plot(fit, type = ty))
})
