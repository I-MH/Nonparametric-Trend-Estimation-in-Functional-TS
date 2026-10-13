test_that("simulated data have the documented structure", {
  sim <- simulate_fts(N = 50, m = 40, trend = "T4", seed = 1)
  expect_equal(dim(sim$Y), c(40L, 50L))
  expect_equal(sim$t, (1:50) / 50)
  expect_equal(sim$Y, sim$trend + sim$X)
  expect_equal(sim$trend[40, 50], 2 * (0.5 + 4)^2)
})

test_that("FAR(1) kernel is scaled to the requested Hilbert-Schmidt norm", {
  hs1 <- sqrt(integrate(function(u) exp(-u^2), 0, 1)$value^2)
  expect_equal(hs1, 0.746824, tolerance = 1e-6)  # constant in the supplementary code
  expect_equal(formals(simulate_fts)$norm, 0.5)
})

test_that("trend functions match the paper", {
  expect_equal(trend_function("T1")(0.5, 0.5), 16)
  expect_equal(trend_function("T3")(0.1, 1), 20)
  expect_equal(trend_function("T5")(0.25, 0), 28 * sin(0.25))
  expect_error(trend_function("T6"))
  f <- trend_function("T2b")
  expect_equal(f(0.25, 0.25), 25)
  # mean over s is zero; the deviation is nonlinear in t
  expect_equal(integrate(function(s) f(s, 0.3), 0, 1)$value, 0, tolerance = 1e-10)
  expect_gt(abs(f(0.25, 0.25) - (f(0.25, 0) + f(0.25, 0.5)) / 2), 1)
})

test_that("linear trend estimator and ISE behave as expected", {
  s <- seq(0, 1, length.out = 20)
  Y <- outer(s, 1:30, function(s, n) sin(s) + 0.1 * n * s)
  lt <- linear_trend(Y)
  expect_equal(lt$fitted, Y, tolerance = 1e-12)
  expect_equal(ise(Y, Y, s, seq(0, 1, length.out = 30)), 0)
  expect_equal(ise(matrix(1, 20, 30), matrix(0, 20, 30), s, seq(0, 1, length.out = 30)), 1)
})
