test_that("Gauss-Legendre rule is exact for polynomials of degree <= 2n - 1", {
  for (n in 1:6) {
    gl <- ftrend:::gauss_legendre(n)
    for (p in 0:(2 * n - 1)) {
      exact <- if (p %% 2 == 1) 0 else 2 / (p + 1)
      expect_equal(sum(gl$weights * gl$nodes^p), exact, tolerance = 1e-12)
    }
  }
})

test_that("B-spline basis has the right size and forms a partition of unity", {
  b <- bspline_basis(k = 10, order = 4)
  x <- seq(0, 1, length.out = 37)
  B <- basis_eval(b, x)
  expect_equal(dim(B), c(37L, 10L))
  expect_equal(rowSums(B), rep(1, 37), tolerance = 1e-12)
  # equidistant interior knots (Assumption 3)
  expect_equal(diff(unique(b$knots)), rep(1 / 7, 7), tolerance = 1e-14)
})

test_that("Gram and penalty matrices agree with adaptive quadrature", {
  b <- bspline_basis(k = 8, order = 4, penalty_order = 2)
  f <- function(i, j, d) function(x)
    basis_eval(b, x, deriv = d)[, i] * basis_eval(b, x, deriv = d)[, j]
  for (ij in list(c(1, 1), c(2, 3), c(4, 4), c(5, 8), c(8, 8))) {
    J_ij <- integrate(f(ij[1], ij[2], 0), 0, 1, rel.tol = 1e-12, subdivisions = 500)$value
    P_ij <- integrate(f(ij[1], ij[2], 2), 0, 1, rel.tol = 1e-12, subdivisions = 500)$value
    expect_equal(b$J[ij[1], ij[2]], J_ij, tolerance = 1e-8)
    expect_equal(b$P[ij[1], ij[2]], P_ij, tolerance = 1e-8)
  }
  expect_true(isSymmetric(b$J))
  expect_true(isSymmetric(b$P))
  # J integrates the partition of unity: 1^T J 1 = int 1 = 1
  expect_equal(sum(b$J), 1, tolerance = 1e-12)
})

test_that("second-derivative penalty annihilates linear functions", {
  b <- bspline_basis(k = 12, order = 4, penalty_order = 2)
  # B-spline coefficients of f(x) = x are the Greville abscissae
  kn <- b$knots
  grev <- sapply(1:12, function(i) mean(kn[(i + 1):(i + 3)]))
  x <- seq(0, 1, length.out = 11)
  expect_equal(as.vector(basis_eval(b, x) %*% grev), x, tolerance = 1e-12)
  expect_equal(as.vector(b$P %*% rep(1, 12)), rep(0, 12), tolerance = 1e-8)
  expect_equal(as.vector(b$P %*% grev), rep(0, 12), tolerance = 1e-8)
  ev <- eigen(b$P, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(sum(ev > 1e-8 * max(ev)), 10L)   # rank k - 2
})

test_that("invalid basis arguments fail informatively", {
  expect_error(bspline_basis(k = 3, order = 4), "k >= order")
  expect_error(bspline_basis(k = 10, order = 4, penalty_order = 4), "penalty_order")
  expect_error(basis_eval(bspline_basis(), c(-0.1, 0.5)), "\\[0, 1\\]")
})
