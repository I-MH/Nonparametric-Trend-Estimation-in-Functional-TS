# Basis component ------------------------------------------------------------
#
# A *basis object* bundles everything the estimator needs from one marginal
# basis: an evaluation method, the Gram matrix J = int b b^T and the penalty
# matrix P = int b^(q) b^(q)^T. The estimator only touches a basis through the
# generics basis_eval(), basis_gram() and basis_penalty(), so that the domain
# basis nu can later be replaced (Phase 2: d-dimensional domains) without
# changing the solver.
#
# All bases live on the rescaled interval [0, 1]; mapping data to [0, 1] is
# done by the caller (see rescale_s() and rescale_t() in data.R).

#' Gauss-Legendre quadrature nodes and weights on [-1, 1]
#'
#' Computed by the Golub-Welsch algorithm (eigen-decomposition of the Jacobi
#' matrix of the Legendre polynomials). A rule with \code{n} nodes integrates
#' polynomials of degree \eqn{\le 2n - 1} exactly.
#'
#' @param n Number of nodes (positive integer).
#' @return A list with components \code{nodes} and \code{weights}.
#' @keywords internal
#' @noRd
gauss_legendre <- function(n) {
  n <- as.integer(n)
  if (n < 1L) stop("'n' must be a positive integer.", call. = FALSE)
  if (n == 1L) return(list(nodes = 0, weights = 2))
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(nodes = e$values[o], weights = 2 * e$vectors[1L, o]^2)
}

#' Quadrature rule over the knot intervals of a spline basis
#'
#' Places an \code{n}-point Gauss-Legendre rule on every non-degenerate knot
#' interval. Because a B-spline is a polynomial of degree \code{order - 1}
#' on each interval, products of two B-splines (or of their derivatives) are
#' integrated exactly when \code{n >= order}.
#' @keywords internal
#' @noRd
knot_quadrature <- function(knots, n) {
  gl <- gauss_legendre(n)
  u <- unique(knots)
  a <- u[-length(u)]
  b <- u[-1L]
  half <- (b - a) / 2
  mid <- (a + b) / 2
  list(
    nodes = as.vector(outer(gl$nodes, half) + rep(mid, each = n)),
    weights = as.vector(outer(gl$weights, half))
  )
}

#' B-spline basis with equidistant knots and a derivative penalty
#'
#' Constructs a B-spline basis on \eqn{[0, 1]} with \code{k} functions of
#' order \code{order} (degree \code{order - 1}) and equidistant interior knots
#' (Assumption 3 of Martínez-Hernández and Genton 2021). The Gram matrix
#' \eqn{J = \int_0^1 b(x) b(x)^\top dx} and the penalty matrix
#' \eqn{P = \int_0^1 b^{(q)}(x) b^{(q)}(x)^\top dx}, with
#' \eqn{q} = \code{penalty_order}, are computed exactly by Gauss-Legendre
#' quadrature on each knot interval.
#'
#' @param k Number of basis functions; must satisfy \code{k >= order}.
#' @param order Spline order (4 = cubic).
#' @param penalty_order Order \eqn{q} of the derivative in the roughness
#'   penalty; must satisfy \code{1 <= penalty_order < order}.
#' @return An object of class \code{c("ft_bspline", "ft_basis")}.
#' @references
#' Martínez-Hernández, Israel, and Marc G. Genton. 2021. "Nonparametric Trend
#' Estimation in Functional Time Series with Application to Annual Mortality
#' Rates." \emph{Biometrics} 77 (3): 866–78.
#' @examples
#' b <- bspline_basis(k = 10)
#' dim(basis_eval(b, seq(0, 1, length.out = 5)))
#' @export
bspline_basis <- function(k = 10, order = 4, penalty_order = 2) {
  k <- as.integer(k)
  order <- as.integer(order)
  penalty_order <- as.integer(penalty_order)
  if (length(order) != 1L || is.na(order) || order < 1L)
    stop("'order' must be a single positive integer.", call. = FALSE)
  if (length(k) != 1L || is.na(k) || k < order)
    stop("'k' must be a single integer with k >= order (", order, ").",
         call. = FALSE)
  if (length(penalty_order) != 1L || is.na(penalty_order) ||
      penalty_order < 1L || penalty_order >= order)
    stop("'penalty_order' must be an integer with 1 <= penalty_order < order.",
         call. = FALSE)
  n_interior <- k - order
  interior <- if (n_interior > 0L) seq_len(n_interior) / (n_interior + 1) else numeric(0)
  knots <- c(rep(0, order), interior, rep(1, order))

  basis <- structure(
    list(type = "bspline", k = k, order = order, penalty_order = penalty_order,
         knots = knots, J = NULL, P = NULL),
    class = c("ft_bspline", "ft_basis")
  )
  q <- knot_quadrature(knots, order)
  B0 <- basis_eval(basis, q$nodes)
  Bq <- basis_eval(basis, q$nodes, deriv = penalty_order)
  basis$J <- symmetrise(crossprod(B0, q$weights * B0))
  basis$P <- symmetrise(crossprod(Bq, q$weights * Bq))
  basis
}

#' Evaluate a basis, or one of its derivatives, at given points
#'
#' @param basis A basis object (class \code{ft_basis}).
#' @param x Numeric vector of points in \eqn{[0, 1]}.
#' @param deriv Order of the derivative (0 = the functions themselves).
#' @return A \code{length(x)} by \code{k} matrix.
#' @export
basis_eval <- function(basis, x, deriv = 0) UseMethod("basis_eval")

#' @export
basis_eval.ft_bspline <- function(basis, x, deriv = 0) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < 0 - 1e-12 | x > 1 + 1e-12))
    stop("Basis evaluation points must lie in [0, 1].", call. = FALSE)
  x <- pmin(pmax(x, 0), 1)
  splines::splineDesign(basis$knots, x, ord = basis$order,
                        derivs = rep(as.integer(deriv), length(x)),
                        outer.ok = FALSE)
}

#' @rdname basis_eval
#' @export
basis_gram <- function(basis) UseMethod("basis_gram")

#' @export
basis_gram.ft_basis <- function(basis) basis$J

#' @rdname basis_eval
#' @export
basis_penalty <- function(basis) UseMethod("basis_penalty")

#' @export
basis_penalty.ft_basis <- function(basis) basis$P

#' @export
print.ft_bspline <- function(x, ...) {
  cat("B-spline basis on [0, 1]: k = ", x$k, ", order = ", x$order,
      ", penalty on derivative ", x$penalty_order, "\n", sep = "")
  invisible(x)
}

symmetrise <- function(A) (A + t(A)) / 2
