# Simulation tools (MHG, Section 4) ---------------------------------------------------

#' Functional trends used in the simulation study
#'
#' The trend surfaces \eqn{T_1}--\eqn{T_5} of Section 4.2 in
#' Martínez-Hernández and Genton (2021), as functions of
#' \eqn{(s, t) \in [0, 1]^2}:
#' \itemize{
#'   \item \code{T1}: \eqn{2s + 30t};
#'   \item \code{T2}: \eqn{25 t \sin(2\pi s)};
#'   \item \code{T3}: \eqn{20t^2 - 5t + 5};
#'   \item \code{T4}: \eqn{2(0.5s + 4t)^2};
#'   \item \code{T5}: \eqn{28\sin(2\pi t + s)};
#'   \item \code{T2b}: \eqn{25 \sin(2\pi t) \sin(2\pi s)} (not in the paper).
#' }
#' The surface \eqn{T_6} of the paper is not part of this simulation design.
#' \code{T2b} is a variant of \code{T2} for which the mean over \eqn{s},
#' \eqn{\int_0^1 T(s, t)\,ds}, is identically zero and the deviation
#' \eqn{T(s, t) - \int_0^1 T(u, t)\,du} is nonlinear in \eqn{t}; in
#' \code{T1}--\code{T4} that deviation is zero or linear in \eqn{t}.
#' @param name One of \code{"T1"}, ..., \code{"T5"}, \code{"T2b"}.
#' @return A vectorised function of \code{(s, t)}.
#' @export
trend_function <- function(name = c("T1", "T2", "T3", "T4", "T5", "T2b")) {
  name <- match.arg(name)
  switch(name,
    T1 = function(s, t) 2 * s + 30 * t,
    T2 = function(s, t) 25 * t * sin(2 * pi * s),
    T3 = function(s, t) 20 * t^2 - 5 * t + 5 + 0 * s,
    T4 = function(s, t) 2 * (0.5 * s + 4 * t)^2,
    T5 = function(s, t) 28 * sin(2 * pi * t + s),
    T2b = function(s, t) 25 * sin(2 * pi * t) * sin(2 * pi * s))
}

#' Simulate a functional time series with a functional trend
#'
#' Simulates \eqn{Y_n(s) = T(s, n/N) + X_n(s)} on an equispaced grid of
#' \eqn{m} points in \eqn{[0, 1]}, where \eqn{X_n} follows the FAR(1) model
#' \eqn{X_n(s) = \int_0^1 \beta(u, s) X_{n-1}(u)\,du + W_n(s)} with
#' \eqn{\beta(u, v) = C \exp\{-(u^2 + v^2)/2\}} and independent standard
#' Brownian motions \eqn{W_n}.
#'
#' @param N Number of curves.
#' @param m Number of grid points.
#' @param trend A trend name (see \code{\link{trend_function}}) or a function
#'   of \code{(s, t)}.
#' @param norm Hilbert-Schmidt norm
#'   \eqn{(\int\!\int \beta^2)^{1/2}} of the FAR(1) kernel; \eqn{C} is chosen
#'   to attain it. The simulation design uses 0.5 (the default), as in the
#'   supplementary code of the paper. The Hilbert-Schmidt norm bounds the operator norm from
#'   above, so the process is stationary when \code{norm < 1}.
#' @param burn Burn-in length.
#' @param seed Optional random seed.
#' @return A list with \code{Y} (\eqn{m \times N} data), \code{trend}
#'   (\eqn{m \times N} true trend), \code{X} (stationary component),
#'   \code{s}, \code{t} (\eqn{= n/N}) and the trend function \code{T}.
#' @details The integral in the FAR(1) recursion is approximated by the
#'   trapezoidal rule on the grid. The supplementary code of the paper uses a
#'   left Riemann sum on \eqn{m + 1} points and times
#'   \eqn{(n-1)/(N-1)}; see \code{inst/comparison} for an exact replication.
#' @export
simulate_fts <- function(N = 100, m = 50, trend = "T1", norm = 0.5, burn = 50,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tfun <- if (is.function(trend)) trend else trend_function(trend)
  s <- seq(0, 1, length.out = m)
  h <- s[2L] - s[1L]
  # (int_0^1 exp(-u^2) du)^2 = int int exp(-(u^2 + v^2)); its square root is
  # the Hilbert-Schmidt norm of exp{-(u^2 + v^2)/2}.
  hs1 <- sqrt(pi) / 2 * (2 * stats::pnorm(sqrt(2)) - 1)
  C <- norm / hs1
  K <- C * exp(-outer(s^2, s^2, `+`) / 2)        # K[u, v] = beta(u, v)
  wq <- quad_weights(s, "trapezoid")
  Op <- t(K * wq)                                 # (Op x)(v) = sum_u w_u beta(u, v) x(u)
  W <- matrix(stats::rnorm((N + burn) * (m - 1L), sd = sqrt(h)), m - 1L)
  W <- rbind(0, apply(W, 2L, cumsum))             # Brownian motions, W(0) = 0
  X <- matrix(0, m, N + burn)
  X[, 1L] <- W[, 1L]
  for (n in 2L:(N + burn)) X[, n] <- Op %*% X[, n - 1L] + W[, n]
  X <- X[, (burn + 1L):(burn + N), drop = FALSE]
  t <- seq_len(N) / N
  TT <- outer(s, t, Tfun)
  list(Y = TT + X, trend = TT, X = X, s = s, t = t, T = Tfun)
}

#' Linear functional trend estimator
#'
#' The least-squares linear trend of Kokoszka and Young (2017),
#' \eqn{Y_n(s) = \mu(s) + n f(s) + X_n(s)}, fitted pointwise in \eqn{s}:
#' \eqn{\hat f(s) = \sum_n (n - \bar n) Y_n(s) / \sum_n (n - \bar n)^2},
#' \eqn{\hat\mu(s) = \bar Y(s) - \hat f(s)\bar n}, \eqn{\bar n = (N+1)/2}.
#'
#' @param Y An \eqn{m \times N} data matrix (common grid).
#' @return A list with \code{mu}, \code{f}, the fitted matrix \code{fitted}
#'   and a function \code{trend(u)} returning \eqn{\hat\mu + N u \hat f} on
#'   the grid for rescaled times \eqn{u = n/N}.
#' @references
#' Kokoszka, Piotr, and Gabriel Young. 2017. "Testing Trend Stationarity of
#' Functional Time Series with Application to Yield and Daily Price Curves."
#' \emph{Statistics and Its Interface} 10 (1): 81–92.
#' @export
linear_trend <- function(Y) {
  N <- ncol(Y)
  n <- seq_len(N)
  nb <- (N + 1) / 2
  f <- as.vector(Y %*% (n - nb)) / sum((n - nb)^2)
  mu <- rowMeans(Y) - f * nb
  list(mu = mu, f = f, fitted = mu + outer(f, n),
       trend = function(u) mu + outer(f, N * u))
}

#' Integrated squared error of a trend estimate
#'
#' \eqn{\mathrm{ISE}_T^2 = \int_0^1\int_0^1 \{T(s, t) - \hat T(s, t)\}^2
#' \,ds\,dt}, computed by the trapezoidal rule in both directions.
#'
#' @param Ttrue,That \eqn{m_s \times m_t} matrices of the true and estimated
#'   trend on the grid \code{s} \eqn{\times} \code{t}.
#' @param s,t Grids in \eqn{[0, 1]}.
#' @return The value of \eqn{\mathrm{ISE}_T^2}.
#' @export
ise <- function(Ttrue, That, s, t) {
  E2 <- (Ttrue - That)^2
  as.numeric(crossprod(quad_weights(s, "trapezoid"), E2 %*% quad_weights(t, "trapezoid")))
}
