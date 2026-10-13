# Estimation of the functional trend -------------------------------------------
#
# Notation (Martínez-Hernández and Genton 2021, hereafter MHG):
#   model (1)    Y_n(s) = T(s, n/N) + X_n(s)
#   eq.   (2)    T(s, t) = nu(s)^T Theta eta(t)
#   eq.   (3)    P(T) = lambda1 int int (d^2 T/ds^2)^2 + lambda2 int int (d^2 T/dt^2)^2
#   eq.   (6)    P(T) = vec(Theta)^T {lambda1 J_eta (x) P_1 + lambda2 P_2 (x) J_nu} vec(Theta)
#   eq.   (4)    Theta_hat = argmin sum_n int_D {Y_n(s) - T(s, n/N)}^2 ds + P(T)
#   eq.   (8)    [(Z (x) V)^T (Z (x) V) + penalty] vec(Theta) = (Z (x) V)^T vec(Y)
#
# Both D and the time axis are rescaled to [0, 1]. The integral over D in (4)
# is approximated by quadrature weights w_nj on the observed points of curve n,
#   L(Theta) = sum_n sum_j w_nj {Y_n(s_nj) - nu(s_nj)^T Theta eta(t_n)}^2,
# so that with w_nj = 1 (weights = "none") and a common grid, L is exactly the
# Frobenius-norm criterion behind (8).

# Prepare per-curve quantities on the rescaled scale.
prepare_fts <- function(x, weights) {
  tm <- time_map(x$t)
  curves <- lapply(x$curves, function(cv) {
    u <- rescale_s(cv$s, x$domain)
    list(u = u, y = cv$y, w = quad_weights(u, weights))
  })
  list(curves = curves, u_t = rescale_t(x$t, tm), time_map = tm,
       common = is_common_grid(x), N = length(curves))
}

# Penalised least squares for given (lambda1, lambda2) ----------------------------

#' @keywords internal
#' @noRd
penalised_fit <- function(prep, bs, bt, lambda) {
  k1 <- ncol(basis_gram(bs))
  k2 <- ncol(basis_gram(bt))
  Z <- basis_eval(bt, prep$u_t)                       # N x k2, Z = {eta(t)}^T
  J_nu <- basis_gram(bs)
  J_eta <- basis_gram(bt)
  P1 <- basis_penalty(bs)
  P2 <- basis_penalty(bt)

  if (prep$common) {
    # Common grid: (Z (x) V)^T W (Z (x) V) = Z^T Z (x) V^T W V and
    # (Z (x) V)^T W vec(Y) = vec(V^T W Y Z), so the mN x k1k2 design matrix
    # of (8) is never formed.
    cv1 <- prep$curves[[1L]]
    V <- basis_eval(bs, cv1$u)                        # m x k1, V = {nu(s)}^T
    G <- crossprod(V, cv1$w * V)                      # V^T W V
    Ymat <- vapply(prep$curves, `[[`, numeric(length(cv1$y)), "y")
    Ymat <- matrix(Ymat, ncol = prep$N)               # m x N
    A_data <- kronecker(crossprod(Z), G)              # eq. (8), first term
    b <- as.vector(crossprod(V, cv1$w * Ymat) %*% Z)  # vec(V^T W Y Z)
  } else {
    # Curve-specific grids (MHG, Remark 2): the design matrix is the row-wise
    # Kronecker (face-splitting) product with rows eta(t_n)^T (x) nu(s_nj)^T.
    # Grouping its rows by curve gives
    #   X^T W X = sum_n eta(t_n) eta(t_n)^T (x) V_n^T W_n V_n,
    #   X^T W y = vec( sum_n V_n^T W_n y_n eta(t_n)^T ),
    # which is accumulated curve by curve without forming X.
    A_data <- matrix(0, k1 * k2, k1 * k2)
    B <- matrix(0, k1, k2)
    for (n in seq_len(prep$N)) {
      cv <- prep$curves[[n]]
      Vn <- basis_eval(bs, cv$u)
      zn <- Z[n, ]
      A_data <- A_data + kronecker(tcrossprod(zn), crossprod(Vn, cv$w * Vn))
      B <- B + tcrossprod(crossprod(Vn, cv$w * cv$y), zn)
    }
    b <- as.vector(B)
  }

  # Marginal penalty, eq. (6). vec() stacks the columns of Theta (k1 x k2),
  # so the Kronecker order is (time) (x) (domain).
  S <- lambda[1L] * kronecker(J_eta, P1) + lambda[2L] * kronecker(P2, J_nu)
  A <- symmetrise(A_data + S)
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(R))
    stop("The penalised normal matrix is not positive definite. Increase ",
         "lambda or reduce the number of basis functions (k).", call. = FALSE)
  theta <- backsolve(R, forwardsolve(t(R), b))
  # Effective degrees of freedom: tr{(A_data + S)^{-1} A_data}
  edf <- sum(diag(backsolve(R, forwardsolve(t(R), A_data))))
  list(Theta = matrix(theta, k1, k2), edf = edf, Z = Z)
}

# Smoothing parameter selection -------------------------------------------------

# REML for a single penalised regression  sum_i w_i (y_i - X_i beta)^2 +
# lambda beta^T P beta, with the penalty used exactly as supplied (paraPen in
# mgcv does not rescale it). The returned lambda is therefore on the scale of
# that criterion.
reml_lambda <- function(y, X, P, w) {
  dat <- list(y = y, X = X, w = w)
  fit <- mgcv::gam(y ~ X - 1, data = dat, weights = w,
                   paraPen = list(X = list(P)), method = "REML")
  unname(fit$sp)
}

# Marginal data used for smoothing parameter selection (MHG, Section 3.3 and
# Web Appendix 4.2).
#   * lambda1: the cross-sectional mean curve  ybar(s) = N^{-1} sum_n Y_n(s).
#     For curve-specific grids each curve is linearly interpolated onto a
#     reference grid covering the range observed by all curves.
#   * lambda2: the scalar series  z_n = int_D Y_n(s) ds / int_D 1 ds,
#     computed with the quadrature weights of curve n.
marginal_data <- function(prep, weights, m_ref = NULL) {
  if (prep$common) {
    cv1 <- prep$curves[[1L]]
    Ymat <- matrix(vapply(prep$curves, `[[`, numeric(length(cv1$y)), "y"),
                   ncol = prep$N)
    u_ref <- cv1$u
    ybar <- rowMeans(Ymat)
    w_ref <- cv1$w
  } else {
    lo <- max(vapply(prep$curves, function(cv) min(cv$u), numeric(1L)))
    hi <- min(vapply(prep$curves, function(cv) max(cv$u), numeric(1L)))
    if (hi - lo <= 0)
      stop("The curves do not share a common range of evaluation points, so ",
           "the mean curve needed for REML selection of lambda1 is not ",
           "defined. Supply lambda1 directly.", call. = FALSE)
    if (hi - lo < 0.5)
      warning("The range observed by all curves covers only ",
              round(100 * (hi - lo)), "% of the domain; lambda1 is selected ",
              "from the mean curve on that range.", call. = FALSE)
    if (is.null(m_ref))
      m_ref <- stats::median(vapply(prep$curves, function(cv) length(cv$u), numeric(1L)))
    m_ref <- max(as.integer(round(m_ref)), 5L)
    u_ref <- seq(lo, hi, length.out = m_ref)
    ybar <- rowMeans(vapply(prep$curves, function(cv)
      stats::approx(cv$u, cv$y, xout = u_ref)$y, numeric(m_ref)))
    w_ref <- quad_weights(u_ref, weights)
  }
  W_n <- vapply(prep$curves, function(cv) sum(cv$w), numeric(1L))
  z <- vapply(prep$curves, function(cv) sum(cv$w * cv$y), numeric(1L)) / W_n
  list(u_ref = u_ref, ybar = ybar, w_ref = w_ref, z = z, W_n = W_n)
}

# Marginal REML on the estimator's own bases (select = "reml-bspline").
#
# Restricting the criterion L(Theta) + P(T) to trends constant in t, T = f(s),
# gives  sum_n sum_j w_j {Y_nj - f(s_j)}^2 = N sum_j w_j {ybar_j - f(s_j)}^2
# + const, with penalty lambda1 int f''^2. REML on the mean curve with prior
# weights w_j and the same basis nu and penalty P_1 estimates lambda_s for the
# criterion  sum_j w_j (ybar_j - f_j)^2 + lambda_s int f''^2, hence
#   lambda1 = N * lambda_s.
# Restricting to trends constant in s, T = g(t), gives
#   sum_n W_n {z_n - g(t_n)}^2 + const,  W_n = sum_j w_nj,
# with penalty lambda2 int g''^2, hence lambda2 = lambda_t from REML on z_n
# with prior weights W_n, basis eta and penalty P_2.
select_lambda_bspline <- function(md, prep, bs, bt) {
  lam_s <- reml_lambda(md$ybar, basis_eval(bs, md$u_ref), basis_penalty(bs), md$w_ref)
  lam_t <- reml_lambda(md$z, basis_eval(bt, prep$u_t), basis_penalty(bt), md$W_n)
  list(lambda = c(prep$N * lam_s, lam_t), marginal = c(lambda_s = lam_s, lambda_t = lam_t))
}

# Marginal REML in the manner of the supplementary code of MHG
# (select = "reml-supplement"): the smoothing parameters of
# gam(ybar ~ te(s, bs = "cr")) and gam(z ~ te(t, bs = "cr")) are used as
# lambda1 and lambda2 of the unweighted criterion of equation (8).
#   k = NULL       te() default basis dimension (5), as in the original code;
#   k = c(k1, k2)  cubic regression splines with the same numbers of basis
#                  functions as nu and eta ("matched").
#   ar = "time"    AR(1) errors for the series z_n (mgcv::gamm with corAR1);
#   ar = "both"    also AR(1) errors along the grid for the mean curve.
# The values are on mgcv's internal scale (te() rescales its penalty), so they
# are not REML estimates for the B-spline bases nu and eta. To use them with
# quadrature weights w_nj the unweighted criterion is matched on average:
# lambda_w = lambda * mean_n(W_n / m_n), W_n = sum_j w_nj (identity for
# weights = "none").
select_lambda_supplement <- function(md, prep, k = NULL, ar = c("none", "time", "both")) {
  ar <- match.arg(ar)
  sp_fit <- function(y, x, kk, ar1) {
    dat <- data.frame(y = y, x = x)
    if (is.null(kk)) {
      f <- y ~ te(x, bs = "cr")
    } else {
      kk <- min(as.integer(kk), length(unique(x)) - 1L)
      f <- y ~ te(x, bs = "cr", k = kk)
    }
    if (ar1) {
      g <- tryCatch(mgcv::gamm(f, data = dat, correlation = nlme::corAR1(),
                               method = "REML"), error = function(e) NULL)
      if (!is.null(g))
        return(c(sp = unname(g$gam$sp),
                 phi = unname(stats::coef(g$lme$modelStruct$corStruct,
                                          unconstrained = FALSE))))
      warning("The AR(1) marginal fit failed; independent errors were used.",
              call. = FALSE)
    }
    c(sp = unname(mgcv::gam(f, data = dat, method = "REML")$sp), phi = NA_real_)
  }
  fs <- sp_fit(md$ybar, md$u_ref, if (is.null(k)) NULL else k[1L], ar == "both")
  ft <- sp_fit(md$z, prep$u_t, if (is.null(k)) NULL else k[2L], ar != "none")
  scale <- mean(vapply(prep$curves, function(cv) sum(cv$w) / length(cv$w), numeric(1L)))
  list(lambda = scale * c(fs[["sp"]], ft[["sp"]]),
       marginal = c(lambda_s = fs[["sp"]], lambda_t = ft[["sp"]]),
       ar_phi = c(s = fs[["phi"]], t = ft[["phi"]]))
}

#' Control parameters for smoothing-parameter selection
#'
#' @param cv_eps Tolerance \eqn{\varepsilon} for the tail fraction
#'   \eqn{\tau(l)} of the long-run trace variance used to choose the gap
#'   \eqn{l} in cross-validation (\code{select = "cv"}).
#' @param cv_one_se Apply the one-standard-error rule along \eqn{\lambda_2}
#'   after cross-validation.
#' @param supp_k For \code{select = "reml-supplement"}: \code{"te-default"}
#'   (the \code{te()} default of 5 basis functions, as in the supplementary
#'   code of the paper) or \code{"matched"} (cubic regression splines with
#'   \eqn{k_1} and \eqn{k_2} functions, the numbers used by the estimator).
#'   Because mgcv rescales the penalty by a constant that depends on the
#'   basis, matching the dimension changes the scale of the selected values;
#'   in simulations \code{"matched"} over-smooths trends whose mean over
#'   \eqn{D} is nonlinear in time.
#' @param supp_ar For \code{select = "reml-supplement"}: error model of the
#'   marginal fits, \code{"none"} (independent), \code{"time"} (AR(1) for the
#'   series \eqn{\int_D Y_n}) or \code{"both"} (also AR(1) along the grid for
#'   the mean curve; this tends to absorb the mean curve into the error and
#'   to select a trend that is almost linear in \eqn{s}).
#' @return A list of class \code{ftrend_control}.
#' @export
ftrend_control <- function(cv_eps = 0.05, cv_one_se = TRUE,
                           supp_k = c("te-default", "matched"),
                           supp_ar = c("none", "time", "both")) {
  structure(list(cv_eps = cv_eps, cv_one_se = isTRUE(cv_one_se),
                 supp_k = match.arg(supp_k), supp_ar = match.arg(supp_ar)),
            class = "ftrend_control")
}

# Main function ---------------------------------------------------------------------

#' Nonparametric trend estimation for functional time series
#'
#' Estimates the deterministic trend \eqn{T(s, t)} in the model
#' \deqn{Y_n(s) = T(s, n/N) + X_n(s), \quad s \in D,\ n = 1, \dots, N,}
#' by a penalised tensor product of B-splines,
#' \eqn{T(s, t) = \nu(s)^\top \Theta \eta(t)}, with the marginal penalty
#' \deqn{\lambda_1 \int_0^1\!\int_D \{\partial^2 T/\partial s^2\}^2 \,ds\,dt +
#'       \lambda_2 \int_D\!\int_0^1 \{\partial^2 T/\partial t^2\}^2 \,dt\,ds}
#' (Martínez-Hernández and Genton 2021).
#'
#' @param Y The data: a numeric matrix (rows = evaluation points, columns =
#'   curves in time order), a long data frame, or an \code{fts_irreg} object.
#' @param s,t,domain,time,arg,value Passed to \code{\link{as_fts_irreg}}.
#' @param k Integer vector \code{c(k1, k2)}: number of B-spline functions on
#'   \eqn{D} and in time. Default \code{c(10, 15)}.
#' @param order Spline order(s), recycled to length 2. Default 4 (cubic).
#' @param penalty_order Derivative order(s) in the penalty, recycled to
#'   length 2. Default 2.
#' @param lambda Numeric vector \code{c(lambda1, lambda2)}. \code{NA} entries
#'   are selected by the method in \code{select}. Values refer to the rescaled problem
#'   (\eqn{D} and time mapped to \eqn{[0, 1]}) and to the quadrature weights in
#'   \code{weights}.
#' @param weights Quadrature weights approximating \eqn{\int_D} in the loss:
#'   \code{"trapezoid"} (default), \code{"equal"} (\eqn{1/m_n}) or
#'   \code{"none"} (unit weights; the unweighted criterion of equation (8) in
#'   the paper).
#' @param select How missing \code{lambda} values are selected:
#'   \code{"cv"} (default; leave-block-out cross-validation over the whole
#'   surface, see Details), \code{"reml-bspline"} (marginal REML on the same
#'   B-spline bases and penalties as the estimator, converted exactly to the
#'   scale of the full criterion) or \code{"reml-supplement"} (marginal REML
#'   with cubic regression splines in \code{mgcv}, as in the supplementary
#'   code of the paper; see \code{\link{ftrend_control}} for the basis
#'   dimension and AR(1) options).
#' @param control Options from \code{\link{ftrend_control}}.
#' @param reml Defunct; use \code{select}.
#' @param m_ref Number of points of the reference grid used to form the mean
#'   curve when curves are observed on different grids. Default: the median
#'   number of points per curve.
#' @param basis_s Optional basis object (class \code{ft_basis}) on \eqn{[0,1]}
#'   for the domain, replacing the default B-spline basis.
#' @param ... Further arguments passed to \code{ftrend_fit}.
#'
#' @details
#' \strong{Rescaling.} \eqn{D = [a, b]} is mapped linearly to \eqn{[0, 1]}.
#' Times are mapped by \eqn{u = (t - t_0)/(t_N - t_0)} with
#' \eqn{t_0 = t_1 - (t_N - t_1)/(N - 1)}, so that equispaced times give
#' \eqn{u_n = n/N} as in model (1).
#'
#' \strong{Estimator.} With quadrature weights \eqn{w_{nj}},
#' \eqn{\hat\Theta} minimises
#' \eqn{\sum_n \sum_j w_{nj}\{Y_n(s_{nj}) - \nu(s_{nj})^\top\Theta\eta(t_n)\}^2
#' + P(T)}. On a common grid this is solved through the identities
#' \eqn{(Z \otimes V)^\top W (Z \otimes V) = Z^\top Z \otimes V^\top W V} and
#' \eqn{(Z \otimes V)^\top W \mathrm{vec}(Y) = \mathrm{vec}(V^\top W Y Z)};
#' for curve-specific grids the normal equations are accumulated curve by
#' curve.
#'
#' \strong{Smoothing parameters: cross-validation (default).} With
#' \code{select = "cv"}, \eqn{(\lambda_1, \lambda_2)} minimise
#' \deqn{\mathrm{CV}_l(\lambda) = \sum_n \sum_j w_{nj}\{Y_n(s_{nj}) -
#' \hat T^{(-n,l)}(s_{nj}, t_n)\}^2,}
#' where \eqn{\hat T^{(-n,l)}} is the same estimator fitted without the curves
#' at times \eqn{|n' - n| \le l}. Leaving out a block removes the bias caused
#' by the serial dependence of \eqn{X_n} (Chu and Marron 1991). The gap is
#' \eqn{l = \max(2, \hat l)}, with \eqn{\hat l} the smallest lag at which the
#' estimated tail fraction of the long-run trace variance is at most
#' \code{cv_eps}, computed from the residuals of a pilot fit with
#' \eqn{l_0 = \max(5, \lceil N^{1/3}\rceil)} and bounded so that
#' \eqn{2l + 1 \le N/(2\,\mathrm{df}_t)}. Finally, the largest
#' \eqn{\lambda_2} whose criterion exceeds the minimum by at most one standard
#' error (allowing for serial dependence) is taken. Each \eqn{\lambda} is
#' capped where the fit stops changing, to avoid ill-conditioning.
#' Only the selection changes; the estimator is that of the paper.
#'
#' \strong{Smoothing parameters: marginal REML.} Following Section 3.3 of the
#' paper, \eqn{\lambda_1} is selected by REML from the mean curve
#' \eqn{\bar Y(s)} and \eqn{\lambda_2} from the series
#' \eqn{\int_D Y_n(s)\,ds}. With \code{select = "reml-bspline"} the
#' marginal fits use the bases \eqn{\nu}, \eqn{\eta} and their unscaled
#' penalties. Writing \eqn{T = \bar T + \tilde T} with
#' \eqn{\bar T(t) = \int_D T(s, t)\,ds}, the identity
#' \eqn{J_\nu = \bar\nu\bar\nu^\top + \int (\nu - \bar\nu)(\nu - \bar\nu)^\top ds}
#' (with \eqn{\bar\nu = \int_D \nu}) splits the criterion exactly into a
#' problem for \eqn{\bar T} and one for \eqn{\tilde T}; the first is the
#' penalised spline of \eqn{\int_D Y_n} with parameter \eqn{\lambda_2}, so
#' \eqn{\lambda_2 = \hat\lambda_t}. The analogous split in time holds up
#' to \eqn{O(N^{-1})} and gives \eqn{\lambda_1 = N\hat\lambda_s}. The
#' marginal series carry no information on \eqn{\tilde T}, which shares the
#' same smoothing parameters, and REML treats their correlated errors as
#' independent; the selected values tend to under-smooth.
#' With \code{select = "reml-supplement"} the selection follows the
#' supplementary code (with the basis dimension and the optional AR(1) error
#' model set in \code{control}); the resulting values refer to a different
#' basis and to a penalty rescaled by mgcv, so they have no direct
#' interpretation in the criterion above.
#'
#' @return An object of class \code{ftrend}; see \code{\link{ftrend-methods}}.
#' @references
#' Martínez-Hernández, Israel, and Marc G. Genton. 2021. "Nonparametric Trend
#' Estimation in Functional Time Series with Application to Annual Mortality
#' Rates." \emph{Biometrics} 77 (3): 866–78.
#'
#' Chu, C. K., and J. S. Marron. 1991. "Comparison of Two Bandwidth Selectors
#' with Dependent Errors." \emph{Annals of Statistics} 19 (4): 1906–18.
#'
#' Wood, Simon N. 2006. "Low-Rank Scale-Invariant Tensor Product Smooths for
#' Generalized Additive Mixed Models." \emph{Biometrics} 62 (4): 1025–36.
#' @examples
#' sim <- simulate_fts(N = 100, m = 50, trend = "T2", seed = 1)
#' fit <- ftrend(sim$Y, s = sim$s)
#' fit
#' @export
ftrend <- function(Y, ...) UseMethod("ftrend")

#' @rdname ftrend
#' @export
ftrend.matrix <- function(Y, s = NULL, t = NULL, domain = NULL, ...) {
  x <- as_fts_irreg(Y, s = s, t = t, domain = domain)
  fit <- ftrend_fit(x, ...)
  fit$input <- list(type = "matrix", dimnames = dimnames(Y))
  fit$fitted <- format_like_input(fit, fit$fitted_internal)
  fit$residuals <- format_like_input(fit, fit$residuals_internal)
  fit$call <- match.call(expand.dots = TRUE)
  fit$call[[1L]] <- as.name("ftrend")
  fit
}

#' @rdname ftrend
#' @export
ftrend.data.frame <- function(Y, time = "t", arg = "s", value = "y",
                              domain = NULL, ...) {
  x <- as_fts_irreg(Y, time = time, arg = arg, value = value, domain = domain)
  fit <- ftrend_fit(x, ...)
  fit$input <- list(type = "data.frame", data = Y, value = value,
                    row_order = attr(x, "row_order"))
  fit$fitted <- format_like_input(fit, fit$fitted_internal)
  fit$residuals <- format_like_input(fit, fit$residuals_internal)
  fit$call <- match.call(expand.dots = TRUE)
  fit$call[[1L]] <- as.name("ftrend")
  fit
}

#' @rdname ftrend
#' @export
ftrend.fts_irreg <- function(Y, ...) {
  x <- validate_fts_irreg(Y)
  fit <- ftrend_fit(x, ...)
  fit$input <- list(type = "fts_irreg")
  fit$fitted <- format_like_input(fit, fit$fitted_internal)
  fit$residuals <- format_like_input(fit, fit$residuals_internal)
  fit$call <- match.call(expand.dots = TRUE)
  fit$call[[1L]] <- as.name("ftrend")
  fit
}

#' @export
ftrend.default <- function(Y, ...) {
  stop("'Y' must be a numeric matrix, a data frame or an 'fts_irreg' object.",
       call. = FALSE)
}

#' @rdname ftrend
#' @param x An \code{fts_irreg} object.
#' @export
ftrend_fit <- function(x, k = c(10, 15), order = 4, penalty_order = 2,
                       lambda = c(NA, NA), weights = c("trapezoid", "equal", "none"),
                       select = c("cv", "reml-bspline", "reml-supplement"),
                       control = ftrend_control(), m_ref = NULL,
                       basis_s = NULL, reml = NULL) {
  if (!is.null(reml))
    stop("Argument 'reml' has been replaced by 'select' (\"cv\", ",
         "\"reml-bspline\" or \"reml-supplement\").", call. = FALSE)
  weights <- match.arg(weights)
  select <- match.arg(select)
  if (!inherits(control, "ftrend_control"))
    stop("'control' must be created by ftrend_control().", call. = FALSE)
  if (length(k) != 2L) stop("'k' must be a vector c(k1, k2).", call. = FALSE)
  order <- rep_len(order, 2L)
  penalty_order <- rep_len(penalty_order, 2L)
  lambda <- as.numeric(rep_len(lambda, 2L))
  if (any(!is.na(lambda) & (lambda < 0 | !is.finite(lambda))))
    stop("'lambda' values must be non-negative and finite (or NA).", call. = FALSE)

  bs <- if (is.null(basis_s)) bspline_basis(k[1L], order[1L], penalty_order[1L]) else basis_s
  if (!inherits(bs, "ft_basis")) stop("'basis_s' must be an 'ft_basis' object.", call. = FALSE)
  bt <- bspline_basis(k[2L], order[2L], penalty_order[2L])

  prep <- prepare_fts(x, weights)

  label <- switch(select,
    cv = "CV (leave-block-out)",
    `reml-bspline` = "REML (bspline)",
    `reml-supplement` = paste0("REML (supplement, k ", control$supp_k,
                               if (control$supp_ar != "none") paste0(", AR(1) ", control$supp_ar), ")"))
  lambda_method <- ifelse(is.na(lambda), label, "supplied")
  marginal <- NULL
  selection <- list(method = select)
  if (anyNA(lambda)) {
    if (select == "cv") {
      sel <- select_lambda_cv_rule(prep, bs, bt, eps = control$cv_eps,
                                   one_se = control$cv_one_se, fixed = lambda)
      selection <- c(selection, list(gap = sel$gap$l, gap_dependence = sel$gap$l_dep,
                                     gap_max = sel$gap$l_max, flat_top_lag = sel$gap$M,
                                     lambda_cv_min = sel$lambda_min,
                                     lambda_pilot = sel$lambda_pilot, l0 = sel$l0,
                                     caps = sel$caps))
    } else {
      md <- marginal_data(prep, weights, m_ref)
      sel <- if (select == "reml-bspline") select_lambda_bspline(md, prep, bs, bt) else
        select_lambda_supplement(md, prep,
          k = if (control$supp_k == "matched") c(ncol(basis_gram(bs)), ncol(basis_gram(bt))),
          ar = control$supp_ar)
      marginal <- sel$marginal
      if (!is.null(sel$ar_phi)) selection$ar_phi <- sel$ar_phi
    }
    lambda[is.na(lambda)] <- sel$lambda[is.na(lambda)]
  }
  names(lambda) <- c("lambda1", "lambda2")
  names(lambda_method) <- names(lambda)

  pf <- penalised_fit(prep, bs, bt, lambda)

  # Fitted trend at the observed points and detrended series
  # X_tilde_n(s) = Y_n(s) - T_hat(s, n/N)  (MHG, Section 2.3).
  fitted <- lapply(seq_len(prep$N), function(n) {
    cv <- prep$curves[[n]]
    as.vector(basis_eval(bs, cv$u) %*% (pf$Theta %*% pf$Z[n, ]))
  })
  resid <- lapply(seq_len(prep$N), function(n) x$curves[[n]]$y - fitted[[n]])

  structure(list(
    coefficients = pf$Theta,
    basis_s = bs, basis_t = bt,
    lambda = lambda, lambda_method = lambda_method, marginal_lambda = marginal,
    weights = weights, select = select, control = control, selection = selection,
    edf = pf$edf,
    domain = x$domain, time_map = prep$time_map,
    data = x,
    fitted_internal = fitted, residuals_internal = resid,
    common_grid = prep$common
  ), class = "ftrend")
}

# Return a list of per-curve vectors in the format of the original input.
format_like_input <- function(object, vals) {
  x <- object$data
  switch(object$input$type,
         matrix = {
           M <- matrix(unlist(vals), ncol = length(vals))
           dimnames(M) <- object$input$dimnames
           M
         },
         data.frame = {
           out <- object$input$data
           v <- numeric(nrow(out))
           v[object$input$row_order] <- unlist(vals)
           out[[object$input$value]] <- v
           out
         },
         fts_irreg = {
           x$curves <- lapply(seq_along(vals), function(n)
             list(s = x$curves[[n]]$s, y = vals[[n]]))
           x
         })
}
