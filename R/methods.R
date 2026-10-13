# Methods for class "ftrend" -------------------------------------------------------

#' Methods for fitted functional trends
#'
#' @param object,x An object of class \code{ftrend}.
#' @param ... Unused (or passed to graphics functions in \code{plot}).
#' @name ftrend-methods
#' @return \code{coef}: the \eqn{k_1 \times k_2} coefficient matrix
#'   \eqn{\hat\Theta}. \code{fitted}: \eqn{\hat T(s, t_n)} at the observed
#'   points, and \code{residuals}: the detrended series
#'   \eqn{\tilde X_n(s) = Y_n(s) - \hat T(s, t_n)}, both in the format of the
#'   input data (matrix, data frame or \code{fts_irreg}).
NULL

#' @rdname ftrend-methods
#' @export
coef.ftrend <- function(object, ...) object$coefficients

#' @rdname ftrend-methods
#' @export
fitted.ftrend <- function(object, ...) object$fitted

#' @rdname ftrend-methods
#' @export
residuals.ftrend <- function(object, ...) object$residuals

#' @rdname ftrend-methods
#' @export
print.ftrend <- function(x, ...) {
  cat("Functional trend estimate (penalised tensor-product B-splines)\n")
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  N <- length(x$data$curves)
  cat("N = ", N, " curves; domain [", format(x$domain[1L]), ", ",
      format(x$domain[2L]), "]; ",
      if (x$common_grid) "common grid" else "curve-specific grids", "\n", sep = "")
  cat("Bases: k1 = ", x$basis_s$k, ", k2 = ", x$basis_t$k,
      " (order ", x$basis_s$order, ", ", x$basis_t$order, ")\n", sep = "")
  cat("lambda1 = ", format(signif(x$lambda[1L], 4)), " [", x$lambda_method[1L],
      "], lambda2 = ", format(signif(x$lambda[2L], 4)), " [",
      x$lambda_method[2L], "]\n", sep = "")
  cat("Effective degrees of freedom: ", format(round(x$edf, 2)), "\n", sep = "")
  invisible(x)
}

#' @rdname ftrend-methods
#' @export
summary.ftrend <- function(object, ...) {
  prep <- prepare_fts(object$data, object$weights)
  w <- unlist(lapply(prep$curves, `[[`, "w"))
  y <- unlist(lapply(prep$curves, `[[`, "y"))
  r <- unlist(object$residuals_internal)
  ybar <- sum(w * y) / sum(w)
  mn <- vapply(object$data$curves, function(cv) length(cv$s), integer(1L))
  out <- list(
    call = object$call, N = length(mn), m = range(mn),
    common_grid = object$common_grid, domain = object$domain,
    t_range = range(object$data$t),
    k = c(k1 = object$basis_s$k, k2 = object$basis_t$k),
    order = c(object$basis_s$order, object$basis_t$order),
    lambda = object$lambda, lambda_method = object$lambda_method,
    marginal_lambda = object$marginal_lambda, weights = object$weights,
    selection = object$selection,
    edf = object$edf,
    # Weighted residual mean square: approximates N^{-1} sum_n int_D X_n^2 / |D|
    residual_ms = sum(w * r^2) / sum(w),
    r_squared = 1 - sum(w * r^2) / sum(w * (y - ybar)^2)
  )
  class(out) <- "summary.ftrend"
  out
}

#' @export
print.summary.ftrend <- function(x, ...) {
  cat("Functional trend estimate\n")
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  cat("\nData: N = ", x$N, " curves, t in [", format(x$t_range[1L]), ", ",
      format(x$t_range[2L]), "], domain [", format(x$domain[1L]), ", ",
      format(x$domain[2L]), "]\n      ",
      if (x$common_grid) paste0("common grid, m = ", x$m[1L])
      else paste0("curve-specific grids, m_n from ", x$m[1L], " to ", x$m[2L]),
      "; quadrature weights: ", x$weights, "\n", sep = "")
  cat("Bases: k1 = ", x$k[1L], ", k2 = ", x$k[2L], ", orders ", x$order[1L],
      " and ", x$order[2L], "\n", sep = "")
  cat("\nSmoothing parameters (rescaled problem):\n")
  print(data.frame(value = signif(x$lambda, 5), method = x$lambda_method))
  sel <- x$selection
  if (!is.null(sel$gap))
    cat("Cross-validation: gap l = ", sel$gap, " (dependence ", sel$gap_dependence,
        ", upper bound ", sel$gap_max, "), pilot gap ", sel$l0,
        "; one-SE lambda2 = ", format(signif(x$lambda[2L], 4)), " (CV minimum ",
        format(signif(sel$lambda_cv_min[2L], 4)), ")\n", sep = "")
  if (!is.null(sel$ar_phi) && any(!is.na(sel$ar_phi)))
    cat("AR(1) coefficients of the marginal fits: ",
        paste(names(sel$ar_phi), "=", signif(sel$ar_phi, 3), collapse = ", "), "\n", sep = "")
  if (!is.null(x$marginal_lambda)) {
    cat("Marginal REML estimates: ")
    cat(paste(names(x$marginal_lambda), "=", signif(x$marginal_lambda, 5)),
        sep = ", ")
    cat("\n")
  }
  cat("\nEffective degrees of freedom: ", format(round(x$edf, 2)), "\n", sep = "")
  cat("Residual mean square (weighted): ", format(signif(x$residual_ms, 5)), "\n",
      sep = "")
  cat("Proportion of variation explained by the trend: ",
      format(round(x$r_squared, 4)), "\n", sep = "")
  invisible(x)
}

# Time-basis rows, with linear extrapolation beyond [0, 1].
# For u > 1 the row is eta(1) + (u - 1) eta'(1), so that
#   T_hat(s, u) = T_hat(s, 1) + (u - 1) dT_hat/dt (s, 1),
# which is the first-order Taylor forecast of equation (5) in MHG (with the
# time step 1/N between consecutive curves). Similarly for u < 0.
time_rows <- function(bt, u) {
  Z <- matrix(0, length(u), ncol(basis_gram(bt)))
  inside <- u >= 0 & u <= 1
  if (any(inside)) Z[inside, ] <- basis_eval(bt, u[inside])
  hi <- u > 1
  if (any(hi)) {
    e1 <- basis_eval(bt, 1)
    d1 <- basis_eval(bt, 1, deriv = 1)
    Z[hi, ] <- matrix(e1, sum(hi), ncol(Z), byrow = TRUE) +
      outer(u[hi] - 1, as.vector(d1))
  }
  lo <- u < 0
  if (any(lo)) {
    e0 <- basis_eval(bt, 0)
    d0 <- basis_eval(bt, 0, deriv = 1)
    Z[lo, ] <- matrix(e0, sum(lo), ncol(Z), byrow = TRUE) +
      outer(u[lo], as.vector(d0))
  }
  Z
}

#' Evaluate the estimated trend at new points
#'
#' @param object An object of class \code{ftrend}.
#' @param s Points of the domain (original scale). Default: 101 equispaced
#'   points on the domain.
#' @param t Times (original scale). Default: the observed times.
#' @param grid If \code{TRUE} (default), evaluate on the grid
#'   \code{s} \eqn{\times} \code{t} and return a \code{length(s)} by
#'   \code{length(t)} matrix. If \code{FALSE}, \code{s} and \code{t} are
#'   paired (recycled to a common length) and a vector is returned.
#' @param ... Unused.
#' @details A warning is given when \code{t} lies outside the observed time
#'   range \eqn{[t_1, t_N]}: a B-spline trend is not reliable there. Beyond the
#'   range of the time basis the trend is extrapolated linearly from the
#'   boundary, which coincides with the first-order Taylor forecast in
#'   equation (5) of the paper.
#' @return A matrix or a vector of trend values.
#' @export
predict.ftrend <- function(object, s = NULL, t = NULL, grid = TRUE, ...) {
  predict_trend(object, s, t, grid, warn = TRUE)
}

predict_trend <- function(object, s = NULL, t = NULL, grid = TRUE, warn = TRUE) {
  d <- object$domain
  if (is.null(s)) s <- seq(d[1L], d[2L], length.out = 101L)
  if (is.null(t)) t <- object$data$t
  s <- as.numeric(s)
  t <- as.numeric(t)
  if (anyNA(s) || anyNA(t)) stop("'s' and 't' must not contain missing values.",
                                 call. = FALSE)
  tol <- 1e-10 * (d[2L] - d[1L])
  if (any(s < d[1L] - tol | s > d[2L] + tol))
    stop("'s' must lie in the domain [", format(d[1L]), ", ", format(d[2L]), "].",
         call. = FALSE)
  tr <- range(object$data$t)
  if (warn && any(t < tr[1L] | t > tr[2L]))
    warning("Some values of 't' lie outside the observed time range [",
            format(tr[1L]), ", ", format(tr[2L]), "]; extrapolation of a ",
            "B-spline trend is unreliable.", call. = FALSE)
  u_s <- pmin(pmax(rescale_s(s, d), 0), 1)
  u_t <- rescale_t(t, object$time_map)
  V <- basis_eval(object$basis_s, u_s)
  Z <- time_rows(object$basis_t, u_t)
  Theta <- object$coefficients
  if (grid) {
    out <- V %*% Theta %*% t(Z)
    dimnames(out) <- list(s = format(s), t = format(t))
    out
  } else {
    n <- max(length(s), length(t))
    V <- V[rep_len(seq_len(nrow(V)), n), , drop = FALSE]
    Z <- Z[rep_len(seq_len(nrow(Z)), n), , drop = FALSE]
    rowSums((V %*% Theta) * Z)
  }
}

#' Forecast the functional trend
#'
#' h-step-ahead forecast of the trend by the first-order Taylor expansion of
#' equation (5) in Martínez-Hernández and Genton (2021), iterated \eqn{h}
#' times:
#' \eqn{\hat T_{N+h}(s) = \hat T(s, 1) + (h/N)\,\partial_t \hat T(s, t)|_{t=1}}.
#'
#' @param object An object of class \code{ftrend}.
#' @param h Forecast horizon (positive integer).
#' @param s Points of the domain. Default: 101 equispaced points.
#' @return A \code{length(s)} by \code{h} matrix; columns correspond to times
#'   \eqn{t_N + j\Delta}, \eqn{j = 1, \dots, h}, where \eqn{\Delta} is the
#'   average time step.
#' @details Equation (5) in the paper uses the step \eqn{1/(N+1)}; since
#'   \eqn{t_{N+1} - t_N = 1/N} on the rescaled time axis, the step used here
#'   is \eqn{1/N}. Iterating the one-step forecast with the derivative
#'   re-estimated from the last two values (as in the supplementary code)
#'   gives the same straight line.
#' @export
forecast_trend <- function(object, h = 1, s = NULL) {
  h <- as.integer(h)
  if (length(h) != 1L || is.na(h) || h < 1L)
    stop("'h' must be a positive integer.", call. = FALSE)
  tm <- object$time_map
  t_new <- tm$tN + tm$step * seq_len(h)
  predict_trend(object, s = s, t = t_new, grid = TRUE, warn = FALSE)
}

#' Plot a fitted functional trend
#'
#' @param x An object of class \code{ftrend}.
#' @param type \code{"surface"} (perspective plot of \eqn{\hat T(s, t)}),
#'   \code{"curves"} (the trend curves \eqn{\hat T(\cdot, t_n)} coloured by
#'   time), \code{"data"} (the observed curves) or \code{"residuals"} (the
#'   detrended curves \eqn{\tilde X_n}).
#' @param n_s Number of points in \eqn{D} at which the trend is evaluated.
#' @param col Colours for the curves (one per time point). Default: a
#'   sequential palette from early (dark) to late (light) times.
#' @param theta,phi Viewing angles for \code{type = "surface"}.
#' @param ... Further graphical parameters.
#' @return \code{x}, invisibly.
#' @export
plot.ftrend <- function(x, type = c("surface", "curves", "data", "residuals"),
                        n_s = 60, col = NULL, theta = -40, phi = 25, ...) {
  type <- match.arg(type)
  N <- length(x$data$curves)
  if (is.null(col)) col <- grDevices::hcl.colors(N, "viridis")
  s <- seq(x$domain[1L], x$domain[2L], length.out = n_s)
  tt <- x$data$t
  if (type == "surface") {
    M <- predict_trend(x, s, tt, warn = FALSE)
    nr <- nrow(M); nc <- ncol(M)
    zf <- (M[-1, -1] + M[-1, -nc] + M[-nr, -1] + M[-nr, -nc]) / 4
    pal <- grDevices::hcl.colors(100, "viridis")
    fcol <- pal[cut(zf, 100, include.lowest = TRUE)]
    graphics::persp(s, tt, M, theta = theta, phi = phi, col = fcol,
                    border = NA, shade = 0.3, ticktype = "detailed",
                    xlab = "s", ylab = "t", zlab = "T(s, t)", ...)
  } else if (type == "curves") {
    M <- predict_trend(x, s, tt, warn = FALSE)
    graphics::matplot(s, M, type = "l", lty = 1, col = col, xlab = "s",
                      ylab = expression(hat(T)(s, t[n])), ...)
  } else {
    obj <- if (type == "data") x$data$curves else
      lapply(seq_len(N), function(n) list(s = x$data$curves[[n]]$s,
                                          y = x$residuals_internal[[n]]))
    yl <- range(unlist(lapply(obj, `[[`, "y")))
    graphics::plot(NA, xlim = x$domain, ylim = yl, xlab = "s",
                   ylab = if (type == "data") expression(Y[n](s))
                          else expression(tilde(X)[n](s)), ...)
    for (n in seq_len(N)) graphics::lines(obj[[n]]$s, obj[[n]]$y, col = col[n])
  }
  invisible(x)
}
