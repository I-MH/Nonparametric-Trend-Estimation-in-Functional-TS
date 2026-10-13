# Functional time series container ---------------------------------------------
#
# Every accepted input form (matrix, long data frame, fts_irreg) is converted
# to one internal representation: an object of class "fts_irreg" holding a list
# of curves, each with its own evaluation points s and values y, together with
# the time vector t and the domain D = [a, b].

#' Functional time series observed on curve-specific grids
#'
#' Constructor for the S3 class \code{fts_irreg}: a functional time series
#' \eqn{Y_1, \dots, Y_N} in which curve \eqn{n} is observed at its own points
#' \eqn{s_{n1} < \dots < s_{n m_n}} of the domain \eqn{D = [a, b]}.
#'
#' @param curves A list of length \eqn{N}; element \eqn{n} is a list with
#'   numeric components \code{s} (strictly increasing evaluation points) and
#'   \code{y} (observed values, same length as \code{s}).
#' @param t Numeric vector of length \eqn{N} with the observation times,
#'   strictly increasing. Defaults to \eqn{t_n = n/N}.
#' @param domain Numeric vector \code{c(a, b)} with the domain \eqn{D}.
#'   Defaults to the range of all evaluation points.
#' @return An object of class \code{fts_irreg} (validated).
#' @seealso \code{\link{as_fts_irreg}}, \code{\link{validate_fts_irreg}}
#' @examples
#' x <- fts_irreg(list(list(s = c(0, 0.5, 1), y = c(1, 2, 1)),
#'                     list(s = c(0, 0.3, 0.8, 1), y = c(2, 3, 3, 2)),
#'                     list(s = c(0, 0.6, 1), y = c(3, 4, 3))),
#'                t = c(2001, 2002, 2003))
#' x
#' @export
fts_irreg <- function(curves, t = NULL, domain = NULL) {
  if (!is.list(curves) || length(curves) == 0L)
    stop("'curves' must be a non-empty list of curves.", call. = FALSE)
  curves <- lapply(curves, function(cv) {
    if (!is.list(cv) || is.null(cv$s) || is.null(cv$y))
      stop("Each curve must be a list with components 's' and 'y'.",
           call. = FALSE)
    list(s = as.numeric(cv$s), y = as.numeric(cv$y))
  })
  N <- length(curves)
  if (is.null(t)) t <- seq_len(N) / N
  if (is.null(domain)) {
    all_s <- unlist(lapply(curves, `[[`, "s"))
    domain <- if (all(is.finite(all_s))) range(all_s) else c(NA_real_, NA_real_)
  }
  x <- structure(list(curves = curves, t = as.numeric(t),
                      domain = as.numeric(domain)),
                 class = "fts_irreg")
  validate_fts_irreg(x)
}

#' Validate an \code{fts_irreg} object
#'
#' Checks for missing or non-finite values, unsorted or duplicated evaluation
#' points, time points that are not strictly increasing, an invalid domain,
#' and evaluation points outside the domain. Fails with an informative message
#' at the first problem found; data are never dropped silently.
#'
#' @param x An object of class \code{fts_irreg}.
#' @return \code{x}, invisibly, if all checks pass.
#' @export
validate_fts_irreg <- function(x) {
  if (!inherits(x, "fts_irreg"))
    stop("Object is not of class 'fts_irreg'.", call. = FALSE)
  N <- length(x$curves)
  t <- x$t
  if (N < 3L)
    stop("At least 3 curves are needed to estimate a trend (found ", N, ").",
         call. = FALSE)
  if (length(t) != N)
    stop("'t' has length ", length(t), " but there are ", N, " curves.",
         call. = FALSE)
  if (anyNA(t)) stop("'t' contains missing values.", call. = FALSE)
  if (any(!is.finite(t))) stop("'t' contains non-finite values.", call. = FALSE)
  if (any(diff(t) <= 0))
    stop("'t' must be strictly increasing (problem at position ",
         which(diff(t) <= 0)[1L] + 1L, ").", call. = FALSE)
  d <- x$domain
  if (length(d) != 2L || anyNA(d) || any(!is.finite(d)) || d[2L] <= d[1L])
    stop("'domain' must be a finite interval c(a, b) with a < b.", call. = FALSE)
  tol <- 1e-10 * (d[2L] - d[1L])
  for (n in seq_len(N)) {
    s <- x$curves[[n]]$s
    y <- x$curves[[n]]$y
    where <- paste0("curve ", n, " (t = ", format(t[n]), ")")
    if (length(s) != length(y))
      stop(where, ": 's' and 'y' have different lengths.", call. = FALSE)
    if (length(s) < 2L)
      stop(where, ": at least 2 evaluation points are needed.", call. = FALSE)
    if (anyNA(s)) stop(where, ": 's' contains missing values.", call. = FALSE)
    if (anyNA(y)) stop(where, ": 'y' contains missing values.", call. = FALSE)
    if (any(!is.finite(s))) stop(where, ": 's' contains non-finite values.", call. = FALSE)
    if (any(!is.finite(y))) stop(where, ": 'y' contains non-finite values.", call. = FALSE)
    ds <- diff(s)
    if (any(ds == 0)) stop(where, ": duplicated evaluation points.", call. = FALSE)
    if (any(ds < 0)) stop(where, ": evaluation points are not sorted increasingly.",
                          call. = FALSE)
    if (s[1L] < d[1L] - tol || s[length(s)] > d[2L] + tol)
      stop(where, ": evaluation points lie outside the domain [",
           format(d[1L]), ", ", format(d[2L]), "].", call. = FALSE)
  }
  invisible(x)
}

#' Coerce data to an \code{fts_irreg} object
#'
#' @param x A numeric matrix (rows = evaluation points, columns = curves in time
#'   order), a long data frame, or an \code{fts_irreg} object.
#' @param s For the matrix method: evaluation points (length = \code{nrow(x)}).
#'   Defaults to an equispaced grid on \eqn{[0, 1]}.
#' @param t Observation times (length = number of curves). Defaults to
#'   \eqn{t_n = n/N}. For the data-frame method the times are taken from the
#'   \code{time} column.
#' @param domain Domain \code{c(a, b)}; defaults to the range of the
#'   evaluation points.
#' @param time,arg,value For the data-frame method: names of the columns with
#'   the curve/time identifier (numeric), the evaluation point and the
#'   observed value.
#' @param ... Unused.
#' @details In the data-frame method, rows are sorted by time and, within each
#'   curve, by evaluation point. Row order carries no information in the long
#'   format, so sorting does not change the data; duplicated evaluation points
#'   within a curve are an error.
#' @return An object of class \code{fts_irreg}.
#' @export
as_fts_irreg <- function(x, ...) UseMethod("as_fts_irreg")

#' @rdname as_fts_irreg
#' @export
as_fts_irreg.fts_irreg <- function(x, ...) validate_fts_irreg(x)

#' @rdname as_fts_irreg
#' @export
as_fts_irreg.matrix <- function(x, s = NULL, t = NULL, domain = NULL, ...) {
  if (!is.numeric(x)) stop("The data matrix must be numeric.", call. = FALSE)
  m <- nrow(x)
  N <- ncol(x)
  if (is.null(s)) s <- seq(0, 1, length.out = m)
  if (length(s) != m)
    stop("'s' has length ", length(s), " but the matrix has ", m,
         " rows (evaluation points).", call. = FALSE)
  if (!is.null(t) && length(t) != N)
    stop("'t' has length ", length(t), " but the matrix has ", N,
         " columns (curves).", call. = FALSE)
  if (anyNA(s) || any(!is.finite(s)))
    stop("'s' contains missing or non-finite values.", call. = FALSE)
  if (any(diff(s) == 0)) stop("'s' contains duplicated points.", call. = FALSE)
  if (any(diff(s) < 0)) stop("'s' is not sorted increasingly.", call. = FALSE)
  if (anyNA(x)) {
    bad <- which(is.na(x), arr.ind = TRUE)[1L, ]
    stop("The data matrix contains missing values (first at row ", bad[1L],
         ", column ", bad[2L], "). Use the data-frame or 'fts_irreg' input ",
         "for curves observed on different grids.", call. = FALSE)
  }
  if (any(!is.finite(x)))
    stop("The data matrix contains non-finite values.", call. = FALSE)
  curves <- lapply(seq_len(N), function(n) list(s = as.numeric(s), y = as.numeric(x[, n])))
  fts_irreg(curves, t = t, domain = domain)
}

#' @rdname as_fts_irreg
#' @export
as_fts_irreg.data.frame <- function(x, time = "t", arg = "s", value = "y",
                                    domain = NULL, ...) {
  for (nm in c(time, arg, value))
    if (!nm %in% names(x))
      stop("Column '", nm, "' not found in the data frame.", call. = FALSE)
  tt <- x[[time]]
  if (!is.numeric(tt))
    stop("The time column '", time, "' must be numeric (e.g. the year).",
         call. = FALSE)
  if (anyNA(tt)) stop("The time column contains missing values.", call. = FALSE)
  ss <- as.numeric(x[[arg]])
  yy <- as.numeric(x[[value]])
  ord <- order(tt, ss)
  times <- sort(unique(tt))
  idx <- split(ord, factor(tt[ord], levels = times))
  curves <- lapply(idx, function(i) list(s = ss[i], y = yy[i]))
  names(curves) <- NULL
  out <- fts_irreg(curves, t = times, domain = domain)
  attr(out, "row_order") <- ord
  out
}

#' @export
as.matrix.fts_irreg <- function(x, ...) {
  if (!is_common_grid(x))
    stop("Curves are observed on different grids; no matrix form.", call. = FALSE)
  M <- vapply(x$curves, `[[`, numeric(length(x$curves[[1L]]$y)), "y")
  M <- matrix(M, ncol = length(x$curves))
  rownames(M) <- format(x$curves[[1L]]$s)
  colnames(M) <- format(x$t)
  M
}

#' @export
as.data.frame.fts_irreg <- function(x, row.names = NULL, optional = FALSE, ...) {
  mn <- vapply(x$curves, function(cv) length(cv$s), integer(1L))
  data.frame(t = rep(x$t, mn),
             s = unlist(lapply(x$curves, `[[`, "s")),
             y = unlist(lapply(x$curves, `[[`, "y")))
}

#' @export
print.fts_irreg <- function(x, ...) {
  mn <- vapply(x$curves, function(cv) length(cv$s), integer(1L))
  cat("Functional time series: N = ", length(x$curves), " curves, t from ",
      format(min(x$t)), " to ", format(max(x$t)), "\n", sep = "")
  cat("Domain: [", format(x$domain[1L]), ", ", format(x$domain[2L]), "]; ",
      if (is_common_grid(x)) paste0("common grid of m = ", mn[1L], " points")
      else paste0("curve-specific grids, ", min(mn), " to ", max(mn), " points"),
      "\n", sep = "")
  invisible(x)
}

#' Do all curves share the same evaluation grid?
#' @param x An \code{fts_irreg} object.
#' @return Logical.
#' @export
is_common_grid <- function(x) {
  s1 <- x$curves[[1L]]$s
  all(vapply(x$curves, function(cv) identical(cv$s, s1), logical(1L)))
}

# Maps between original and rescaled coordinates ---------------------------------

# Domain D = [a, b] is mapped linearly to [0, 1].
rescale_s <- function(s, domain) (s - domain[1L]) / (domain[2L] - domain[1L])

# Time map. With t_0 = t_1 - (t_N - t_1)/(N - 1) (one average time step before
# the first observation), t is mapped to u = (t - t_0)/(t_N - t_0). For
# equispaced times this gives u_n = n/N exactly, as in model (1) of the paper.
time_map <- function(t) {
  N <- length(t)
  step <- (t[N] - t[1L]) / (N - 1)
  list(t0 = t[1L] - step, tN = t[N], step = step)
}
rescale_t <- function(t, tm) (t - tm$t0) / (tm$tN - tm$t0)

# Quadrature weights approximating int_D g(s) ds on the rescaled domain [0, 1].
#  "trapezoid": trapezoidal rule over the observed points of each curve;
#  "equal":     1/m_n for each point (sums to one);
#  "none":      1 for each point, i.e. the unweighted sum of squares of
#               equation (8) in the paper (lambda then on the paper's scale).
quad_weights <- function(u, type = c("trapezoid", "equal", "none")) {
  type <- match.arg(type)
  m <- length(u)
  switch(type,
         trapezoid = {
           du <- diff(u)
           (c(du, 0) + c(0, du)) / 2
         },
         equal = rep(1 / m, m),
         none = rep(1, m))
}
