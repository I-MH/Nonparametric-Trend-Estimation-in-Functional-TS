#' ftrend: Nonparametric Trend Estimation in Functional Time Series
#'
#' Implements the penalised tensor-product spline estimator of the functional
#' trend \eqn{T(s, t)} in \eqn{Y_n(s) = T(s, n/N) + X_n(s)} proposed by
#' Martínez-Hernández and Genton (2021). The main function is
#' \code{\link{ftrend}}.
#'
#' @references
#' Martínez-Hernández, Israel, and Marc G. Genton. 2021. "Nonparametric Trend
#' Estimation in Functional Time Series with Application to Annual Mortality
#' Rates." \emph{Biometrics} 77 (3): 866–78.
#' @importFrom stats coef fitted predict residuals
#' @keywords internal
"_PACKAGE"
