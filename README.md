# ftrend

R implementation of the nonparametric functional trend estimator of

> Martínez-Hernández, Israel, and Marc G. Genton. 2021. "Nonparametric Trend
> Estimation in Functional Time Series with Application to Annual Mortality
> Rates." *Biometrics* 77 (3): 866–78.

The model is `Y_n(s) = T(s, n/N) + X_n(s)`, with `T(s, t) = ν(s)ᵀ Θ η(t)` a
tensor product of cubic B-splines and marginal second-derivative penalties.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Martinez-Hernandez/Nonparametric-Trend-Estimation-in-Functional-TS",
                        build_vignettes = TRUE)
```

Building the vignette needs `knitr`, `rmarkdown` and pandoc (included with
RStudio). Omit `build_vignettes = TRUE` to install without it.

The vignette shows a simulated example and the analysis of annual mortality
rates in France:

```r
browseVignettes("ftrend")
vignette("ftrend", package = "ftrend")
```

Its mortality section runs only when the data are available locally: either
the `demography` package is installed (it contains `fr.mort`), or the
environment variable `FTREND_HMD_MX` gives the path to a copy of the Human
Mortality Database file `Mx_1x1.txt` for France. The vignette never downloads
data.

## Original supplementary code

The R scripts published with the paper (supplementary material of
Martínez-Hernández and Genton 2021) are preserved at the tag
[`biometrics-2021`](https://github.com/Martinez-Hernandez/Nonparametric-Trend-Estimation-in-Functional-TS/tree/biometrics-2021).
`ftrend` reproduces their estimate exactly with
`select = "reml-supplement"`, `weights = "none"` and the supplement's time grid
(see `inst/comparison/compare_supplement.R`).

## Usage

```r
library(ftrend)
sim <- simulate_fts(N = 200, m = 50, trend = "T2", seed = 1)

fit <- ftrend(sim$Y, s = sim$s)          # matrix: rows = points, columns = curves
summary(fit)
plot(fit, type = "surface")
Xtilde <- residuals(fit)                 # detrended series, same format as input
predict(fit, s = c(0.25, 0.75), t = c(0.5, 1))
forecast_trend(fit, h = 4)               # equation (5) of the paper
```

Curves observed on their own grids can be given as a long data frame
(`ftrend(df, time = "year", arg = "age", value = "rate")`) or as an
`fts_irreg` object (`fts_irreg(list(list(s = ..., y = ...), ...), t = ...)`).

## Main arguments of `ftrend()`

| argument | default | meaning |
|---|---|---|
| `k` | `c(10, 15)` | numbers of B-splines in `s` and `t` |
| `order`, `penalty_order` | `4`, `2` | spline order and penalised derivative |
| `lambda` | `c(NA, NA)` | smoothing parameters; `NA` = selected by `select` |
| `select` | `"cv"` | `"cv"`: leave-block-out cross-validation with a data-driven gap and the one-standard-error rule; `"reml-bspline"`: marginal REML on the estimator's own bases, converted exactly to the full criterion; `"reml-supplement"`: marginal REML as in the paper's supplementary code |
| `control` | `ftrend_control()` | `cv_eps`, `cv_one_se`; for `"reml-supplement"`: `supp_k` (`"te-default"` or `"matched"`) and `supp_ar` (`"none"`, `"time"`, `"both"`: AR(1) errors in the marginal fits) |
| `weights` | `"trapezoid"` | quadrature for `∫_D`; `"none"` gives the unweighted criterion of eq. (8) |

## Contents

* `R/basis.R`: replaceable basis component (`bspline_basis()`, `basis_eval()`,
  Gram matrix `J`, penalty matrix `P`), exact Gauss–Legendre quadrature.
* `R/data.R`: input classes and validation (`fts_irreg`, `as_fts_irreg()`).
* `R/fit.R`: penalised estimator (eq. 8 via Kronecker identities; curve-wise
  accumulation for curve-specific grids) and marginal REML.
* `R/select_cv.R`: leave-block-out cross-validation for the smoothing
  parameters (the default selection).
* `R/methods.R`: `print`, `summary`, `plot`, `fitted`, `residuals`, `coef`,
  `predict`, `forecast_trend`.
* `R/simulate.R`: simulation design of Section 4, linear trend of Kokoszka and
  Young (2017), integrated squared error.
* `inst/comparison/compare_supplement.R`: reproduces the supplementary code's
  example to rounding error.
* `inst/simulation/run_simulation.R`: simulation study (T1–T5,
  N ∈ {100, 300, 500, 1000}).
* `inst/prototype/`: the studies used to choose the default selection rule.
* `vignettes/ftrend.Rmd`: simulated example and French mortality analysis
  (uses `demography::fr.mort` or a local HMD file; never downloads data).

## Tests

```r
testthat::test_local()
```

## Licence

GPL (>= 2); see `LICENSE.md`.
