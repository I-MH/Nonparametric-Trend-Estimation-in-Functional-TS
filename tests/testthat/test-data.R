test_that("matrix input is converted with the documented defaults", {
  Y <- matrix(rnorm(20 * 6), 20, 6)
  x <- as_fts_irreg(Y)
  expect_s3_class(x, "fts_irreg")
  expect_equal(x$t, (1:6) / 6)
  expect_equal(x$curves[[1]]$s, seq(0, 1, length.out = 20))
  expect_equal(x$curves[[3]]$y, Y[, 3])
  expect_true(is_common_grid(x))
  expect_equal(unname(as.matrix(x)), Y)
})

test_that("long data frame input handles curve-specific grids and row order", {
  set.seed(2)
  df <- do.call(rbind, lapply(1:5, function(n) {
    s <- sort(runif(5 + n))
    data.frame(year = 2000 + n, age = s, rate = sin(s) + n)
  }))
  df <- df[sample(nrow(df)), ]
  x <- as_fts_irreg(df, time = "year", arg = "age", value = "rate")
  expect_equal(x$t, 2001:2005)
  expect_false(is_common_grid(x))
  expect_equal(vapply(x$curves, function(cv) length(cv$s), 1L), 6:10)
  expect_true(all(vapply(x$curves, function(cv) all(diff(cv$s) > 0), TRUE)))
  back <- as.data.frame(x)
  expect_equal(nrow(back), nrow(df))
})

test_that("fts_irreg constructor and coercions round-trip", {
  cv <- lapply(1:4, function(n) list(s = c(0, 0.2 * n, 1), y = c(1, n, 2)))
  x <- fts_irreg(cv, t = 1:4, domain = c(0, 1))
  expect_identical(as_fts_irreg(x), x)
  df <- as.data.frame(x)
  x2 <- as_fts_irreg(df)
  expect_equal(x2$curves, x$curves)
  expect_output(print(x), "curve-specific grids")
})

test_that("validation catches bad input with informative messages", {
  Y <- matrix(rnorm(30), 10, 3)
  Yna <- Y; Yna[4, 2] <- NA
  expect_error(as_fts_irreg(Yna), "missing values \\(first at row 4, column 2\\)")
  Yinf <- Y; Yinf[1, 1] <- Inf
  expect_error(as_fts_irreg(Yinf), "non-finite")
  expect_error(as_fts_irreg(Y, s = c(1:5, 7, 6, 8:10)), "not sorted")
  expect_error(as_fts_irreg(Y, s = c(1:5, 5, 7:10)), "duplicated")
  expect_error(as_fts_irreg(Y, s = 1:9), "length 9")
  expect_error(as_fts_irreg(Y, t = c(1, 3, 2)), "strictly increasing")
  expect_error(as_fts_irreg(Y, t = 1:4), "length 4")
  expect_error(as_fts_irreg(Y, domain = c(0, 0.5)), "outside the domain")
  expect_error(as_fts_irreg(Y[, 1:2]), "At least 3 curves")
  cv <- list(list(s = c(0, 1), y = 1:2), list(s = c(0, 1), y = 1:2),
             list(s = c(0, 1), y = c(1, NA)))
  expect_error(fts_irreg(cv), "curve 3 .* 'y' contains missing")
  cv[[3]] <- list(s = 0.5, y = 1)
  expect_error(fts_irreg(cv), "at least 2 evaluation points")
  df <- data.frame(t = c(1, 1, 2, 2, 3, 3), s = c(0, 0, 0, 1, 0, 1), y = 1)
  expect_error(as_fts_irreg(df), "duplicated evaluation points")
  expect_error(as_fts_irreg(data.frame(a = 1)), "Column 't' not found")
  expect_error(ftrend(list(1, 2)), "must be a numeric matrix")
})
