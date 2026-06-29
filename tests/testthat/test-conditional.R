make_df <- function(n = 1500) {
  set.seed(7); z <- rnorm(n)
  age <- round(30 + 8 * z + rnorm(n, 0, 3)); age[age < 0] <- 0
  data.frame(age = age, z = z)
}
test_that(".draw_cond (ranger) respects per-row bounds and falls back when sparse", {
  skip_if_not_installed("ranger")
  df <- make_df(); tr <- df[1:1400, ]; nw <- df[1401:1410, ]
  fit <- .fit_cond(tr, age ~ z, engine = "ranger", weight_col = NULL, case.weights = NULL)
  fm  <- .fit_marginal(tr$age, "lnorm", NULL)
  d <- .draw_cond(fit, nw, llow = rep(28, 10), lup = rep(32, 10), fit_marginal = fm)
  expect_length(d, 10)
  expect_true(all(d >= 28 & d <= 32))
})
test_that(".draw_cond (lm) draws truncated normal around the fitted value", {
  df <- make_df(); tr <- df[1:1400, ]; nw <- df[1401:1410, ]
  fit <- .fit_cond(tr, age ~ z, engine = "lm", weight_col = NULL, case.weights = NULL)
  fm  <- .fit_marginal(tr$age, "norm", NULL)
  d <- .draw_cond(fit, nw, llow = rep(20, 10), lup = rep(40, 10), fit_marginal = fm)
  expect_true(all(d >= 20 & d <= 40)); expect_length(d, 10)
})
