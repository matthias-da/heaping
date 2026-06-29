test_that(".draw_marginal stays within bounds for every method", {
  set.seed(2)
  trusted <- round(rlnorm(2000, 2.47, 1.0)); trusted <- trusted[trusted < 93]
  for (m in c("lnorm", "norm", "kernel", "unif")) {
    fit <- .fit_marginal(trusted, method = m, sd = NULL)
    d <- .draw_marginal(50, fit, llow = 28, lup = 32)
    expect_true(all(d >= 28 & d <= 32), info = m)
    expect_length(d, 50)
  }
})
