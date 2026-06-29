test_that("correctHeapsMI returns m reproducible, distinct imputations", {
  set.seed(1); age <- round(rlnorm(3000, 2.47, 1.4)); age <- age[age < 93]
  y5 <- seq(0, max(age), 5); a <- sample(c(age, age[age %in% y5]))
  mi1 <- correctHeapsMI(a, m = 5, heaps = "5year", seed = 99)
  mi2 <- correctHeapsMI(a, m = 5, heaps = "5year", seed = 99)
  expect_s3_class(mi1, "heapingMI")
  expect_length(mi1$imputations, 5)
  expect_identical(mi1$imputations[[1]], mi2$imputations[[1]])        # reproducible
  expect_false(identical(mi1$imputations[[1]], mi1$imputations[[2]])) # distinct draws
})

test_that("correctHeapsMI validates m and prints a summary", {
  set.seed(2); age <- round(rlnorm(2000, 2.47, 1.4)); age <- age[age < 93]
  y5 <- seq(0, max(age), 5); a <- sample(c(age, age[age %in% y5]))
  expect_error(correctHeapsMI(a, m = 0), "positive")
  mi <- correctHeapsMI(a, m = 3, heaps = "5year", seed = 7)
  expect_length(mi$seeds, 3)
  expect_output(print(mi), "heapingMI")
})
