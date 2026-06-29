test_that(".heap_ratios reduces to count ratio when unweighted", {
  x <- c(rep(30, 20), rep(29, 5), rep(31, 5), 28, 32)
  r <- .heap_ratios(x, s = 30, weight = NULL)
  # count(30)=20, neighbours mean = (5+5)/2 = 5 -> ratio 4
  expect_equal(unname(r$ratio["30"]), 4)
  expect_equal(unname(r$expected["30"]), 5)
})

test_that(".heap_ratios uses weighted counts when weights given", {
  x <- c(rep(30, 10), rep(29, 10), rep(31, 10))
  w <- c(rep(2, 10), rep(1, 10), rep(1, 10))   # heap mass 20, neighbours 10 & 10
  r <- .heap_ratios(x, s = 30, weight = w)
  expect_equal(unname(r$ratio["30"]), 2)        # 20 / mean(10,10)
})
