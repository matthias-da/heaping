test_that("equal-weight selection moves ceiling(excess) records (legacy count)", {
  set.seed(1)
  idx <- 1:20
  # old ssize = ceiling(n - n/ratio); ratio 4, n 20 -> ceiling(15)=15; excess=15
  sel <- .select_to_correct(idx, weight = NULL, excess = 15)
  expect_length(sel, 15)
  expect_true(all(sel %in% idx))
})

test_that("weighted selection removes at least the excess mass, minimally", {
  set.seed(1)
  idx <- 1:10
  w <- rep(2, 10)                 # each record carries mass 2
  sel <- .select_to_correct(idx, weight = w, excess = 7)
  expect_equal(length(sel), 4)    # cumsum 2,4,6,8 -> first >=7 is the 4th
})

test_that("selection caps at available mass", {
  sel <- .select_to_correct(1:3, weight = c(1, 1, 1), excess = 99)
  expect_length(sel, 3)
})

test_that("zero or NA excess selects nothing", {
  expect_length(.select_to_correct(1:5, NULL, 0), 0)
  expect_length(.select_to_correct(1:5, NULL, NA_real_), 0)
})
