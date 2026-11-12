# One-sided rule for a heap at the lower boundary (boundary = "lower").

make_heaped <- function(n = 6000, intensity = 0.3, interval = 5, seed = 1) {
  set.seed(seed)
  a <- round(rlnorm(n, 2.466869, 1.652772)); a <- a[a < 93]
  k <- sample(length(a), round(intensity * length(a)))
  a[k] <- round(a[k] / interval) * interval
  a
}

test_that(".heap_ratios skips the boundary heap by default and flags it with 'lower'", {
  x <- c(rep(0, 40), rep(1, 10), rep(2, 8), rep(5, 20), rep(4, 5), rep(6, 5))
  r0 <- .heap_ratios(x, s = c(0, 5))
  expect_true(is.na(r0$ratio["0"]))
  expect_false(any(r0$onesided))
  r1 <- .heap_ratios(x, s = c(0, 5), boundary = "lower")
  expect_equal(unname(r1$expected["0"]), 10)
  expect_equal(unname(r1$ratio["0"]), 4)
  expect_equal(unname(r1$onesided), c(TRUE, FALSE))
  expect_equal(r1$ratio["5"], r0$ratio["5"])          # interior heap unchanged
  r2 <- .heap_ratios(x, s = c(0, 5), boundary = "lower", bin = 0.5)
  expect_equal(unname(r2$expected["0"]), 5)
  expect_equal(unname(r2$ratio["0"]), 8)
})

test_that("default leaves the heap at age 0 untouched and equals boundary = 'none'", {
  h <- make_heaped()
  a <- correctHeaps(h, "5year", seed = 3)
  b <- correctHeaps(h, "5year", seed = 3, boundary = "none")
  expect_identical(a, b)
  expect_equal(sum(a == 0), sum(h == 0))
})

test_that("boundary = 'lower' corrects the heap at age 0 with a one-sided window", {
  h <- make_heaped()
  c0 <- sum(h == 0); c1 <- sum(h == 1)
  for (bin in c(1, 0.5)) {
    r <- correctHeaps(h, "5year", seed = 3, boundary = "lower", boundary.bin = bin,
                      verbose = TRUE)
    out <- r$corrected
    n_corr <- ceiling(c0 - bin * c1)
    expect_equal(unname(r$changes_by_heap["0"]), n_corr)
    expect_equal(sum(out == 0), c0 - n_corr)
    moved <- out[h == 0 & out != 0]
    expect_length(moved, n_corr)
    expect_true(all(moved %in% 1:2))                  # window 1..width, heap age excluded
    expect_true(all(out >= 0))
    expect_equal(length(out), length(h))
  }
})

test_that("one-sided window holds for every method and the model arm", {
  h <- make_heaped(n = 3000)
  for (mt in c("lnorm", "norm", "unif", "kernel")) {
    out <- correctHeaps(h, "10year", method = mt, seed = 5, boundary = "lower")
    moved <- out[h == 0 & out != 0]
    expect_gt(length(moved), 0)
    expect_true(all(moved %in% 1:5), info = mt)
    expect_true(all(out == round(out)), info = mt)
  }
  set.seed(9); z <- as.numeric(scale(h)) + rnorm(length(h), 0, 0.6)
  dm <- data.frame(age = h, z = z)
  engines <- c("lm", if (requireNamespace("ranger", quietly = TRUE)) "ranger")
  for (eng in engines) {
    out <- correctHeaps(h, "5year", model = age ~ z, dataModel = dm,
                        model.engine = eng, seed = 5, boundary = "lower")
    moved <- out[h == 0 & out != 0]
    expect_gt(length(moved), 0)
    expect_true(all(moved %in% 1:2), info = eng)
  }
})

test_that("a boundary heap with an empty neighbour stays uncorrected", {
  x <- c(rep(0, 50), rep(2, 10), rep(3, 10), rep(4, 5), rep(5, 30), rep(6, 5))
  out <- correctHeaps(x, "5year", seed = 1, boundary = "lower")
  expect_equal(sum(out == 0), 50)
})

test_that("the rule applies to a heap at the smallest observed age, not only age 0", {
  set.seed(2)
  x <- c(rep(15, 120), sample(16:60, 1500, replace = TRUE))
  c15 <- sum(x == 15); c16 <- sum(x == 16)
  r <- correctHeaps(x, heaps = 15, seed = 1, boundary = "lower", verbose = TRUE)
  expect_equal(unname(r$changes_by_heap["15"]), ceiling(c15 - c16))
  moved <- r$corrected[x == 15 & r$corrected != 15]
  expect_true(all(moved %in% 16:17))
  expect_true(all(r$corrected >= 15))
})

test_that("weights and the MI wrapper work with boundary = 'lower'", {
  h <- make_heaped(n = 3000)
  set.seed(4); w <- runif(length(h), 0.5, 2)
  out <- correctHeaps(h, "5year", weight = w, seed = 1, boundary = "lower",
                      boundary.bin = 0.5)
  expect_lt(sum(w[out == 0]), sum(w[h == 0]))
  expect_true(all(out[h == 0 & out != 0] %in% 1:2))
  mi <- correctHeapsMI(h, m = 3, heaps = "5year", seed = 1, boundary = "lower")
  expect_true(all(sapply(mi$imputations, function(i) sum(i == 0) < sum(h == 0))))
})

test_that("boundary arguments are validated", {
  h <- make_heaped(n = 1000)
  expect_error(correctHeaps(h, "5year", boundary = "upper"))
  expect_error(correctHeaps(h, "5year", boundary = "lower", boundary.bin = 0),
               "boundary.bin")
  expect_error(correctHeaps(h, "5year", boundary = "lower", boundary.bin = c(1, 2)),
               "boundary.bin")
})
