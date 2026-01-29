# Tests for whipple function

library(testthat)

# =============================================================================
# Standard Whipple Index Tests
# =============================================================================

test_that("whipple returns 100 for uniform data", {
  set.seed(42)
  # Create uniform age distribution (no heaping)
  age <- sample(23:62, 10000, replace = TRUE)

  result <- whipple(age, method = "standard")

  # Should be close to 100 for uniform distribution
  expect_true(result > 90 && result < 110)
})

test_that("whipple returns 500 for complete heaping", {
  # Ages ending in 0 or 5 only
  age <- sample(seq(25, 60, by = 5), 1000, replace = TRUE)

  result <- whipple(age, method = "standard")

  expect_equal(result, 500)
})

test_that("whipple handles weighted data correctly", {
  set.seed(42)
  age <- sample(23:62, 1000, replace = TRUE)
  weights <- runif(1000, 0.5, 1.5)

  result_weighted <- whipple(age, method = "standard", weight = weights)
  result_unweighted <- whipple(age, method = "standard")

  # Both should be close to 100 for uniform data
  expect_true(result_weighted > 80 && result_weighted < 120)
  expect_true(result_unweighted > 80 && result_unweighted < 120)
})

# =============================================================================
# Modified Whipple Index Tests
# =============================================================================

test_that("modified whipple returns close to 0 for uniform data", {
  set.seed(42)
  age <- sample(0:100, 10000, replace = TRUE)

  result <- whipple(age, method = "modified")

  # Should be close to 0 for uniform distribution
  expect_true(result < 0.1)
})

test_that("modified whipple returns close to 1 for extreme heaping", {
  # All ages ending in 0 only
  age <- sample(seq(0, 100, by = 10), 1000, replace = TRUE)

  result <- whipple(age, method = "modified")

  # Should be close to 1 for complete single-digit heaping
  expect_true(result > 0.9)
})

test_that("modified whipple handles 5-year heaping", {
  # Ages ending in 0 or 5 only
  age <- sample(seq(0, 100, by = 5), 1000, replace = TRUE)

  result <- whipple(age, method = "modified")

  # Should be between 0 and 1, indicating heaping
  expect_true(result > 0.5 && result < 1)
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("whipple validates method parameter", {
  age <- 20:50

  expect_error(
    whipple(age, method = "invalid"),
    "Unsupported value"
  )
})

test_that("whipple validates numeric input", {
  expect_error(
    whipple(c("a", "b", "c")),
    "must be a numeric"
  )
})

test_that("whipple validates weight length", {
  age <- 20:50
  weights <- 1:10

  expect_error(
    whipple(age, weight = weights),
    "same length"
  )
})

test_that("whipple validates weight type", {
  age <- 20:50
  weights <- as.character(1:31)

  expect_error(
    whipple(age, weight = weights),
    "must be a numeric"
  )
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("whipple handles small datasets", {
  age <- c(23, 25, 30, 35, 40)

  expect_no_error(whipple(age, method = "standard"))
  expect_no_error(whipple(age, method = "modified"))
})

test_that("whipple handles data outside 23-62 range", {
  # Standard Whipple only considers ages 23-62
  age <- c(10, 15, 20, 70, 80, 90)

  # Should handle gracefully (even if all data is filtered out)
  expect_true(is.numeric(whipple(age, method = "standard")) ||
                is.na(whipple(age, method = "standard")))
})
