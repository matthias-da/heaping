# Tests for sprague function

library(testthat)

# =============================================================================
# Basic Functionality Tests
# =============================================================================

test_that("sprague returns correct length", {
  # Example population data for 17 age groups
  x <- c(
    1971990, 2095820, 2157190, 2094110, 2116580,
    2003840, 1785690, 1502990, 1214170, 796934,
    627551, 530305, 488014, 364498, 259029, 158047, 125941
  )

  result <- sprague(x)

  # Should return 81 values (ages 0-79 plus 80+)
  expect_equal(length(result), 81)
})

test_that("sprague preserves total population", {
  x <- c(
    1971990, 2095820, 2157190, 2094110, 2116580,
    2003840, 1785690, 1502990, 1214170, 796934,
    627551, 530305, 488014, 364498, 259029, 158047, 125941
  )

  result <- sprague(x)

  # Sum should be preserved
  expect_equal(sum(result), sum(x), tolerance = 1e-6)
})

test_that("sprague returns named vector", {
  x <- c(
    1971990, 2095820, 2157190, 2094110, 2116580,
    2003840, 1785690, 1502990, 1214170, 796934,
    627551, 530305, 488014, 364498, 259029, 158047, 125941
  )

  result <- sprague(x)

  # Should have names 0, 1, 2, ..., 79, 80+
  expect_equal(names(result)[1], "0")
  expect_equal(names(result)[80], "79")
  expect_equal(names(result)[81], "80+")
})

test_that("sprague 80+ group is unchanged", {
  x <- c(
    1971990, 2095820, 2157190, 2094110, 2116580,
    2003840, 1785690, 1502990, 1214170, 796934,
    627551, 530305, 488014, 364498, 259029, 158047, 125941
  )

  result <- sprague(x)

  # 80+ should be unchanged
  expect_equal(result["80+"], x[17], ignore_attr = TRUE)
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("sprague validates input length", {
  x <- c(100, 200, 300)

  expect_error(
    sprague(x),
    "exactly 17"
  )
})

test_that("sprague validates numeric input", {
  x <- letters[1:17]

  expect_error(
    sprague(x),
    "must be a numeric"
  )
})

# =============================================================================
# Distribution Properties Tests
# =============================================================================

test_that("sprague produces reasonable distribution", {
  # Flat population distribution
  x <- rep(1000, 17)

  result <- sprague(x)

  # Within each 5-year group, values should sum to approximately 1000
  # (except for boundary effects)
  group_10_14 <- sum(result[11:15])  # ages 10-14
  expect_true(abs(group_10_14 - 1000) < 100)
})

test_that("sprague handles zero population groups", {
  x <- c(
    1000, 1000, 1000, 1000, 0,
    0, 1000, 1000, 1000, 1000,
    1000, 1000, 1000, 1000, 1000, 1000, 1000
  )

  # Should not error with zero values
  expect_no_error(sprague(x))

  result <- sprague(x)
  expect_equal(sum(result), sum(x), tolerance = 1e-6)
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("sprague handles very small populations", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)

  result <- sprague(x)

  expect_equal(length(result), 81)
  expect_equal(sum(result), sum(x), tolerance = 1e-6)
})

test_that("sprague handles very large populations", {
  x <- rep(1e9, 17)

  result <- sprague(x)

  expect_equal(length(result), 81)
  expect_equal(sum(result), sum(x), tolerance = 1)
})
