# Tests for correctHeaps function
# Compares new implementation against expected behavior

library(testthat)

# Helper function to create test data
create_test_data <- function(seed = 123) {
  set.seed(seed)
  age <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
  age <- round(age[age < 93])

  # Introduce 5-year heaping
  year5 <- seq(0, max(age), 5)
  age5 <- sample(c(age, age[age %in% year5]))

  list(original = age, heaped = age5)
}

# =============================================================================
# Basic Functionality Tests
# =============================================================================

test_that("correctHeaps returns vector of same length", {
  data <- create_test_data()
  result <- correctHeaps(data$heaped, heaps = "5year", seed = 42)
  expect_equal(length(result), length(data$heaped))
})

test_that("correctHeaps with seed is reproducible", {
  data <- create_test_data()

  result1 <- correctHeaps(data$heaped, heaps = "5year", seed = 42)
  result2 <- correctHeaps(data$heaped, heaps = "5year", seed = 42)

expect_identical(result1, result2)
})

test_that("correctHeaps without seed produces different results", {
  data <- create_test_data()

  result1 <- correctHeaps(data$heaped, heaps = "5year")
  result2 <- correctHeaps(data$heaped, heaps = "5year")

  # Results should differ (with very high probability)
  expect_false(identical(result1, result2))
})

test_that("correctHeaps reduces heaping", {
  data <- create_test_data()
  year5 <- seq(0, max(data$heaped), 5)

  # Count heaped values before
  heaped_count_before <- sum(data$heaped %in% year5)

  result <- correctHeaps(data$heaped, heaps = "5year", seed = 42)

  # Count heaped values after
  heaped_count_after <- sum(result %in% year5)

  expect_lt(heaped_count_after, heaped_count_before)
})

# =============================================================================
# Method Tests
# =============================================================================

test_that("all methods work without error", {
  data <- create_test_data()

  expect_no_error(correctHeaps(data$heaped, heaps = "5year", method = "lnorm", seed = 42))
  expect_no_error(correctHeaps(data$heaped, heaps = "5year", method = "norm", seed = 42))
  expect_no_error(correctHeaps(data$heaped, heaps = "5year", method = "unif", seed = 42))
  expect_no_error(correctHeaps(data$heaped, heaps = "5year", method = "kernel", seed = 42))
})

test_that("norm method uses estimated sd when not provided", {
  data <- create_test_data()

  result <- correctHeaps(data$heaped, heaps = "5year", method = "norm",
                         verbose = TRUE, seed = 42)

  expect_true(!is.null(result$fit_params$sd))
  expect_true(result$fit_params$sd > 0)
})

test_that("norm method uses provided sd", {
  data <- create_test_data()

  result <- correctHeaps(data$heaped, heaps = "5year", method = "norm",
                         sd = 5, verbose = TRUE, seed = 42)

  expect_equal(result$fit_params$sd, 5)
})

# =============================================================================
# NA Handling Tests
# =============================================================================

test_that("NA values are handled with na.action='omit'", {
  data <- create_test_data()
  data_with_na <- data$heaped
  data_with_na[c(1, 50, 100)] <- NA

  result <- correctHeaps(data_with_na, heaps = "5year", na.action = "omit", seed = 42)

  # Length should be preserved
  expect_equal(length(result), length(data_with_na))

  # NA positions should be preserved
  expect_true(all(is.na(result[c(1, 50, 100)])))

  # Non-NA values should be numeric
  expect_true(all(!is.na(result[-c(1, 50, 100)])))
})

test_that("NA values cause error with na.action='fail'", {
  data <- create_test_data()
  data_with_na <- data$heaped
  data_with_na[1] <- NA

  expect_error(
    correctHeaps(data_with_na, heaps = "5year", na.action = "fail"),
    "NA values found"
  )
})

# =============================================================================
# Verbose Output Tests
# =============================================================================

test_that("verbose=TRUE returns list with diagnostics", {
  data <- create_test_data()

  result <- correctHeaps(data$heaped, heaps = "5year", verbose = TRUE, seed = 42)

  expect_type(result, "list")
  expect_true("corrected" %in% names(result))
  expect_true("n_changed" %in% names(result))
  expect_true("changes_by_heap" %in% names(result))
  expect_true("ratios" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("seed" %in% names(result))
})

test_that("verbose output has correct values", {
  data <- create_test_data()

  result <- correctHeaps(data$heaped, heaps = "5year", verbose = TRUE, seed = 42)

  expect_equal(length(result$corrected), length(data$heaped))
  expect_true(result$n_changed > 0)
  expect_equal(result$method, "lnorm")
  expect_equal(result$seed, 42)
})

# =============================================================================
# Custom Heaps Tests
# =============================================================================

test_that("custom heap positions work", {
  data <- create_test_data()

  # Create custom heaps at specific ages
  custom <- c(20, 30, 40, 50)

  result <- correctHeaps(data$heaped, heaps = custom, seed = 42)

  expect_equal(length(result), length(data$heaped))
})

test_that("custom heaps only affect specified positions", {
  data <- create_test_data()

  custom <- c(25)  # Only correct age 25
  result_custom <- correctHeaps(data$heaped, heaps = custom, verbose = TRUE, seed = 42)

  # Only age 25 should be in the changes
  expect_true(all(names(result_custom$changes_by_heap) == "25"))
})

# =============================================================================
# Fixed Observations Tests
# =============================================================================

test_that("fixed observations are not changed", {
  data <- create_test_data()
  fixed_idx <- 1:100

  result <- correctHeaps(data$heaped, heaps = "5year", fixed = fixed_idx, seed = 42)

  # First 100 observations should be unchanged
  expect_equal(result[fixed_idx], data$heaped[fixed_idx])
})

# =============================================================================
# Heap Type Tests
# =============================================================================

test_that("10year heaps work correctly", {
  data <- create_test_data()

  result <- correctHeaps(data$heaped, heaps = "10year", seed = 42)

  expect_equal(length(result), length(data$heaped))
})

test_that("start parameter shifts heap positions", {
  data <- create_test_data()

  result_0 <- correctHeaps(data$heaped, heaps = "5year", start = 0, verbose = TRUE, seed = 42)
  result_2 <- correctHeaps(data$heaped, heaps = "5year", start = 2, verbose = TRUE, seed = 42)

  # Different starting points should affect different ages
  expect_false(identical(names(result_0$ratios), names(result_2$ratios)))
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("invalid method causes error", {
  data <- create_test_data()

  expect_error(
    correctHeaps(data$heaped, method = "invalid"),
    "Unsupported value"
  )
})

test_that("invalid heaps causes error", {
  data <- create_test_data()

  expect_error(
    correctHeaps(data$heaped, heaps = "invalid"),
    "Unsupported value"
  )
})

test_that("non-numeric input causes error", {
  expect_error(
    correctHeaps(c("a", "b", "c")),
    "must be a numeric"
  )
})

# =============================================================================
# correctSingleHeap Tests
# =============================================================================

test_that("correctSingleHeap works correctly", {
  set.seed(123)
  age <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
  age <- round(age[age < 93])
  age23 <- c(age, rep(23, length = sum(age == 23)))

  result <- correctSingleHeap(age23, heap = 23, before = 5, after = 5, seed = 42)

  expect_equal(length(result), length(age23))

  # Count at age 23 should decrease
  expect_lt(sum(result == 23), sum(age23 == 23))
})

test_that("correctSingleHeap verbose mode works", {
  set.seed(123)
  age <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
  age <- round(age[age < 93])
  age23 <- c(age, rep(23, length = sum(age == 23)))

  result <- correctSingleHeap(age23, heap = 23, before = 5, after = 5,
                              verbose = TRUE, seed = 42)

  expect_type(result, "list")
  expect_true("n_changed" %in% names(result))
  expect_true(result$n_changed > 0)
})

test_that("correctSingleHeap respects fixed observations", {
  set.seed(123)
  age <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
  age <- round(age[age < 93])
  age23 <- c(age, rep(23, length = sum(age == 23)))

  fixed_idx <- which(age23 == 23)[1:10]

  result <- correctSingleHeap(age23, heap = 23, before = 5, after = 5,
                              fixed = fixed_idx, seed = 42)

  expect_equal(result[fixed_idx], age23[fixed_idx])
})

# =============================================================================
# Edge Case Tests
# =============================================================================

test_that("handles data with no heaping gracefully", {
  # Create uniform data with no heaping
  set.seed(42)
  uniform_ages <- sample(1:80, 1000, replace = TRUE)

  result <- correctHeaps(uniform_ages, heaps = "5year", verbose = TRUE, seed = 42)

  # Should still return valid output
  expect_equal(length(result$corrected), length(uniform_ages))
})

test_that("handles small datasets", {
  small_data <- c(20, 25, 30, 25, 25, 35, 40)

  expect_no_error(suppressWarnings(correctHeaps(small_data, heaps = "5year", seed = 42)))
})

test_that("handles single value datasets", {
  single_val <- rep(30, 100)

  # Should not error but may warn
  expect_no_error(suppressWarnings(correctHeaps(single_val, heaps = "5year", seed = 42)))
})
