# Tests for heaping indices functions

library(testthat)

# =============================================================================
# Test Data Setup
# =============================================================================

# Helper function to create test data
create_test_ages <- function(seed = 42) {
  set.seed(seed)
  list(
    # Uniform distribution (no heaping)
    uniform = sample(20:70, 10000, replace = TRUE),
    # Heaping on 0 and 5 (5-year)
    heap5 = sample(seq(20, 70, by = 5), 5000, replace = TRUE),
    # Heaping on 0 only (10-year)
    heap10 = sample(seq(20, 70, by = 10), 5000, replace = TRUE),
    # Old ages for Kannisto/Jdanov tests
    old_uniform = sample(85:105, 3000, replace = TRUE),
    old_heaped = c(sample(85:105, 2000, replace = TRUE),
                   rep(c(90, 100), each = 300))
  )
}

# =============================================================================
# Myers' Index Tests
# =============================================================================

test_that("myers returns close to 0 for uniform data", {
  data <- create_test_ages()
  result <- myers(data$uniform)
  expect_true(result < 5)  # Should be very close to 0
})

test_that("myers returns high value for heaped data", {
  data <- create_test_ages()
  result <- myers(data$heap5)
  expect_true(result > 30)  # Should be high
})

test_that("myers validates input", {
  expect_error(myers(c("a", "b", "c")), "must be a numeric")
})

test_that("myers handles weights", {
  data <- create_test_ages()
  weights <- runif(length(data$uniform))
  expect_no_error(myers(data$uniform, weight = weights))
})

# =============================================================================
# Bachi's Index Tests
# =============================================================================

test_that("bachi returns close to 0 for uniform data", {
  data <- create_test_ages()
  result <- bachi(data$uniform)
  expect_true(result < 5)  # Should be very close to 0
})

test_that("bachi returns high value for heaped data", {
  data <- create_test_ages()
  result <- bachi(data$heap5)
  expect_true(result > 20)  # Should be high
})

test_that("bachi validates input", {
  expect_error(bachi(c("a", "b", "c")), "must be a numeric")
})

test_that("bachi handles weights", {
  data <- create_test_ages()
  weights <- runif(length(data$uniform))
  expect_no_error(bachi(data$uniform, weight = weights))
})

# =============================================================================
# Noumbissi's Index Tests
# =============================================================================

test_that("noumbissi returns close to 1 for uniform data", {
  data <- create_test_ages()
  result0 <- noumbissi(data$uniform, digit = 0)
  result5 <- noumbissi(data$uniform, digit = 5)

  expect_true(result0 > 0.8 && result0 < 1.2)
  expect_true(result5 > 0.8 && result5 < 1.2)
})

test_that("noumbissi returns high value for heaped digit", {
  data <- create_test_ages()
  result <- noumbissi(data$heap10, digit = 0)
  expect_true(result > 3)  # Should be much greater than 1
})

test_that("noumbissi validates digit parameter", {
  data <- create_test_ages()
  expect_error(noumbissi(data$uniform, digit = 10), "must be an integer from 0 to 9")
})

test_that("noumbissi validates input", {
  expect_error(noumbissi(c("a", "b", "c")), "must be a numeric")
})

# =============================================================================
# Spoorenberg's Total Modified Whipple Index Tests
# =============================================================================

test_that("spoorenberg returns close to 0 for uniform data", {
  data <- create_test_ages()
  result <- spoorenberg(data$uniform)
  expect_true(result < 2)  # Should be close to 0
})

test_that("spoorenberg returns high value for heaped data", {
  data <- create_test_ages()
  result <- spoorenberg(data$heap5)
  expect_true(result > 5)  # Should be high
})

test_that("spoorenberg validates input", {
  expect_error(spoorenberg(c("a", "b", "c")), "must be a numeric")
})

# =============================================================================
# Coale-Li Index Tests
# =============================================================================

test_that("coale_li returns close to 1 for uniform data", {
  set.seed(42)
  age <- sample(60:90, 5000, replace = TRUE)
  result <- coale_li(age, digit = 0)
  expect_true(result > 0.7 && result < 1.3)
})

test_that("coale_li returns high value for heaped data", {
  set.seed(42)
  age <- c(sample(60:90, 4000, replace = TRUE),
           rep(seq(60, 90, by = 10), each = 200))
  result <- coale_li(age, digit = 0)
  expect_true(result > 1.2)
})

test_that("coale_li validates digit parameter", {
  age <- sample(60:90, 1000, replace = TRUE)
  expect_error(coale_li(age, digit = 10), "must be an integer from 0 to 9")
})

test_that("coale_li validates input", {
  expect_error(coale_li(c("a", "b", "c")), "must be a numeric")
})

# =============================================================================
# Jdanov's Index Tests
# =============================================================================

test_that("jdanov returns close to 100 for uniform data", {
  data <- create_test_ages()
  result <- jdanov(data$old_uniform)
  expect_true(result > 80 && result < 120)
})

test_that("jdanov returns high value for heaped data", {
  data <- create_test_ages()
  result <- jdanov(data$old_heaped, Agei = c(90, 100))
  expect_true(result > 120)
})

test_that("jdanov validates input", {
  expect_error(jdanov(c("a", "b", "c")), "must be a numeric")
})

test_that("jdanov handles weights", {
  data <- create_test_ages()
  weights <- runif(length(data$old_uniform))
  expect_no_error(jdanov(data$old_uniform, weight = weights))
})

# =============================================================================
# Kannisto's Index Tests
# =============================================================================

test_that("kannisto returns close to 1 for uniform data", {
  data <- create_test_ages()
  result <- kannisto(data$old_uniform, Agei = 95)
  expect_true(result > 0.5 && result < 1.5)
})

test_that("kannisto returns high value for heaped data", {
  data <- create_test_ages()
  result <- kannisto(data$old_heaped, Agei = 90)
  expect_true(result > 1.2)
})

test_that("kannisto validates Agei parameter", {
  data <- create_test_ages()
  expect_error(kannisto(data$old_uniform, Agei = c(90, 95)), "must be a single age")
})

test_that("kannisto validates input", {
  expect_error(kannisto(c("a", "b", "c")), "must be a numeric")
})

# =============================================================================
# heaping_indices Convenience Function Tests
# =============================================================================

test_that("heaping_indices returns all indices", {
  data <- create_test_ages()
  result <- heaping_indices(data$uniform)

  expect_type(result, "list")
  expect_true("whipple_standard" %in% names(result))
  expect_true("whipple_modified" %in% names(result))
  expect_true("myers" %in% names(result))
  expect_true("bachi" %in% names(result))
  expect_true("spoorenberg" %in% names(result))
  expect_true("noumbissi_0" %in% names(result))
  expect_true("noumbissi_5" %in% names(result))
})

test_that("heaping_indices shows low values for uniform data", {
  data <- create_test_ages()
  result <- heaping_indices(data$uniform)

  # Standard Whipple should be close to 100
  expect_true(result$whipple_standard > 80 && result$whipple_standard < 120)

  # Modified Whipple should be close to 0
  expect_true(result$whipple_modified < 0.1)

  # Myers should be close to 0
  expect_true(result$myers < 5)

  # Noumbissi should be close to 1
  expect_true(result$noumbissi_0 > 0.8 && result$noumbissi_0 < 1.2)
})

test_that("heaping_indices shows high values for heaped data", {
  data <- create_test_ages()
  result <- heaping_indices(data$heap5)

  # Standard Whipple should be 500
  expect_equal(result$whipple_standard, 500)

  # Modified Whipple should be high
  expect_true(result$whipple_modified > 0.5)

  # Myers should be high
  expect_true(result$myers > 30)

  # Spoorenberg should be high
  expect_true(result$spoorenberg > 5)
})

test_that("heaping_indices handles weights", {
  data <- create_test_ages()
  weights <- runif(length(data$uniform))

  expect_no_error(heaping_indices(data$uniform, weight = weights))
})

# =============================================================================
# Weight Validation Tests (for all functions)
# =============================================================================

test_that("all functions validate weight length", {
  data <- create_test_ages()
  bad_weights <- 1:10  # Wrong length

  expect_error(myers(data$uniform, weight = bad_weights), "same length")
  expect_error(bachi(data$uniform, weight = bad_weights), "same length")
  expect_error(noumbissi(data$uniform, weight = bad_weights), "same length")
  expect_error(spoorenberg(data$uniform, weight = bad_weights), "same length")
  expect_error(coale_li(data$uniform, ageMin = 20, weight = bad_weights), "same length")
  expect_error(jdanov(data$old_uniform, weight = bad_weights), "same length")
  expect_error(kannisto(data$old_uniform, weight = bad_weights), "same length")
})

test_that("all functions validate weight type", {
  data <- create_test_ages()
  bad_weights <- as.character(1:length(data$uniform))

  expect_error(myers(data$uniform, weight = bad_weights), "must be a numeric")
  expect_error(bachi(data$uniform, weight = bad_weights), "must be a numeric")
  expect_error(noumbissi(data$uniform, weight = bad_weights), "must be a numeric")
})
