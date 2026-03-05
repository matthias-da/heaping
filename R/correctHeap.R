#' Correct Age Heaping
#'
#' Correct for age heaping at regular intervals using truncated distributions.
#'
#' @description
#' Age heaping can cause substantial bias in important demographic measures
#' and thus should be corrected. This function corrects heaping at regular
#' intervals (every 5 or 10 years) by replacing a proportion of heaped
#' observations with draws from fitted truncated distributions.
#'
#' @details
#' For method \dQuote{lnorm}, a truncated log-normal distribution is fit to
#' the whole age distribution. Then for each age heap (at 0, 5, 10, 15, ...
#' or 0, 10, 20, ...) random numbers from a truncated log-normal distribution
#' (with lower and upper bounds) are drawn.
#'
#' The correction range depends on the heap type:
#' \itemize{
#'   \item For 5-year heaps: values are drawn from \eqn{\pm 2} years around the heap
#'   \item For 10-year heaps: values are drawn in two groups, \eqn{\pm 4} and
#'     \eqn{\pm 5} years around the heap
#' }
#'
#' The ratio of observations to replace is calculated by comparing the count
#' at each heap age to the arithmetic mean of the two neighboring ages. For
#' example, for age heap 5, the ratio is: count(age=5) / mean(count(age=4), count(age=6)).
#'
#' Method \dQuote{norm} uses truncated normal distributions instead. The choice
#' between \dQuote{lnorm} and \dQuote{norm} depends on whether the age
#' distribution is right-skewed (use \dQuote{lnorm}) or more symmetric
#' (use \dQuote{norm}). Many distributions with heaping problems are right-skewed.
#'
#' Method \dQuote{unif} draws from truncated uniform distributions around the
#' age heaps, providing a simpler baseline approach.
#'
#' Method \dQuote{kernel} uses kernel density estimation to sample replacement
#' values, providing a nonparametric alternative that adapts to the local
#' data distribution.
#'
#' Repeated calls of this function mimic multiple imputation, i.e., repeating
#' this procedure m times provides m corrected datasets that properly reflect
#' the uncertainty from the correction process. Use the \code{seed} parameter
#' to ensure reproducibility.
#'
#' @param x numeric vector of ages (typically integers).
#' @param heaps character string specifying the heaping pattern:
#'   \describe{
#'     \item{\code{"5year"}}{heaps are assumed every 5 years (0, 5, 10, 15, ...)}
#'     \item{\code{"10year"}}{heaps are assumed every 10 years (0, 10, 20, ...)}
#'   }
#'   Alternatively, a numeric vector specifying custom heap positions.
#' @param method character string specifying the distribution used for correction:
#'   \describe{
#'     \item{\code{"lnorm"}}{truncated log-normal distribution (default).
#'       Parameters are estimated from the input data.}
#'     \item{\code{"norm"}}{truncated normal distribution.
#'       Parameters are estimated from the input data.}
#'     \item{\code{"unif"}}{uniform distribution within the truncation bounds.}
#'     \item{\code{"kernel"}}{kernel density estimation for nonparametric sampling.}
#'   }
#' @param start numeric value for the starting point of the heap sequence
#'   (default 0). Use 5 if heaps occur at 5, 15, 25, ... instead of 0, 10, 20, ...
#'   Ignored if \code{heaps} is a numeric vector.
#' @param fixed numeric vector of indices indicating observations that should
#'   not be changed. Useful for preserving known accurate values.
#' @param model optional formula for model-based correction. When provided,
#'   a random forest model is fit to predict age from other variables, and
#'   the correction direction is adjusted to be consistent with this prediction.
#'   Requires packages \pkg{ranger} and \pkg{VIM}.
#' @param dataModel data frame containing variables for the model formula.
#'   Required when \code{model} is specified. Missing values are imputed
#'   using k-nearest neighbors via \code{\link[VIM]{kNN}}.
#' @param seed optional integer for random seed to ensure reproducibility.
#'   If \code{NULL} (default
#' ), no seed is set.
#' @param na.action character string specifying how to handle \code{NA} values:
#'   \describe{
#'     \item{\code{"omit"}}{remove NA values before processing, then restore positions (default)}
#'     \item{\code{"fail"}}{stop with an error if NA values are present}
#'   }
#' @param verbose logical. If \code{TRUE}, return a list with corrected values
#'   and diagnostic information. If \code{FALSE} (default), return only the
#'   corrected vector.
#' @param sd optional numeric value for standard deviation when \code{method = "norm"}.
#'   If \code{NULL} (default), estimated from the data using MAD (median absolute deviation)
#'   of non-heap ages, which is robust to the heaping.
#'
#' @return If \code{verbose = FALSE}, a numeric vector of the same length as
#'   \code{x} with heaping corrected. If \code{verbose = TRUE}, a list with:
#'   \describe{
#'     \item{corrected}{the corrected numeric vector}
#'     \item{n_changed}{total number of values changed}
#'     \item{changes_by_heap}{named vector of changes per heap age}
#'     \item{ratios}{named vector of heaping ratios per heap age}
#'     \item{method}{method used}
#'     \item{seed}{seed used (if any)}
#'   }
#'
#' @seealso \code{\link{correctSingleHeap}} for correcting a single specific heap.
#'
#' @family heaping correction
#'
#' @author Matthias Templ, Bernhard Meindl
#'
#' @references
#' Templ, M. (2026). Correction of heaping on individual level.
#' \emph{Journal TBD}.
#'
#' Templ, M., Meindl, B., Kowarik, A., Alfons, A., Dupriez, O. (2017).
#' Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Software}, \strong{79}(10), 1-38.
#' \doi{10.18637/jss.v079.i10}
#'
#' @export
#'
#' @examples
#' # Create artificial age data with log-normal distribution
#' set.seed(123)
#' age <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
#' age <- round(age[age < 93])
#'
#' # Artificially introduce 5-year heaping
#' year5 <- seq(0, max(age), 5)
#' age5 <- sample(c(age, age[age %in% year5]))
#'
#' # Correct with reproducible results
#' age5_corrected <- correctHeaps(age5, heaps = "5year", method = "lnorm", seed = 42)
#'
#' # Get diagnostic information
#' result <- correctHeaps(age5, heaps = "5year", verbose = TRUE, seed = 42)
#' print(result$n_changed)
#' print(result$ratios)
#'
#' # Use kernel method for nonparametric correction
#' age5_kernel <- correctHeaps(age5, heaps = "5year", method = "kernel", seed = 42)
#'
#' # Custom heap positions (e.g., heaping at 12, 18, 21)
#' custom_heaps <- c(12, 18, 21)
#' age_custom <- correctHeaps(age5, heaps = custom_heaps, method = "lnorm", seed = 42)
correctHeaps <- function(x, heaps = "10year", method = "lnorm", start = 0,
                         fixed = NULL, model = NULL, dataModel = NULL,
                         seed = NULL, na.action = "omit", verbose = FALSE,
                         sd = NULL) {


  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!method %in% c("lnorm", "norm", "unif", "kernel")) {
    stop("Unsupported value in argument 'method'. ",
         "Must be one of: 'lnorm', 'norm', 'unif', 'kernel'")
  }

  if (!na.action %in% c("omit", "fail")) {
    stop("'na.action' must be one of: 'omit', 'fail'")
  }

  # Handle NA values
  na_idx <- which(is.na(x))
  if (length(na_idx) > 0) {
    if (na.action == "fail") {
      stop("NA values found in 'x'. Set na.action = 'omit' to handle them.")
    }
    # Store original positions and remove NAs
    x_complete <- x[!is.na(x)]
    # Adjust fixed indices if needed
    if (!is.null(fixed)) {
      fixed <- fixed[!fixed %in% na_idx]
      # Recalculate indices after NA removal
      fixed <- sapply(fixed, function(i) i - sum(na_idx < i))
    }
  } else {
    x_complete <- x
  }

  # Handle heaps parameter - can be character or numeric vector

 if (is.character(heaps)) {
    if (!heaps %in% c("5year", "10year")) {
      stop("Unsupported value in argument 'heaps'. ",
           "Must be one of: '5year', '10year', or a numeric vector of heap positions")
    }
    heap_interval <- if (heaps == "10year") 10 else 5
    s <- seq(start, max(x_complete), by = heap_interval)
  } else if (is.numeric(heaps)) {
    # Custom heap positions
    s <- sort(unique(heaps))
    s <- s[s >= 0 & s <= max(x_complete)]
    heap_interval <- NA
  } else {
    stop("'heaps' must be a character ('5year', '10year') or numeric vector")
  }

  # Warn if heap positions cover most unique values
  n_unique <- length(unique(x_complete))
  heap_coverage <- length(s[s %in% unique(x_complete)]) / n_unique
  if (heap_coverage > 0.5) {
    warning("More than 50% of unique values in 'x' are declared as heaps (",
            round(heap_coverage * 100), "%). ",
            "Heaping correction assumes heaps occur at sparse positions ",
            "(e.g., multiples of 5 or 10). Declaring most values as heaps ",
            "is likely a misspecification.")
  }

  if (start > max(x_complete, na.rm = TRUE)) {
    stop("Starting year is greater than the maximum age in the data.")
  }

  # Model-based correction setup
  pred <- NULL
  if (!is.null(model)) {
    if (!inherits(model, "formula")) {
      stop("'model' must be a valid formula.")
    }
    if (is.null(dataModel)) {
      stop("'dataModel' must be provided when 'model' is specified.")
    }
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Package 'ranger' is required for model-based correction.")
    }
    if (!requireNamespace("VIM", quietly = TRUE)) {
      stop("Package 'VIM' is required for model-based correction.")
    }
    if (any(is.na(dataModel))) {
      dataModel <- VIM::kNN(dataModel, imp_var = FALSE)
    }
    rf <- ranger::ranger(formula = model, data = dataModel)
    pred <- predict(object = rf, data = dataModel, type = "response")
  }

  xorig <- x_complete

  # Build complete age frequency table
  tab <- table(x_complete)
  complete_ages <- seq(0, max(x_complete, na.rm = TRUE))
  complete_tab <- setNames(rep(0, length(complete_ages)), as.character(complete_ages))
  complete_tab[names(tab)] <- as.numeric(tab)

  # Calculate replacement ratios using named indexing for clarity
  calc_ratio <- function(heap_age, freq_table) {
    heap_char <- as.character(heap_age)
    before_char <- as.character(heap_age - 1)
    after_char <- as.character(heap_age + 1)

    if (!before_char %in% names(freq_table) || !after_char %in% names(freq_table)) {
      return(NA)
    }

    neighbor_mean <- mean(c(freq_table[before_char], freq_table[after_char]))
    if (neighbor_mean == 0) return(NA)

    freq_table[heap_char] / neighbor_mean
  }

  ratios <- sapply(s, calc_ratio, freq_table = complete_tab)
  names(ratios) <- as.character(s)

  # Fit distribution parameters
  fit_params <- list()
  if (method == "lnorm") {
    age0 <- as.numeric(x_complete)
    age0[age0 == 0] <- 0.01
    fit_params$logn <- fitdistrplus::fitdist(age0, "lnorm")
  } else if (method == "norm") {
    # Estimate sd from non-heap ages using MAD (robust to heaping)
    non_heap_ages <- x_complete[!x_complete %in% s]
    if (is.null(sd)) {
      if (length(non_heap_ages) > 10) {
        fit_params$sd <- stats::mad(non_heap_ages, constant = 1.4826)
      } else {
        fit_params$sd <- stats::sd(x_complete)
      }
    } else {
      fit_params$sd <- sd
    }
  } else if (method == "kernel") {
    # For kernel method, estimate density from non-heap ages
    non_heap_ages <- x_complete[!x_complete %in% s]
    if (length(non_heap_ages) > 10) {
      fit_params$density <- stats::density(non_heap_ages,
                                           from = 0,
                                           to = max(x_complete),
                                           n = 512)
    } else {
      # Fallback to full data if not enough non-heap observations
      fit_params$density <- stats::density(x_complete,
                                           from = 0,
                                           to = max(x_complete),
                                           n = 512)
    }
  }

  # Track changes for verbose output
  changes_by_heap <- setNames(rep(0L, length(s)), as.character(s))

  # Determine bounds based on heap type
  get_bounds <- function(heap_age, heap_interval, max_val) {
    if (is.na(heap_interval)) {
      # Custom heaps: use ±2 by default
      list(
        list(lower = max(heap_age - 2, 0), upper = min(heap_age + 2, max_val))
      )
    } else if (heap_interval == 10) {
      # 10-year heaps: two groups
      list(
        list(lower = max(heap_age - 4, 0), upper = min(heap_age + 4, max_val)),
        list(lower = max(heap_age - 5, 0), upper = min(heap_age + 5, max_val))
      )
    } else {
      # 5-year heaps: single group
      list(
        list(lower = max(heap_age - 2, 0), upper = min(heap_age + 2, max_val))
      )
    }
  }

  # Process each heap
  for (j in seq_along(s)) {
    i <- s[j]
    ratio <- ratios[j]
    index <- which(xorig == i)

    # Skip if no heaping detected (ratio <= 1) or no observations
    if (is.na(ratio) || ratio <= 1 || length(index) == 0) {
      next
    }

    # Calculate sample size
    ssize <- ceiling(length(index) - length(index) / ratio)
    if (ssize <= 0) next

    # Get bounds for this heap
    bounds_list <- get_bounds(i, heap_interval, max(x_complete))

    # Handle fixed observations
    available_idx <- if (is.null(fixed)) {
      index
    } else {
      index[!index %in% fixed]
    }

    if (length(available_idx) == 0) {
      warning("No suitable observations to change at heap ", i)
      next
    }

    # Adjust sample size if not enough available
    ssize <- min(ssize, length(available_idx))

    if (length(bounds_list) == 1) {
      # Single bound group (5-year heaps or custom)
      r <- available_idx[sample.int(length(available_idx), size = ssize)]
      x_complete[r] <- .draw_replacements_v2(
        n = length(r),
        method = method,
        fit_params = fit_params,
        center = i,
        llow = bounds_list[[1]]$lower,
        lup = bounds_list[[1]]$upper
      )
      changes_by_heap[as.character(i)] <- length(r)

    } else {
      # Two bound groups (10-year heaps)
      size1 <- ceiling(ssize / 2)
      size2 <- ssize - size1

      r1 <- available_idx[sample.int(length(available_idx), size = min(size1, length(available_idx)))]

      x_complete[r1] <- .draw_replacements_v2(
        n = length(r1),
        method = method,
        fit_params = fit_params,
        center = i,
        llow = bounds_list[[1]]$lower,
        lup = bounds_list[[1]]$upper
      )

      remaining_idx <- setdiff(available_idx, r1)
      if (length(remaining_idx) > 0 && size2 > 0) {
        r2 <- remaining_idx[sample.int(length(remaining_idx), size = min(size2, length(remaining_idx)))]

        x_complete[r2] <- .draw_replacements_v2(
          n = length(r2),
          method = method,
          fit_params = fit_params,
          center = i,
          llow = bounds_list[[2]]$lower,
          lup = bounds_list[[2]]$upper
        )

        changes_by_heap[as.character(i)] <- length(r1) + length(r2)
      } else {
        changes_by_heap[as.character(i)] <- length(r1)
      }
    }
  }

  # Apply model-based sign adjustment if requested
  if (!is.null(model) && !is.null(pred)) {
    x_complete <- .adjust_signs(xorig, x_complete, pred$predictions)
  }

  # Restore NA values to original positions
  if (length(na_idx) > 0) {
    result <- rep(NA_real_, length(x))
    result[-na_idx] <- x_complete
  } else {
    result <- x_complete
  }

  # Return result
  if (verbose) {
    list(
      corrected = result,
      n_changed = sum(changes_by_heap),
      changes_by_heap = changes_by_heap[changes_by_heap > 0],
      ratios = ratios,
      method = method,
      seed = seed,
      fit_params = if (method == "norm") list(sd = fit_params$sd) else NULL
    )
  } else {
    result
  }
}

#' @rdname correctHeaps
#' @export
correctHeaps2 <- correctHeaps


#' Correct a Single Age Heap
#'
#' Correct a specific age heap in a vector containing ages.
#'
#' @description
#' While \code{\link{correctHeaps}} corrects regular heaping patterns,
#' this function allows correction of a single specific heap value.
#' This is useful when heaping occurs at irregular intervals or when
#' only a particular age shows excessive heaping.
#'
#' @param x numeric vector representing ages (typically integers).
#' @param heap numeric value specifying the age for which heaping should
#'   be corrected. Must be present in \code{x}.
#' @param before numeric value specifying the number of years before
#'   the heap to use as the lower bound for replacement values.
#'   Will be rounded to an integer. Default is 2.
#' @param after numeric value specifying the number of years after
#'   the heap to use as the upper bound for replacement values.
#'   Will be rounded to an integer. Default is 2.
#' @param method character string specifying the distribution used for correction:
#'   \describe{
#'     \item{\code{"lnorm"}}{truncated log-normal distribution (default).
#'       Parameters are estimated from the input data.}
#'     \item{\code{"norm"}}{truncated normal distribution.
#'       Parameters are estimated from the input data.}
#'     \item{\code{"unif"}}{uniform distribution within the truncation bounds.}
#'     \item{\code{"kernel"}}{kernel density estimation for nonparametric sampling.}
#'   }
#' @param fixed numeric vector of indices indicating observations that should
#'   not be changed. Useful for preserving known accurate values.
#' @param seed optional integer for random seed to ensure reproducibility.
#' @param na.action character string specifying how to handle \code{NA} values:
#'   \code{"omit"} (default) or \code{"fail"}.
#' @param verbose logical. If \code{TRUE}, return diagnostic information.
#' @param sd optional numeric value for standard deviation when \code{method = "norm"}.
#'
#' @return A numeric vector of the same length as \code{x} with the specified
#'   heap corrected, or a list with diagnostics if \code{verbose = TRUE}.
#'
#' @seealso \code{\link{correctHeaps}} for correcting regular heaping patterns.
#'
#' @family heaping correction
#'
#' @author Matthias Templ, Bernhard Meindl
#'
#' @export
#'
#' @examples
#' # Create artificial age data
#' set.seed(123)
#' age <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
#' age <- round(age[age < 93])
#'
#' # Artificially introduce a heap at age 23
#' age23 <- c(age, rep(23, length = sum(age == 23)))
#'
#' # Correct with reproducible results
#' age23_corrected <- correctSingleHeap(age23, heap = 23, before = 5, after = 5,
#'                                      method = "lnorm", seed = 42)
#'
#' # Get diagnostic information
#' result <- correctSingleHeap(age23, heap = 23, before = 5, after = 5,
#'                             verbose = TRUE, seed = 42)
#' print(result$n_changed)
correctSingleHeap <- function(x, heap, before = 2, after = 2,
                              method = "lnorm", fixed = NULL,
                              seed = NULL, na.action = "omit",
                              verbose = FALSE, sd = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!method %in% c("lnorm", "norm", "unif", "kernel")) {
    stop("Unsupported value in argument 'method'. ",
         "Must be one of: 'lnorm', 'norm', 'unif', 'kernel'")
  }

  if (!na.action %in% c("omit", "fail")) {
    stop("'na.action' must be one of: 'omit', 'fail'")
  }

  # Handle NA values
  na_idx <- which(is.na(x))
  if (length(na_idx) > 0) {
    if (na.action == "fail") {
      stop("NA values found in 'x'. Set na.action = 'omit' to handle them.")
    }
    x_complete <- x[!is.na(x)]
    if (!is.null(fixed)) {
      fixed <- fixed[!fixed %in% na_idx]
      fixed <- sapply(fixed, function(i) i - sum(na_idx < i))
    }
  } else {
    x_complete <- x
  }

  if (!heap %in% unique(x_complete)) {
    stop("Specified heap value is not present in 'x'.")
  }
  if (length(heap) != 1) {
    stop("Only one heap value can be specified at a time.")
  }
  if (!is.numeric(before) || !is.numeric(after)) {
    stop("Arguments 'before' and 'after' must be numeric.")
  }
  if (length(before) != 1 || length(after) != 1) {
    stop("Arguments 'before' and 'after' must be single values.")
  }
  if (before < 0 || after < 0) {
    stop("Arguments 'before' and 'after' must be non-negative.")
  }

  before <- round(before)
  after <- round(after)

  llow <- heap - before
  lup <- heap + after

  if (llow < 0 || lup > max(x_complete)) {
    stop("Parameters 'before' or 'after' result in bounds outside the data range.")
  }

  # Calculate ratio
  tab <- table(x_complete)
  complete_ages <- seq(0, max(x_complete))
  complete_tab <- setNames(rep(0, length(complete_ages)), as.character(complete_ages))
  complete_tab[names(tab)] <- as.numeric(tab)

  heap_char <- as.character(heap)
  before_char <- as.character(heap - 1)
  after_char <- as.character(heap + 1)

  neighbor_mean <- mean(c(complete_tab[before_char], complete_tab[after_char]))
  ratio <- if (neighbor_mean > 0) complete_tab[heap_char] / neighbor_mean else NA

  # Fit distribution parameters
  fit_params <- list()
  if (method == "lnorm") {
    age0 <- as.numeric(x_complete)
    age0[age0 == 0] <- 0.01
    fit_params$logn <- fitdistrplus::fitdist(age0, "lnorm")
  } else if (method == "norm") {
    non_heap_ages <- x_complete[x_complete != heap]
    if (is.null(sd)) {
      if (length(non_heap_ages) > 10) {
        fit_params$sd <- stats::mad(non_heap_ages, constant = 1.4826)
      } else {
        fit_params$sd <- stats::sd(x_complete)
      }
    } else {
      fit_params$sd <- sd
    }
  } else if (method == "kernel") {
    non_heap_ages <- x_complete[x_complete != heap]
    if (length(non_heap_ages) > 10) {
      fit_params$density <- stats::density(non_heap_ages, from = 0, to = max(x_complete), n = 512)
    } else {
      fit_params$density <- stats::density(x_complete, from = 0, to = max(x_complete), n = 512)
    }
  }

  # Find indices to potentially modify
  index <- which(x_complete == heap)

  if (is.na(ratio) || ratio <= 1) {
    ssize <- 0
  } else {
    ssize <- ceiling(length(index) - length(index) / ratio)
  }

  n_changed <- 0
  if (ssize > 0) {
    available_idx <- if (is.null(fixed)) {
      index
    } else {
      index[!index %in% fixed]
    }

    if (length(available_idx) == 0) {
      warning("No suitable observations to change at heap ", heap)
    } else {
      ssize <- min(ssize, length(available_idx))
      r <- available_idx[sample.int(length(available_idx), size = ssize)]

      x_complete[r] <- .draw_replacements_v2(
        n = length(r),
        method = method,
        fit_params = fit_params,
        center = heap,
        llow = llow,
        lup = lup
      )
      n_changed <- length(r)
    }
  }

  # Restore NA values
  if (length(na_idx) > 0) {
    result <- rep(NA_real_, length(x))
    result[-na_idx] <- x_complete
  } else {
    result <- x_complete
  }

  if (verbose) {
    list(
      corrected = result,
      n_changed = n_changed,
      ratio = ratio,
      method = method,
      seed = seed
    )
  } else {
    result
  }
}


# Internal helper function to draw replacement values (version 2)
#
# @param n number of values to draw
# @param method distribution method
# @param fit_params list of fitted parameters
# @param center center value for normal distribution
# @param llow lower bound
# @param lup upper bound
# @return numeric vector of rounded replacement values
# @keywords internal
.draw_replacements_v2 <- function(n, method, fit_params, center, llow, lup) {
  if (n == 0) return(numeric(0))

  if (method == "lnorm") {
    round(EnvStats::rlnormTrunc(n,
                                meanlog = fit_params$logn$estimate[1],
                                sdlog = as.numeric(fit_params$logn$estimate[2]),
                                min = llow, max = lup))
  } else if (method == "norm") {
    round(EnvStats::rnormTrunc(n,
                               mean = center,
                               sd = fit_params$sd,
                               min = llow, max = lup))
  } else if (method == "kernel") {
    # Sample from kernel density within bounds
    .sample_from_density(n, fit_params$density, llow, lup)
  } else {
    # uniform
    sample(llow:lup, n, replace = TRUE)
  }
}


# Internal helper function to sample from kernel density
#
# @param n number of values to draw
# @param dens density object from stats::density
# @param llow lower bound
# @param lup upper bound
# @return numeric vector of rounded values sampled from the density
# @keywords internal
.sample_from_density <- function(n, dens, llow, lup) {
  # Subset density to the bounds
  in_range <- dens$x >= llow & dens$x <= lup
  x_sub <- dens$x[in_range]
  y_sub <- dens$y[in_range]

  if (length(x_sub) < 2) {
    # Fallback to uniform if density doesn't cover range
    return(sample(llow:lup, n, replace = TRUE))
  }

  # Normalize density in range
  y_sub <- y_sub / sum(y_sub)

  # Sample from the density using inverse CDF
  sampled <- sample(x_sub, size = n, replace = TRUE, prob = y_sub)

  # Add small noise and round to integers
  sampled <- sampled + stats::runif(n, -0.5, 0.5)
  sampled <- pmax(llow, pmin(lup, round(sampled)))

  as.integer(sampled)
}


# Internal helper function for model-based sign adjustment
#
# @param xorig original values before correction
# @param x corrected values
# @param predictions model predictions
# @return adjusted numeric vector
# @keywords internal
.adjust_signs <- function(xorig, x, predictions) {
  w <- !(xorig == x)
  if (sum(w) == 0) return(x)

  signs <- (xorig - x)[w] > 0
  signsModel <- (xorig - predictions)[w] > 0
  changesigns <- signs != signsModel

  difference <- x[w] - xorig[w]
  difference[changesigns] <- -difference[changesigns]

  adjusted_x <- x
  adjusted_x[w] <- xorig[w] + difference

  adjusted_x
}


# Legacy function for backward compatibility
# Uses old algorithm to draw replacements
.draw_replacements <- function(n, method, logn, center, llow, lup) {
  if (method == "lnorm") {
    round(EnvStats::rlnormTrunc(n,
                                meanlog = logn$estimate[1],
                                sdlog = as.numeric(logn$estimate[2]),
                                min = llow, max = lup))
  } else if (method == "norm") {
    round(EnvStats::rnormTrunc(n,
                               mean = center, sd = 1,
                               min = llow, max = lup))
  } else {
    sample(llow:lup, n, replace = TRUE)
  }
}
