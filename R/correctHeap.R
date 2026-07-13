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
#' For method \dQuote{lnorm}, a truncated log-normal distribution is fit to the
#' \emph{trusted} (non-selected) records. Then for each age heap (at 0, 5, 10,
#' 15, ... or 0, 10, 20, ...) random numbers from a truncated log-normal
#' distribution (with lower and upper bounds) are drawn for the selected records.
#'
#' The correction range is controlled by \code{width}: replacement values are
#' drawn from a truncated distribution on \eqn{[\,\mathrm{heap}-\mathrm{width},\,
#' \mathrm{heap}+\mathrm{width}\,]}. The default half-width is 2 years for 5-year
#' and custom heaps and 5 years for 10-year heaps (a single symmetric window,
#' replacing the earlier two-stage \eqn{\pm 4} / \eqn{\pm 5} correction).
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
#' @param model optional formula (e.g. \code{age ~ x1 + x2}) for model-based
#'   correction. When provided, an imputation model for age given the covariates
#'   is fit on the \emph{trusted} (non-selected) records and used to draw
#'   covariate-conditional replacements for the selected heaped records. This
#'   replaces the earlier random-forest sign-adjustment heuristic. Requires
#'   \pkg{ranger} (for \code{model.engine = "ranger"}); \pkg{VIM} is used only
#'   when \code{dataModel} contains missing values.
#' @param dataModel data frame containing variables for the model formula.
#'   Required when \code{model} is specified. Missing values are imputed
#'   using k-nearest neighbors via \code{\link[VIM]{kNN}}.
#' @param model.engine character string selecting the conditional model engine
#'   when \code{model} is supplied: \code{"ranger"} (quantile regression forest,
#'   the default) or \code{"lm"} (linear model with a truncated-normal residual
#'   draw). Ignored when \code{model} is \code{NULL}.
#' @param weight optional numeric vector of survey weights, the same length as
#'   \code{x}. When supplied, heaping ratios and the selection of records to
#'   correct use weighted counts: records are drawn uniformly at random and
#'   accumulated until their cumulative weight covers the excess mass, and the
#'   weight is added as a predictor in the model-based arm
#'   (Reiter-Raghunathan-Zanutto). Defaults to \code{NULL} (unweighted).
#' @param width optional half-width (in years) of the truncation window around
#'   each heap. Defaults to 2 for 5-year and custom heaps and 5 for 10-year
#'   heaps.
#' @param case.weights logical; if \code{TRUE} and \code{weight} is supplied,
#'   the weights are additionally passed as case weights to the model engine.
#'   Defaults to \code{FALSE}.
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
#'     \item{n_fallback}{number of model-based draws that fell back to the
#'       marginal distribution because no fitted conditional quantile lay within
#'       the heap window (0 for the simple, non-model arm)}
#'     \item{weighted}{logical; whether weighted counts were used}
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
#' # Artificially introduce 5-year heaping by relabelling (misreporting)
#' # a share of records to the nearest multiple of five
#' k <- sample(length(age), round(0.3 * length(age)))
#' age5 <- age
#' age5[k] <- round(age5[k] / 5) * 5
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
                         model.engine = "ranger", weight = NULL, width = NULL,
                         case.weights = FALSE, seed = NULL, na.action = "omit",
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

  if (!is.null(weight)) {
    if (!is.numeric(weight) || length(weight) != length(x)) {
      stop("'weight' must be numeric and the same length as 'x'.")
    }
  }

  # --- NA handling: drop NAs, remember positions, align weight/dataModel/fixed ---
  na_idx <- which(is.na(x))
  if (length(na_idx) > 0 && na.action == "fail") {
    stop("NA values found in 'x'. Set na.action = 'omit' to handle them.")
  }
  keep    <- if (length(na_idx)) setdiff(seq_along(x), na_idx) else seq_along(x)
  xc      <- x[keep]
  wc      <- if (is.null(weight))    NULL else weight[keep]
  dm      <- if (is.null(dataModel)) NULL else dataModel[keep, , drop = FALSE]
  fixed_c <- if (is.null(fixed))     NULL else match(intersect(fixed, keep), keep)

  # --- heap sequence + truncation half-width ---
  if (is.character(heaps)) {
    if (!heaps %in% c("5year", "10year")) {
      stop("Unsupported value in argument 'heaps'. ",
           "Must be one of: '5year', '10year', or a numeric vector of heap positions")
    }
    heap_interval <- if (heaps == "10year") 10 else 5
    s <- seq(start, max(xc), by = heap_interval)
    if (is.null(width)) width <- if (heap_interval == 10) 5 else 2
  } else if (is.numeric(heaps)) {
    s <- sort(unique(heaps))
    s <- s[s >= 0 & s <= max(xc)]
    if (is.null(width)) width <- 2
  } else {
    stop("'heaps' must be a character ('5year', '10year') or numeric vector")
  }

  if (start > max(xc, na.rm = TRUE)) {
    stop("Starting year is greater than the maximum age in the data.")
  }

  # Warn if heap positions cover most unique values (likely misspecification)
  n_unique <- length(unique(xc))
  heap_coverage <- length(s[s %in% unique(xc)]) / n_unique
  if (heap_coverage > 0.5) {
    warning("More than 50% of unique values in 'x' are declared as heaps (",
            round(heap_coverage * 100), "%). ",
            "Heaping correction assumes heaps occur at sparse positions ",
            "(e.g., multiples of 5 or 10). Declaring most values as heaps ",
            "is likely a misspecification.")
  }

  # --- weighted heaping ratios + expected (local-smooth) counts ---
  rr     <- .heap_ratios(xc, s, wc)
  ratios <- rr$ratio
  maxage <- max(xc)

  # --- PASS 1: select records to correct (uniform draw, weight-based stopping) ---
  # excess = (weighted) count at the heap minus the expected count; with unit
  # weights the number selected reduces to ceiling(n - n/ratio), the legacy size.
  changes_by_heap <- setNames(rep(0L, length(s)), as.character(s))
  sel_idx <- integer(0); sel_low <- numeric(0); sel_up <- numeric(0)
  for (a in s) {
    hc <- as.character(a)
    rt <- ratios[hc]
    if (is.na(rt) || rt <= 1) next
    at   <- which(xc == a)
    cand <- if (is.null(fixed_c)) at else setdiff(at, fixed_c)
    if (!length(cand)) next
    excess <- rr$ratio[hc] * rr$expected[hc] - rr$expected[hc]
    wsel   <- if (is.null(wc)) NULL else wc[cand]
    chosen <- .select_to_correct(cand, wsel, excess)
    if (!length(chosen)) next
    sel_idx <- c(sel_idx, chosen)
    sel_low <- c(sel_low, rep(max(a - width, 0),      length(chosen)))
    sel_up  <- c(sel_up,  rep(min(a + width, maxage), length(chosen)))
    changes_by_heap[hc] <- length(chosen)
  }

  # --- PASS 2: fit on the trusted complement, then draw replacements ---
  result_c   <- xc
  n_fallback <- 0L
  fit_marg   <- NULL
  if (length(sel_idx)) {
    trusted  <- setdiff(seq_along(xc), sel_idx)
    fit_marg <- .fit_marginal(xc[trusted], method, sd)

    if (is.null(model)) {
      # simple (marginal) arm: draw per distinct window, vectorised
      key <- paste(sel_low, sel_up, sep = ":")
      for (g in unique(key)) {
        ix <- which(key == g)
        result_c[sel_idx[ix]] <- .draw_marginal(length(ix), fit_marg,
                                                 sel_low[ix[1]], sel_up[ix[1]])
      }
    } else {
      # model arm: covariate-conditional predictive draw fit on the trusted set
      if (!inherits(model, "formula")) {
        stop("'model' must be a valid formula.")
      }
      if (is.null(dm)) {
        stop("'dataModel' must be provided when 'model' is specified.")
      }
      if (any(is.na(dm))) {
        if (!requireNamespace("VIM", quietly = TRUE)) {
          stop("Package 'VIM' is required to impute missing covariates.")
        }
        dm <- VIM::kNN(dm, imp_var = FALSE)
      }
      wcol <- NULL
      if (!is.null(wc)) {
        dm[[".weight"]] <- wc
        wcol <- ".weight"
        message("Adding sampling weight as a predictor in the imputation model ",
                "(Reiter-Raghunathan-Zanutto).")
      }
      cw    <- if (isTRUE(case.weights) && !is.null(wc)) wc[trusted] else NULL
      fit_c <- .fit_cond(dm[trusted, , drop = FALSE], model, model.engine, wcol, cw)
      drawn <- .draw_cond(fit_c, dm[sel_idx, , drop = FALSE], sel_low, sel_up, fit_marg)
      n_fallback <- attr(drawn, "n_fallback")
      if (is.null(n_fallback)) n_fallback <- 0L
      result_c[sel_idx] <- as.numeric(drawn)
    }
  }

  # Restore NA values to original positions
  if (length(na_idx) > 0) {
    result <- rep(NA_real_, length(x))
    result[keep] <- result_c
  } else {
    result <- result_c
  }

  # Return result
  if (verbose) {
    list(
      corrected = result,
      n_changed = length(sel_idx),
      changes_by_heap = changes_by_heap[changes_by_heap > 0],
      ratios = ratios,
      method = method,
      seed = seed,
      fit_params = if (method == "norm" && !is.null(fit_marg)) list(sd = fit_marg$sd) else NULL,
      n_fallback = n_fallback,
      weighted = !is.null(weight)
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
#' @param weight optional numeric vector of survey weights, the same length as
#'   \code{x}. When supplied, the heaping ratio and the selection of records to
#'   correct use weighted counts. Defaults to \code{NULL} (unweighted).
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
#' # Artificially introduce a heap at age 23 by relabelling (misreporting)
#' # half of the nearby ages 21-25 to 23
#' k <- which(age %in% 21:25)
#' k <- sample(k, round(0.5 * length(k)))
#' age23 <- age
#' age23[k] <- 23
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
                              method = "lnorm", fixed = NULL, weight = NULL,
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
  if (!is.null(weight)) {
    if (!is.numeric(weight) || length(weight) != length(x)) {
      stop("'weight' must be numeric and the same length as 'x'.")
    }
  }

  # --- NA handling: drop NAs, remember positions, align weight + fixed ---
  na_idx <- which(is.na(x))
  if (length(na_idx) > 0 && na.action == "fail") {
    stop("NA values found in 'x'. Set na.action = 'omit' to handle them.")
  }
  keep    <- if (length(na_idx)) setdiff(seq_along(x), na_idx) else seq_along(x)
  xc      <- x[keep]
  wc      <- if (is.null(weight)) NULL else weight[keep]
  fixed_c <- if (is.null(fixed)) NULL else match(intersect(fixed, keep), keep)

  if (!heap %in% unique(xc)) {
    stop("Specified heap value is not present in 'x'.")
  }

  before <- round(before)
  after <- round(after)
  llow <- heap - before
  lup <- heap + after

  if (llow < 0 || lup > max(xc)) {
    stop("Parameters 'before' or 'after' result in bounds outside the data range.")
  }

  # --- weighted heaping ratio at this single heap ---
  rr    <- .heap_ratios(xc, heap, wc)
  hc    <- as.character(heap)
  ratio <- rr$ratio[hc]

  # --- PASS 1: select records at the heap to correct (uniform draw, weight stop) ---
  at     <- which(xc == heap)
  cand   <- if (is.null(fixed_c)) at else setdiff(at, fixed_c)
  chosen <- integer(0)
  if (!is.na(ratio) && ratio > 1 && length(cand)) {
    excess <- rr$ratio[hc] * rr$expected[hc] - rr$expected[hc]
    wsel   <- if (is.null(wc)) NULL else wc[cand]
    chosen <- .select_to_correct(cand, wsel, excess)
  }

  # --- PASS 2: fit on the trusted complement, then draw replacements ---
  result_c <- xc
  if (length(chosen)) {
    trusted  <- setdiff(seq_along(xc), chosen)
    fit_marg <- .fit_marginal(xc[trusted], method, sd)
    result_c[chosen] <- .draw_marginal(length(chosen), fit_marg, llow, lup)
  }
  n_changed <- length(chosen)

  # Restore NA values
  if (length(na_idx) > 0) {
    result <- rep(NA_real_, length(x))
    result[keep] <- result_c
  } else {
    result <- result_c
  }

  if (verbose) {
    list(
      corrected = result,
      n_changed = n_changed,
      ratio = ratio,
      method = method,
      seed = seed,
      weighted = !is.null(weight)
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
