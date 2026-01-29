#' Whipple Index (Original and Modified)
#'
#' Calculate the original or modified Whipple index to evaluate age heaping.
#'
#' @description
#' The Whipple index is a demographic measure used to detect and quantify age
#' heaping (digit preference) in population data. This function implements both
#' the original (standard) and modified versions of the index.
#'
#' @details
#' The original Whipple index is obtained by summing the number of persons in
#' the age range between 23 and 62, and calculating the ratio of reported ages
#' ending in 0 or 5 to one-fifth of the total sample. A linear decrease in the
#' number of persons of each age within the age range is assumed. Therefore,
#' low ages (0-22 years) and high ages (63 years and above) are excluded from
#' analysis since this assumption is not plausible.
#'
#' The original Whipple index ranges from:
#' \itemize{
#'   \item 0: when digits 0 and 5 are never reported
#'   \item 100: no preference for 0 or 5 (perfect data)
#'   \item 500: when only digits 0 and 5 are reported (maximum heaping)
#' }
#'
#' For the modified Whipple index, age heaping is calculated for all ten digits
#' (0-9). For each digit, the degree of preference or avoidance is determined,
#' and the modified Whipple index is given by the absolute sum of these
#' (indices - 1), scaled between 0 and 1:
#' \itemize{
#'   \item 0: ages are distributed perfectly equally across all digits
#'   \item 1: all age values end with the same digit
#' }
#'
#' @param x numeric vector holding the ages of persons.
#' @param method character string specifying which index to calculate:
#'   \describe{
#'     \item{\code{"standard"}}{Original Whipple index (default). Ranges 0-500,
#'       with 100 indicating no heaping.
#'     }
#'     \item{\code{"modified"}}{Modified Whipple index. Ranges 0-1, with 0
#'       indicating no heaping.
#'     }
#'   }
#' @param weight optional numeric vector holding the sampling weights of each
#'   person. Must be the same length as \code{x}. If \code{NULL} (default),
#'   unweighted counts are used.
#'
#' @return A single numeric value representing the Whipple index.
#'
#' @author Matthias Templ, Alexander Kowarik
#'
#' @seealso \code{\link{sprague}} for disaggregating 5-year age groups.
#'
#' @references
#' Shryock, H. S. and Siegel, J. S. (1976). \emph{The Methods and Materials of
#' Demography}. New York: Academic Press.
#'
#' Spoorenberg, T. and Dutreuilh, C. (2007). Quality of age reporting:
#' extension and application of the modified Whipple's index.
#' \emph{Population}, \strong{62}(4), 729-741.
#'
#' @family heaping indices
#'
#' @export
#'
#' @examples
#' # Equally distributed ages (no heaping)
#' set.seed(42)
#' age_uniform <- sample(1:100, 5000, replace = TRUE)
#' whipple(age_uniform)                    # Should be close to 100
#' whipple(age_uniform, method = "modified")  # Should be close to 0
#'
#' # Strong heaping on 5 and 10 (ages ending in 0 or 5 only)
#' age_5year <- sample(seq(0, 100, by = 5), 5000, replace = TRUE)
#' whipple(age_5year)                     # Should be 500
#' whipple(age_5year, method = "modified")   # Should be close to 0.8
#'
#' # Extreme heaping on 10 only (ages ending in 0 only)
#' age_10year <- sample(seq(0, 100, by = 10), 5000, replace = TRUE)
#' whipple(age_10year)                    # Should be 500
#' whipple(age_10year, method = "modified")  # Should be close to 1
#'
#' # Using weights
#' weights <- runif(5000)
#' whipple(age_uniform, weight = weights)
#'
whipple <- function(x, method = "standard", weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!method %in% c("standard", "modified")) {
    stop("Unsupported value in argument 'method'. ",
         "Must be one of: 'standard', 'modified'")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  if (method == "standard") {
    if (is.null(weight)) {
      x <- x[x >= 23 & x <= 62]
      xm <- x %% 5
      return((length(xm[xm == 0]) / length(x)) * 500)
    } else {
      weight <- weight[x >= 23 & x <= 62]
      x <- x[x >= 23 & x <= 62]
      xm <- x %% 5
      return((sum(weight[xm == 0]) / sum(weight)) * 500)
    }

  } else if (method == "modified") {
    if (is.null(weight)) {
      tab <- table(x)
    } else {
      tab <- .tableWt(x, weight)
    }

    W <- numeric(10)
    for (i in 1:10) {
      W[i] <- sum(tab[as.numeric(names(tab)) %in% seq(i - 10, 200, by = 10)]) /
        (sum(tab) / 10)
    }
    return(sum(abs(W - 1), na.rm = TRUE) / 18)
  }
}


#' Weighted Cross Tabulation (Internal)
#'
#' Compute weighted frequency table for a single variable.
#'
#' @param x numeric vector to tabulate.
#' @param weights numeric vector of sample weights.
#' @return A named numeric vector of weighted counts.
#' @keywords internal
#' @noRd
.tableWt <- function(x, weights) {
  if (is.null(weights)) {
    return(table(x))
  }

  # Simple weighted tabulation
  unique_vals <- sort(unique(x))
  tab <- sapply(unique_vals, function(v) sum(weights[x == v]))
  names(tab) <- as.character(unique_vals)
  tab <- round(tab)
  class(tab) <- "table"
  return(tab)
}
