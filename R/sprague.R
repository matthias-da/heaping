#' Sprague Index (Multipliers)
#'
#' Disaggregate 5-year age group counts into single-year ages using Sprague
#' multipliers.
#'
#' @description
#' The Sprague method uses multipliers to estimate population counts for each
#' single year of age from 5-year interval data. This is useful for creating
#' smooth single-year age distributions from grouped census data.
#'
#' @details
#' The input must be population counts for 17 five-year age groups:
#' 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49,
#' 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, and 80+.
#'
#' The Sprague multipliers are applied differently depending on the position
#' of the age group:
#' \itemize{
#'   \item \strong{Lowest groups} (0-4): Uses only following age groups
#'   \item \strong{Low groups} (5-9): Uses mostly following age groups
#'   \item \strong{Normal groups} (10-74): Uses symmetric weighting
#'   \item \strong{High groups} (75-79): Uses mostly preceding age groups
#'   \item \strong{Highest groups} (80+): Returned as-is (open-ended)
#' }
#'
#' The total population is preserved: sum of output equals sum of input.
#'
#' @param x numeric vector of population counts in five-year age intervals.
#'   Must have exactly 17 elements corresponding to age groups 0-4, 5-9, ...,
#'   75-79, 80+.
#'
#' @return A named numeric vector with 81 elements: single-year population
#'   counts for ages 0, 1, 2, ..., 79, and the 80+ group.
#'
#' @author Matthias Templ
#'
#' @seealso \code{\link{whipple}} for measuring age heaping.
#'
#' @references
#' Calot, G. and Sardon, J.-P. (1998). \emph{Methodology for the calculation of
#' Eurostat's demographic indicators}. Detailed report by the European
#' Demographic Observatory.
#'
#' Sprague, T. B. (1880). Explanation of a new formula for interpolation.
#' \emph{Journal of the Institute of Actuaries}, \strong{22}, 270-285.
#'
#' @family graduation methods
#'
#' @export
#'
#' @examples
#' # Example from World Bank data
#' x <- data.frame(
#'   age = as.factor(c(
#'     "0-4", "5-9", "10-14", "15-19", "20-24",
#'     "25-29", "30-34", "35-39", "40-44", "45-49",
#'     "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+"
#'   )),
#'   pop = c(
#'     1971990, 2095820, 2157190, 2094110, 2116580,
#'     2003840, 1785690, 1502990, 1214170, 796934,
#'     627551, 530305, 488014, 364498, 259029, 158047, 125941
#'   )
#' )
#'
#' # Apply Sprague multipliers
#' s <- sprague(x$pop)
#' head(s, 20)  # First 20 single-year ages
#'
#' # Verify population is preserved
#' all.equal(sum(s), sum(x$pop))
#'
sprague <- function(x) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (length(x) != 17) {
    stop("Input must have exactly 17 elements for age groups: ",
         "0-4, 5-9, 10-14, ..., 70-74, 75-79, 80+")
  }

  breaks <- c(seq(0, 80, 5), 150)
  Ns <- x

  plus80 <- Ns[length(Ns)]

  # Sprague multipliers matrix
  multipliers <- data.frame(
    G1 = c(0.3616, 0.264, 0.184, 0.12, 0.0704,
           0.0336, 0.008, -0.008, -0.016, -0.0176,
           -0.0128, -0.0016, 0.0064, 0.0064, 0.0016,
           rep(0, 10)),
    G2 = c(-0.2768, -0.096, 0.04, 0.136, 0.1968,
           0.2272, 0.232, 0.216, 0.184, 0.1408,
           0.0848, 0.0144, -0.0336, -0.0416, -0.024,
           -0.0144, -0.008, 0, 0.008, 0.0144,
           0.0176, 0.016, 0.008, -0.008, -0.0336),
    G3 = c(0.1488, 0.04, -0.032, -0.072, -0.0848,
           -0.0752, -0.048, -0.008, 0.04, 0.0912,
           0.1504, 0.2224, 0.2544, 0.2224, 0.1504,
           0.0912, 0.04, -0.008, -0.048, -0.0752,
           -0.0848, -0.072, -0.032, 0.04, 0.1488),
    G4 = c(-0.0336, -0.008, 0.008, 0.016, 0.0176,
           0.0144, 0.008, 0, -0.008, -0.0144,
           -0.024, -0.0416, -0.0336, 0.0144, 0.0848,
           0.1408, 0.184, 0.216, 0.232, 0.2272,
           0.1968, 0.136, 0.04, -0.096, -0.2768),
    G5 = c(0, 0, 0, 0, 0,
           0, 0, 0, 0, 0,
           0.0016, 0.0064, 0.0064, -0.0016, -0.0128,
           -0.0176, -0.016, -0.008, 0.008, 0.0336,
           0.0704, 0.12, 0.184, 0.264, 0.3616)
  )

  multipliers <- cbind(
    groups = rep(c("lowest", "low", "normal", "high", "highest"), each = 5),
    multipliers
  )

  # Internal function to calculate population for a single year
  infoGroup <- function(n, mult = multipliers, mybreaks = breaks, popN = Ns) {
    # Which of the five years within a 5-year group
    groups <- n %% 5

    # Determine which multiplier set to use based on age
    if (n < 5) {
      tab <- subset(mult, subset = groups == "lowest")
    } else if (n >= 5 & n < 10) {
      tab <- subset(mult, subset = groups == "low")
    } else if (n >= 75 & n < 80) {
      tab <- subset(mult, subset = groups == "high")
    } else if (n >= 79) {
      tab <- subset(mult, subset = groups == "highest")
    } else {
      tab <- subset(mult, subset = groups == "normal")
    }

    # Determine which 5-year group this age belongs to
    ng <- cut(n, mybreaks, right = FALSE)
    mygroup <- which(levels(ng) %in% ng)

    # Get appropriate multiplier row
    rowsm <- tab[groups + 1, 2:ncol(tab)]

    # Determine which 5 adjacent groups to use
    if (mygroup == 1) {
      s <- seq(mygroup, mygroup + 4, 1)
    } else if (mygroup == 2) {
      s <- seq(mygroup - 1, mygroup + 3, 1)
    } else if (mygroup == 17) {
      s <- seq(mygroup - 4, mygroup, 1)
    } else if (mygroup == 16) {
      s <- seq(mygroup - 3, mygroup + 1, 1)
    } else {
      s <- seq(mygroup - 2, mygroup + 2, 1)
    }

    # Calculate interpolated population
    cohort <- sum(rowsm * popN[s])
    return(cohort)
  }

  # Apply to all ages 0-79

  cohorts <- sapply(0:79, infoGroup)

  # Add 80+ group
  cohorts <- c(cohorts, plus80)
  names(cohorts) <- c(0:79, "80+")

  return(cohorts)
}
