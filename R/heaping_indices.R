#' Myers' Blended Index of Age Heaping
#'
#' Calculate Myers' blended index to measure digit preference in age data.
#'
#' @description
#' Myers' index measures preferences for each of the ten possible terminal
#' digits (0-9) as a blended index. It is based on the principle that in the
#' absence of age heaping, the aggregate population of each age ending in one
#' of the digits 0 to 9 should represent 10 percent of the total population.
#'
#' @details
#' The index uses a blending technique that weights earlier ages more for
#' digit preference calculation and later ages more for avoidance, creating
#' a balanced measure across the age range.
#'
#' The theoretical range is 0 to 90:
#' \itemize{
#'   \item 0: no digit preference (perfect data)
#'   \item 90: all ages reported with same terminal digit (maximum heaping)
#' }
#'
#' @param x numeric vector of individual ages.
#' @param ageMin minimum age to include (default 23).
#' @param ageMax maximum age to include (default 82).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing Myers' blended index.
#'
#' @author Matthias Templ
#'
#' @references
#' Myers, R. J. (1940). Errors and bias in the reporting of ages in census
#' data. \emph{Transactions of the Actuarial Society of America}, \strong{41},
#' 395-415.
#'
#' Myers, R. J. (1954). Accuracy of age reporting in the 1950 United States
#' Census. \emph{Journal of the American Statistical Association},
#' \strong{49}(268), 826-831.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{bachi}} for Bachi's index,
#'   \code{\link{whipple}} for Whipple's index.
#'
#' @export
#'
#' @examples
#' # No heaping (uniform ages)
#' set.seed(42)
#' age_uniform <- sample(23:82, 10000, replace = TRUE)
#' myers(age_uniform)  # Should be close to 0
#'
#' # Strong heaping on ages ending in 0 or 5
#' age_heaped <- sample(seq(25, 80, by = 5), 5000, replace = TRUE)
#' myers(age_heaped)  # Should be high
#'
myers <- function(x, ageMin = 23, ageMax = 82, weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Ensure age range is a multiple of 10
  diff <- ageMax - ageMin
  age_interval <- diff - (diff + 1) %% 10
  ageMax_used <- ageMin + age_interval

  # Filter to age range
  ind <- x >= ageMin & x <= ageMax_used
  if (sum(ind) == 0) {
    warning("No ages in specified range.")
    return(NA)
  }

  ages <- x[ind]
  if (!is.null(weight)) {
    wts <- weight[ind]
  } else {
    wts <- rep(1, length(ages))
  }

  # Create weighted frequency table
  digits <- ages %% 10
  period <- 10

  # Blending weights
  w <- 1:10
  names(w) <- as.character(0:9)

  # Calculate blended sums for each digit
  tab1 <- numeric(10)
  tab2 <- numeric(10)
  names(tab1) <- names(tab2) <- as.character(0:9)

  # Group ages into decades
  decades <- (ages - ageMin) %/% 10

  for (d in 0:9) {
    digit_mask <- digits == d
    for (dec in unique(decades)) {
      dec_mask <- decades == dec
      combined <- digit_mask & dec_mask

      # tab1: all but last decade for this digit
      if (dec < max(decades)) {
        tab1[as.character(d)] <- tab1[as.character(d)] + sum(wts[combined])
      }
      # tab2: all but first decade for this digit
      if (dec > 0) {
        tab2[as.character(d)] <- tab2[as.character(d)] + sum(wts[combined])
      }
    }
  }

  # Apply blending
  TAB <- tab1 * w + tab2 * (10 - w)

  fractions <- TAB / sum(TAB)

  # Myers' index
  out <- sum(abs(fractions - 0.1)) * 50

  return(out)
}


#' Bachi's Index of Age Heaping
#'
#' Calculate Bachi's index to measure digit preference in age data.
#'
#' @description
#' Bachi's index involves applying the Whipple method repeatedly to determine
#' the extent of preference for each terminal digit (0-9). It equals the sum
#' of positive deviations from 10 percent.
#'
#' @details
#' The theoretical range is 0 to 90:
#' \itemize{
#'   \item 0: no digit preference (each digit represents 10% of population)
#'   \item 90: maximum heaping (all ages end in same digit)
#' }
#'
#' For populations with no age heaping, each digit should appear in
#' approximately 10% of reported ages.
#'
#' @param x numeric vector of individual ages.
#' @param ageMin minimum age to include (default 23).
#' @param ageMax maximum age to include (default 77, adjusted to fit decades).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing Bachi's index.
#'
#' @author Matthias Templ
#'
#' @references
#' Bachi, R. (1951). The tendency to round off age returns: measurement and
#' correction. \emph{Bulletin of the International Statistical Institute},
#' \strong{33}(4), 195-222.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{myers}} for Myers' index,
#'   \code{\link{whipple}} for Whipple's index.
#'
#' @export
#'
#' @examples
#' # No heaping
#' set.seed(42)
#' age_uniform <- sample(23:77, 10000, replace = TRUE)
#' bachi(age_uniform)  # Should be close to 0
#'
#' # Strong heaping on 0 and 5
#' age_heaped <- sample(seq(25, 75, by = 5), 5000, replace = TRUE)
#' bachi(age_heaped)  # Should be high
#'
bachi <- function(x, ageMin = 23, ageMax = 77, weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Adjust ageMax to fit complete decades
  maxA <- max(x)
  if (ageMax > maxA) {
    ageMax <- ageMin + 4 + 10 * floor((maxA - ageMin - 4) / 10)
  }

  diff <- ageMax - ageMin
  age_interval <- diff - diff %% 10 - 1

  # Two overlapping sets of age ranges for blending
  lower_agesI <- ageMin + 0:4
  upper_agesI <- lower_agesI + age_interval
  lower_agesII <- lower_agesI + 1
  upper_agesII <- upper_agesI + 1

  if (any(upper_agesII > ageMax)) {
    upper_agesII <- upper_agesII - 10
    upper_agesI <- upper_agesI - 10
  }

  numeratorI <- numeric(10)
  numeratorII <- numeric(10)
  denominatorI <- numeric(10)
  denominatorII <- numeric(10)

  trailing_digits <- x %% 10

  for (dig in 0:9) {
    # Find which age range set to use for this digit
    ind <- (lower_agesI %% 5) == (dig %% 5)

    # Set I
    maskI <- x >= lower_agesI[ind] & x <= upper_agesI[ind]
    if (!is.null(weight)) {
      numeratorI[dig + 1] <- sum(weight[maskI & trailing_digits == dig])
      denominatorI[dig + 1] <- sum(weight[maskI])
    } else {
      numeratorI[dig + 1] <- sum(maskI & trailing_digits == dig)
      denominatorI[dig + 1] <- sum(maskI)
    }

    # Set II
    maskII <- x >= lower_agesII[ind] & x <= upper_agesII[ind]
    if (!is.null(weight)) {
      numeratorII[dig + 1] <- sum(weight[maskII & trailing_digits == dig])
      denominatorII[dig + 1] <- sum(weight[maskII])
    } else {
      numeratorII[dig + 1] <- sum(maskII & trailing_digits == dig)
      denominatorII[dig + 1] <- sum(maskII)
    }
  }

  # Average the two sets and compute deviations
  fractions <- ((numeratorI / denominatorI) + (numeratorII / denominatorII)) / 2
  fractions[is.nan(fractions)] <- 0

  out <- 100 * sum(abs(fractions - 0.1)) / 2

  return(out)
}


#' Noumbissi's Digit Heaping Index
#'
#' Calculate Noumbissi's index for a specific terminal digit.
#'
#' @description
#' Noumbissi's method improves on Whipple's method by extending its basic
#' principle to all ten digits. It compares the count of ages ending in a
#' specific digit to the count in 5-year age groups centered on that digit.
#'
#' @details
#' The index compares the number of persons reporting ages ending in a
#' specific digit to one-fifth of the population in the 5-year age groups
#' centered on those ages.
#'
#' Interpretation:
#' \itemize{
#'   \item 1.0: no preference for the digit
#'   \item >1.0: preference (attraction) to the digit
#'   \item <1.0: avoidance of the digit
#' }
#'
#' @param x numeric vector of individual ages.
#' @param digit integer (0-9) specifying which terminal digit to evaluate
#'   (default 0).
#' @param ageMin minimum age to include (default 20 + digit).
#' @param ageMax maximum age to include (default ageMin + 30).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing Noumbissi's index for the
#'   specified digit.
#'
#' @author Matthias Templ
#'
#' @references
#' Noumbissi, A. (1992). L'indice de Whipple modifie: une application aux
#' donnees du Cameroun, de la Suede et de la Belgique. \emph{Population},
#' \strong{47}(4), 1038-1041.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{spoorenberg}} for Total Modified Whipple index,
#'   \code{\link{whipple}} for original Whipple's index.
#'
#' @export
#'
#' @examples
#' # No heaping
#' set.seed(42)
#' age_uniform <- sample(20:70, 10000, replace = TRUE)
#' noumbissi(age_uniform, digit = 0)  # Should be close to 1
#' noumbissi(age_uniform, digit = 5)  # Should be close to 1
#'
#' # Heaping on digit 0
#' age_heap0 <- sample(seq(20, 70, by = 10), 5000, replace = TRUE)
#' noumbissi(age_heap0, digit = 0)  # Should be > 1
#'
noumbissi <- function(x, digit = 0, ageMin = 20 + digit, ageMax = ageMin + 30,
                      weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!digit %in% 0:9) {
    stop("'digit' must be an integer from 0 to 9.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Filter to valid age range
  ind <- x >= (ageMin - 2) & x <= (ageMax + 2)
  if (sum(ind) == 0) {
    warning("No ages in specified range.")
    return(NA)
  }

  ages <- x[ind]
  if (!is.null(weight)) {
    wts <- weight[ind]
  } else {
    wts <- rep(1, length(ages))
  }

  # Numerator: ages ending in specified digit within range
  numer_mask <- ages >= ageMin & ages <= ageMax & (ages %% 10) == digit

  # Denominator: 5-year groups centered on ages ending in digit
  # For each age ending in digit, include ages within ±2
  denom_mask <- rep(FALSE, length(ages))
  target_ages <- ages[numer_mask]
  for (ta in unique(target_ages)) {
    denom_mask <- denom_mask | (ages >= (ta - 2) & ages <= (ta + 2))
  }

  numerator <- sum(wts[numer_mask])
  denominator <- sum(wts[denom_mask])

  if (denominator == 0) {
    return(NA)
  }

  return(5 * numerator / denominator)
}


#' Spoorenberg's Total Modified Whipple Index
#'
#' Calculate the Total Modified Whipple Index (Wtot) proposed by Spoorenberg.
#'
#' @description
#' The Total Modified Whipple Index extends Noumbissi's approach by summing
#' the absolute deviations from 1 for all ten digits, providing an overall
#' measure of age heaping across all terminal digits.
#'
#' @details
#' The index is calculated as:
#' \deqn{W_{tot} = \sum_{i=0}^{9} |1 - W_i|}
#' where \eqn{W_i} is Noumbissi's index for digit \eqn{i}.
#'
#' Interpretation:
#' \itemize{
#'   \item 0: no heaping (perfect data)
#'   \item Higher values indicate more heaping
#'   \item Maximum theoretical value is 16 (if all ages end in one digit)
#' }
#'
#' @param x numeric vector of individual ages.
#' @param ageMin minimum age to include (default 20).
#' @param ageMax maximum age to include (default 64).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing the Total Modified Whipple Index.
#'
#' @author Matthias Templ
#'
#' @references
#' Spoorenberg, T. and Dutreuilh, C. (2007). Quality of age reporting:
#' extension and application of the modified Whipple's index.
#' \emph{Population}, \strong{62}(4), 729-741.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{noumbissi}} for single-digit index,
#'   \code{\link{whipple}} for original Whipple's index.
#'
#' @export
#'
#' @examples
#' # No heaping
#' set.seed(42)
#' age_uniform <- sample(20:64, 10000, replace = TRUE)
#' spoorenberg(age_uniform)  # Should be close to 0
#'
#' # Strong heaping on 0 and 5
#' age_heaped <- sample(seq(20, 60, by = 5), 5000, replace = TRUE)
#' spoorenberg(age_heaped)  # Should be high
#'
spoorenberg <- function(x, ageMin = 20, ageMax = 64, weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Calculate Noumbissi index for each digit
  Wi <- sapply(0:9, function(d) {
    noumbissi(x, digit = d,
              ageMin = ageMin + d,
              ageMax = min(ageMax, ageMin + d + 30),
              weight = weight)
  })

  Wi[is.nan(Wi) | is.na(Wi)] <- 1  # Treat missing as no heaping

  # Total Modified Whipple Index
  Wtot <- sum(abs(1 - Wi))

  return(Wtot)
}


#' Coale-Li Age Heaping Index
#'
#' Calculate the Coale-Li index for detecting age heaping at older ages.
#'
#' @description
#' The Coale-Li index was developed to detect age heaping in populations
#' with high proportions of elderly persons. It compares actual counts at
#' specific ages to smoothed reference values using moving averages.
#'
#' @details
#' The method applies double moving averages to create a smooth reference
#' distribution, then calculates the ratio of observed to expected counts
#' for ages ending in a specified digit.
#'
#' Interpretation:
#' \itemize{
#'   \item 1.0: no preference for the digit
#'   \item >1.0: attraction to the digit (heaping)
#'   \item <1.0: avoidance of the digit
#' }
#'
#' This index is particularly useful for evaluating data quality at older ages
#' (60+) where heaping on round numbers is common.
#'
#' @param x numeric vector of individual ages.
#' @param digit integer (0-9) specifying which terminal digit to evaluate
#'   (default 0).
#' @param ageMin minimum age to include (default 60).
#' @param ageMax maximum age to include (default max(x)).
#' @param terms number of terms for moving average smoothing (default 5).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing the Coale-Li index.
#'
#' @author Matthias Templ
#'
#' @references
#' Coale, A. J. and Li, S. (1991). The effect of age misreporting in China
#' on the calculation of mortality rates at very high ages.
#' \emph{Demography}, \strong{28}(2), 293-301.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{kannisto}} for Kannisto's index,
#'   \code{\link{jdanov}} for Jdanov's index.
#'
#' @export
#'
#' @examples
#' # Create age data with heaping at older ages
#' set.seed(42)
#' age <- c(sample(60:99, 5000, replace = TRUE),
#'          rep(seq(60, 90, by = 10), each = 200))  # Add heaping on 0s
#' coale_li(age, digit = 0)  # Should be > 1
#' coale_li(age, digit = 5)  # Should be closer to 1
#'
coale_li <- function(x, digit = 0, ageMin = 60, ageMax = max(x),
                     terms = 5, weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!digit %in% 0:9) {
    stop("'digit' must be an integer from 0 to 9.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Create age frequency table
  ages <- ageMin:ageMax
  if (!is.null(weight)) {
    tab <- sapply(ages, function(a) sum(weight[x == a]))
  } else {
    tab <- sapply(ages, function(a) sum(x == a))
  }
  names(tab) <- ages

  if (length(tab) < 2 * terms) {
    warning("Age range too narrow for smoothing.")
    return(NA)
  }

  # Double moving average for smoothing
  ma <- function(v, n) {
    stats::filter(v, rep(1/n, n), sides = 2)
  }
  reference <- ma(ma(tab, terms), terms)

  # Calculate ratios
  ratio <- tab / reference
  ratio[is.nan(ratio) | is.infinite(ratio)] <- NA

  # Average ratio for specified digit
  digit_ages <- ages[ages %% 10 == digit]
  digit_idx <- which(ages %in% digit_ages)

  if (length(digit_idx) == 0) {
    return(NA)
  }

  mean(ratio[digit_idx], na.rm = TRUE)
}


#' Jdanov's Old-Age Heaping Index
#'
#' Calculate Jdanov's index for detecting heaping at very old ages.
#'
#' @description
#' Jdanov's index is designed to detect age heaping at very old ages
#' (typically 95+), where data quality is often poorest. It applies
#' the Whipple principle to specific old-age values.
#'
#' @details
#' The index compares counts at specified old ages to the surrounding
#' 5-year age groups, similar to the standard Whipple approach but
#' focused on the oldest ages where heaping is most problematic.
#'
#' Interpretation:
#' \itemize{
#'   \item 100: no heaping
#'   \item >100: preference for the specified ages
#'   \item 500: maximum heaping (all ages at specified values)
#' }
#'
#' @param x numeric vector of individual ages.
#' @param Agei numeric vector of specific ages to evaluate (default
#'   c(95, 100, 105)).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing Jdanov's index.
#'
#' @author Matthias Templ
#'
#' @references
#' Jdanov, D. A., Scholz, R. D., and Shkolnikov, V. M. (2008).
#' Official population statistics and the Human Mortality Database
#' estimates of populations aged 80+ in Germany and nine other
#' European countries. \emph{Demographic Research}, \strong{19},
#' 1169-1196.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{kannisto}} for Kannisto's index,
#'   \code{\link{coale_li}} for Coale-Li index.
#'
#' @export
#'
#' @examples
#' # Create old-age data with heaping
#' set.seed(42)
#' age <- c(sample(90:110, 2000, replace = TRUE),
#'          rep(c(95, 100, 105), each = 100))  # Add heaping
#' jdanov(age)  # Should be > 100
#'
#' # No heaping
#' age_uniform <- sample(90:110, 2000, replace = TRUE)
#' jdanov(age_uniform)  # Should be close to 100
#'
jdanov <- function(x, Agei = c(95, 100, 105), weight = NULL) {


  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Numerator: ages at specified values
  numer_mask <- x %in% Agei

  # Denominator: 5-year groups centered on specified ages
  denom_mask <- rep(FALSE, length(x))
  for (a in Agei) {
    denom_mask <- denom_mask | (x >= (a - 2) & x <= (a + 2))
  }

  if (!is.null(weight)) {
    numerator <- sum(weight[numer_mask])
    denominator <- sum(weight[denom_mask])
  } else {
    numerator <- sum(numer_mask)
    denominator <- sum(denom_mask)
  }

  if (denominator == 0) {
    return(NA)
  }

  return(500 * numerator / denominator)
}


#' Kannisto's Age Heaping Index
#'
#' Calculate Kannisto's index for detecting heaping at a specific old age.
#'
#' @description
#' Kannisto's index compares the count at a specific age to a geometric
#' mean of surrounding ages, providing a measure of heaping that is
#' robust to exponentially declining populations at old ages.
#'
#' @details
#' Unlike other indices that use arithmetic means, Kannisto's index uses
#' geometric means of neighboring ages, which is more appropriate for
#' old-age populations where counts decline exponentially.
#'
#' The index is calculated as the ratio of the count at age \code{Agei}
#' to the geometric mean of counts at ages \code{Agei-2} through
#' \code{Agei+2}.
#'
#' Interpretation:
#' \itemize{
#'   \item 1.0: no heaping at the specified age
#'   \item >1.0: heaping (attraction to the age)
#'   \item <1.0: avoidance of the age
#' }
#'
#' @param x numeric vector of individual ages.
#' @param Agei single age value to evaluate (default 90).
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A single numeric value representing Kannisto's index.
#'
#' @author Matthias Templ
#'
#' @references
#' Kannisto, V. (1999). Assessing the information on age at death of
#' old persons in national vital statistics. \emph{Validation of
#' Exceptional Longevity, Odense Monographs on Population Aging},
#' \strong{6}, 235-249.
#'
#' @family heaping indices
#'
#' @seealso \code{\link{jdanov}} for Jdanov's index,
#'   \code{\link{coale_li}} for Coale-Li index.
#'
#' @export
#'
#' @examples
#' # Create old-age data with heaping at 90
#' set.seed(42)
#' age <- c(sample(85:95, 2000, replace = TRUE),
#'          rep(90, 200))  # Add heaping at 90
#' kannisto(age, Agei = 90)  # Should be > 1
#'
#' # No heaping
#' age_uniform <- sample(85:95, 2000, replace = TRUE)
#' kannisto(age_uniform, Agei = 90)  # Should be close to 1
#'
kannisto <- function(x, Agei = 90, weight = NULL) {

  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (length(Agei) != 1) {
    stop("'Agei' must be a single age value.")
  }

  if (!is.null(weight)) {
    if (!is.numeric(weight)) {
      stop("'weight' must be a numeric vector.")
    }
    if (length(weight) != length(x)) {
      stop("'weight' must have the same length as 'x'.")
    }
  }

  # Get counts for Agei and surrounding ages
  ages_needed <- (Agei - 2):(Agei + 2)

  if (!is.null(weight)) {
    counts <- sapply(ages_needed, function(a) sum(weight[x == a]))
  } else {
    counts <- sapply(ages_needed, function(a) sum(x == a))
  }

  # Check for zeros (can't compute geometric mean with zeros)
  if (any(counts == 0)) {
    # Use arithmetic mean as fallback
    denom <- mean(counts)
  } else {
    # Geometric mean
    denom <- exp(mean(log(counts)))
  }

  if (denom == 0) {
    return(NA)
  }

  # Numerator is count at Agei
  if (!is.null(weight)) {
    numer <- sum(weight[x == Agei])
  } else {
    numer <- sum(x == Agei)
  }

  return(numer / denom)
}


#' Calculate All Heaping Indices
#'
#' Convenience function to calculate multiple heaping indices at once.
#'
#' @description
#' This function calculates all available heaping indices for a given age
#' vector, providing a comprehensive assessment of data quality.
#'
#' @param x numeric vector of individual ages.
#' @param weight optional numeric vector of sampling weights.
#'
#' @return A named list with all heaping indices:
#'   \describe{
#'     \item{whipple_standard}{Standard Whipple index (100 = no heaping)}
#'     \item{whipple_modified}{Modified Whipple index (0 = no heaping)}
#'     \item{myers}{Myers' blended index (0 = no heaping)}
#'     \item{bachi}{Bachi's index (0 = no heaping)}
#'     \item{spoorenberg}{Total Modified Whipple index (0 = no heaping)}
#'     \item{noumbissi_0}{Noumbissi's index for digit 0 (1 = no heaping)}
#'     \item{noumbissi_5}{Noumbissi's index for digit 5 (1 = no heaping)}
#'   }
#'
#' @author Matthias Templ
#'
#' @family heaping indices
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' # Uniform ages (no heaping)
#' age_uniform <- sample(20:70, 10000, replace = TRUE)
#' heaping_indices(age_uniform)
#'
#' # Heaped ages
#' age_heaped <- sample(seq(20, 70, by = 5), 5000, replace = TRUE)
#' heaping_indices(age_heaped)
#'
heaping_indices <- function(x, weight = NULL) {

  list(
    whipple_standard = whipple(x, method = "standard", weight = weight),
    whipple_modified = whipple(x, method = "modified", weight = weight),
    myers = myers(x, weight = weight),
    bachi = bachi(x, weight = weight),
    spoorenberg = spoorenberg(x, weight = weight),
    noumbissi_0 = noumbissi(x, digit = 0, weight = weight),
    noumbissi_5 = noumbissi(x, digit = 5, weight = weight)
  )
}
