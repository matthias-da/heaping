#' Multiple-Imputation Heaping Correction
#'
#' @description
#' Run \code{\link{correctHeaps}} \code{m} times with distinct, deterministically
#' derived seeds to produce \code{m} corrected datasets for multiple-imputation
#' inference. Because each corrected value is a random draw from a (truncated)
#' predictive distribution, a single correction understates the uncertainty it
#' introduces. Computing an estimate on each of the \code{m} datasets and pooling
#' with Rubin's rules propagates that correction uncertainty into downstream
#' standard errors and confidence intervals.
#'
#' @details
#' Child seeds are drawn once from the (optionally user-supplied) \code{seed}, so
#' the whole set of \code{m} imputations is reproducible while the individual
#' imputations differ from one another. Any further arguments in \code{\dots}
#' (e.g. \code{heaps}, \code{method}, \code{weight}, \code{model},
#' \code{dataModel}, \code{model.engine}) are passed unchanged to
#' \code{\link{correctHeaps}}.
#'
#' @param x numeric vector of ages (typically integers).
#' @param m integer number of imputations (default 50).
#' @param ... further arguments passed to \code{\link{correctHeaps}}.
#' @param seed optional integer seed. The \code{m} per-imputation seeds are
#'   derived deterministically from it, so the result is reproducible.
#'
#' @return An object of class \code{heapingMI}: a list with
#'   \describe{
#'     \item{imputations}{a list of \code{m} corrected numeric vectors}
#'     \item{m}{the number of imputations}
#'     \item{seeds}{the integer seeds used for the individual imputations}
#'   }
#'
#' @seealso \code{\link{correctHeaps}} for a single correction.
#'
#' @family heaping correction
#'
#' @author Matthias Templ
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' age <- round(rlnorm(2000, meanlog = 2.466869, sdlog = 1.652772))
#' age <- age[age < 93]
#' # introduce 5-year heaping by relabelling a share of records
#' k <- sample(length(age), round(0.3 * length(age)))
#' age5 <- age
#' age5[k] <- round(age5[k] / 5) * 5
#'
#' # Five corrected datasets for multiple-imputation inference
#' mi <- correctHeapsMI(age5, m = 5, heaps = "5year", seed = 42)
#' mi
#'
#' # Pool a simple estimate (the mean) across imputations
#' means <- sapply(mi$imputations, mean)
#' mean(means)        # point estimate
#' var(means)         # between-imputation variance (cost of correction)
correctHeapsMI <- function(x, m = 50, ..., seed = NULL) {
  if (!is.numeric(m) || length(m) != 1 || m < 1) {
    stop("'m' must be a positive integer.")
  }
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, m)
  imps <- lapply(seeds, function(sd_k) correctHeaps(x, ..., seed = sd_k))
  structure(list(imputations = imps, m = m, seeds = seeds), class = "heapingMI")
}

#' @rdname correctHeapsMI
#' @export
print.heapingMI <- function(x, ...) {
  cat("<heapingMI>", x$m, "imputed datasets of length",
      length(x$imputations[[1]]), "\n")
  invisible(x)
}
