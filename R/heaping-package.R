#' heaping: Correction of Heaping on Individual Level
#'
#' Provides methods for correcting heaping (digit preference) in survey data
#' at the individual record level. Age heaping, where respondents
#' disproportionately report ages ending in 0 or 5, is a common phenomenon that
#' can distort demographic analyses.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{correctHeaps}}}{Correct regular age heaping patterns
#'     (5-year or 10-year intervals)}
#'   \item{\code{\link{correctSingleHeap}}}{Correct a specific single age heap}
#'   \item{\code{\link{correctHeapsMI}}}{Produce multiple corrected datasets for
#'     multiple-imputation inference}
#' }
#'
#' @section Methodology:
#' Unlike traditional smoothing methods that only correct aggregated statistics,
#' this package corrects individual values by replacing a calculated proportion
#' of heaped observations with draws from fitted truncated distributions
#' (log-normal, normal, or uniform).
#'
#' The correction ratio is determined by comparing the (optionally weighted)
#' count at each heap to the mean of neighboring ages. Records exceeding this
#' expected count are selected and treated as untrusted; the replacement
#' distribution is then fitted to the remaining (trusted) records and used to
#' draw truncated replacements. Survey weights are supported via the
#' \code{weight} argument.
#'
#' @section Model-Based Correction:
#' When a \code{model} formula and \code{dataModel} are supplied, an imputation
#' model for age given the covariates is fitted on the trusted records and used
#' to draw covariate-conditional replacements for the heaped records, preserving
#' relationships with other variables. The engine is selectable via
#' \code{model.engine}: a quantile regression forest (\pkg{ranger}, the default)
#' or a linear model. \pkg{VIM} is used only to impute missing covariates.
#'
#' @section Multiple Imputation:
#' \code{\link{correctHeapsMI}} repeats the correction \code{m} times with
#' distinct seeds to produce \code{m} corrected datasets. Pooling an estimate
#' across them with Rubin's rules propagates the correction uncertainty into
#' downstream standard errors and confidence intervals.
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
#' @author Matthias Templ \email{matthias.templ@@fhnw.ch}
#'
#' @docType package
#' @name heaping-package
#' @aliases heaping
#'
#' @importFrom stats setNames rnorm plogis predict
#' @importFrom fitdistrplus fitdist
#' @importFrom EnvStats rlnormTrunc rnormTrunc
"_PACKAGE"
