# Internal engine for trusted-set, weight-aware heaping correction.
# Used by correctHeaps() and correctSingleHeap(). Not exported.
# (.sample_from_density lives in correctHeap.R and is reused by .draw_marginal.)

# Weighted (or unweighted) frequency over ages 0:maxage.
# @keywords internal
.wcount <- function(x, maxage, weight = NULL) {
  ages <- seq.int(0, maxage)
  out <- setNames(rep(0, length(ages)), as.character(ages))
  if (is.null(weight)) {
    tab <- table(x)
    out[names(tab)] <- as.numeric(tab)
  } else {
    agg <- tapply(weight, factor(x, levels = ages), sum)
    agg[is.na(agg)] <- 0
    out[] <- as.numeric(agg)
  }
  out
}

# Heaping ratio and expected (local-smooth) count at each heap age in s.
# Returns list(ratio, expected), both named numeric vectors over as.character(s).
# @keywords internal
.heap_ratios <- function(x, s, weight = NULL) {
  maxage <- max(x, na.rm = TRUE)
  wc <- .wcount(x, maxage, weight)
  ratio <- expected <- setNames(rep(NA_real_, length(s)), as.character(s))
  for (a in s) {
    bc <- as.character(a - 1); ac <- as.character(a + 1); hc <- as.character(a)
    if (!bc %in% names(wc) || !ac %in% names(wc)) next
    e <- mean(c(wc[bc], wc[ac]))
    expected[hc] <- e
    ratio[hc] <- if (e > 0) wc[hc] / e else NA_real_
  }
  list(ratio = ratio, expected = expected)
}

# Uniformly shuffle candidate indices, then take the smallest prefix whose
# cumulative weight reaches `excess`. Unit weights -> ceiling(excess) records,
# matching the legacy ssize = ceiling(n - n/ratio).
# @keywords internal
.select_to_correct <- function(idx, weight, excess) {
  n <- length(idx)
  if (n == 0 || is.na(excess) || excess <= 0) return(integer(0))
  ord <- sample.int(n)                      # uniform selection order
  w <- if (is.null(weight)) rep(1, n) else weight[ord]
  cw <- cumsum(w)
  k <- which(cw >= excess)[1]
  if (is.na(k)) k <- n                      # not enough mass: take all
  idx[ord[seq_len(k)]]
}

# Fit the marginal age distribution on a trusted subset of ages.
# @keywords internal
.fit_marginal <- function(trusted, method = "lnorm", sd = NULL) {
  fit <- list(method = method)
  if (method == "lnorm") {
    a0 <- as.numeric(trusted); a0[a0 == 0] <- 0.01
    fit$logn <- fitdistrplus::fitdist(a0, "lnorm")
  } else if (method == "norm") {
    fit$sd <- if (!is.null(sd)) sd else if (length(trusted) > 10)
      stats::mad(trusted, constant = 1.4826) else stats::sd(trusted)
  } else if (method == "kernel") {
    fit$density <- stats::density(trusted, from = 0, to = max(trusted), n = 512)
  }
  fit
}

# Draw n truncated replacements from a fitted marginal.
# @keywords internal
.draw_marginal <- function(n, fit, llow, lup) {
  if (n == 0) return(numeric(0))
  center <- (llow + lup) / 2
  if (fit$method == "lnorm") {
    round(EnvStats::rlnormTrunc(n, meanlog = fit$logn$estimate[1],
                                sdlog = as.numeric(fit$logn$estimate[2]),
                                min = llow, max = lup))
  } else if (fit$method == "norm") {
    round(EnvStats::rnormTrunc(n, mean = center, sd = fit$sd, min = llow, max = lup))
  } else if (fit$method == "kernel") {
    .sample_from_density(n, fit$density, llow, lup)
  } else {
    sample(llow:lup, n, replace = TRUE)
  }
}

# Fit age | covariates on the trusted set. engine: "ranger" (quantile forest) or "lm".
# weight_col, if non-NULL, is already a column of trusted_df (RRZ predictor).
# @keywords internal
.fit_cond <- function(trusted_df, formula, engine = "ranger",
                      weight_col = NULL, case.weights = NULL) {
  if (!is.null(weight_col)) formula <- stats::update(formula, paste(". ~ . +", weight_col))
  if (engine == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE))
      stop("Package 'ranger' is required for model.engine = 'ranger'.")
    m <- ranger::ranger(formula, data = trusted_df, quantreg = TRUE,
                        case.weights = case.weights)
  } else if (engine == "lm") {
    # Build the call with do.call: lm() evaluates `weights` non-standardly, so a bare
    # `weights = case.weights` fails to resolve inside a function (and errors when NULL).
    args <- list(formula = formula, data = trusted_df)
    if (!is.null(case.weights)) args$weights <- case.weights
    m <- do.call(stats::lm, args)
  } else stop("Unsupported model.engine: ", engine)
  attr(m, "engine") <- engine
  m
}

# Draw one truncated value per row of newdata from the conditional predictive.
# Falls back to the marginal fit when no fitted quantile lands in the row's window.
# @keywords internal
.draw_cond <- function(model, newdata, llow, lup, fit_marginal) {
  engine <- attr(model, "engine")
  n <- nrow(newdata); out <- numeric(n); nfb <- 0L
  if (engine == "ranger") {
    qs <- seq(0.005, 0.995, by = 0.01)
    pq <- predict(model, newdata, type = "quantiles", quantiles = qs)$predictions
    for (i in seq_len(n)) {
      cand <- pq[i, ]; cand <- cand[cand >= llow[i] & cand <= lup[i]]
      if (length(cand) >= 2) {
        out[i] <- round(sample(cand, 1))
      } else {                               # no fitted quantile in window -> marginal
        out[i] <- .draw_marginal(1, fit_marginal, llow[i], lup[i]); nfb <- nfb + 1L
      }
    }
  } else {                                   # lm
    mu <- as.numeric(predict(model, newdata))
    sdv <- stats::sd(stats::residuals(model))
    for (i in seq_len(n))
      out[i] <- round(EnvStats::rnormTrunc(1, mu[i], sdv, min = llow[i], max = lup[i]))
  }
  attr(out, "n_fallback") <- nfb
  out
}
