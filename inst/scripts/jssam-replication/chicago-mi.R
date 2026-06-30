# JSSAM-2026-017 R1 -- Workstream 2b: Chicago wage-gap under the NEW method + MI
# -----------------------------------------------------------------------------
# The paper's headline: model-based correction stays within <2% of the true
# coefficients while naive correction deviates >10%. Here we re-run that claim
# under the REDESIGNED correctHeaps() (conditional-draw model arm) with proper
# multiple imputation (Rubin-pooled), and we settle the leakage question the
# old code raised: the old model arm imputed age USING ln.real.wage (the
# outcome). We report the model arm BOTH with and without the outcome.
#
# Estimand: the age coefficient in the wage regression
#   ln.real.wage ~ age + female + LTHS + some.college + college +
#                  advanced.degree + foreign.born
# Deviation = |beta_age(arm) - beta_age(true)| / |beta_age(true)|.
# Heaping: 27% to nearest 5 (the paper's setting). m = 50 imputations.

suppressMessages({ library(heaping)
                   library(oaxaca) })
data("chicago")

# heaping mechanism, matching inst/scripts/cmd-functions.R::introduce_heaping
introduce_heaping <- function(df, age_var, heap_ratio, interval = 5) {
  df <- as.data.frame(df); n <- nrow(df)
  idx <- sample(n, round(heap_ratio * n))
  df[idx, age_var] <- round(df[idx, age_var] / interval) * interval
  df
}

vars    <- c("ln.real.wage", "age", "female", "LTHS", "some.college",
             "college", "advanced.degree", "foreign.born")
chic    <- na.omit(chicago[, vars])
form    <- ln.real.wage ~ age + female + LTHS + some.college + college +
                          advanced.degree + foreign.born
beta_t  <- unname(coef(lm(form, chic))["age"])               # TRUE age coefficient

rubin <- function(qs, us) {
  m <- length(qs); Qbar <- mean(qs); Ubar <- mean(us)
  B  <- if (m > 1) stats::var(qs) else 0
  Tt <- Ubar + (1 + 1/m) * B
  r  <- (1 + 1/m) * B / Ubar
  df <- if (r > 0) (m - 1) * (1 + 1/r)^2 else Inf
  tc <- stats::qt(0.975, df)
  c(est = Qbar, se = sqrt(Tt), lo = Qbar - tc*sqrt(Tt), hi = Qbar + tc*sqrt(Tt),
    share_between = (1 + 1/m) * B / Tt)
}
beta_age <- function(d) { s <- unname(summary(lm(form, d))$coefficients["age", ]); c(q = s[1], u = s[2]^2) }
pool_mi  <- function(imps, base) {
  e <- sapply(imps, function(a) { d <- base; d$age <- a; beta_age(d) })
  rubin(e["q", ], e["u", ])
}

set.seed(20260629)
heap <- introduce_heaping(chic, "age", 0.27, interval = 5)

# single-estimate arms
ci_h <- { s <- unname(summary(lm(form, heap))$coefficients["age", ])
          c(est = s[1], se = s[2], lo = s[1]-1.96*s[2], hi = s[1]+1.96*s[2], share_between = 0) }

m <- 50
mi_simple   <- correctHeapsMI(heap$age, m = m, heaps = "5year", method = "lnorm", seed = 1)
mi_mod_wage <- correctHeapsMI(heap$age, m = m, heaps = "5year",
                 model = age ~ ln.real.wage + female + LTHS + some.college +
                               college + advanced.degree + foreign.born,
                 dataModel = heap, model.engine = "ranger", seed = 1)
mi_mod_no   <- correctHeapsMI(heap$age, m = m, heaps = "5year",
                 model = age ~ female + LTHS + some.college +
                               college + advanced.degree + foreign.born,
                 dataModel = heap, model.engine = "ranger", seed = 1)

tab <- rbind(
  heaped               = ci_h,
  `simple+MI`          = pool_mi(mi_simple$imputations,   heap),
  `model+MI (w/ wage)` = pool_mi(mi_mod_wage$imputations, heap),
  `model+MI (no wage)` = pool_mi(mi_mod_no$imputations,   heap))

cat(sprintf("True age coefficient (beta_age) = %.4f   [n = %d, 27%% heaped to nearest 5, m = %d]\n\n",
            beta_t, nrow(chic), m))
print(data.frame(
  arm           = rownames(tab),
  beta_age      = round(tab[, "est"], 4),
  dev_from_true = sprintf("%.1f%%", 100 * abs(tab[, "est"] - beta_t) / abs(beta_t)),
  CI95          = sprintf("[%.4f, %.4f]", tab[, "lo"], tab[, "hi"]),
  covers_true   = ifelse(tab[, "lo"] <= beta_t & tab[, "hi"] >= beta_t, "yes", "NO"),
  MI_unc_share  = sprintf("%.0f%%", 100 * tab[, "share_between"])
), row.names = FALSE, right = FALSE)
