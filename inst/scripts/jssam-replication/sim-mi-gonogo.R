# JSSAM-2026-017 R1 -- Workstream 2 go/no-go ---------------------------------
# Does the covariate-conditional (model) arm keep its edge once multiple
# imputation makes the correction uncertainty honest -- AND is naive (marginal)
# correction actually worse than leaving the data heaped (the spine, R2.21)?
#
# Estimand: the AGE-COVARIATE LINK, slope of lm(z ~ age), which is what marginal
# correction is supposed to break. Heaping is COARSE (10-year: round to nearest
# 10, error up to +/-5) so that correction genuinely moves ages. z is the signal
# covariate (correlated with true age); z_null is an independent decoy.
# Imputation uses z (+ z_null); it never sees the estimand's response role.
#
# Arms: oracle (true age) | heaped | simple+MI (marginal) | model+MI (age | z).
# MI arms pooled by Rubin's rules. Truth = population slope cov(z,age)/var(age).

suppressMessages(library(heaping))

# age = 45 + 9 z + N(0, 6^2)  =>  slope of z ~ age = cov/var = 9 / (81 + 36)
TRUTH <- 9 / (81 + 36)

gen <- function(n) {
  z      <- rnorm(n)
  z_null <- rnorm(n)
  age    <- round(45 + 9 * z + rnorm(n, 0, 6))   # corr(z, age) ~ 0.83
  age    <- pmin(pmax(age, 15), 92)
  data.frame(age = age, z = z, z_null = z_null)
}

relabel10 <- function(age, p = 0.50) {           # misreport to nearest 10 (coarse)
  k <- sample(length(age), round(p * length(age)))
  age[k] <- round(age[k] / 10) * 10
  age
}

est <- function(age, target) {                   # slope of target ~ age, + variance
  s <- summary(lm(target ~ age))$coefficients["age", ]
  c(q = unname(s[1]), u = unname(s[2])^2)
}

rubin <- function(qs, us) {
  m <- length(qs); Qbar <- mean(qs); Ubar <- mean(us)
  B  <- if (m > 1) stats::var(qs) else 0
  Tt <- Ubar + (1 + 1/m) * B
  r  <- (1 + 1/m) * B / Ubar
  df <- if (r > 0) (m - 1) * (1 + 1/r)^2 else Inf
  tc <- stats::qt(0.975, df)
  c(est = Qbar, se = sqrt(Tt), lo = Qbar - tc * sqrt(Tt), hi = Qbar + tc * sqrt(Tt),
    within = Ubar, between = B)
}

ci_single <- function(e) {
  q <- unname(e["q"]); u <- unname(e["u"])
  c(est = q, se = sqrt(u), lo = q - 1.96 * sqrt(u), hi = q + 1.96 * sqrt(u),
    within = u, between = 0)
}

one_run <- function(n = 8000, m = 50, seed = 1) {
  set.seed(seed)
  d <- gen(n); d$heaped <- relabel10(d$age)

  ci_o <- ci_single(est(d$age,    d$z))
  ci_h <- ci_single(est(d$heaped, d$z))

  mi_s <- correctHeapsMI(d$heaped, m = m, heaps = "10year", method = "lnorm", seed = seed)
  es   <- sapply(mi_s$imputations, function(a) est(a, d$z))
  ci_s <- rubin(es["q", ], es["u", ])

  dm   <- data.frame(age = d$heaped, z = d$z, z_null = d$z_null)
  mi_m <- correctHeapsMI(d$heaped, m = m, heaps = "10year",
                         model = age ~ z + z_null, dataModel = dm,
                         model.engine = "ranger", seed = seed)
  em   <- sapply(mi_m$imputations, function(a) est(a, d$z))
  ci_m <- rubin(em["q", ], em["u", ])

  corr <- c(oracle = cor(d$age, d$z), heaped = cor(d$heaped, d$z),
            `simple+MI` = mean(sapply(mi_s$imputations, function(a) cor(a, d$z))),
            `model+MI`  = mean(sapply(mi_m$imputations, function(a) cor(a, d$z))))

  list(truth = TRUTH,
       tab = rbind(oracle = ci_o, heaped = ci_h, `simple+MI` = ci_s, `model+MI` = ci_m),
       corr = corr,
       fab = c(true = cor(d$age, d$z_null),
               model = mean(sapply(mi_m$imputations, function(a) cor(a, d$z_null)))))
}

report <- function(n, m = 50) {
  r <- one_run(n = n, m = m); tab <- r$tab
  cov <- tab[, "lo"] <= r$truth & tab[, "hi"] >= r$truth
  share <- 100 * (1 + 1/m) * tab["model+MI", "between"] /
           (tab["model+MI", "within"] + (1 + 1/m) * tab["model+MI", "between"])
  cat(sprintf("\n=== n = %d, m = %d, 50%% heaped to nearest 10   (truth slope = %.4f) ===\n",
              n, m, r$truth))
  print(data.frame(
    arm        = rownames(tab),
    est        = round(tab[, "est"], 4),
    bias_pct   = round(100 * (tab[, "est"] - r$truth) / r$truth, 1),
    corr_age_z = round(r$corr[rownames(tab)], 3),
    CI95       = sprintf("[%.4f, %.4f]", tab[, "lo"], tab[, "hi"]),
    width      = round(tab[, "hi"] - tab[, "lo"], 4),
    covers     = ifelse(cov, "yes", "NO")
  ), row.names = FALSE, right = FALSE)
  cat(sprintf("model+MI correction uncertainty (between-imputation) share of total variance: %.1f%%\n",
              share))
  cat(sprintf("fabrication corr(age, z_null): true %.3f | model %.3f\n", r$fab["true"], r$fab["model"]))
}
report(8000)
report(400)
