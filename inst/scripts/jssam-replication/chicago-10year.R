# JSSAM-2026-017 R2 -- Reviewer 3, Section 3.1: "What about 10year heaping or commonly
# used heaping level? The results from a real-life data set are not convincing ..."
#
# Chicago wage-gap application under five- AND ten-year heaping, at two shares of
# misreported records (0.27 = the paper's setting, 0.50 = severe). For each cell:
#   (a) twofold Blinder-Oaxaca: unexplained age coefficient A/B   (m = 20, as Table 2)
#   (b) wage regression: pooled age slope with 95% CI             (m = 50, as Table 3)
#   (c) Whipple index of the age variable: original / heaped / simple / model
# The (5-year, 0.27) cell reproduces Tables 2 and 3 of the R1 manuscript exactly
# (same seeds), which serves as the built-in check of this script.

suppressMessages({ library(heaping)
                   library(oaxaca) })
data("chicago")
vars <- c("ln.real.wage","age","female","LTHS","some.college","college","advanced.degree","foreign.born")
chic <- na.omit(chicago[, vars])
bo_form <- ln.real.wage ~ age + female + LTHS + some.college + college +
             advanced.degree | foreign.born | LTHS + some.college + college + advanced.degree
form    <- ln.real.wage ~ age + female + LTHS + some.college + college +
             advanced.degree + foreign.born
f_nowage <- age ~ female + LTHS + some.college + college + advanced.degree + foreign.born

introduce_heaping <- function(df, age_var, p, interval = 5) {
  df <- as.data.frame(df); n <- nrow(df); k <- sample(n, round(p * n))
  df[k, age_var] <- round(df[k, age_var] / interval) * interval; df
}
unexpl_age <- function(d) {
  r <- suppressMessages(oaxaca(bo_form, data = d, R = 2))
  v <- r$twofold$variables[[5]]["age", c("coef(unexplained A)", "coef(unexplained B)")]
  c(A = unname(v[1]), B = unname(v[2]))
}
rubin <- function(qs, us) {
  m <- length(qs); Qbar <- mean(qs); Ubar <- mean(us)
  B  <- if (m > 1) stats::var(qs) else 0
  Tt <- Ubar + (1 + 1/m) * B
  r  <- (1 + 1/m) * B / Ubar
  df <- if (r > 0) (m - 1) * (1 + 1/r)^2 else Inf
  tc <- stats::qt(0.975, df)
  c(est = Qbar, lo = Qbar - tc*sqrt(Tt), hi = Qbar + tc*sqrt(Tt), share = (1 + 1/m) * B / Tt)
}
beta_age <- function(d) { s <- unname(summary(lm(form, d))$coefficients["age", ]); c(q = s[1], u = s[2]^2) }
pct <- function(a, b) 100 * abs(a - b) / abs(b)

orig_bo <- unexpl_age(chic)
beta_t  <- unname(coef(lm(form, chic))["age"])
whip_t  <- whipple(chic$age)
cat(sprintf("n = %d; age range %d-%d\n", nrow(chic), min(chic$age), max(chic$age)))
cat(sprintf("R2 of age on the imputation covariates: %.3f without the wage, %.3f with the wage\n",
            summary(lm(f_nowage, chic))$r.squared,
            summary(lm(update(f_nowage, . ~ . + ln.real.wage), chic))$r.squared))
cat(sprintf("TRUTH: B-O unexplained age A = %.4f, B = %.4f; beta_age = %.4f; Whipple = %.1f\n",
            orig_bo["A"], orig_bo["B"], beta_t, whip_t))

cells <- expand.grid(interval = c(5, 10), share = c(0.27, 0.50))
for (j in seq_len(nrow(cells))) {
  iv <- cells$interval[j]; p <- cells$share[j]
  hp <- if (iv == 5) "5year" else "10year"
  set.seed(20260629); heap <- introduce_heaping(chic, "age", p, iv)
  cat(sprintf("\n================ %d-year heaping, share %.2f (records changed: %d) ================\n",
              iv, p, sum(heap$age != chic$age)))

  ## (a) Blinder-Oaxaca, m = 20
  s20 <- correctHeapsMI(heap$age, m = 20, heaps = hp, method = "lnorm", seed = 1)
  m20 <- correctHeapsMI(heap$age, m = 20, heaps = hp, model = f_nowage, dataModel = heap,
                        model.engine = "ranger", seed = 1)
  bo  <- function(imps) rowMeans(sapply(imps, function(a) { d <- heap; d$age <- a; unexpl_age(d) }))
  H <- unexpl_age(heap); S <- bo(s20$imputations); M <- bo(m20$imputations)
  cat("(a) twofold B-O, unexplained age coefficient, deviation from truth in % (m = 20)\n")
  print(data.frame(arm = c("heaped", "simple", "model"),
                   coefA = round(c(H["A"], S["A"], M["A"]), 4),
                   coefB = round(c(H["B"], S["B"], M["B"]), 4),
                   diffA = round(pct(c(H["A"], S["A"], M["A"]), orig_bo["A"]), 2),
                   diffB = round(pct(c(H["B"], S["B"], M["B"]), orig_bo["B"]), 2)),
        row.names = FALSE)

  ## (b) wage regression slope, m = 50
  s50 <- correctHeapsMI(heap$age, m = 50, heaps = hp, method = "lnorm", seed = 1)
  m50 <- correctHeapsMI(heap$age, m = 50, heaps = hp, model = f_nowage, dataModel = heap,
                        model.engine = "ranger", seed = 1)
  pool <- function(imps) { e <- sapply(imps, function(a) { d <- heap; d$age <- a; beta_age(d) })
                           rubin(e["q", ], e["u", ]) }
  hs <- unname(summary(lm(form, heap))$coefficients["age", ])
  tab <- rbind(heaped = c(est = hs[1], lo = hs[1] - 1.96*hs[2], hi = hs[1] + 1.96*hs[2], share = 0),
               simple = pool(s50$imputations), model = pool(m50$imputations))
  cat("(b) wage regression, age slope (m = 50, Rubin)\n")
  print(data.frame(arm = rownames(tab), beta = sprintf("%.5f", tab[, "est"]),
                   dev_pct = round(pct(tab[, "est"], beta_t), 2),
                   CI95 = sprintf("[%.4f, %.4f]", tab[, "lo"], tab[, "hi"]),
                   covers = ifelse(tab[, "lo"] <= beta_t & tab[, "hi"] >= beta_t, "yes", "NO"),
                   between_share_pct = round(100 * tab[, "share"], 1)), row.names = FALSE)

  ## (c) Whipple index of the age variable
  wm <- function(imps) mean(sapply(imps, whipple))
  cat(sprintf("(c) Whipple: original %.1f | heaped %.1f | simple %.1f | model %.1f  (mean over m = 50)\n",
              whip_t, whipple(heap$age), wm(s50$imputations), wm(m50$imputations)))
}
