# JSSAM-2026-017 R2 -- Reviewer 1: computational cost for very large survey data sets.
# One run each of correctHeaps() at n = 1e5 and n = 1e6, marginal (log-normal) and
# model-based (quantile regression forest and linear model), five-year heaping at 30%.

suppressMessages(library(heaping))
cat(R.version.string, "|", Sys.info()[["machine"]], "| cores:", parallel::detectCores(), "\n")

for (n in c(1e5, 1e6)) {
  set.seed(2026)
  age <- round(rlnorm(1.15 * n, 2.466869, 1.652772)); age <- age[age < 93][seq_len(n)]
  z   <- as.numeric(scale(age)) + rnorm(n, 0, 0.6)
  k   <- sample(n, round(0.30 * n)); h <- age; h[k] <- round(h[k] / 5) * 5
  d   <- data.frame(age = h, z = z)
  t_s <- system.time(correctHeaps(h, "5year", method = "lnorm", seed = 1))[["elapsed"]]
  t_l <- system.time(correctHeaps(h, "5year", model = age ~ z, dataModel = d,
                                  model.engine = "lm", seed = 1))[["elapsed"]]
  cat(sprintf("n = %.0e: marginal %.1f s | linear model %.1f s\n", n, t_s, t_l)); flush.console()
  t_r <- system.time(correctHeaps(h, "5year", model = age ~ z, dataModel = d,
                                  model.engine = "ranger", seed = 1))[["elapsed"]]
  cat(sprintf("n = %.0e: quantile regression forest %.1f s\n", n, t_r)); flush.console()
}
