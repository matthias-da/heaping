# JSSAM-2026-017 R2 -- the heap at age 0 (follow-up to Figures 4-5).
# (a) How much of the residual Jensen-Shannon divergence of the corrected data comes from the
#     boundary heap at age 0 when it is left uncorrected (package default, boundary = "none")?
# (b) What does the one-sided rule (boundary = "lower", boundary.bin = 0.5) change?
# Same base population and heaping mechanism as sim-distribution.R; 10 replications.
suppressMessages({ library(heaping); library(mgcv) })
set.seed(2026)
true_age <- round(rlnorm(12000, 2.466869, 1.652772)); true_age <- true_age[true_age < 93]
ages <- 0:max(true_age); counts <- function(x) tabulate(factor(x, levels = ages))
relabel <- function(age, it, iv) { k <- sample(length(age), round(it * length(age))); age[k] <- round(age[k] / iv) * iv; age }
gam_dist <- function(age) { d <- data.frame(a = ages, n = counts(age))
  sm <- as.numeric(predict(mgcv::gam(n ~ s(a), family = poisson, data = d), type = "response")); sm[sm < 0] <- 0; sm }
jsd_terms <- function(p, q) { eps <- 1e-12; p <- (p+eps)/sum(p+eps); q <- (q+eps)/sum(q+eps); m <- (p+q)/2
  0.5 * p * log(p/m) + 0.5 * q * log(q/m) }                      # sums to the JSD
mape <- function(cnt, ct) { p <- cnt / sum(cnt); pt <- ct / sum(ct); i <- pt > 0; 100 * mean(abs(p[i] - pt[i]) / pt[i]) }
p_true <- counts(true_age)
cat(sprintf("share of records at ages 0-2: %.1f%%; at ages 0-5: %.1f%%; true count at age 0: %d\n",
            100*mean(true_age <= 2), 100*mean(true_age <= 5), p_true[1]))
for (cfg in list(c(5, .1), c(5, .3), c(5, .5), c(10, .3))) {
  iv <- cfg[1]; it <- cfg[2]; low <- if (iv == 5) 0:2 else 0:5; hp <- if (iv == 5) "5year" else "10year"
  r <- sapply(1:10, function(rp) {
    set.seed(as.integer(1000 * it + rp)); h <- relabel(true_age, it, iv)
    s <- correctHeaps(h, hp, method = "lnorm", seed = rp)
    b <- correctHeaps(h, hp, method = "lnorm", seed = rp, boundary = "lower", boundary.bin = 0.5)
    ts <- jsd_terms(counts(s), p_true); tb <- jsd_terms(counts(b), p_true)
    tg <- jsd_terms(gam_dist(h), p_true); th <- jsd_terms(counts(h), p_true)
    c(heaped = sum(th), corr = sum(ts), corr_low = sum(ts[ages %in% low]), corr_rest = sum(ts[!ages %in% low]),
      gam = sum(tg), gam_rest = sum(tg[!ages %in% low]),
      bnd = sum(tb), bnd_low = sum(tb[ages %in% low]), bnd_rest = sum(tb[!ages %in% low]),
      n0_heaped = sum(h == 0), n0_corr = sum(s == 0), n0_bnd = sum(b == 0), min_bnd = min(b),
      mape_corr = mape(counts(s), p_true), mape_bnd = mape(counts(b), p_true))
  })
  m <- rowMeans(r)
  cat(sprintf("%2d-year, intensity %.1f | JSD heaped %.4f | corrected %.4f, of which ages %d-%d: %.0f%% (rest %.4f) | GAM %.4f (rest %.4f)\n",
              iv, it, m["heaped"], m["corr"], min(low), max(low), 100*m["corr_low"]/m["corr"], m["corr_rest"], m["gam"], m["gam_rest"]))
  cat(sprintf("   one-sided rule       | corrected %.4f, of which ages %d-%d: %.0f%% (rest %.4f) | count at age 0: heaped %.0f, default %.0f, one-sided %.0f | MAPE default %.1f, one-sided %.1f | minimum age %.0f\n",
              m["bnd"], min(low), max(low), 100*m["bnd_low"]/m["bnd"], m["bnd_rest"], m["n0_heaped"], m["n0_corr"], m["n0_bnd"],
              m["mape_corr"], m["mape_bnd"], min(r["min_bnd", ])))
}
