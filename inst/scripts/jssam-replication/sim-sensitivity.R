# JSSAM-2026-017 R1 -- Workstream 4: sensitivity analysis (Reviewer 1.3).
# Are the results robust to the design choices R1 named -- the truncated
# distribution, the truncation width, the heaping interval, and the model
# strategy? We report distributional metrics (Whipple/MAPE/JSD) across the first
# three, and the age-covariate-link recovery across model strategies.

suppressMessages(library(heaping))

set.seed(2026)
true_age <- round(rlnorm(12000, 2.466869, 1.652772)); true_age <- true_age[true_age < 93]
z <- as.numeric(scale(true_age)) + rnorm(length(true_age), 0, 0.6)
maxage <- max(true_age)

distn <- function(x) { p <- tabulate(factor(x, levels = 0:maxage)); p / sum(p) }
mape  <- function(pa, pt) { i <- pt > 0; 100 * mean(abs(pa[i] - pt[i]) / pt[i]) }
jsd   <- function(p, q) { eps <- 1e-12; p <- (p+eps)/sum(p+eps); q <- (q+eps)/sum(q+eps)
  m <- (p+q)/2; kl <- function(a,b) sum(a*log(a/b)); 0.5*kl(p,m) + 0.5*kl(q,m) }
p_true     <- distn(true_age)
link_truth <- unname(coef(lm(z ~ true_age))[2])

relabel <- function(age, intensity, interval) {
  k <- sample(length(age), round(intensity * length(age)))
  age[k] <- round(age[k] / interval) * interval; age
}
dmet <- function(a) c(Whipple = round(unname(whipple(a)), 1),
                      MAPE = round(mape(distn(a), p_true), 1),
                      JSD = round(jsd(distn(a), p_true), 4))
link <- function(a) round(unname(coef(lm(z ~ a))[2]), 4)

set.seed(30); h5 <- relabel(true_age, 0.3, 5)
dm5 <- data.frame(age = h5, z = z)

cat(sprintf("References: true Whipple=%.1f (MAPE 0, JSD 0) | heaped 5yr/30%%: Whipple=%.1f MAPE=%.1f JSD=%.4f\n",
            whipple(true_age), unname(whipple(h5)), mape(distn(h5), p_true), jsd(distn(h5), p_true)))

cat("\nA) Distribution method (simple arm, 5-year, width 2, intensity 0.3):\n")
print(as.data.frame(t(sapply(c("lnorm","norm","unif","kernel"),
      function(mt) dmet(correctHeaps(h5, "5year", method = mt, width = 2, seed = 1))))))

cat("\nB) Truncation half-width (lnorm, simple arm, 5-year, intensity 0.3):\n")
B <- t(sapply(1:4, function(w) dmet(correctHeaps(h5, "5year", method = "lnorm", width = w, seed = 1))))
rownames(B) <- paste0("width=", 1:4); print(as.data.frame(B))

cat("\nD) Heaping interval (lnorm, simple arm, default width, intensity 0.3):\n")
set.seed(31); h10 <- relabel(true_age, 0.3, 10)
print(as.data.frame(rbind(`5-year`  = dmet(correctHeaps(h5,  "5year",  "lnorm", seed = 1)),
                          `10-year` = dmet(correctHeaps(h10, "10year", "lnorm", seed = 1)))))

cat(sprintf("\nC) Model strategy: age-covariate-link recovery (10-year, intensity 0.4; truth slope=%.4f):\n", link_truth))
set.seed(40); h10b <- relabel(true_age, 0.4, 10); dm10 <- data.frame(age = h10b, z = z)
Cc <- c(heaped = link(h10b),
        none   = link(correctHeaps(h10b, "10year", "lnorm", seed = 1)),
        lm     = link(correctHeaps(h10b, "10year", "lnorm", model = age ~ z, dataModel = dm10, model.engine = "lm",     seed = 1)),
        ranger = link(correctHeaps(h10b, "10year", "lnorm", model = age ~ z, dataModel = dm10, model.engine = "ranger", seed = 1)))
print(data.frame(arm = names(Cc), link_slope = unname(Cc),
                 bias_pct = round(100 * (unname(Cc) - link_truth) / link_truth, 1)),
      row.names = FALSE, right = FALSE)
