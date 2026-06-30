# JSSAM-2026-017 R1 -- Workstream 3: distributional recovery + manuscript figures.
# Extends the earlier version: adds the GAM (aggregated-smoothing) comparison arm,
# a fine heaping-intensity grid with replication, the relabel-not-duplicate
# mechanism (finding A1), and writes the three manuscript figures
# (whipplesim/jsdsim/mapesim) with COLOUR-BLIND-SAFE styling (Okabe-Ito palette +
# line type + point shape, so the arms are distinguishable without colour).
#
# Metrics vs the TRUE age distribution: Whipple (100 = no heaping), MAPE of the
# age-frequency distribution, Jensen-Shannon divergence. Arms: heaped / simple /
# model / GAM-smoothed, across heaping intensities, averaged over replications.

suppressMessages({
  library(heaping)
  library(ggplot2); library(mgcv)
})

set.seed(2026)
true_age <- round(rlnorm(12000, 2.466869, 1.652772))
true_age <- true_age[true_age < 93]
z <- as.numeric(scale(true_age)) + rnorm(length(true_age), 0, 0.6)   # auxiliary covariate
maxage <- max(true_age); ages <- 0:maxage

relabel5 <- function(age, intensity) {
  k <- sample(length(age), round(intensity * length(age)))
  age[k] <- round(age[k] / 5) * 5; age
}
counts <- function(x) tabulate(factor(x, levels = ages))
mape   <- function(pa, pt) { i <- pt > 0; 100 * mean(abs(pa[i] - pt[i]) / pt[i]) }
jsd    <- function(p, q) { eps <- 1e-12; p <- (p+eps)/sum(p+eps); q <- (q+eps)/sum(q+eps)
  m <- (p+q)/2; kl <- function(a,b) sum(a*log(a/b)); 0.5*kl(p,m) + 0.5*kl(q,m) }
whip_c <- function(cnt) {            # Whipple from counts (ages 25..60 over 23..62)
  rng <- ages >= 23 & ages <= 62; mult <- ages %in% seq(25, 60, 5)
  sum(cnt[mult]) / (sum(cnt[rng]) / 5) * 100
}
gam_dist <- function(age) {          # aggregated GAM smoothing of the age histogram
  d <- data.frame(a = ages, n = counts(age))
  sm <- as.numeric(predict(mgcv::gam(n ~ s(a), family = poisson, data = d), type = "response"))
  sm[sm < 0] <- 0; sm
}
p_true <- counts(true_age) / sum(counts(true_age))

intensities <- seq(0, 0.5, by = 0.025)
nrep <- 25
grid <- expand.grid(intensity = intensities, rep = seq_len(nrep))

res <- do.call(rbind, Map(function(it, rp) {
  set.seed(as.integer(1000 * it + rp))
  h  <- relabel5(true_age, it)
  s  <- correctHeaps(h, "5year", method = "lnorm", seed = rp)
  mo <- correctHeaps(h, "5year", model = age ~ z,
                     dataModel = data.frame(age = h, z = z),
                     model.engine = "ranger", seed = rp)
  arms <- list(Heaped = counts(h), `Corrected (simple)` = counts(s),
               `Corrected (model)` = counts(mo), `Smoothed (GAM)` = gam_dist(h))
  do.call(rbind, lapply(names(arms), function(nm) {
    cnt <- arms[[nm]]; p <- cnt / sum(cnt)
    data.frame(intensity = it, rep = rp, arm = nm,
               Whipple = whip_c(cnt), MAPE = mape(p, p_true), JSD = jsd(p, p_true))
  }))
}, grid$intensity, grid$rep))

agg <- aggregate(cbind(Whipple, MAPE, JSD) ~ intensity + arm, res, mean)
lev <- c("Heaped", "Corrected (simple)", "Corrected (model)", "Smoothed (GAM)")
agg$arm <- factor(agg$arm, levels = lev)

# Okabe-Ito (colour-blind safe) + redundant line type + shape
ok <- c("#D55E00", "#0072B2", "#009E73", "#E69F00"); names(ok) <- lev
lt <- c("solid", "dashed", "dotdash", "dotted");      names(lt) <- lev
sh <- c(16, 17, 15, 4);                               names(sh) <- lev

mk <- function(yv, ylab, hline = NA) {
  g <- ggplot(agg, aes(intensity, .data[[yv]], colour = arm, linetype = arm, shape = arm))
  if (!is.na(hline)) g <- g + geom_hline(yintercept = hline, colour = "grey55", linewidth = 0.3)
  g + geom_line(linewidth = 0.8) + geom_point(size = 1.7) +
    scale_colour_manual(values = ok) + scale_linetype_manual(values = lt) +
    scale_shape_manual(values = sh) +
    labs(x = "Heaping intensity", y = ylab, colour = NULL, linetype = NULL, shape = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

ggsave("whipplesim.pdf", mk("Whipple", "Whipple index", 100), width = 6.5, height = 4.2)
ggsave("jsdsim.pdf",     mk("JSD", "Jensen-Shannon divergence"), width = 6.5, height = 4.2)
ggsave("mapesim.pdf",    mk("MAPE", "Mean absolute percentage error (%)"), width = 6.5, height = 4.2)

cat(sprintf("true-age Whipple = %.1f; n = %d; nrep = %d; %d intensities (0 to 0.5 by 0.025)\n",
            whip_c(counts(true_age)), length(true_age), nrep, length(intensities)))
print(agg[agg$intensity %in% c(0, 0.1, 0.2, 0.3, 0.4, 0.5), ], row.names = FALSE)
cat("\nFigures written: whipplesim.pdf jsdsim.pdf mapesim.pdf\n")
