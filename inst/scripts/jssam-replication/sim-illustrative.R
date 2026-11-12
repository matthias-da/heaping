# JSSAM-2026-017 R2 -- Illustrative figures (Figures 4 and 5), Reviewer 3 on Section 3.3:
#   * state the heaping intensity shown, and show different intensities,
#   * visualise the data-generating process,
#   * show the GAM result so that its limitation is visible,
#   * regular axis labels.
#
# Base population = the simulation's base population (sim-distribution.R), so that one
# description of the data-generating process (Section 3.2.1) covers the simulation and
# the illustration. The TRUE (unheaped) distribution is drawn as a black outline in
# every panel. The correction uses the one-sided rule for the heap at age 0
# (boundary = "lower", boundary.bin = 0.5, as in sim-distribution.R; needs heaping >= 0.2.1).
#
# Figure 4 (illustration5year.pdf): five-year heaping at intensities 0.1, 0.3, 0.5;
#   columns: heaped | corrected (individual level, log-normal) | smoothed (GAM).
# Figure 5 (illustration10year.pdf): ten-year heaping at intensity 0.3;
#   columns: heaped | corrected (individual level, log-normal).
# Writes illustration5year.pdf and illustration10year.pdf to the working directory.

suppressMessages({
  library(heaping)
  library(ggplot2); library(mgcv)
})

## ---- base population: identical to sim-distribution.R ------------------------
set.seed(2026)
true_age <- round(rlnorm(12000, 2.466869, 1.652772))
true_age <- true_age[true_age < 93]
maxage <- max(true_age); ages <- 0:maxage
counts <- function(x) tabulate(factor(x, levels = ages))
relabel <- function(age, intensity, interval) {      # misreport to the nearest multiple
  k <- sample(length(age), round(intensity * length(age)))
  age[k] <- round(age[k] / interval) * interval; age
}
gam_dist <- function(age) {                          # as in sim-distribution.R
  d <- data.frame(a = ages, n = counts(age))
  sm <- as.numeric(predict(mgcv::gam(n ~ s(a), family = poisson, data = d), type = "response"))
  sm[sm < 0] <- 0; sm
}
whip_c <- function(cnt) {
  rng <- ages >= 23 & ages <= 62; mult <- ages %in% seq(25, 60, 5)
  sum(cnt[mult]) / (sum(cnt[rng]) / 5) * 100
}
truth <- data.frame(age = ages, count = counts(true_age))
cat(sprintf("base population: n = %d, max age = %d, true Whipple = %.1f, count at age 0 = %d\n",
            length(true_age), maxage, whip_c(truth$count), truth$count[1]))

## ---- panel data ---------------------------------------------------------------
build <- function(intensity, interval, heaps, gam = TRUE) {
  set.seed(as.integer(1000 * intensity + 1))         # as replication 1 of sim-distribution.R
  h <- relabel(true_age, intensity, interval)
  s <- correctHeaps(h, heaps, method = "lnorm", seed = 1,
                    boundary = "lower", boundary.bin = 0.5)
  stopifnot(length(h) == length(true_age), length(s) == length(true_age))
  arms <- list(Heaped = counts(h), `Corrected (individual level)` = counts(s))
  if (gam) arms$`Smoothed (GAM)` <- gam_dist(h)
  cat(sprintf("%2d-year, intensity %.1f | Whipple: %s | count at age 0: %s\n", interval, intensity,
              paste(sprintf("%s %.0f", names(arms), sapply(arms, whip_c)), collapse = ", "),
              paste(sprintf("%.0f", sapply(arms, `[`, 1)), collapse = " / ")))
  do.call(rbind, lapply(names(arms), function(nm)
    data.frame(age = ages, count = arms[[nm]], arm = nm,
               intensity = sprintf("Heaping intensity %.1f", intensity),
               heap = ifelse(ages %% interval == 0, "heap age", "other age"))))
}

plot_panels <- function(d, xbreaks) {
  d$arm <- factor(d$arm, levels = unique(d$arm))
  ggplot(d, aes(age, count)) +
    geom_col(aes(fill = heap), width = 0.9) +
    geom_step(data = truth, aes(age, count), direction = "mid", colour = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(`heap age` = "gold2", `other age` = "grey60"), guide = "none") +
    scale_x_continuous(breaks = xbreaks, expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Age", y = "Frequency") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92"))
}

## ---- Figure 4: five-year heaping, three intensities ----------------------------
d5 <- do.call(rbind, lapply(c(0.1, 0.3, 0.5), build, interval = 5, heaps = "5year"))
g5 <- plot_panels(d5, seq(0, 80, 20)) + facet_grid(intensity ~ arm)
ggsave("illustration5year.pdf", g5, width = 7.5, height = 6.6)

## ---- Figure 5: ten-year heaping, intensity 0.3 ----------------------------------
d10 <- build(0.3, interval = 10, heaps = "10year", gam = FALSE)
g10 <- plot_panels(d10, seq(0, 90, 10)) + facet_grid(. ~ arm)
ggsave("illustration10year.pdf", g10, width = 7.5, height = 3.3)

cat("Figures written: illustration5year.pdf illustration10year.pdf\n")
