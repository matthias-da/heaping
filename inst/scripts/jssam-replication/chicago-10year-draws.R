# JSSAM-2026-017 R2 -- robustness of the ten-year heaping check in the Chicago application.
# chicago-10year.R uses ONE heaping draw (the draw behind Tables 2 and 3). Here the heaping is
# redrawn 20 times, and the model-based correction is run both without and with the wage in
# the imputation model. Share of misreported records 0.27; five- and ten-year heaping.
# Per draw: twofold Blinder-Oaxaca unexplained age coefficient (A, B) and wage-regression age
# slope, each as the mean over m = 20 corrected data sets; Whipple index of the age variable.
suppressMessages({ library(heaping); library(oaxaca); library(parallel) })
data("chicago")
vars <- c("ln.real.wage","age","female","LTHS","some.college","college","advanced.degree","foreign.born")
chic <- na.omit(chicago[, vars])
bo_form <- ln.real.wage ~ age + female + LTHS + some.college + college +
             advanced.degree | foreign.born | LTHS + some.college + college + advanced.degree
form     <- ln.real.wage ~ age + female + LTHS + some.college + college + advanced.degree + foreign.born
f_nowage <- age ~ female + LTHS + some.college + college + advanced.degree + foreign.born
f_wage   <- update(f_nowage, . ~ . + ln.real.wage)
unexpl_age <- function(d) {
  r <- suppressMessages(oaxaca(bo_form, data = d, R = 2))
  v <- r$twofold$variables[[5]]["age", c("coef(unexplained A)", "coef(unexplained B)")]
  c(A = unname(v[1]), B = unname(v[2]))
}
slope <- function(d) unname(coef(lm(form, d))["age"])
truth <- c(unexpl_age(chic), slope = slope(chic)); m <- 20; ndraw <- 20

one_draw <- function(draw, iv) {
  set.seed(20261100 + draw)
  heap <- chic; k <- sample(nrow(chic), round(0.27 * nrow(chic)))
  heap$age[k] <- round(heap$age[k] / iv) * iv
  hp <- if (iv == 5) "5year" else "10year"
  est <- function(imps) rowMeans(sapply(imps, function(a) { d <- heap; d$age <- a
    c(unexpl_age(d), slope = slope(d), whipple = whipple(a)) }))
  arms <- rbind(
    heaped = c(unexpl_age(heap), slope = slope(heap), whipple = whipple(heap$age)),
    simple = est(correctHeapsMI(heap$age, m = m, heaps = hp, method = "lnorm", seed = draw)$imputations),
    `model (no wage)` = est(correctHeapsMI(heap$age, m = m, heaps = hp, model = f_nowage, dataModel = heap,
                                           model.engine = "ranger", seed = draw)$imputations),
    `model (with wage)` = est(correctHeapsMI(heap$age, m = m, heaps = hp, model = f_wage, dataModel = heap,
                                             model.engine = "ranger", seed = draw)$imputations))
  dev <- 100 * abs(sweep(arms[, c("A", "B", "slope")], 2, truth, "-")) / abs(matrix(truth, 4, 3, byrow = TRUE))
  data.frame(draw = draw, arm = rownames(arms), devA = dev[, "A"], devB = dev[, "B"],
             devSlope = dev[, "slope"], whipple = arms[, "whipple"], row.names = NULL)
}
cat(sprintf("n = %d; truth: A = %.4f, B = %.4f, slope = %.5f; Whipple original = %.1f\n",
            nrow(chic), truth["A"], truth["B"], truth["slope"], whipple(chic$age)))
cat(sprintf("%d heaping draws, share 0.27, m = %d imputations per arm and draw\n", ndraw, m))
for (iv in c(5, 10)) {
  res <- do.call(rbind, mclapply(seq_len(ndraw), one_draw, iv = iv, mc.cores = 8))
  res$arm <- factor(res$arm, levels = c("heaped", "simple", "model (no wage)", "model (with wage)"))
  agg <- do.call(rbind, lapply(split(res, res$arm), function(d) data.frame(
    arm = d$arm[1], devA = mean(d$devA), devB = mean(d$devB), devSlope = mean(d$devSlope),
    sdA = sd(d$devA), sdB = sd(d$devB), whipple = mean(d$whipple))))
  cat(sprintf("\n===== %d-year heaping: mean absolute deviation from the unheaped value in %% (SD over draws) =====\n", iv))
  print(data.frame(arm = agg$arm, devA = sprintf("%.2f (%.2f)", agg$devA, agg$sdA),
                   devB = sprintf("%.2f (%.2f)", agg$devB, agg$sdB),
                   devSlope = sprintf("%.2f", agg$devSlope), Whipple = sprintf("%.1f", agg$whipple)), row.names = FALSE)
  h <- res[res$arm == "heaped", ]
  for (a in levels(res$arm)[-1]) { d <- res[res$arm == a, ]
    cat(sprintf("draws in which '%s' is closer to the truth than the heaped data: A %d/%d, B %d/%d, slope %d/%d\n",
                a, sum(d$devA < h$devA), ndraw, sum(d$devB < h$devB), ndraw, sum(d$devSlope < h$devSlope), ndraw)) }
}
