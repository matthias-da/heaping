# JSSAM-2026-017 R1 -- honest recomputation of Table 2 (twofold B-O, age's
# unexplained component A/B) under the redesigned method + MI. Fixes E1 (numbers
# recomputed at full precision so the diffs reproduce) and E5 (the model arm's
# imputation model EXCLUDES the outcome ln.real.wage -- no leakage).

suppressMessages({ library(heaping)
                   library(oaxaca) })
data("chicago")
vars <- c("ln.real.wage","age","female","LTHS","some.college","college","advanced.degree","foreign.born")
chic <- na.omit(chicago[, vars])
bo_form <- ln.real.wage ~ age + female + LTHS + some.college + college +
             advanced.degree | foreign.born | LTHS + some.college + college + advanced.degree

introduce_heaping <- function(df, age_var, p, interval = 5) {
  df <- as.data.frame(df); n <- nrow(df); k <- sample(n, round(p * n))
  df[k, age_var] <- round(df[k, age_var] / interval) * interval; df
}
unexpl_age <- function(d) {
  r <- suppressMessages(oaxaca(bo_form, data = d, R = 2))
  v <- r$twofold$variables[[5]]["age", c("coef(unexplained A)", "coef(unexplained B)")]
  c(A = unname(v[1]), B = unname(v[2]))
}
pool <- function(imps) {
  M <- sapply(imps, function(a) { d <- heap; d$age <- a; unexpl_age(d) })   # 2 x m (A,B)
  c(A = mean(M["A", ]), B = mean(M["B", ]), A_sd = sd(M["A", ]), B_sd = sd(M["B", ]))
}

orig <- unexpl_age(chic)
set.seed(20260629); heap <- introduce_heaping(chic, "age", 0.27, 5)
hpd  <- unexpl_age(heap)

m <- 20
f_nowage <- as.formula("age ~ female + LTHS + some.college + college + advanced.degree + foreign.born")
f_wage   <- as.formula("age ~ ln.real.wage + female + LTHS + some.college + college + advanced.degree + foreign.born")
mi_s  <- correctHeapsMI(heap$age, m = m, heaps = "5year", method = "lnorm", seed = 1)
mi_mn <- correctHeapsMI(heap$age, m = m, heaps = "5year", model = f_nowage, dataModel = heap, model.engine = "ranger", seed = 1)
mi_mw <- correctHeapsMI(heap$age, m = m, heaps = "5year", model = f_wage,   dataModel = heap, model.engine = "ranger", seed = 1)
ps <- pool(mi_s$imputations); pmn <- pool(mi_mn$imputations); pmw <- pool(mi_mw$imputations)

A <- c(orig["A"], hpd["A"], ps["A"], pmn["A"], pmw["A"])
B <- c(orig["B"], hpd["B"], ps["B"], pmn["B"], pmw["B"])
lab <- c("original", "heaping", "after (simple+MI)", "after model+MI (no wage)", "after model+MI (w/ wage)")
cat(sprintf("Twofold B-O: age's unexplained coefficient (group.weight = -1). m = %d. n = %d.\n", m, nrow(chic)))
cat("Original (truth): A = 0.181, B = 0.104  [matches submitted Table 2]\n\n")
print(data.frame(
  data  = lab,
  coefA = round(A, 4),
  coefB = round(B, 4),
  diffA = sprintf("%.2f%%", 100 * abs(A - A[1]) / abs(A[1])),
  diffB = sprintf("%.2f%%", 100 * abs(B - B[1]) / abs(B[1]))
), row.names = FALSE, right = FALSE)
cat(sprintf("\nMI between-imputation SD (simple A=%.4f, model-no-wage A=%.4f)\n", ps["A_sd"], pmn["A_sd"]))
