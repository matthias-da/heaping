# JSSAM-2026-017 R1 -- Illustrative before/after age-heaping figures (Figures 4-5).
#
# Regenerates age5heaps1/2.pdf and age10heaps1/2.pdf using the RELABELLING heaping
# mechanism (a proportion of existing records is misreported to the nearest round
# age), NOT duplication -- so length(age) is preserved. Left panel = heaped,
# right panel = after correctHeaps() (log-normal). Heap ages highlighted yellow.

suppressMessages(library(heaping))

set.seed(2026)
age <- round(rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772))
age <- age[age < 93]

# relabel: misreport round(intensity * n) records to the nearest multiple of
# 'interval'. length(age) is preserved (no duplication).
relabel <- function(age, intensity, interval) {
  k <- sample(length(age), round(intensity * length(age)))
  age[k] <- round(age[k] / interval) * interval
  age
}

intensity <- 0.35

# tabulate heaped and corrected over a common level set so the two panels share an
# identical x-axis; colour a bar yellow iff its integer age is a multiple of 'interval'.
make_panel <- function(x_heaped, x_corr, interval, file_heaped, file_corr) {
  maxa <- max(c(x_heaped, x_corr)); lev <- 0:maxa
  t1 <- table(factor(x_heaped, levels = lev))
  t2 <- table(factor(x_corr,   levels = lev))
  last <- max(which((t1 + t2) > 0))
  t1 <- t1[seq_len(last)]; t2 <- t2[seq_len(last)]
  cc <- ifelse(as.integer(names(t1)) %% interval == 0, "yellow", "darkgrey")
  pdf(file_heaped); barplot(t1, col = cc, xlab = "age", ylab = "frequency"); dev.off()
  pdf(file_corr);   barplot(t2, col = cc, xlab = "age", ylab = "frequency"); dev.off()
}

## 5-year heaping (Figure 4)
set.seed(101)
age5     <- relabel(age, intensity, 5)
age5corr <- correctHeaps(age5, heaps = "5year", method = "lnorm", seed = 1)
make_panel(age5, age5corr, 5, "age5heaps1.pdf", "age5heaps2.pdf")

## 10-year heaping (Figure 5)
set.seed(102)
age10     <- relabel(age, intensity, 10)
age10corr <- correctHeaps(age10, heaps = "10year", method = "lnorm", seed = 1)
make_panel(age10, age10corr, 10, "age10heaps1.pdf", "age10heaps2.pdf")

stopifnot(length(age5) == length(age), length(age10) == length(age),   # n preserved
          sum(age5 == 25)  > sum(age5corr == 25),    # *heaps1 = heaped (before), *heaps2 = corrected
          sum(age10 == 30) > sum(age10corr == 30))
cat("Figures written: age5heaps1.pdf age5heaps2.pdf age10heaps1.pdf age10heaps2.pdf\n")
