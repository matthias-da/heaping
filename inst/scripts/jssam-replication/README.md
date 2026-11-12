# Replication scripts — "Correction of Heaping at the Individual Level" (JSSAM-2026-017)

These scripts reproduce every table and figure in the paper using the `heaping` package
(version 0.2.0 or later). The demographic simulation scripts (`sim-distribution.R`,
`sim-sensitivity.R`, `sim-illustrative.R`, `boundary-check.R`) use the one-sided rule for the heap at
age 0 (`boundary = "lower"`) and need version 0.2.1, which is on GitHub until its release on CRAN
(`remotes::install_github("matthias-da/heaping")`). Each script is self-contained and sets its own
random seed.

## Requirements

```r
install.packages("heaping")   # the methods (CRAN)
install.packages(c("oaxaca", "ranger", "EnvStats", "fitdistrplus",
                   "VIM", "ggplot2", "mgcv"))   # used by the scripts below
```

Run any script with `Rscript <name>.R` (or `source()` it in an R session).

## Scripts

| Script | Reproduces |
|---|---|
| `chicago-twofold.R` | **Table 2** — twofold Blinder–Oaxaca decomposition of the wage gap (original / heaped / simple / model corrections). |
| `chicago-mi.R` | **§3.1 MI table** — wage-regression age coefficient under multiple imputation, pooled by Rubin's rules with 95% confidence intervals. |
| `sim-mi-gonogo.R` | **§3.2 simulation MI table** — recovery of the age–covariate slope under coarse 10-year heaping, with Rubin-pooled CIs at two sample sizes and a fabrication check. |
| `sim-sensitivity.R` | **§3.2 sensitivity table** — robustness to the truncated distribution, the truncation half-width, the heaping interval, and the model engine. |
| `sim-distribution.R` | **Figures 1–3 (`whipplesim`, `jsdsim`, `mapesim`)** — distributional recovery (Whipple, MAPE, JSD) across heaping intensities for the heaped data, the two individual-level corrections, and aggregated GAM smoothing. Runs in parallel; `BOUNDARY=none Rscript sim-distribution.R` reruns the grid without the one-sided rule and writes only `out-figs-R2-none.csv`. |
| `sim-illustrative.R` | **Figures 4–5 (`illustration5year`, `illustration10year`)** — illustrative bar plots on the simulation's base population with the true distribution as outline: five-year heaping at intensities 0.1, 0.3, 0.5 (heaped / corrected / GAM-smoothed) and ten-year heaping at intensity 0.3 (heaped / corrected). |
| `chicago-10year.R` | **§3.1, ten-year heaping check** — the Chicago application under five- and ten-year heaping at shares 0.27 and 0.50: Blinder–Oaxaca coefficients, pooled age slope, and the Whipple index of the age variable before and after correction. |
| `chicago-10year-draws.R` | **§3.1, ten-year heaping over 20 heaping draws** — mean deviation of the Blinder–Oaxaca coefficients and the age slope, with the imputation model without and with the wage; these are the numbers quoted in §3.1. |
| `boundary-check.R` | **§3.2.3 and §3.2.6, residual divergence** — share of the Jensen–Shannon divergence that remains after correction and sits around age 0, without and with the one-sided rule, and the count at age 0. |
| `timing.R` | **Conclusion, computational cost** — run times of `correctHeaps()` at n = 100,000 and 1,000,000 for the marginal, linear-model and quantile-forest variants. |

The `chicago` dataset is the 2013 CPS-ORG Chicago sample shipped with the `oaxaca` package; the
demographic simulations generate synthetic age data from the parameters described in the paper.
