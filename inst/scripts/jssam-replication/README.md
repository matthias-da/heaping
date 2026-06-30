# Replication scripts — "Correction of Heaping at the Individual Level" (JSSAM-2026-017)

These scripts reproduce every table and figure in the paper using the `heaping` package
(version 0.2.0 or later). Each is self-contained and sets its own random seed.

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
| `sim-distribution.R` | **Figures (`whipplesim`, `jsdsim`, `mapesim`)** — distributional recovery (Whipple, MAPE, JSD) across heaping intensities for the heaped data, the two individual-level corrections, and aggregated GAM smoothing. |

The `chicago` dataset is the 2013 CPS-ORG Chicago sample shipped with the `oaxaca` package; the
demographic simulations generate synthetic age data from the parameters described in the paper.
