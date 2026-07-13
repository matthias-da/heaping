# heaping 0.2.0 — submission comments

Update from 0.1.0. Major changes (details in NEWS.md):

* Redesigned model-based correction: `correctHeaps(model=, dataModel=)` now fits an
  imputation model on the retained ("trusted") records and draws covariate-conditional
  replacements; engine selectable via `model.engine = c("ranger", "lm")`. This is a
  documented breaking change relative to 0.1.0.
* New `correctHeapsMI()`: repeats the stochastic correction to produce multiple corrected
  data sets for multiple-imputation inference.
* Survey-weight support (`weight=`) for the correction functions and all heaping indices.
* Help-page examples now introduce artificial heaping by relabelling (misreporting)
  existing records rather than by duplication, matching the method's assumptions.

## Test environments

* Local: macOS (Apple Silicon), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 1 note

* "checking for future file timestamps ... unable to verify current time" — local
  time-server artifact; not reproducible on CRAN machines.

## Reverse dependencies

There are no reverse dependencies on CRAN.
