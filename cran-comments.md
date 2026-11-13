# heaping 0.2.1 — submission comments

Minor update from 0.2.0 (details in NEWS.md):

* `correctHeaps()` gains `boundary = c("none", "lower")` and `boundary.bin`. With
  `boundary = "lower"`, a heap at the smallest observed age (typically age 0), which has no
  neighbour below it and was left uncorrected so far, is corrected by a one-sided rule.
  `correctHeapsMI()` passes both arguments through.
* The default `boundary = "none"` reproduces the results of 0.2.0, so existing code is
  unaffected.
* Fix: the uniform method now handles a window that contains a single age.

## Test environments

* Local: macOS (Apple Silicon), R 4.5.2, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
