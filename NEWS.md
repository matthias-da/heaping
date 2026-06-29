# heaping 0.2.0

## Major changes

* `correctHeaps()` and `correctSingleHeap()` gain survey-weight support via a
  new `weight=` argument. Heaping ratios and the selection of records to correct
  then use weighted counts (records are drawn uniformly and accumulated until
  their cumulative weight covers the excess mass).

* **Breaking:** model-based correction is redesigned. Instead of a random-forest
  sign-adjustment heuristic, `correctHeaps(model=, dataModel=)` now fits an
  imputation model for age given the covariates on the *retained* ("trusted")
  records and draws covariate-conditional replacements for the selected heaped
  records. The engine is selectable via `model.engine = c("ranger", "lm")`
  (a `ranger` quantile forest by default). Corrected values therefore differ
  from 0.1.x, even with the same seed.

* New `width=` argument controls the truncation half-window; 10-year heaps now
  use a single symmetric window instead of the previous two-stage +/-4 / +/-5
  correction.

* New `correctHeapsMI()` produces `m` corrected datasets (with deterministically
  derived seeds) for multiple-imputation inference, with a `print()` method for
  the returned `heapingMI` object.

## Internal

* Removed the internal helpers `.adjust_signs()` and the legacy
  `.draw_replacements()`; the correction engine now lives in `R/impute-model.R`.

# heaping 0.1.1

## Bug fixes

* Fixed cascading drift in `correctHeaps()` when using custom heap positions.
  When heaps were specified at consecutive integers (e.g.,
  `heaps = seq(2, max(x), by = 1)`), observations corrected for one heap could
  be picked up and re-corrected at subsequent heaps, causing values to drift
  far from their original position (reported by Saskia Schirmer).

* Fixed R's `sample()` single-value trap in both `correctHeaps()` and
  `correctSingleHeap()`. When only one observation was available for correction
  at a heap, `sample(n, size = 1)` would sample from `1:n` instead of returning
  `n`, potentially writing replacement values to wrong indices.

* Added a warning when more than 50% of unique values in the data are declared
  as heaps, indicating likely misspecification of the `heaps` argument. Heaping
  correction is designed for sparse heap positions (e.g., multiples of 5 or 10),
  not for every value in the data.

# heaping 0.1.0

* Initial release.
* `correctHeaps()` and `correctSingleHeap()` for individual-level heaping
  correction using truncated log-normal, normal, uniform, or kernel density
  distributions.
* Heaping indices: `whipple()`, `myers()`, `bachi()`, `noumbissi()`,
  `spoorenberg()`, `coale_li()`, `jdanov()`, `kannisto()`, and
  `heaping_indices()`.
* `sprague()` for disaggregating 5-year age groups using Sprague multipliers.
* Support for sampling weights in all heaping indices.
* Optional model-based correction using random forest predictions.
* Vignette with comprehensive examples.
