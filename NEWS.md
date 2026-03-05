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
