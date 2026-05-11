# ri2 0.5.0

* Fixed multi-arm trial bug (#33): arms 2+ now each get a fresh conditional permutation matrix instead of reusing the first arm's matrix, giving correct p-values for all arms.
* Fixed silent outcome transformation dropping (#20): `log(Y) ~ Z` and other formula transformations now apply correctly via `model.frame()`.
* Fixed S3 class name conflict with `bit::ri` (#28): class is now `"ri2"`.
* Fixed histogram bar stacking (#6): use `fill` aesthetic instead of `alpha` in `plot.ri2()`.
* Replaced floating-point rounding hack with tolerance-based comparison in `summary.ri2()` and `plot.ri2()`.
* Implemented `sampling_weights` argument (#16): pass a column name containing sampling weights; multiplied with the auto-computed IPW weights when both are present.
* Removed `IPW_weights` argument: pre-specified IPW weights are incoherent for RI because weights must vary with each permuted assignment. IPW is always computed from the declaration via `obtain_condition_probabilities()`.
* Implemented clustered standard errors for `studentize = TRUE`: add a `clusters` argument; CR2 SEs are used automatically when clusters are supplied.
* Added `ri_ci()` for randomization inference confidence intervals by test inversion (#31).
* Added explicit NA checks with informative messages (#2).
* Moved `randomizr` and `estimatr` from `Depends` to `Imports`.
* Bumped `randomizr` requirement to >= 2.0.0.

# ri2 0.4.1

* Minor patch release fixing package documentation and link anchors.

# ri2 0.4.0

* disabled tests that check against the `ri` package because it is no longer on CRAN.

# ri2 0.2.0

* added support for user-supplied permutation matrix
* fixed a non-zero null bug
* Added a `NEWS.md` file to track changes to the package.
