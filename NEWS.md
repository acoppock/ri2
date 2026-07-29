# ri2 0.5.0

* Fixed wrong null distributions when the assignment variable appears in more than one term of the formula, such as `Y ~ Z * X`. Only the assignment columns of the design matrix were updated for each permutation, so interaction columns kept their observed values and the design matrix disagreed with itself about what the assignment was. The resulting null distribution was too narrow and p-values were anti-conservative. Such formulas now rebuild the design matrix from the permuted assignment; the single-term case keeps the faster column-patching path.
* Added a `potential_outcomes` argument to `conduct_ri()` (#19). Naming one column of `data` per treatment condition states the sharp hypothesis as a full schedule of potential outcomes, which allows hypothesized effects that vary across units; `sharp_hypothesis` can only express a constant shift per arm. Build the columns however you like, holding each unit's observed outcome fixed for the condition it was assigned to. Supplied columns are checked against the observed outcomes, which catches columns given in the wrong order and schedules that do not preserve what was observed.
* Supplying a non-zero `sharp_hypothesis` to `conduct_ri()` alongside `test_function` without also naming `outcome` is now an error. `outcome` defaults to `NULL`, in which case the switching equation never ran and the hypothesis was silently ignored.
* `conduct_ri()` now errors when the assignment variable does not appear in the formula, or in either model of a nested-model comparison. Only the assignment column is permuted, so columns derived from it are never recomputed; previously such a call ran and returned a degenerate null distribution with a p-value of 1 regardless of the data.
* Supplying `sampling_weights` alongside `test_function` is now an error rather than being silently ignored. ri2 cannot weight an arbitrary scalar statistic, so weighting belongs inside the test function; the error message says how.
* Vignette gains a "Factorial designs" section (#27) showing that a factorial is a multi-arm trial over its cells, and how to test any factorial estimand by deriving the factors inside a `test_function`.
* Fixed multi-arm trial bug (#33): arms 2+ now each get a fresh conditional permutation matrix instead of reusing the first arm's matrix, giving correct p-values for all arms.
* Fixed silent outcome transformation dropping (#20): `log(Y) ~ Z` and other formula transformations now apply correctly via `model.frame()`.
* Fixed S3 class name conflict with `bit::ri` (#28): class is now `"ri2"`.
* Fixed histogram bar stacking (#6): use `fill` aesthetic instead of `alpha` in `plot.ri2()`.
* Fixed `plot()` failing to draw any bars when the number of simulations is not a multiple of 20 (#32, thanks @mreece13): the bin count is now floored to a whole number, as `geom_histogram()` requires.
* `plot()` no longer over-bins multi-arm results. The bin count came from the total number of simulated estimates, but the plot facets by term, so a four-arm trial drew three times as many bins per panel as a two-arm trial over the same number of simulations per panel. Bins are now chosen per panel.
* Replaced floating-point rounding hack with tolerance-based comparison in `summary.ri2()` and `plot.ri2()`.
* Implemented `sampling_weights` argument (#16): pass a column name containing sampling weights; multiplied with the auto-computed IPW weights when both are present.
* Removed `IPW_weights` argument: pre-specified IPW weights are incoherent for RI because weights must vary with each permuted assignment. IPW is always computed from the declaration via `obtain_condition_probabilities()`.
* Implemented clustered standard errors for `studentize = TRUE`: add a `clusters` argument; CR2 SEs are used automatically when clusters are supplied.
* Added `ri_ci()` for randomization inference confidence intervals by test inversion (#31).
* Added explicit NA checks with informative messages (#2).
* Moved `randomizr` and `estimatr` from `Depends` to `Imports`. To keep `library(ri2)` self-sufficient, `declare_ra()`, `conduct_ra()`, the `*_ra()` assignment functions, and the `obtain_*()` helpers are now re-exported from ri2, so existing scripts that relied on `randomizr` being attached continue to work.

# ri2 0.4.1

* Minor patch release fixing package documentation and link anchors.

# ri2 0.4.0

* disabled tests that check against the `ri` package because it is no longer on CRAN.

# ri2 0.2.0

* added support for user-supplied permutation matrix
* fixed a non-zero null bug
* Added a `NEWS.md` file to track changes to the package.
