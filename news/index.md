# Changelog

## ri2 0.5.1

- Fixed a silent failure with custom declarations.
  [`obtain_num_permutations()`](https://declaredesign.org/r/randomizr/reference/obtain_num_permutations.html)
  on a declaration made from a `permutation_matrix` is the width of that
  matrix, so a matrix with more columns than `sims` skipped the
  exact-test branch and fell through to the conditional re-randomization
  path, which has no method for that class. Every draw came back as the
  observed assignment, giving a null distribution with a single value
  and a p-value of 1 whatever the data. A custom declaration’s
  permutation matrix is the complete reference set, so it is now always
  used, however `sims` compares to its width.
- `conduct_conditional_ra()` now errors on a declaration class it does
  not handle. Its branches were a chain of `if` statements with no
  `else`, over a variable initialized to the observed assignment, so any
  unrecognized class returned the observed assignment for every draw
  rather than raising anything. That is the mechanism behind the bug
  above, and it would have silently repeated for any design type added
  to randomizr in future.

## ri2 0.5.0

CRAN release: 2026-07-29

- Fixed wrong null distributions when the assignment variable appears in
  more than one term of the formula, such as `Y ~ Z * X`. Only the
  assignment columns of the design matrix were updated for each
  permutation, so interaction columns kept their observed values and the
  design matrix disagreed with itself about what the assignment was. The
  resulting null distribution was too narrow and p-values were
  anti-conservative. Such formulas now rebuild the design matrix from
  the permuted assignment; the single-term case keeps the faster
  column-patching path.
- Added a `potential_outcomes` argument to
  [`conduct_ri()`](https://alexandercoppock.com/ri2/reference/conduct_ri.md)
  ([\#19](https://github.com/acoppock/ri2/issues/19)). Naming one column
  of `data` per treatment condition states the sharp hypothesis as a
  full schedule of potential outcomes, which allows hypothesized effects
  that vary across units; `sharp_hypothesis` can only express a constant
  shift per arm. Build the columns however you like, holding each unit’s
  observed outcome fixed for the condition it was assigned to. Supplied
  columns are checked against the observed outcomes, which catches
  columns given in the wrong order and schedules that do not preserve
  what was observed.
- Supplying a non-zero `sharp_hypothesis` to
  [`conduct_ri()`](https://alexandercoppock.com/ri2/reference/conduct_ri.md)
  alongside `test_function` without also naming `outcome` is now an
  error. `outcome` defaults to `NULL`, in which case the switching
  equation never ran and the hypothesis was silently ignored.
- [`conduct_ri()`](https://alexandercoppock.com/ri2/reference/conduct_ri.md)
  now errors when the assignment variable does not appear in the
  formula, or in either model of a nested-model comparison. Only the
  assignment column is permuted, so columns derived from it are never
  recomputed; previously such a call ran and returned a degenerate null
  distribution with a p-value of 1 regardless of the data.
- Supplying `sampling_weights` alongside `test_function` is now an error
  rather than being silently ignored. ri2 cannot weight an arbitrary
  scalar statistic, so weighting belongs inside the test function; the
  error message says how.
- Vignette gains a “Factorial designs” section
  ([\#27](https://github.com/acoppock/ri2/issues/27)) showing that a
  factorial is a multi-arm trial over its cells, and how to test any
  factorial estimand by deriving the factors inside a `test_function`.
- Fixed multi-arm trial bug
  ([\#33](https://github.com/acoppock/ri2/issues/33)): arms 2+ now each
  get a fresh conditional permutation matrix instead of reusing the
  first arm’s matrix, giving correct p-values for all arms.
- Fixed silent outcome transformation dropping
  ([\#20](https://github.com/acoppock/ri2/issues/20)): `log(Y) ~ Z` and
  other formula transformations now apply correctly via
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html).
- Fixed S3 class name conflict with `bit::ri`
  ([\#28](https://github.com/acoppock/ri2/issues/28)): class is now
  `"ri2"`.
- Fixed histogram bar stacking
  ([\#6](https://github.com/acoppock/ri2/issues/6)): use `fill`
  aesthetic instead of `alpha` in `plot.ri2()`.
- Fixed [`plot()`](https://rdrr.io/r/graphics/plot.default.html) failing
  to draw any bars when the number of simulations is not a multiple of
  20 ([\#32](https://github.com/acoppock/ri2/issues/32), thanks
  [@mreece13](https://github.com/mreece13)): the bin count is now
  floored to a whole number, as
  [`geom_histogram()`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
  requires.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) no longer
  over-bins multi-arm results. The bin count came from the total number
  of simulated estimates, but the plot facets by term, so a four-arm
  trial drew three times as many bins per panel as a two-arm trial over
  the same number of simulations per panel. Bins are now chosen per
  panel.
- Replaced floating-point rounding hack with tolerance-based comparison
  in `summary.ri2()` and `plot.ri2()`.
- Implemented `sampling_weights` argument
  ([\#16](https://github.com/acoppock/ri2/issues/16)): pass a column
  name containing sampling weights; multiplied with the auto-computed
  IPW weights when both are present.
- Removed `IPW_weights` argument: pre-specified IPW weights are
  incoherent for RI because weights must vary with each permuted
  assignment. IPW is always computed from the declaration via
  [`obtain_condition_probabilities()`](https://declaredesign.org/r/randomizr/reference/obtain_condition_probabilities.html).
- Implemented clustered standard errors for `studentize = TRUE`: add a
  `clusters` argument; CR2 SEs are used automatically when clusters are
  supplied.
- Added [`ri_ci()`](https://alexandercoppock.com/ri2/reference/ri_ci.md)
  for randomization inference confidence intervals by test inversion
  ([\#31](https://github.com/acoppock/ri2/issues/31)).
- Added explicit NA checks with informative messages
  ([\#2](https://github.com/acoppock/ri2/issues/2)).
- Moved `randomizr` and `estimatr` from `Depends` to `Imports`. To keep
  [`library(ri2)`](https://alexandercoppock.com/ri2/) self-sufficient,
  [`declare_ra()`](https://declaredesign.org/r/randomizr/reference/declare_ra.html),
  [`conduct_ra()`](https://declaredesign.org/r/randomizr/reference/conduct_ra.html),
  the `*_ra()` assignment functions, and the `obtain_*()` helpers are
  now re-exported from ri2, so existing scripts that relied on
  `randomizr` being attached continue to work.

## ri2 0.4.1

CRAN release: 2025-10-14

- Minor patch release fixing package documentation and link anchors.

## ri2 0.4.0

CRAN release: 2022-05-26

- disabled tests that check against the `ri` package because it is no
  longer on CRAN.

## ri2 0.2.0

CRAN release: 2020-12-07

- added support for user-supplied permutation matrix
- fixed a non-zero null bug
- Added a `NEWS.md` file to track changes to the package.
