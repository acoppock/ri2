## Release version 0.5.0 of ri2

This release fixes several correctness bugs, adds features, and moves two
packages from Depends to Imports.

Bug fixes affecting results:

* Multi-arm trials reused the first arm's conditional permutation matrix for
  every subsequent arm, so all comparisons beyond the first were tested against
  the wrong null distribution.
* Formulas in which the assignment variable appeared in more than one term
  (`Y ~ Z * X`) updated only the assignment columns of the design matrix for
  each permutation, leaving interaction columns at their observed values. The
  resulting null distribution was too narrow and p-values were anti-conservative.
* Outcome transformations in a formula (`log(Y) ~ Z`) were silently dropped and
  the untransformed outcome was used.
* `plot()` drew no bars when the number of simulations was not a multiple of 20,
  because a fractional bin count was passed to `geom_histogram()`.

Interface changes:

* `randomizr` and `estimatr` moved from Depends to Imports. The random
  assignment functions users need (`declare_ra()`, `conduct_ra()`, the `*_ra()`
  family, and the `obtain_*()` helpers) are re-exported, so existing scripts
  calling them after `library(ri2)` continue to work.
* The `IPW_weights` argument has been removed. It was documented but silently
  ignored, and a fixed column of inverse probability weights is not coherent
  for randomization inference, since the weights must be recomputed for each
  permuted assignment.
* Several arguments that were previously accepted and ignored now raise errors
  explaining what to do instead.

New features: a `clusters` argument for CR2 standard errors, a
`sampling_weights` argument, a `potential_outcomes` argument allowing sharp
hypotheses that vary across units, and `ri_ci()` for confidence intervals by
test inversion.

## Test environments

* local: aarch64-apple-darwin, R 4.6.0
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release,
  oldrel-1)

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs on any of the above.

## Downstream dependencies

There are currently no reverse dependencies for this package on CRAN
(checked against the current package database).
