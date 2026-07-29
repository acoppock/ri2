# Compute a randomization inference confidence interval

Inverts the RI test over a grid of sharp null hypotheses to find the set
of hypotheses that cannot be rejected at level `alpha`. The bounds of
that set form the confidence interval.

## Usage

``` r
ri_ci(..., alpha = 0.05, n_grid = 40)
```

## Arguments

- ...:

  Arguments passed to
  [`conduct_ri`](https://alexandercoppock.com/ri2/reference/conduct_ri.md).

- alpha:

  Significance level. Defaults to 0.05.

- n_grid:

  Number of candidate sharp hypotheses to evaluate. Defaults to 40.
  Increase for a finer grid and more precise bounds.

## Value

A data frame with columns `term`, `ci_lower`, `ci_upper`, and `alpha`.

## Details

The permutation matrix is generated once and reused across all grid
points, so the cost is roughly `n_grid` times the cost of a single
`conduct_ri` call (without the permutation matrix generation step).

Currently only supported for two-arm trials (single-term formulas). For
multi-arm designs, call `ri_ci` separately for each pairwise comparison
by setting `condition1` and `condition2` in the formula.

## Examples

``` r
declaration <- randomizr::declare_ra(N = 40, m = 20)
Z <- randomizr::conduct_ra(declaration)
Y <- 0.5 * Z + rnorm(40)
dat <- data.frame(Y, Z)
ri_ci(Y ~ Z, declaration = declaration, assignment = "Z", data = dat,
      sims = 200)
#>   term ci_lower ci_upper alpha
#> 1    Z 0.128207 1.226221  0.05
```
