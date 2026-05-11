# ri2 logfile

## 2026-05-11 — 0.5.0 overhaul

### Goal

Fix all 11 open GitHub issues, clean up the codebase, and add multi-arm trial documentation to the vignette.

### Bugs fixed

**#33 — Multi-arm permutation matrix reuse (correctness bug)**
In `conduct_ri_ATE`, the loop over arms generated a conditional RA permutation matrix for arm i=2 (T1/T2) and then cached it. On i=3+, `is.null(permutation_matrix)` was FALSE, so the T1/T2 matrix was reused for the T1/T3 comparison — giving wrong null distributions and wrong p-values for all arms beyond the first. Fix: restructure the loop to generate a fresh conditional RA matrix per arm (or reuse the user-supplied / exact-enumeration matrix, which is valid for all arms).

**#28 — S3 class "ri" conflicts with bit::ri**
When `readr` >= 2.0.1 is loaded, the `bit` package (a dependency) defines its own `ri` class. Printing/summarizing ri2 results dispatched to bit's method and errored. Fix: rename class from `"ri"` to `"ri2"` throughout all three handler files and all S3 methods.

**#20 — Outcome transformation silently dropped**
`log(Y) ~ Z` used `all.vars(formula[[2]])` to extract the outcome name ("Y") and then `data[["Y"]]` — ignoring the `log()`. Fix: extract outcome via `model.frame()` + `model.response()`, so any transformation is applied correctly.

**#6 — Histogram bars stacking in plot.ri**
`aes(alpha = extreme)` caused ggplot to stack histogram bars by alpha level, distorting the plot. Fix: use `aes(fill = extreme)` + `scale_fill_manual()` instead.

### Removed

**`IPW_weights` argument removed**
The parameter was documented but silently ignored in the original code. On reflection, a pre-specified column of IPW weights is incoherent for RI: weights are `1/P(Z_sim_i = z_sim_i)` and must vary with each permuted assignment. Any fixed column is specific to the observed assignment. Removed entirely; `IPW = TRUE` always computes weights fresh from the declaration via `obtain_condition_probabilities()`.

### New features

**`sampling_weights`** — pass a column name containing sampling weights. These are fixed across permutations (reflect sampling design, not assignment) and are multiplied with auto-computed IPW weights when both are present.

**`clusters`** — pass a cluster variable column name. When `studentize = TRUE`, CR2 clustered SEs are used instead of HC2. Previously documented as "CLUSTERING NOT YET IMPLEMENTED."

**`ri_ci()`** — randomization inference confidence interval by test inversion. Generates the permutation matrix once, evaluates the test over a grid of `n_grid` sharp hypotheses, and returns the interval where p ≥ alpha. Single-term (two-arm) formulas only.

**NA checks** — explicit `anyNA()` checks on outcome and assignment variables at the top of each handler, with informative error messages. Closes #2.

### Infrastructure

- `randomizr` and `estimatr` moved from `Depends` to `Imports`; `randomizr` requirement bumped to >= 2.0.0
- Floating-point rounding hack (`round(est_sim, 10)`) replaced with adaptive tolerance in `ri_extreme()`: `abs(est_obs - H0) * sqrt(.Machine$double.eps)` as tolerance
- 21 new tests in `test-fixes.R` covering all fixed bugs and new features
- Version bumped to 0.5.0

### Vignette additions

Added a "Multi-arm trials" section covering:
- The formula interface for 3+ arm designs (one result per pairwise arm vs. baseline)
- Vector `sharp_hypothesis` for arm-specific nulls
- **Why conditioning matters**: a worked comparison of naive full permutation vs. conditional RA for the T2 vs T1 test statistic
- **Leveling treatment Monte Carlo**: demonstrates that when T3 has low-variance outcomes (a leveling treatment), naive full permutation is anti-conservative for T2 (~10% false rejection rate vs 5% nominal) AND reduces power for T3 (~21% vs ~29% correct). Explains the mechanism: low-variance T3 units narrow the naive T2 null; high-variance T1/T2 units inflate the naive T3 null.

### Key technical insight documented

The direction of naive's distortion depends on T3 outcome variance relative to T1/T2:
- T3 high-variance (e.g., large constant effect): naive null is WIDER → conservative (under-rejects T2, misses T3 effects)
- T3 low-variance (leveling treatment): naive null is NARROWER → anti-conservative (false-rejects T2, still misses T3 effects)
Both are wrong. Conditional permutation is always correct.

### Resume this session

Session: /Users/alexandercoppock/.claude/projects/-Users-alexandercoppock-git-projects-metaprep/0e4e8439-5cbc-4e27-9ad8-17c3a95c0efd.jsonl
