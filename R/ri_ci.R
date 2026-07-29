#' Compute a randomization inference confidence interval
#'
#' Inverts the RI test over a grid of sharp null hypotheses to find the set of
#' hypotheses that cannot be rejected at level \code{alpha}. The bounds of that
#' set form the confidence interval.
#'
#' The permutation matrix is generated once and reused across all grid points,
#' so the cost is roughly \code{n_grid} times the cost of a single
#' \code{conduct_ri} call (without the permutation matrix generation step).
#'
#' Currently only supported for two-arm trials (single-term formulas). For
#' multi-arm designs, call \code{ri_ci} separately for each pairwise comparison
#' by setting \code{condition1} and \code{condition2} in the formula.
#'
#' @param ... Arguments passed to \code{\link{conduct_ri}}.
#' @param alpha Significance level. Defaults to 0.05.
#' @param n_grid Number of candidate sharp hypotheses to evaluate. Defaults to
#'   40. Increase for a finer grid and more precise bounds.
#'
#' @return A data frame with columns \code{term}, \code{ci_lower},
#'   \code{ci_upper}, and \code{alpha}.
#' @export
#'
#' @examples
#' declaration <- randomizr::declare_ra(N = 40, m = 20)
#' Z <- randomizr::conduct_ra(declaration)
#' Y <- 0.5 * Z + rnorm(40)
#' dat <- data.frame(Y, Z)
#' # sims and n_grid are kept small here so the example runs quickly;
#' # use larger values in practice for a finer, less noisy interval.
#' ri_ci(Y ~ Z, declaration = declaration, assignment = "Z", data = dat,
#'       sims = 100, n_grid = 20)
ri_ci <- function(..., alpha = 0.05, n_grid = 40) {
  dots <- list(...)

  if (is.null(dots[["declaration"]])) {
    stop("ri_ci() requires a `declaration` argument.")
  }
  if (!is.null(dots[["model_1"]]) || !is.null(dots[["test_function"]])) {
    stop("ri_ci() only supports the formula interface (single-term ATE tests).")
  }

  # First run to get the observed estimate and the scale of the null distribution
  out_0 <- do.call(conduct_ri, c(dots, list(sharp_hypothesis = 0)))

  terms <- unique(out_0$sims_df$term)
  if (length(terms) > 1) {
    stop(
      "ri_ci() currently supports only two-arm trials (one term). ",
      "For multi-arm designs, conduct pairwise comparisons separately."
    )
  }

  beta_hat <- unique(out_0$sims_df$est_obs)
  perm_sd  <- stats::sd(out_0$sims_df$est_sim)
  if (perm_sd == 0) perm_sd <- abs(beta_hat) + 1

  # Generate permutation matrix once and reuse across the grid
  declaration <- dots[["declaration"]]
  sims_val    <- if (!is.null(dots[["sims"]])) dots[["sims"]] else 1000L
  pm <- randomizr::obtain_permutation_matrix(declaration, maximum_permutations = sims_val)

  tau_grid <- seq(
    beta_hat - 5 * perm_sd,
    beta_hat + 5 * perm_sd,
    length.out = n_grid
  )

  p_vals <- vapply(tau_grid, function(tau) {
    out_i <- do.call(conduct_ri, c(
      dots,
      list(sharp_hypothesis = tau, permutation_matrix = pm)
    ))
    tidy(out_i)$p.value[1]
  }, numeric(1))

  accept <- tau_grid[p_vals >= alpha]

  if (length(accept) == 0L) {
    warning(
      "No grid points fell in the acceptance region. ",
      "The true CI may be outside the search range. ",
      "Consider increasing n_grid or sims."
    )
    return(data.frame(
      term      = terms,
      ci_lower  = NA_real_,
      ci_upper  = NA_real_,
      alpha     = alpha
    ))
  }

  data.frame(
    term     = terms,
    ci_lower = min(accept),
    ci_upper = max(accept),
    alpha    = alpha
  )
}
