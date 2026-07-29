conduct_ri_test_function <- function(test_function,
                                     assignment = "Z",
                                     outcome = "Y",
                                     declaration,
                                     sharp_hypothesis = 0,
                                     potential_outcomes = NULL,
                                     sampling_weights = NULL,
                                     permutation_matrix = NULL,
                                     data,
                                     sims = 1000,
                                     progress_bar = FALSE) {
  # ri2 cannot weight an arbitrary scalar statistic on the user's behalf, and
  # this argument was previously accepted and silently ignored. Weighting is the
  # test function's job: the sampling weights column is passed through untouched,
  # and inverse probability weights can be recomputed from the declaration, which
  # the function captures by closure.
  if (!is.null(sampling_weights)) {
    stop(
      "sampling_weights is not used with the test_function interface.\n",
      "Apply weights inside your test function instead. The '", sampling_weights,
      "' column is passed through unchanged, and inverse probability weights can ",
      "be recomputed for each permuted assignment with\n",
      "  1 / randomizr::obtain_condition_probabilities(declaration, assignment = data[[\"",
      assignment, "\"]])"
    )
  }

  test_stat_obs <- test_function(data)
  assignment_vec <- data[[assignment]]

  # Without `outcome`, the switching equation never runs and the outcome stays
  # at its observed values. That is correct only under the sharp null of no
  # effect, so any other hypothesis would otherwise be silently ignored.
  if (is.null(outcome) && (any(sharp_hypothesis != 0) || !is.null(potential_outcomes))) {
    stop(
      "A sharp hypothesis was supplied without naming the outcome variable.\n",
      "Set outcome = <name> so ri2 can impute each unit's outcome under the ",
      "permuted assignment. Without it the outcome stays at its observed values, ",
      "which is only correct under the sharp null of no effect."
    )
  }

  if (!is.null(outcome)) {
    if (length(sharp_hypothesis) == 1) {
      sharp_hypothesis <-
        rep(sharp_hypothesis, length(unique(assignment_vec))-1)
    }

    pos_mat <- resolve_pos(
      Y = data[[outcome]],
      assignment_vec = assignment_vec,
      sharp_hypothesis = sharp_hypothesis,
      potential_outcomes = potential_outcomes,
      data = data
    )
  }

  if (is.null(permutation_matrix)) {
    permutation_matrix <- obtain_permutation_matrix(declaration,
      maximum_permutations = sims
    )
  }

  ri_function <- function(Z_sim) {
    if (is.factor(assignment_vec)) {
      Z_sim <- factor(Z_sim, levels = levels(assignment_vec))
    }

    data[, assignment] <- Z_sim

    if (!is.null(outcome)) {
      data[, outcome] <-
        switching_equation(pos_mat = pos_mat, assignment_vec = data[[assignment]])
    }

    test_function(data)
  }

  if (progress_bar) {
    test_stat_sim <- pbapply::pbapply(permutation_matrix, 2, ri_function)
  } else {
    test_stat_sim <- apply(permutation_matrix, 2, ri_function)
  }

  sims_df <-
    data.frame(
      est_sim = test_stat_sim,
      est_obs = test_stat_obs,
      term = "Custom Test Statistic"
    )

  return(structure(list(sims_df = sims_df),
    class = "ri2"
  ))
}
