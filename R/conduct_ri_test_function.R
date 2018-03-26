conduct_ri_test_function <- function(test_function,
                                     assignment = "Z",
                                     outcome = "Y",
                                     declaration,
                                     sharp_hypothesis = 0,
                                     IPW_weights = NULL,
                                     sampling_weights = NULL,
                                     permutation_matrix = NULL,
                                     data,
                                     sims = 1000,
                                     progress_bar = FALSE) {
  test_stat_obs <- test_function(data)
  assignment_vec <- data[[assignment]]


  if (!is.null(outcome)) {
    if (length(sharp_hypothesis) == 1) {
      sharp_hypothesis <-
        rep(sharp_hypothesis, length(unique(assignment_vec))-1)
    }

    pos_mat <- generate_pos(
      Y = data[[outcome]],
      assignment_vec = assignment_vec,
      sharp_hypothesis = sharp_hypothesis
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

    if (!is.null(IPW_weights)) {
      data[, IPW_weights] <-
        1 / obtain_condition_probabilities(declaration, assignment = data[[assignment]])
    }

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
      coefficient = "Custom Test Statistic"
    )

  return(structure(list(sims_df = sims_df),
    class = "ri"
  ))
}
