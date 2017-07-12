conduct_ri_test_function <- function(test_function,
                                     assignment = "Z",
                                     outcome = "Y",
                                     declaration,
                                     sharp_hypothesis = 0,
                                     IPW_weights = NULL,
                                     sampling_weights = NULL,
                                     data,
                                     sims = 1000) {
  test_stat_obs <- test_function(data)

  if (!is.null(outcome)) {
    pos_mat <- generate_pos(
      Y = outcome_vec,
      assignment_vec = assignment_vec,
      sharp_hypothesis = sharp_hypothesis
    )
  }

  ri_function <- function() {
    data[, assignment] <- conduct_ra(declaration)

    if (!is.null(IPW_weights)) {
      data[, IPW_weights] <-
        1 / obtain_condition_probabilities(declaration, assignment = data[, assignment])
    }

    if (!is.null(outcome)) {
      data[, outcome] <-
        switching_equation(pos_mat = pos_mat, assignment_vec = data[, assignment])

    }

    test_function(data)
  }

  test_stat_sim <- pbapply::pbreplicate(sims, ri_function())

  sims_df <-
    data.frame(est_sim = test_stat_sim,
               est_obs = test_stat_obs,
               coefficient = "Custom Test Statistic")

  return(structure(list(sims_df = sims_df),
                   class = "ri"))
}
