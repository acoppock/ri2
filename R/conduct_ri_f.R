
conduct_ri_f <- function(model_1,
                         model_2,
                         assignment = "Z",
                         declaration,
                         sharp_hypothesis = 0,
                         IPW = TRUE,
                         IPW_weights = NULL,
                         sampling_weights = NULL,
                         permutation_matrix = NULL,
                         data = data,
                         sims = 1000) {
  # setup

  design_matrix_1 <- model.matrix.default(model_1, data = data)
  design_matrix_2 <- model.matrix.default(model_2, data = data)

  if (nrow(design_matrix_1) != nrow(design_matrix_2)) {
    stop("Missigness!")
  }
  if (all.vars(model_1[[2]]) != all.vars(model_2[[2]])) {
    stop("Outcome Match!")
  }

  assignment_vec <- data[, assignment]
  outcome_vec <- data[, all.vars(model_1[[2]])]

  # The observed value ------------------------------------------------------

  if (IPW) {
    weights_vec <-
      1 / obtain_condition_probabilities(declaration, assignment = assignment_vec)
    design_matrix_1 <- sqrt(weights_vec) * design_matrix_1
    design_matrix_2 <- sqrt(weights_vec) * design_matrix_2
    outcome_vec <- sqrt(weights_vec) * outcome_vec
  }

  coefs_obs_1 <-
    quick_lm(y = outcome_vec, X = design_matrix_1)$coefficients
  coefs_obs_2 <-
    quick_lm(y = outcome_vec, X = design_matrix_2)$coefficients

  ssr_1 <- sum((outcome_vec - design_matrix_1 %*% coefs_obs_1) ^ 2)
  ssr_2 <- sum((outcome_vec - design_matrix_2 %*% coefs_obs_2) ^ 2)

  f_obs = (ssr_1 - ssr_2) / (ncol(design_matrix_2) - ncol(design_matrix_1)) /
    (ssr_2 / (length(outcome_vec) - ncol(design_matrix_2)))


  # Obtain Hypothesized POs -------------------------------------------------

  pos_mat <- generate_pos(Y = outcome_vec,
                          assignment_vec = assignment_vec,
                          sharp_hypothesis = sharp_hypothesis)

  if (is.null(permutation_matrix)) {
    permutation_matrix <- obtain_permutation_matrix(declaration,
                                                    maximum_permutations = sims)
  }


  ri_function <- function(Z_sim) {
    if (is.factor(assignment_vec)) {
      Z_sim <- factor(Z_sim, levels = levels(assignment_vec))
    }

    data[, assignment] <- Z_sim

    design_matrix_sim_1 <-
      model.matrix.default(model_1, data = data)
    design_matrix_sim_2 <-
      model.matrix.default(model_2, data = data)

    outcome_vec_sim <-
      switching_equation(pos_mat = pos_mat, assignment_vec = data[, assignment])

    if (IPW) {
      weights_vec_sim <-
        1 / obtain_condition_probabilities(declaration, assignment = data[, assignment])
      design_matrix_sim_1 <-
        sqrt(weights_vec_sim) * design_matrix_sim_1
      design_matrix_sim_2 <-
        sqrt(weights_vec_sim) * design_matrix_sim_2
      outcome_vec_sim <- sqrt(weights_vec_sim) * outcome_vec_sim
    }

    coefs_sim_1 <-
      quick_lm(y = outcome_vec_sim, X = design_matrix_sim_1)$coefficients
    coefs_sim_2 <-
      quick_lm(y = outcome_vec_sim, X = design_matrix_sim_2)$coefficients

    ssr_sim_1 <-
      sum((outcome_vec_sim - design_matrix_sim_1 %*% coefs_sim_1) ^ 2)
    ssr_sim_2 <-
      sum((outcome_vec_sim - design_matrix_sim_2 %*% coefs_sim_2) ^ 2)

    f_sim = (ssr_sim_1 - ssr_sim_2) / (ncol(design_matrix_sim_2) - ncol(design_matrix_sim_1)) /
      (ssr_sim_2 / (length(outcome_vec_sim) - ncol(design_matrix_sim_2)))


    return(f_sim)
  }

  null_distribution <- pbapply::pbapply(permutation_matrix, 2, ri_function)

  sims_df <-
    data.frame(est_sim = null_distribution,
               est_obs = f_obs,
               coefficient = "F-statistic")

  return(structure(list(sims_df = sims_df),
                   class = "ri"))
}
