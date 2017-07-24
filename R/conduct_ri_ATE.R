conduct_ri_ATE <- function(formula,
                           assignment = "Z",
                           declaration,
                           sharp_hypothesis = 0,
                           studentize = FALSE,
                           IPW = TRUE,
                           IPW_weights = NULL,
                           sampling_weights = NULL,
                           permutation_matrix = NULL,
                           data,
                           sims = 1000) {
  # setup

  assignment_vec <- data[, assignment]
  design_matrix <- model.matrix.default(formula, data = data)
  outcome_vec <- data[, all.vars(formula[[2]])]
  condition_names <- sort(unique(assignment_vec))

  # Determine coefficient names

  if (is.numeric(assignment_vec)) {
    coefficient_names <- assignment
  } else{
    coefficient_names <- paste0(assignment, condition_names[-1])
  }

  if (length(sharp_hypothesis) != 1 &
      length(sharp_hypothesis) != length(coefficient_names)) {
    stop(
      "If you supply multiple sharp hypotheses, you must supply a number of sharp hypotheses equal to the number of treatment conditions minus 1."
    )
  }

  if (length(sharp_hypothesis) == 1) {
    sharp_hypothesis <-
      rep(sharp_hypothesis, length(coefficient_names))
  }

  pos_mat <- generate_pos(Y = outcome_vec,
                          assignment_vec = assignment_vec,
                          sharp_hypothesis = sharp_hypothesis)


  if(studentize){
    se_type <- "HC2"
  } else {
    se_type <- "none"
  }

  # The observed value ------------------------------------------------------

  if (IPW) {
    weights_vec <-
      1 / obtain_condition_probabilities(declaration, assignment = assignment_vec)
    design_matrix_w <- sqrt(weights_vec) * design_matrix
    outcome_vec_w <- sqrt(weights_vec) * outcome_vec
    fit_obs <-
      estimatr:::lm_robust_helper(y = outcome_vec_w, X = design_matrix_w, type = se_type)
  } else {
    fit_obs <-
      estimatr:::lm_robust_helper(y = outcome_vec, X = design_matrix, type = se_type)
  }

  if(se_type == "none"){
    coefs_obs <- fit_obs$beta_hat
  } else {
    coefs_obs <- fit_obs$beta_hat / sqrt(diag(fit_obs$Vcov_hat))
  }

  rownames(coefs_obs) <- colnames(design_matrix)
  coefs_obs <- as.list(coefs_obs[coefficient_names, ])


  # set up functions --------------------------------------------------------


  null_distributions <- vector("list",
                               length = length(condition_names) - 1)

  names(null_distributions) <- coefficient_names


  if (is.null(permutation_matrix) &
      sims >= obtain_num_permutations(declaration)) {
    permutation_matrix <- obtain_permutation_matrix(declaration,
                                                    maximum_permutations = sims)
  }


  for (i in 2:length(condition_names)) {
    if (is.null(permutation_matrix)) {
      permutation_matrix <-
        replicate(
          sims,
          conduct_conditional_ra(
            declaration,
            assignment_vec = assignment_vec,
            conditions = condition_names[c(1, i)]
          )
        )
    }

    ri_function <- function(Z_sim) {
      if (is.factor(assignment_vec)) {
        Z_sim <- factor(Z_sim, levels = levels(assignment_vec))
      }

      design_matrix[, coefficient_names] <-
        model.matrix.default(~ Z_sim)[,-1]

      if (sharp_hypothesis[i - 1] == 0) {
        outcome_vec_sim <- outcome_vec
      } else{
        outcome_vec_sim <-
          switching_equation(pos_mat = pos_mat, assignment_vec = Z_sim)
      }

      if (IPW) {
        weights_vec <-
          1 / obtain_condition_probabilities(declaration, assignment = Z_sim)
        design_matrix_w <- sqrt(weights_vec) * design_matrix
        outcome_vec_sim_w <- sqrt(weights_vec) * outcome_vec_sim
        fit_sim <-
        estimatr:::lm_robust_helper(y = outcome_vec_sim_w, X = design_matrix_w, type = se_type)
      } else {
        fit_sim <-
          estimatr:::lm_robust_helper(y = outcome_vec_sim, X = design_matrix, type = se_type)
      }

      if(se_type == "none"){
        coefs_sim <- fit_sim$beta_hat
      } else {
        coefs_sim <- fit_sim$beta_hat / sqrt(diag(fit_sim$Vcov_hat))
      }

      rownames(coefs_sim) <- colnames(design_matrix)
      coefs_sim[coefficient_names[i - 1], ]
    }

    null_distributions[[i - 1]] <-
      pbapply::pbapply(permutation_matrix, 2, ri_function)
  }


  sharp_hypothesis <- as.list(sharp_hypothesis)
  names(sharp_hypothesis) <- coefficient_names

  sims_df <-
    mapply(
      FUN = data.frame,
      est_sim = null_distributions,
      est_obs = coefs_obs,
      SIMPLIFY = FALSE
    ) %>%
    bind_rows(.id = "coefficient")

  if(studentize){
    sims_df$coefficient <- paste0(sims_df$coefficient, " (studentized)")
  }

  return(structure(list(sims_df = sims_df),
                   class = "ri"))
}
