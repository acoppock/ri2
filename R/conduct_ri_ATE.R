#' @importFrom randomizr obtain_permutation_matrix obtain_num_permutations
#' @importFrom estimatr lm_fit
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

  if (studentize) {
    se_type <- "HC2"
  } else {
    se_type <- "none"
  }

  # The observed value ------------------------------------------------------

  if (IPW) {
    weights_vec <-
      1 / obtain_condition_probabilities(declaration, assignment = assignment_vec)
  } else {
    weights_vec <- NULL
  }

    fit_obs <- estimatr:::lm_fit(y = outcome_vec,
                                 design_matrix = design_matrix,
                                 weights = weights_vec,
                                 ci = FALSE,
                                 coefficient_name = coefficient_names,
                                 cluster = NULL,
                                 alpha = 0.05,
                                 se_type = se_type)

    if (studentize) {
      coefs_obs <- fit_obs$est / fit_obs$se
    } else {
      coefs_obs <- fit_obs$est
    }

  #rownames(coefs_obs) <- colnames(design_matrix[coefficient_names,])
  #coefs_obs <- as.list(coefs_obs[coefficient_names, ])

    names(coefs_obs) <- coefficient_names
    coefs_obs <- as.list(coefs_obs)

  # set up functions --------------------------------------------------------


  null_distributions <- vector("list", length = length(condition_names) - 1)

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
    } else {
      weights_vec <- NULL
    }

    fit_sim <- estimatr:::lm_fit(y = outcome_vec_sim,
                                 design_matrix = design_matrix,
                                 weights = weights_vec,
                                 ci = FALSE,
                                 coefficient_name = coefficient_names[i - 1],
                                 cluster = NULL,
                                 alpha = 0.05,
                                 se_type = se_type)

    if (studentize) {
      coefs_sim <- fit_sim$est / fit_sim$se
    } else {
      coefs_sim <- fit_sim$est
    }


    names(coefs_sim) <- coefficient_names[i - 1]
    return(coefs_sim)
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

  if (studentize) {
    sims_df$coefficient <- paste0(sims_df$coefficient, " (studentized)")
  }

  return(structure(list(sims_df = sims_df),
                   class = "ri"))
}
