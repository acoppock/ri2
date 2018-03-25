
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
                         sims = 1000,
                         progress_bar = FALSE) {
  # setup

  model_1 <- as.formula(model_1)
  model_2 <- as.formula(model_2)

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
  } else {
    weights_vec <- NULL
  }

  coefs_obs_1 <-
    estimatr::tidy.lm_robust(lm_robust_fit(
      y = outcome_vec,
      X = design_matrix_1,
      weights = weights_vec,
      ci = FALSE,
      cluster = NULL,
      alpha = 0.05,
      se_type = "none",
      return_vcov = FALSE,
      try_cholesky = FALSE,
      has_int = TRUE
    ))

  jx <- intersect(c("term", "coefficient_name"),  names(coefs_obs_1))[1]
  beta_ix <- intersect(c("estimate", "coefficients"),  names(coefs_obs_1))[1]
  se_ix <- intersect(c("se", "std.error"),  names(coefs_obs_1))[1]


  coefs_obs_1 <- coefs_obs_1[coefs_obs_1[[jx]] %in% colnames(design_matrix_1), beta_ix, drop = TRUE]

  coefs_obs_2 <-
    estimatr::tidy.lm_robust(lm_robust_fit(
      y = outcome_vec,
      X = design_matrix_2,
      weights = weights_vec,
      ci = FALSE,
      cluster = NULL,
      alpha = 0.05,
      se_type = "none",
      return_vcov = FALSE,
      try_cholesky = FALSE,
      has_int = TRUE
    ))

  coefs_obs_2 <- coefs_obs_2[coefs_obs_2[[jx]] %in% colnames(design_matrix_2), beta_ix, drop = TRUE]

  ssr_1 <- sum((outcome_vec - design_matrix_1 %*% coefs_obs_1) ^ 2)
  ssr_2 <- sum((outcome_vec - design_matrix_2 %*% coefs_obs_2) ^ 2)

  f_obs <- (ssr_1 - ssr_2) / (ncol(design_matrix_2) - ncol(design_matrix_1)) /
    (ssr_2 / (length(outcome_vec) - ncol(design_matrix_2)))


  # Obtain Hypothesized POs -------------------------------------------------

  if (length(sharp_hypothesis) == 1) {
    sharp_hypothesis <-
      rep(sharp_hypothesis, length(unique(assignment_vec))-1)
  }

  pos_mat <- generate_pos(
    Y = outcome_vec,
    assignment_vec = assignment_vec,
    sharp_hypothesis = sharp_hypothesis
  )

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

    design_matrix_sim_1 <-
      model.matrix.default(model_1, data = data)
    design_matrix_sim_2 <-
      model.matrix.default(model_2, data = data)

    outcome_vec_sim <-
      switching_equation(pos_mat = pos_mat, assignment_vec = data[, assignment])

    if (IPW) {
      weights_vec <-
        1 / obtain_condition_probabilities(declaration, assignment = Z_sim)
    } else {
      weights_vec <- NULL
    }

    coefs_sim_1 <-
      estimatr::tidy.lm_robust(lm_robust_fit(
        y = outcome_vec_sim,
        X = design_matrix_sim_1,
        weights = weights_vec,
        ci = FALSE,
        cluster = NULL,
        alpha = 0.05,
        se_type = "none",
        return_vcov = FALSE,
        try_cholesky = FALSE,
        has_int = TRUE
      ))

    coefs_sim_1 <- coefs_sim_1[coefs_sim_1[[jx]] %in% colnames(design_matrix_sim_1), beta_ix, drop = TRUE]

    coefs_sim_2 <-
      estimatr::tidy.lm_robust(lm_robust_fit(
        y = outcome_vec_sim,
        X = design_matrix_sim_2,
        weights = weights_vec,
        ci = FALSE,
        cluster = NULL,
        alpha = 0.05,
        se_type = "none",
        return_vcov = FALSE,
        try_cholesky = FALSE,
        has_int = TRUE
      ))

    coefs_sim_2 <- coefs_sim_2[coefs_sim_2[[jx]] %in% colnames(design_matrix_sim_2), beta_ix , drop = TRUE]

    ssr_sim_1 <-
      sum((outcome_vec_sim - design_matrix_sim_1 %*% coefs_sim_1) ^ 2)
    ssr_sim_2 <-
      sum((outcome_vec_sim - design_matrix_sim_2 %*% coefs_sim_2) ^ 2)

    f_sim <- (ssr_sim_1 - ssr_sim_2) / (ncol(design_matrix_sim_2) - ncol(design_matrix_sim_1)) /
      (ssr_sim_2 / (length(outcome_vec_sim) - ncol(design_matrix_sim_2)))

    return(f_sim)
  }

  if (progress_bar) {
    null_distribution <- pbapply::pbapply(permutation_matrix, 2, ri_function)
  } else {
    null_distribution <- apply(permutation_matrix, 2, ri_function)
  }

  sims_df <-
    data.frame(
      est_sim = null_distribution,
      est_obs = f_obs,
      coefficient = "F-statistic"
    )

  return(structure(list(sims_df = sims_df),
    class = "ri"
  ))
}
