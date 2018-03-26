#' @importFrom randomizr obtain_permutation_matrix obtain_num_permutations
#' @importFrom estimatr lm_robust_fit
#' @importFrom stats model.matrix.default as.formula
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
                           sims = 1000,
                           progress_bar = FALSE) {
  # setup

  formula <- as.formula(formula)
  assignment_vec <- data[[assignment]]
  design_matrix <- model.matrix.default(formula, data = data)
  outcome_vec <- data[[all.vars(formula[[2]])]]
  condition_names <- sort(unique(assignment_vec))

  # Determine coefficient names

  if (is.numeric(assignment_vec)) {
    coefficient_names <- assignment
  } else {
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

  pos_mat <- generate_pos(
    Y = outcome_vec,
    assignment_vec = assignment_vec,
    sharp_hypothesis = sharp_hypothesis
  )

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

  fit_obs <- lm_robust_fit(
    y = outcome_vec,
    X = design_matrix,
    weights = weights_vec,
    ci = FALSE,
    cluster = NULL,
    alpha = 0.05,
    se_type = se_type,
    return_vcov = FALSE,
    try_cholesky = FALSE,
    has_int = TRUE
  )



  fit_obs <- estimatr::tidy.lm_robust(fit_obs)

  jx <- intersect(c("term", "coefficient_name"),  names(fit_obs))[1]
  beta_ix <- intersect(c("estimate", "coefficients"),  names(fit_obs))[1]
  se_ix <- intersect(c("se", "std.error"),  names(fit_obs))[1]


  fit_obs <- fit_obs[fit_obs[[jx]] %in% coefficient_names, , drop = FALSE]

  if (studentize) {
    coefs_obs <- fit_obs[[beta_ix]] / fit_obs[[se_ix]]
  } else {
    coefs_obs <- fit_obs[[beta_ix]]
  }

  names(coefs_obs) <- coefficient_names
  coefs_obs <- as.list(coefs_obs)

  # set up functions --------------------------------------------------------

  null_distributions <-
    vector("list", length = length(condition_names) - 1)

  names(null_distributions) <- coefficient_names


  if (is.null(permutation_matrix) &
    sims >= obtain_num_permutations(declaration)) {
    permutation_matrix <- obtain_permutation_matrix(declaration,
      maximum_permutations = sims
    )
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
        model.matrix.default(~ Z_sim)[, -1]

      if (sharp_hypothesis[i - 1] == 0) {
        outcome_vec_sim <- outcome_vec
      } else {
        outcome_vec_sim <-
          switching_equation(pos_mat = pos_mat, assignment_vec = Z_sim)
      }

      if (IPW) {
        weights_vec <-
          1 / obtain_condition_probabilities(declaration, assignment = Z_sim)
      } else {
        weights_vec <- NULL
      }

      fit_sim <- lm_robust_fit(
        y = outcome_vec_sim,
        X = design_matrix,
        weights = weights_vec,
        ci = FALSE,
        cluster = NULL,
        alpha = 0.05,
        se_type = se_type,
        return_vcov = FALSE,
        try_cholesky = FALSE,
        has_int = TRUE
      )

      fit_sim <- estimatr::tidy.lm_robust(fit_sim)
      fit_sim <- fit_sim[fit_sim[[jx]] %in% coefficient_names[i - 1], , drop = FALSE]

      if (studentize) {
        coefs_sim <- fit_sim[[beta_ix]] / fit_sim[[se_ix]]
      } else {
        coefs_sim <- fit_sim[[beta_ix]]
      }

      names(coefs_sim) <- coefficient_names[i - 1]
      return(coefs_sim)
    }

    if (progress_bar) {
      null_distributions[[i - 1]] <-
        pbapply::pbapply(permutation_matrix, 2, ri_function)
    } else {
      null_distributions[[i - 1]] <-
        apply(permutation_matrix, 2, ri_function)
    }
  }

  sharp_hypothesis <- as.list(sharp_hypothesis)
  names(sharp_hypothesis) <- coefficient_names

  sims_list <-
    mapply(
      FUN = data.frame,
      est_sim = null_distributions,
      est_obs = coefs_obs,
      SIMPLIFY = FALSE
    )

  sims_df <- do.call("rbind", sims_list)
  sims_df$coefficient <- rep(names(sims_list), sapply(sims_list, nrow))

  if (studentize) {
    sims_df$coefficient <- paste0(sims_df$coefficient, " (studentized)")
  }

  return(structure(list(sims_df = sims_df),
    class = "ri"
  ))
}
