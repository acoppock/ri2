#' @importFrom randomizr obtain_permutation_matrix
#' @importFrom estimatr lm_robust_fit
#' @importFrom stats model.matrix.default as.formula model.frame model.response
conduct_ri_f <- function(model_1,
                         model_2,
                         assignment = "Z",
                         declaration,
                         sharp_hypothesis = 0,
                         IPW = TRUE,
                         sampling_weights = NULL,
                         permutation_matrix = NULL,
                         data = data,
                         sims = 1000,
                         progress_bar = FALSE) {

  model_1 <- as.formula(model_1)
  model_2 <- as.formula(model_2)

  # Outcome variable must match across models ----
  if (!identical(all.vars(model_1[[2]]), all.vars(model_2[[2]]))) {
    stop("The outcome variable must be the same in both models.")
  }

  assignment_vec <- data[[assignment]]

  # Input validation ----
  # Only the assignment column is permuted, so columns derived from it (say
  # z1 and z2 built from a factorial cell variable) keep their observed values.
  # If neither model mentions the assignment, every permutation reproduces the
  # observed F-statistic and the null distribution collapses to a point mass,
  # silently returning p = 1.
  if (!assignment %in% c(all.vars(model_1), all.vars(model_2))) {
    stop(
      "The assignment variable '", assignment, "' does not appear in either model.\n",
      "ri2 permutes only the assignment column, so any column derived from it is ",
      "not recomputed and the null distribution would be degenerate.\n",
      "If your models use variables derived from the assignment (for example ",
      "factors derived from a factorial cell variable), derive them inside a ",
      "test_function instead, which sees the permuted assignment on every draw."
    )
  }
  if (anyNA(assignment_vec)) {
    stop("Missing values in assignment variable '", assignment,
         "'. Please remove NA rows before calling conduct_ri().")
  }
  for (v in all.vars(model_1[[2]])) {
    if (!v %in% names(data)) stop("Variable '", v, "' not found in data.")
    if (anyNA(data[[v]])) {
      stop("Missing values in '", v,
           "'. Please remove NA rows before calling conduct_ri().")
    }
  }

  # Extract (possibly transformed) outcome and design matrices ----
  mf <- model.frame(model_1, data = data)
  outcome_vec <- model.response(mf)
  outcome_name <- names(mf)[1]

  design_matrix_1 <- model.matrix.default(model_1, data = data)
  design_matrix_2 <- model.matrix.default(model_2, data = data)

  if (nrow(design_matrix_1) != nrow(design_matrix_2)) {
    stop("The number of complete observations must be the same in both models.")
  }

  # Weights function — recomputed for every assignment vector ----
  make_weights <- function(Z_vec) {
    w <- if (IPW) {
      1 / obtain_condition_probabilities(declaration, assignment = Z_vec)
    } else {
      NULL
    }
    if (!is.null(sampling_weights)) {
      sw <- data[[sampling_weights]]
      w <- if (!is.null(w)) w * sw else sw
    }
    w
  }

  # Helper: fit a model and return coefficients ----
  fit_coefs <- function(y_vec, dm, weights) {
    fit <- tidy(lm_robust_fit(
      y            = matrix(y_vec, dimnames = list(NULL, outcome_name)),
      X            = dm,
      weights      = weights,
      ci           = FALSE,
      cluster      = NULL,
      alpha        = 0.05,
      se_type      = "none",
      return_vcov  = FALSE,
      try_cholesky = FALSE,
      has_int      = TRUE
    ))
    jx <- intersect(c("term", "coefficient_name"), names(fit))[1]
    beta_ix <- intersect(c("estimate", "coefficients"), names(fit))[1]
    fit[fit[[jx]] %in% colnames(dm), beta_ix, drop = TRUE]
  }

  # Observed F-statistic ----
  coefs_1 <- fit_coefs(outcome_vec, design_matrix_1, make_weights(assignment_vec))
  coefs_2 <- fit_coefs(outcome_vec, design_matrix_2, make_weights(assignment_vec))
  ssr_1 <- sum((outcome_vec - design_matrix_1 %*% coefs_1) ^ 2)
  ssr_2 <- sum((outcome_vec - design_matrix_2 %*% coefs_2) ^ 2)
  f_obs <- (ssr_1 - ssr_2) / (ncol(design_matrix_2) - ncol(design_matrix_1)) /
    (ssr_2 / (length(outcome_vec) - ncol(design_matrix_2)))

  # Potential outcomes under sharp null ----
  if (length(sharp_hypothesis) == 1) {
    sharp_hypothesis <- rep(sharp_hypothesis, length(unique(assignment_vec)) - 1)
  }
  pos_mat <- generate_pos(
    Y              = outcome_vec,
    assignment_vec = assignment_vec,
    sharp_hypothesis = sharp_hypothesis
  )

  # Permutation matrix ----
  if (is.null(permutation_matrix)) {
    permutation_matrix <- obtain_permutation_matrix(declaration,
      maximum_permutations = sims
    )
  }

  # Permuted F-statistics ----
  ri_function <- function(Z_sim) {
    if (is.factor(assignment_vec)) {
      Z_sim <- factor(Z_sim, levels = levels(assignment_vec))
    }

    data_sim <- data
    data_sim[[assignment]] <- Z_sim

    dm_sim_1 <- model.matrix.default(model_1, data = data_sim)
    dm_sim_2 <- model.matrix.default(model_2, data = data_sim)

    y_sim <- switching_equation(pos_mat = pos_mat, assignment_vec = Z_sim)
    w_sim <- make_weights(Z_sim)

    c1 <- fit_coefs(y_sim, dm_sim_1, w_sim)
    c2 <- fit_coefs(y_sim, dm_sim_2, w_sim)
    ssr_s1 <- sum((y_sim - dm_sim_1 %*% c1) ^ 2)
    ssr_s2 <- sum((y_sim - dm_sim_2 %*% c2) ^ 2)

    (ssr_s1 - ssr_s2) / (ncol(dm_sim_2) - ncol(dm_sim_1)) /
      (ssr_s2 / (length(y_sim) - ncol(dm_sim_2)))
  }

  null_distribution <- if (progress_bar) {
    pbapply::pbapply(permutation_matrix, 2, ri_function)
  } else {
    apply(permutation_matrix, 2, ri_function)
  }

  sims_df <- data.frame(
    est_sim = null_distribution,
    est_obs = f_obs,
    term    = "F-statistic"
  )

  structure(list(sims_df = sims_df), class = "ri2")
}
