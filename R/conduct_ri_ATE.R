#' @importFrom randomizr obtain_permutation_matrix obtain_num_permutations
#' @importFrom estimatr lm_robust_fit
#' @importFrom stats model.matrix.default as.formula model.frame model.response
conduct_ri_ATE <- function(formula,
                           assignment = "Z",
                           declaration,
                           sharp_hypothesis = 0,
                           studentize = FALSE,
                           IPW = TRUE,
                           sampling_weights = NULL,
                           clusters = NULL,
                           permutation_matrix = NULL,
                           data,
                           sims = 1000,
                           progress_bar = FALSE) {

  formula <- as.formula(formula)
  assignment_vec <- data[[assignment]]

  # Input validation ----
  if (anyNA(assignment_vec)) {
    stop("Missing values in assignment variable '", assignment,
         "'. Please remove NA rows before calling conduct_ri().")
  }
  for (v in all.vars(formula[[2]])) {
    if (!v %in% names(data)) stop("Variable '", v, "' not found in data.")
    if (anyNA(data[[v]])) {
      stop("Missing values in '", v,
           "'. Please remove NA rows before calling conduct_ri().")
    }
  }

  # Extract outcome (handles log(Y), sqrt(Y), etc.) and design matrix ----
  mf <- model.frame(formula, data = data)
  outcome_vec <- model.response(mf)
  outcome_name <- names(mf)[1]
  design_matrix <- model.matrix.default(formula, data = data)

  # When the assignment variable enters more than one term (an interaction such
  # as Y ~ Z * X, or Y ~ Z + I(Z^2)), overwriting only the assignment columns
  # would leave the other terms at their observed values, so the design matrix
  # would disagree with itself about what the assignment was. Rebuild it from
  # the permuted data in that case; patch columns in the common single-term
  # case, which is much faster.
  term_labels <- attr(stats::terms(formula, data = data), "term.labels")
  assignment_terms <- vapply(
    term_labels,
    function(label) assignment %in% all.vars(str2lang(label)),
    logical(1)
  )
  rebuild_design <- sum(assignment_terms) > 1

  condition_names <- sort(unique(assignment_vec))

  # Coefficient names ----
  if (is.numeric(assignment_vec)) {
    coefficient_names <- assignment
  } else {
    coefficient_names <- paste0(assignment, condition_names[-1])
  }

  if (length(sharp_hypothesis) != 1 &&
      length(sharp_hypothesis) != length(coefficient_names)) {
    stop(
      "If you supply multiple sharp hypotheses, supply one per treatment condition minus 1 (",
      length(coefficient_names), " required)."
    )
  }
  if (length(sharp_hypothesis) == 1) {
    sharp_hypothesis <- rep(sharp_hypothesis, length(coefficient_names))
  }

  # Potential outcomes under sharp null ----
  pos_mat <- generate_pos(
    Y = outcome_vec,
    assignment_vec = assignment_vec,
    sharp_hypothesis = sharp_hypothesis
  )

  # SE type — clustering supported when studentize = TRUE ----
  clusters_vec <- if (!is.null(clusters)) data[[clusters]] else NULL
  if (studentize) {
    se_type <- if (!is.null(clusters_vec)) "CR2" else "HC2"
  } else {
    se_type <- "none"
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

  # Observed fit ----
  fit_obs <- lm_robust_fit(
    y            = matrix(outcome_vec, dimnames = list(NULL, outcome_name)),
    X            = design_matrix,
    weights      = make_weights(assignment_vec),
    ci           = FALSE,
    cluster      = clusters_vec,
    alpha        = 0.05,
    se_type      = se_type,
    return_vcov  = FALSE,
    try_cholesky = FALSE,
    has_int      = TRUE
  )
  fit_obs <- tidy(fit_obs)

  jx      <- intersect(c("term", "coefficient_name"), names(fit_obs))[1]
  beta_ix <- intersect(c("estimate", "coefficients"),  names(fit_obs))[1]
  se_ix   <- intersect(c("se", "std.error"),            names(fit_obs))[1]

  fit_obs_sub <- fit_obs[fit_obs[[jx]] %in% coefficient_names, , drop = FALSE]
  coefs_obs <- if (studentize) {
    fit_obs_sub[[beta_ix]] / fit_obs_sub[[se_ix]]
  } else {
    fit_obs_sub[[beta_ix]]
  }
  names(coefs_obs) <- coefficient_names
  coefs_obs <- as.list(coefs_obs)

  # Permutation matrix setup ----
  # For the exact test, generate once — valid for all pairwise arm comparisons.
  # For the simulated case, generate a fresh conditional-RA matrix per arm so
  # each arm's null distribution is drawn from the correct conditional design.
  null_distributions <- vector("list", length = length(condition_names) - 1)
  names(null_distributions) <- coefficient_names

  exact_pm <- NULL
  if (is.null(permutation_matrix) &&
      sims >= obtain_num_permutations(declaration)) {
    exact_pm <- obtain_permutation_matrix(declaration, maximum_permutations = sims)
  }

  # Loop over treatment arms vs. baseline ----
  for (i in 2:length(condition_names)) {

    pm_i <- if (!is.null(permutation_matrix)) {
      permutation_matrix          # user-supplied: shared across arms
    } else if (!is.null(exact_pm)) {
      exact_pm                    # exact enumeration: same set for all arms
    } else {
      replicate(                  # simulated: fresh draws for each arm
        sims,
        conduct_conditional_ra(
          declaration,
          assignment_vec = assignment_vec,
          conditions     = condition_names[c(1, i)]
        )
      )
    }

    ri_function <- function(Z_sim) {
      if (is.factor(assignment_vec)) {
        Z_sim <- factor(Z_sim, levels = levels(assignment_vec))
      }

      dm_sim <- if (rebuild_design) {
        data_sim <- data
        data_sim[[assignment]] <- Z_sim
        model.matrix.default(formula, data = data_sim)
      } else {
        dm <- design_matrix
        dm[, coefficient_names] <- model.matrix.default(~ Z_sim)[, -1]
        dm
      }

      outcome_vec_sim <- if (all(sharp_hypothesis == 0)) {
        outcome_vec
      } else {
        switching_equation(pos_mat = pos_mat, assignment_vec = Z_sim)
      }

      fit_sim <- lm_robust_fit(
        y            = matrix(outcome_vec_sim, dimnames = list(NULL, outcome_name)),
        X            = dm_sim,
        weights      = make_weights(Z_sim),
        ci           = FALSE,
        cluster      = clusters_vec,
        alpha        = 0.05,
        se_type      = se_type,
        return_vcov  = FALSE,
        try_cholesky = FALSE,
        has_int      = TRUE
      )
      fit_sim <- tidy(fit_sim)
      fit_sim <- fit_sim[fit_sim[[jx]] %in% coefficient_names[i - 1], , drop = FALSE]

      coefs_sim <- if (studentize) {
        fit_sim[[beta_ix]] / fit_sim[[se_ix]]
      } else {
        fit_sim[[beta_ix]]
      }
      names(coefs_sim) <- coefficient_names[i - 1]
      coefs_sim
    }

    null_distributions[[i - 1]] <- if (progress_bar) {
      pbapply::pbapply(pm_i, 2, ri_function)
    } else {
      apply(pm_i, 2, ri_function)
    }
  }

  sims_list <- mapply(
    FUN      = data.frame,
    est_sim  = null_distributions,
    est_obs  = coefs_obs,
    SIMPLIFY = FALSE
  )

  sims_df <- do.call("rbind", sims_list)
  sims_df$term <- rep(names(sims_list), vapply(sims_list, nrow, integer(1)))

  if (studentize) {
    sims_df$term <- paste0(sims_df$term, " (studentized)")
  }

  structure(list(sims_df = sims_df), class = "ri2")
}
