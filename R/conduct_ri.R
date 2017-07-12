


#' Conduct Randomization Inference
#'
#' @param formula an object of class formula, as in \code{\link{lm}}.
#' @param data A data.frame.
#' @param assignment a character string that indicates which variable is randomly assigned. Defaults to "Z".
#' @param declaration A random assignment declaration, created by \code{\link{declare_ra}}.
#'
#' @export
#'
#' @import dplyr
#' @importFrom randomizr conduct_ra obtain_condition_probabilities
#'
conduct_ri <- function(formula,
                       assignment = "Z",
                       declaration,
                       data,
                       IPW = TRUE,
                       sharp_hypothesis = 0,
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


  # The observed value ------------------------------------------------------

  if (IPW) {
    weights_vec <-
      1 / obtain_condition_probabilities(declaration, assignment = assignment_vec)
    design_matrix <- sqrt(weights_vec) * design_matrix
    outcome_vec <- sqrt(weights_vec) * outcome_vec
  }

  coefs_obs <-
    quick_lm(y = outcome_vec, X = design_matrix)$coefficients
  rownames(coefs_obs) <- colnames(design_matrix)
  coefs_obs <- as.list(coefs_obs[coefficient_names, ])


  # Obtain Hypothesized POs -------------------------------------------------

  pos_mat <- generate_pos(Y = outcome_vec,
                          assignment_vec = assignment_vec,
                          sharp_hypothesis = sharp_hypothesis)


  null_distributions <- vector("list",
                               length = length(condition_names) - 1)

  names(null_distributions) <- coefficient_names

  for (i in 2:length(condition_names)) {
    ri_function <- function() {
      Z_sim <- conduct_conditional_ra(declaration,
                                      assignment_vec = assignment_vec,
                                      conditions = as.character(condition_names[c(1, i)]))

      design_matrix[, coefficient_names] <- model.matrix.default(~ Z_sim)[,-1]

      if (sharp_hypothesis[i - 1] == 0) {
        outcome_vec_sim <- data[, all.vars(formula[[2]])]
      } else{
        outcome_vec_sim <-
          switching_equation(pos_mat = pos_mat, assignment_vec = Z_sim)
      }

      if (IPW) {
        weights_vec_sim <-
          1 / obtain_condition_probabilities(declaration, assignment = Z_sim)
        design_matrix <- sqrt(weights_vec_sim) * design_matrix
        outcome_vec_sim <- sqrt(weights_vec_sim) * outcome_vec_sim
      }


      coefs_sim <-
        quick_lm(y = outcome_vec_sim, X = design_matrix)$coefficients
      rownames(coefs_sim) <- colnames(design_matrix)
      coefs_sim[coefficient_names[i-1], ]
    }

    null_distributions[[i - 1]] <-
      pbapply::pbreplicate(sims, ri_function())
    #null_distribution <- replicate(sims, ri_function())
  }

  sharp_hypothesis <- as.list(sharp_hypothesis)
  names(sharp_hypothesis) <- coefficient_names

  sims_df <-
    mapply(FUN = data.frame,
           est_sim = null_distributions,
           est_obs = coefs_obs,
           SIMPLIFY = FALSE) %>%
    bind_rows(.id = "coefficient")

  return(structure(
    list(
      sims_df = sims_df
    ),
    class = "ri"
  ))
}


#' @export
#' @import ggplot2
#' @import dplyr
#'
#'
plot.ri <- function(x, p = "two-tailed", ...) {

  if (p == "two-tailed") {
    x$sims_df <-
      x$sims_df %>%
      mutate(extreme = abs(est_sim) >= abs(est_obs))

  } else if (p == "lower") {

    x$sims_df <-
      x$sims_df %>%
      mutate(extreme = est_sim <= est_obs)

  } else if (p == "upper") {
    x$sims_df <-
      x$sims_df %>%
      mutate(extreme = est_sim >= est_obs)

  } else {
    stop('p must be either "two-tailed" (the default), "lower", or "upper".')
  }

  summary_df <-
    x$sims_df %>%
    group_by(coefficient) %>%
    summarize(est_obs = unique(est_obs),
              Estimate = "Observed Value")

  ggplot(x$sims_df, aes(x = est_sim, alpha = extreme)) +
    geom_histogram(bins = nrow(x$sims_df) / 20) +
    geom_vline(
      data = summary_df,
      aes(
        xintercept = est_obs,
        linetype = Estimate,
        colour = Estimate
      ),
      show.legend = TRUE
    ) +
    scale_alpha_manual(values = c(0.5, 1), guide = FALSE) +
    xlab("Simulated Estimates") +
    ggtitle("Randomization Inference") +
    facet_wrap(~coefficient) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.y = element_blank())
}

#' @export
print.ri <- function(x, p = "two-tailed", ...) {
  print(summary(x, p))
  invisible(summary(x, p))
}

#' @export
#' @import dplyr
summary.ri <- function(object, p = "two-tailed", ...) {
  if (p == "two-tailed") {
    object$sims_df <-
    object$sims_df %>%
      mutate(extreme = abs(est_sim) >= abs(est_obs))

  } else if (p == "lower") {

    object$sims_df <-
      object$sims_df %>%
      mutate(extreme = est_sim <= est_obs)

  } else if (p == "upper") {
    object$sims_df <-
      object$sims_df %>%
      mutate(extreme = est_sim >= est_obs)

  } else {
    stop('p must be either "two-tailed" (the default), "lower", or "upper".')
  }

  return_df <-
    object$sims_df %>%
    group_by(coefficient) %>%
    summarize(estimate = unique(est_obs),
              p_value = mean(extreme),
              null_ci_lower = quantile(est_sim, 0.025),
              null_ci_upper = quantile(est_sim, 0.975))

  return(return_df)
}
