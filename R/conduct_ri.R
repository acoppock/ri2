



#' Conduct Randomization Inference
#'
#' @param formula an object of class formula, as in \code{\link{lm}}.
#'
#' @param data A data.frame.
#' @param alpha The significance level, 0.05 by default.
#' @param assignment a character string that indicates which variable is randomly assigned. Defaults to "Z".
#' @param declaration A random assignment declaration, created by \code{\link{declare_ra}}.
#'
#' @export
#'
#' @importFrom randomizr conduct_ra obtain_condition_probabilities
#'
conduct_ri <- function(formula,
                       alpha = .05,
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

  pos_mat <- generate_pos(Y = outcome_vec,
                          Z = assignment_vec,
                          sharp_hypothesis = sharp_hypothesis)

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
  est_obs <- coefs_obs[assignment,]


  # The ri function ------------------------------------------

  ri_function <- function() {
    Z_sim <- conduct_ra(declaration)
    design_matrix[, assignment] <- Z_sim
    outcome_vec_sim <-
      switching_equation(pos_mat = pos_mat, Z = Z_sim)

    if (IPW) {
      weights_vec_sim <-
        1 / obtain_condition_probabilities(declaration, assignment = Z_sim)
      design_matrix <- sqrt(weights_vec_sim) * design_matrix
      outcome_vec_sim <- sqrt(weights_vec_sim) * outcome_vec_sim
    }


    coefs_sim <-
      quick_lm(y = outcome_vec_sim, X = design_matrix)$coefficients
    rownames(coefs_sim) <- colnames(design_matrix)
    coefs_sim[assignment,]
  }


  null_distribution <- replicate(sims, ri_function())


  return(structure(
    list(
      null_distribution = null_distribution,
      est_obs = est_obs,
      pos_mat = pos_mat,
      ri_function = ri_function,
      sharp_hypothesis = sharp_hypothesis
    ),
    class = "ri"
  ))


}


#' @export
#'
#'
plot.ri <- function(x, ...) {
  ests_df <- data.frame(ests = x$null_distribution,
                        extreme = as.factor(abs(x$null_distribution) >= abs(x$est_obs)))

  results_df <- data.frame(est_obs = x$est_obs,
                           Estimate = "Observed Value")

  ggplot(ests_df, aes(x = ests, alpha = extreme)) +
    geom_histogram(bins = nrow(ests_df)/20) +
    geom_vline(data = results_df, aes(xintercept = est_obs,
                                      linetype = Estimate,
                                      colour = Estimate),
               show.legend = TRUE) +
    scale_alpha_manual(values = c(0.5, 1), guide = FALSE) +
    xlab("Simulated Estimates") +
    ggtitle("Randomization Inference Under Sharp Null Hypothesis",
            paste0("Hypothesized Value of Treatment Effect = ",
                   x$sharp_hypothesis,
                   "; Two-tailed p value = ",
                   round(mean(abs(x$null_distribution) >= abs(x$est_obs)), 3))
            ) +
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
summary.ri <- function(object, p = "two-tailed", ...) {
  if (p == "two-tailed") {
    return_vec <- c(estimate = object$est,
                    p_value = mean(abs(object$null_distribution) >= abs(object$est)),
                    ci_lower = quantile(object$null_distribution, 0.025),
                    ci_upper = quantile(object$null_distribution, 0.975))
  } else if (p == "lower") {
    return_vec <- c(
      estimate = object$est,
      p_value = mean(object$null_distribution <= object$est),
      ci_lower = quantile(object$null_distribution, 0.025),
      ci_upper = quantile(object$null_distribution, 0.975))
  } else {
    return_vec <- c(
      estimate = object$est,
      p_value = mean(object$null_distribution >= object$est),
      ci_lower = quantile(object$null_distribution, 0.025),
      ci_upper = quantile(object$null_distribution, 0.975))
  }
  names(return_vec) <- c("estimate", "p_value", "2.5th Percentile of Null Distribution", "97.5th Percentile of Null Distribution")
  return(return_vec)
}
