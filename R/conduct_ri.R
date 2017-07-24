#' Conduct Randomization Inference
#'
#' This function makes it easy to conduct three kinds of randomization inference.
#'
#' 1. Conduct hypothesis tests under the sharp null when the test statistic is the difference-in-means or covariate-adjusted average treatment effect estimate.
#' 2. Conduct "ANOVA" style hypothesis tests, where the f-statistic from two nested models is the test statistic. This procedure is especially helpful when testing interaction terms under null of constant effects.
#' 3. Arbitrary (scalar) test statistics
#'
#' @param formula an object of class formula, as in \code{\link{lm}}. Use formula when conducting significance tests of an Average Treatment Effect estimate under a sharp null hypothesis. For the difference-in-means estimate, do not include covariates. For the OLS covariate-adjusted estimate, include covariates.
#' @param model_1 an object of class formula, as in \code{\link{lm}}. Models 1 and 2 must be "nested." model_1 should be the "restricted" model and model_2 should be the "unrestricted" model.
#' @param model_2 an object of class formula, as in \code{\link{lm}}. Models 1 and 2 must be "nested." model_1 should be the "restricted" model and model_2 should be the "unrestricted" model.
#' @param test_function A function that takes data and returns a scalar test statistic.
#' @param assignment a character string that indicates which variable is randomly assigned. Defaults to "Z".
#' @param outcome a character string that indicates which variable is the outcome variable. Defaults to NULL.
#' @param declaration A random assignment declaration, created by \code{\link{declare_ra}}.
#' @param sharp_hypothesis either a numeric scalar or a numeric vector of length k - 1, where k is the number of treatment conditions. In a two-arm trial, this number is the *hypothesized* difference between the treated and untreated potential potential outcomes for each unit.. In a multi-arm trial, each number in the vector is the hypothesized difference in potential outcomes between the baseline condition and each successive treatment condition.
#' @param studentize logical, defaults to FALSE. Should the test statistic be the t-ratio rather than the estimated ATE? T-ratios will be calculated using HC2 robust standard errors or their clustered equivalent. CLUSTERING NOT YET IMPLEMENTED.
#' @param IPW logical, defaults to TRUE. Should inverse probability weights be calculated?
#' @param IPW_weights a character string that indicates which variable is the existing inverse probability weights vector. Usually unnecessary, as IPW weights will be incorporated automatically if IPW = TRUE. Defaults to NULL.
#' @param sampling_weights a character string that indicates which variable is the sampling weights vector. Optional, defaults to NULL. NOT YET IMPLEMENTED
#' @param permutation_matrix An optional matrix of random assignmnets, typically created by \code{\link{obtain_permutation_matrix}}.
#' @param data A data.frame.
#' @param sims the number of simulations. Defaults to 1000.
#'
#' @export
#'
#' @import dplyr
#' @importFrom randomizr conduct_ra obtain_condition_probabilities
#'
conduct_ri <- function(formula = NULL,
                       model_1 = NULL,
                       model_2 = NULL,
                       test_function = NULL,
                       assignment = "Z",
                       outcome = NULL,
                       declaration,
                       sharp_hypothesis = 0,
                       studentize = FALSE,
                       IPW = TRUE,
                       IPW_weights = NULL,
                       sampling_weights = NULL,
                       permutation_matrix = NULL,
                       data,
                       sims = 1000) {

# Case 1: ATE -------------------------------------------------------------

  if(!is.null(formula)){

ri_out <- conduct_ri_ATE(formula = formula,
                         assignment = assignment,
                         declaration = declaration,
                         sharp_hypothesis = sharp_hypothesis,
                         studentize = studentize,
                         IPW = IPW,
                         IPW_weights = IPW_weights,
                         sampling_weights = sampling_weights,
                         permutation_matrix = permutation_matrix,
                         data = data,
                         sims = sims)
  }

# Case 2: F-test ----------------------------------------------------------

  if(!is.null(model_1) & !is.null(model_2)){

ri_out <- conduct_ri_f(model_1 = model_1,
                       model_2 = model_2,
                       assignment = assignment,
                       declaration = declaration,
                       sharp_hypothesis = sharp_hypothesis,
                       IPW = IPW,
                       IPW_weights = IPW_weights,
                       sampling_weights = sampling_weights,
                       permutation_matrix = permutation_matrix,
                       data = data,
                       sims = sims)
  }

# Case 3: Arbitrary Function ----------------------------------------------

  if(!is.null(test_function)){

ri_out <- conduct_ri_test_function(test_function = test_function,
                                   assignment = assignment,
                                   outcome = outcome,
                                   declaration = declaration,
                                   sharp_hypothesis = sharp_hypothesis,
                                   IPW_weights = IPW_weights,
                                   sampling_weights = sampling_weights,
                                   permutation_matrix = permutation_matrix,
                                   data = data,
                                   sims = sims)
  }

  if(is.null(formula) & is.null(model_1) & is.null(model_2) & is.null(test_function)){
   stop("You must specify either a formula, models 1 and 2, or a test function.")
  }

return(ri_out)

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
    geom_histogram(bins = max(30, nrow(x$sims_df) / 20)) +
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
