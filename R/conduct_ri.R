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
#' @param permutation_matrix An optional matrix of random assignments, typically created by \code{\link{obtain_permutation_matrix}}.
#' @param data A data.frame.
#' @param sims the number of simulations. Defaults to 1000.
#' @param progress_bar logical, defaults to FALSE.  Should a progress bar be displayed in the console?
#'
#' @export
#'
#' @importFrom randomizr declare_ra conduct_ra obtain_condition_probabilities
#'
#' @examples
#'
#' # Data from Gerber and Green Table 2.2
#'
#'
#' # Randomization Inference for the Average Treatment Effect
#'
#' table_2.2 <-
#'     data.frame(d = c(1, 0, 0, 0, 0, 0, 1),
#'                y = c(15, 15, 20, 20, 10, 15, 30))
#'
#' ## Declare randomization procedure
#' declaration <- declare_ra(N = 7, m = 2)
#'
#' ## Conduct Randomization Inference
#' out <- conduct_ri(y ~ d,
#'                       declaration = declaration,
#'                       assignment = "d",
#'                       sharp_hypothesis = 0,
#'                       data = table_2.2)
#'
#' summary(out)
#' plot(out)
#'
#' # Randomization Inference for an Interaction
#'
#'
#' N <- 100
#' declaration <- randomizr::declare_ra(N = N, m = 50)
#'
#' Z <- randomizr::conduct_ra(declaration)
#' X <- rnorm(N)
#' Y <- .9 * X + .2 * Z + 1 * X * Z + rnorm(N)
#' dat <- data.frame(Y, X, Z)
#'
#' ate_obs <- coef(lm(Y ~ Z, data = dat))[2]
#'
#' out <-
#'   conduct_ri(
#'     model_1 = Y ~ Z + X,
#'     model_2 = Y ~ Z + X + Z * X,
#'     declaration = declaration,
#'     assignment = "Z",
#'     sharp_hypothesis = ate_obs,
#'     data = dat, sims = 100
#'   )
#'
#' plot(out)
#' summary(out)
#'
#' summary(out, p = "two-tailed")
#' summary(out, p = "upper")
#' summary(out, p = "lower")
#'
#' # Randomization Inference for arbitrary test statistics
#'
#' ## In this example we're conducting a randomization check (in this case, a balance test).
#'
#' N <- 100
#' declaration <- randomizr::declare_ra(N = N, m = 50)
#'
#' Z <- randomizr::conduct_ra(declaration)
#' X <- rnorm(N)
#' Y <- .9 * X + .2 * Z + rnorm(N)
#' dat <- data.frame(Y, X, Z)
#'
#' balance_fun <- function(data) {
#'     f_stat <- summary(lm(Z ~ X, data = data))$f[1]
#'     names(f_stat) <- NULL
#'     return(f_stat)
#' }
#'
#' ## confirm function works as expected
#' balance_fun(dat)
#'
#' ## conduct randomization inference
#'
#' out <-
#'   conduct_ri(
#'     test_function = balance_fun,
#'     declaration = declaration,
#'     assignment = "Z",
#'     sharp_hypothesis = 0,
#'     data = dat, sims = 100
#'   )
#'
#' plot(out)
#' summary(out)
#'
conduct_ri <- function(formula = NULL,
                       model_1 = NULL,
                       model_2 = NULL,
                       test_function = NULL,
                       assignment = "Z",
                       outcome = NULL,
                       declaration = NULL,
                       sharp_hypothesis = 0,
                       studentize = FALSE,
                       IPW = TRUE,
                       IPW_weights = NULL,
                       sampling_weights = NULL,
                       permutation_matrix = NULL,
                       data,
                       sims = 1000,
                       progress_bar = FALSE) {
  # some error checking -----------------------------------------------------

  if (is.null(declaration) &
    is.null(permutation_matrix)) {
    stop("Please supply either a random assignment declaration or a permutation matrix")
  }
  if (is.null(declaration) & !is.null(permutation_matrix)) {
    declaration <- randomizr::declare_ra(permutation_matrix = permutation_matrix)
    permutation_matrix <- NULL
  }

  # Case 1: ATE -------------------------------------------------------------

  if (!is.null(formula)) {
    ri_out <- conduct_ri_ATE(
      formula = formula,
      assignment = assignment,
      declaration = declaration,
      sharp_hypothesis = sharp_hypothesis,
      studentize = studentize,
      IPW = IPW,
      IPW_weights = IPW_weights,
      sampling_weights = sampling_weights,
      permutation_matrix = permutation_matrix,
      data = data,
      sims = sims,
      progress_bar = progress_bar
    )
  }

  # Case 2: F-test ----------------------------------------------------------

  if (!is.null(model_1) & !is.null(model_2)) {
    ri_out <- conduct_ri_f(
      model_1 = model_1,
      model_2 = model_2,
      assignment = assignment,
      declaration = declaration,
      sharp_hypothesis = sharp_hypothesis,
      IPW = IPW,
      IPW_weights = IPW_weights,
      sampling_weights = sampling_weights,
      permutation_matrix = permutation_matrix,
      data = data,
      sims = sims,
      progress_bar = progress_bar
    )
  }

  # Case 3: Arbitrary Function ----------------------------------------------

  if (!is.null(test_function)) {
    ri_out <- conduct_ri_test_function(
      test_function = test_function,
      assignment = assignment,
      outcome = outcome,
      declaration = declaration,
      sharp_hypothesis = sharp_hypothesis,
      IPW_weights = IPW_weights,
      sampling_weights = sampling_weights,
      permutation_matrix = permutation_matrix,
      data = data,
      sims = sims,
      progress_bar = progress_bar
    )
  }

  if (is.null(formula) &
    is.null(model_1) & is.null(model_2) & is.null(test_function)) {
    stop("You must specify either a formula, models 1 and 2, or a test function.")
  }

  # deal with numerical instability
  ri_out$sims_df <-
    within(ri_out$sims_df, {
      est_sim <- round(est_sim, 10)
      est_obs <- round(est_obs, 10)
    })

  return(ri_out)
}

#' @export
#' @import ggplot2
#'
#'
plot.ri <- function(x, p = "two-tailed", ...) {
  if (p == "two-tailed") {
    x$sims_df <-
      within(
        x$sims_df,
        extreme <- abs(est_sim) >= abs(est_obs)
      )
  } else if (p == "lower") {
    x$sims_df <-
      within(
        x$sims_df,
        extreme <- est_sim <= est_obs
      )
  } else if (p == "upper") {
    x$sims_df <-
      within(
        x$sims_df,
        extreme <- est_sim >= est_obs
      )
  } else {
    stop('p must be either "two-tailed" (the default), "lower", or "upper".')
  }

  summary_fun <-
    function(dat) {
      with(
        dat,
        data.frame(
          est_obs = unique(est_obs),
          Estimate = "Observed Value"
        )
      )
    }


  summary_df <- split(x$sims_df, x$sims_df$coefficient)
  summary_df <- lapply(summary_df[lapply(summary_df, nrow) != 0], FUN = summary_fun)
  summary_df <- do.call(rbind, summary_df)

  summary_df$coefficient <- rownames(summary_df)

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
    facet_wrap(~ coefficient) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title.y = element_blank()
    )
}

#' @export
print.ri <- function(x, p = "two-tailed", ...) {
  print(summary(x, p))
  invisible(summary(x, p))
}

#' @export
#' @importFrom stats quantile
summary.ri <- function(object, p = "two-tailed", ...) {
  if (p == "two-tailed") {
    object$sims_df <-
      within(
        object$sims_df,
        extreme <- abs(est_sim) >= abs(est_obs)
      )
  } else if (p == "lower") {
    object$sims_df <-
      within(
        object$sims_df,
        extreme <- est_sim <= est_obs
      )
  } else if (p == "upper") {
    object$sims_df <-
      within(
        object$sims_df,
        extreme <- est_sim >= est_obs
      )
  } else {
    stop('p must be either "two-tailed" (the default), "lower", or "upper".')
  }


  summary_fun <-
    function(dat) {
      with(
        dat,
        data.frame(
          estimate = unique(est_obs),
          p_value = mean(extreme),
          null_ci_lower = quantile(est_sim, 0.025),
          null_ci_upper = quantile(est_sim, 0.975)
        )
      )
    }

  return_df <- split(object$sims_df, object$sims_df$coefficient)
  return_df <- lapply(return_df[lapply(return_df, nrow) != 0], FUN = summary_fun)
  return_df <- do.call(rbind, return_df)
  return_df$coefficient <- rownames(return_df)
  rownames(return_df) <- NULL
  return_df <-
    return_df[, c(
      "coefficient",
      "estimate",
      "p_value",
      "null_ci_lower",
      "null_ci_upper"
    )]


  if (p == "two-tailed") {
    colnames(return_df)[3] <- "two_tailed_p_value"
  } else if (p == "lower") {
    colnames(return_df)[3] <- "lower_p_value"
  } else if (p == "upper") {
    colnames(return_df)[3] <- "upper_p_value"
  }
  return(return_df)
}
