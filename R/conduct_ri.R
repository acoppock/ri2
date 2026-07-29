#' Conduct Randomization Inference
#'
#' This function makes it easy to conduct three kinds of randomization inference.
#'
#' 1. Conduct hypothesis tests under the sharp null when the test statistic is the difference-in-means or covariate-adjusted average treatment effect estimate.
#' 2. Conduct "ANOVA" style hypothesis tests, where the f-statistic from two nested models is the test statistic. This procedure is especially helpful when testing interaction terms under null of constant effects.
#' 3. Arbitrary (scalar) test statistics
#'
#' @param formula an object of class formula, as in \code{\link{lm}}. Use formula when conducting significance tests of an Average Treatment Effect estimate under a sharp null hypothesis. For the difference-in-means estimate, do not include covariates. For the OLS covariate-adjusted estimate, include covariates. Transformations of the outcome variable such as \code{log(Y) ~ Z} are supported.
#' @param model_1 an object of class formula, as in \code{\link{lm}}. Models 1 and 2 must be "nested." model_1 should be the "restricted" model and model_2 should be the "unrestricted" model.
#' @param model_2 an object of class formula, as in \code{\link{lm}}. Models 1 and 2 must be "nested." model_1 should be the "restricted" model and model_2 should be the "unrestricted" model.
#' @param test_function A function that takes data and returns a scalar test statistic.
#' @param assignment a character string that indicates which variable is randomly assigned. Defaults to "Z".
#' @param outcome a character string that indicates which variable is the outcome variable. Defaults to NULL.
#' @param declaration A random assignment declaration, created by \code{\link[randomizr]{declare_ra}}.
#' @param sharp_hypothesis either a numeric scalar or a numeric vector of length k - 1, where k is the number of treatment conditions. In a two-arm trial, this number is the hypothesized difference between the treated and untreated potential outcomes for each unit. In a multi-arm trial, each number in the vector is the hypothesized difference in potential outcomes between the baseline condition and each successive treatment condition.
#' @param studentize logical, defaults to FALSE. Should the test statistic be the t-ratio rather than the estimated ATE? T-ratios will be calculated using HC2 robust standard errors, or CR2 clustered standard errors when \code{clusters} is specified.
#' @param IPW logical, defaults to TRUE. Should inverse probability weights be calculated?
#' @param sampling_weights a character string indicating which variable in \code{data} contains sampling weights. Sampling weights are fixed across permutations (they reflect the sampling design, not the assignment). When combined with \code{IPW = TRUE}, sampling and inverse probability weights are multiplied together.
#' @param clusters a character string indicating which variable in \code{data} contains the cluster IDs. When supplied with \code{studentize = TRUE}, CR2 clustered standard errors are used instead of HC2.
#' @param permutation_matrix An optional matrix of random assignments, typically created by \code{\link[randomizr]{obtain_permutation_matrix}}.
#' @param data A data.frame.
#' @param sims the number of simulations. Defaults to 1000.
#' @param progress_bar logical, defaults to FALSE. Should a progress bar be displayed in the console?
#' @param p Should "two-tailed", "upper", or "lower" p-values be reported? Defaults to "two-tailed". For two-tailed p-values, whether or not a simulated value is as large or larger than the observed value is determined with respect to the distance to the sharp null.
#'
#' @export
#'
#' @importFrom randomizr declare_ra conduct_ra obtain_condition_probabilities
#' @importFrom estimatr lm_robust_fit tidy
#'
#' @examples
#'
#' # Data from Gerber and Green Table 2.2
#'
#' table_2.2 <-
#'     data.frame(d = c(1, 0, 0, 0, 0, 0, 1),
#'                y = c(15, 15, 20, 20, 10, 15, 30))
#'
#' ## Declare randomization procedure
#' declaration <- randomizr::declare_ra(N = 7, m = 2)
#'
#' ## Conduct Randomization Inference
#' out <- conduct_ri(y ~ d,
#'                   declaration = declaration,
#'                   assignment = "d",
#'                   sharp_hypothesis = 0,
#'                   data = table_2.2)
#'
#' summary(out)
#' plot(out)
#' tidy(out)
#'
#' # Using a custom permutation matrix
#'
#' permutation_matrix <-
#'  matrix(c(0, 0, 0, 0, 0, 0, 1,
#'           0, 0, 0, 0, 0, 1, 0,
#'           0, 0, 0, 0, 1, 0, 0,
#'           0, 0, 0, 1, 0, 0, 0,
#'           0, 0, 1, 0, 0, 0, 0,
#'           0, 1, 0, 0, 0, 0, 0,
#'           1, 0, 0, 0, 0, 0, 0),
#'         ncol = 7)
#'
#' conduct_ri(y ~ d, assignment = "d", data = table_2.2,
#'            permutation_matrix = permutation_matrix)
#'
#'
#' # Randomization Inference for an Interaction
#'
#' N <- 100
#' declaration <- randomizr::declare_ra(N = N, m = 50)
#'
#' Z <- randomizr::conduct_ra(declaration)
#' X <- rnorm(N)
#' Y <- .9 * X + .2 * Z + 1 * X * Z + rnorm(N)
#' dat <- data.frame(Y, X, Z)
#'
#' ate_obs <- coef(lm(Y ~ Z, data = dat))[[2]]
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
#' # Randomization Inference for arbitrary test statistics
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
#' tidy(out)
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
                       sampling_weights = NULL,
                       clusters = NULL,
                       permutation_matrix = NULL,
                       data,
                       sims = 1000,
                       progress_bar = FALSE,
                       p = "two-tailed") {

  if (is.null(declaration) && is.null(permutation_matrix)) {
    stop("Please supply either a random assignment declaration or a permutation matrix.")
  }
  if (is.null(declaration) && !is.null(permutation_matrix)) {
    declaration <- randomizr::declare_ra(permutation_matrix = permutation_matrix)
    permutation_matrix <- NULL
  }
  if (!assignment %in% names(data)) {
    stop("Assignment variable '", assignment, "' not found in data.")
  }

  # Case 1: ATE ----
  if (!is.null(formula)) {
    ri_out <- conduct_ri_ATE(
      formula           = formula,
      assignment        = assignment,
      declaration       = declaration,
      sharp_hypothesis  = sharp_hypothesis,
      studentize        = studentize,
      IPW               = IPW,
      sampling_weights  = sampling_weights,
      clusters          = clusters,
      permutation_matrix = permutation_matrix,
      data              = data,
      sims              = sims,
      progress_bar      = progress_bar
    )
  }

  # Case 2: F-test ----
  if (!is.null(model_1) && !is.null(model_2)) {
    ri_out <- conduct_ri_f(
      model_1           = model_1,
      model_2           = model_2,
      assignment        = assignment,
      declaration       = declaration,
      sharp_hypothesis  = sharp_hypothesis,
      IPW               = IPW,
      sampling_weights  = sampling_weights,
      permutation_matrix = permutation_matrix,
      data              = data,
      sims              = sims,
      progress_bar      = progress_bar
    )
  }

  # Case 3: Arbitrary function ----
  if (!is.null(test_function)) {
    ri_out <- conduct_ri_test_function(
      test_function     = test_function,
      assignment        = assignment,
      outcome           = outcome,
      declaration       = declaration,
      sharp_hypothesis  = sharp_hypothesis,
      sampling_weights  = sampling_weights,
      permutation_matrix = permutation_matrix,
      data              = data,
      sims              = sims,
      progress_bar      = progress_bar
    )
  }

  if (is.null(formula) &&
      is.null(model_1) && is.null(model_2) &&
      is.null(test_function)) {
    stop("You must specify either a formula, models 1 and 2, or a test function.")
  }

  ri_out$p <- p
  ri_out$sharp_hypothesis <- sharp_hypothesis

  ri_out
}

# Compute the "extreme" flag for a given p-value direction.
# A small tolerance handles floating-point near-ties (the observed statistic
# is itself one of the permuted values and should always count as extreme).
ri_extreme <- function(est_sim, est_obs, sharp_hypothesis, p) {
  tol <- sqrt(.Machine$double.eps)
  if (p == "two-tailed") {
    abs(est_sim - sharp_hypothesis) >=
      abs(est_obs - sharp_hypothesis) - tol * (abs(est_obs - sharp_hypothesis) + 1)
  } else if (p == "lower") {
    est_sim <= est_obs + tol * (abs(est_obs) + 1)
  } else if (p == "upper") {
    est_sim >= est_obs - tol * (abs(est_obs) + 1)
  } else {
    stop('p must be "two-tailed" (the default), "lower", or "upper".')
  }
}

#' @export
#' @import ggplot2
plot.ri2 <- function(x, p = NULL, ...) {
  if (is.null(p)) p <- x$p
  x$sims_df$extreme <- ri_extreme(
    x$sims_df$est_sim, x$sims_df$est_obs, x$sharp_hypothesis, p
  )

  summary_df <- do.call(rbind, lapply(
    split(x$sims_df, x$sims_df$term),
    function(dat) data.frame(est_obs = unique(dat$est_obs), Estimate = "Observed Value")
  ))
  summary_df$term <- rownames(summary_df)

  # geom_histogram() requires a whole number of bins, and nrow / 20 is
  # fractional whenever the row count is not a multiple of 20.
  n_bins <- max(30, floor(nrow(x$sims_df) / 20))

  ggplot(x$sims_df, aes(x = est_sim, fill = extreme)) +
    geom_histogram(bins = n_bins) +
    geom_vline(
      data = summary_df,
      aes(xintercept = est_obs, linetype = Estimate, colour = Estimate),
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = c("FALSE" = "grey80", "TRUE" = "steelblue3"),
      guide  = "none"
    ) +
    xlab("Simulated Estimates") +
    ggtitle("Randomization Inference") +
    facet_wrap(~ term) +
    theme_bw() +
    theme(legend.position = "bottom", axis.title.y = element_blank())
}

#' @export
print.ri2 <- function(x, p = NULL, ...) {
  print(summary(x, p))
  invisible(summary(x, p))
}

#' @export
#' @importFrom stats quantile
summary.ri2 <- function(object, p = NULL, ...) {
  if (is.null(p)) p <- object$p
  object$sims_df$extreme <- ri_extreme(
    object$sims_df$est_sim, object$sims_df$est_obs, object$sharp_hypothesis, p
  )

  return_df <- do.call(rbind, lapply(
    split(object$sims_df, object$sims_df$term),
    function(dat) data.frame(estimate = unique(dat$est_obs), p_value = mean(dat$extreme))
  ))
  return_df$term <- rownames(return_df)
  rownames(return_df) <- NULL
  return_df <- return_df[, c("term", "estimate", "p_value")]

  colnames(return_df)[3] <- switch(p,
    "two-tailed" = "two_tailed_p_value",
    "lower"      = "lower_p_value",
    "upper"      = "upper_p_value"
  )
  return_df
}

#' @export
tidy.ri2 <- function(x, ...) {
  ret <- summary(x)
  colnames(ret)[3] <- "p.value"
  ret
}
