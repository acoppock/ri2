#' Add potential outcome columns under a sharp null hypothesis
#'
#' Materializes the schedule of potential outcomes that \code{\link{conduct_ri}}
#' would otherwise build internally, as columns of \code{data}. This is useful
#' as a starting point for a sharp hypothesis that is not a constant shift: add
#' the columns, edit them, then pass their names back to \code{conduct_ri} via
#' \code{potential_outcomes}.
#'
#' One column is added per treatment condition, named
#' \code{<outcome>_<assignment>_<condition>} by default. Under the sharp null,
#' each unit's outcome in the column matching its observed assignment is its
#' observed outcome; the other columns are the hypothesized counterfactuals.
#'
#' @param data A data.frame.
#' @param outcome a character string naming the outcome variable.
#' @param assignment a character string naming the treatment assignment
#'   variable. Defaults to "Z".
#' @param sharp_hypothesis the hypothesized effect of each treatment condition
#'   relative to baseline. A scalar or a vector of length k - 1 behaves as in
#'   \code{\link{conduct_ri}}. To state an effect that varies across units,
#'   supply a vector of length \code{nrow(data)} (two-arm trials) or an
#'   \code{nrow(data)} by k - 1 matrix.
#' @param prefix a character string prepended to the condition names to form the
#'   new column names. Defaults to \code{<outcome>_<assignment>_}.
#'
#' @return \code{data} with one potential outcome column added per condition.
#' @export
#'
#' @examples
#'
#' declaration <- randomizr::declare_ra(N = 100, m = 50)
#' Z <- randomizr::conduct_ra(declaration)
#' X <- rnorm(100)
#' dat <- data.frame(Y = 0.5 * Z + rnorm(100), Z, X)
#'
#' # The constant-effects null of a 0.2 shift, written out explicitly
#' head(add_potential_outcomes(dat, outcome = "Y", assignment = "Z",
#'                             sharp_hypothesis = 0.2))
#'
#' # A null that sharp_hypothesis cannot express: an effect of 0.5 * X per unit
#' dat <- add_potential_outcomes(dat, outcome = "Y", assignment = "Z",
#'                               sharp_hypothesis = 0.5 * dat$X)
#'
#' conduct_ri(Y ~ Z, declaration = declaration,
#'            potential_outcomes = c("Y_Z_0", "Y_Z_1"),
#'            data = dat, sims = 100)
add_potential_outcomes <- function(data,
                                   outcome,
                                   assignment = "Z",
                                   sharp_hypothesis = 0,
                                   prefix = NULL) {

  if (!outcome %in% names(data)) {
    stop("Outcome variable '", outcome, "' not found in data.")
  }
  if (!assignment %in% names(data)) {
    stop("Assignment variable '", assignment, "' not found in data.")
  }

  assignment_vec <- data[[assignment]]
  condition_names <- sort(unique(assignment_vec))
  n <- nrow(data)
  k <- length(condition_names)

  # A per-unit effect is the whole point of writing the schedule out, so accept
  # one: an n by k - 1 matrix, or a length-n vector in a two-arm trial. Scalars
  # and one-value-per-arm vectors behave as in conduct_ri().
  tau <- if (is.matrix(sharp_hypothesis)) {
    if (!identical(dim(sharp_hypothesis), c(n, k - 1L))) {
      stop("A matrix sharp_hypothesis must be ", n, " by ", k - 1,
           " (one row per unit, one column per treatment condition).")
    }
    sharp_hypothesis
  } else if (length(sharp_hypothesis) == 1) {
    matrix(sharp_hypothesis, nrow = n, ncol = k - 1)
  } else if (k == 2 && length(sharp_hypothesis) == n) {
    matrix(sharp_hypothesis, ncol = 1)
  } else if (length(sharp_hypothesis) == k - 1) {
    matrix(sharp_hypothesis, nrow = n, ncol = k - 1, byrow = TRUE)
  } else {
    stop(
      "sharp_hypothesis must be a scalar, a vector of length ", k - 1,
      " (one per treatment condition minus 1), a vector of length ", n,
      " (one per unit, two-arm trials only), or an ", n, " by ", k - 1, " matrix."
    )
  }

  # Hold each unit's observed outcome fixed for the condition it received, and
  # build the counterfactuals around it.
  Y <- data[[outcome]]
  observed_col <- match(assignment_vec, condition_names)
  baseline <- Y - cbind(0, tau)[cbind(seq_len(n), observed_col)]

  pos_mat <- cbind(baseline, baseline + tau)
  colnames(pos_mat) <- condition_names

  if (is.null(prefix)) prefix <- paste0(outcome, "_", assignment, "_")
  colnames(pos_mat) <- paste0(prefix, condition_names)

  data[colnames(pos_mat)] <- as.data.frame(pos_mat)
  data
}

# Build the potential outcome matrix, either from user-supplied columns or from
# a constant-shift sharp hypothesis. Column order is matched to the sorted
# condition names, which is the order conduct_ri reports results in.
resolve_pos <- function(Y,
                        assignment_vec,
                        sharp_hypothesis,
                        potential_outcomes,
                        data) {

  if (is.null(potential_outcomes)) {
    return(generate_pos(
      Y = Y,
      assignment_vec = assignment_vec,
      sharp_hypothesis = sharp_hypothesis
    ))
  }

  condition_names <- sort(unique(assignment_vec))

  if (length(potential_outcomes) != length(condition_names)) {
    stop(
      "potential_outcomes must name one column per treatment condition (",
      length(condition_names), " required: ",
      paste(condition_names, collapse = ", "),
      "), in that order. Received ", length(potential_outcomes), "."
    )
  }

  absent <- setdiff(potential_outcomes, names(data))
  if (length(absent) > 0) {
    stop("Potential outcome column(s) not found in data: ",
         paste(absent, collapse = ", "), ".")
  }

  pos_mat <- as.matrix(data[, potential_outcomes, drop = FALSE])
  if (anyNA(pos_mat)) {
    stop("Missing values in the potential outcome columns. ",
         "Every unit needs a hypothesized outcome under every condition.")
  }
  colnames(pos_mat) <- condition_names

  # Under any sharp hypothesis, the schedule must reproduce what was actually
  # observed for the assignment each unit actually received. A mismatch almost
  # always means the columns were given in the wrong order.
  implied <- switching_equation(pos_mat = pos_mat, assignment_vec = assignment_vec)
  if (!isTRUE(all.equal(unname(implied), unname(Y)))) {
    stop(
      "The potential outcome columns do not reproduce the observed outcome.\n",
      "For each unit, the column matching its observed assignment must equal ",
      "its observed outcome.\n",
      "Columns are matched to conditions in sorted order (",
      paste(condition_names, collapse = ", "),
      "), so check that potential_outcomes is given in that order."
    )
  }

  pos_mat
}
