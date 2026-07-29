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
  # observed for the assignment each unit received. A mismatch usually means the
  # columns were given in the wrong order, or that the observed outcome was not
  # held fixed for the condition each unit received.
  implied <- switching_equation(pos_mat = pos_mat, assignment_vec = assignment_vec)
  if (!isTRUE(all.equal(unname(implied), unname(Y)))) {
    stop(
      "The potential outcome columns do not reproduce the observed outcome.\n",
      "For each unit, the column matching its observed assignment must equal ",
      "its observed outcome; only the counterfactual columns are hypothesized.\n",
      "Columns are matched to conditions in sorted order (",
      paste(condition_names, collapse = ", "),
      "), so check that potential_outcomes is given in that order."
    )
  }

  pos_mat
}
