switching_equation <- function(pos_mat, assignment_vec) {
  Y <- rep(NA, length(assignment_vec))

  condition_names <- sort(unique(assignment_vec))

  for (cond in condition_names) {
    Y[assignment_vec == cond] <-
      pos_mat[assignment_vec == cond, paste(cond)]
  }

  return(Y)
}
