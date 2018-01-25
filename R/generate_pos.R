
generate_pos <- function(Y,
                         assignment_vec,
                         sharp_hypothesis) {
  condition_names <- sort(unique(assignment_vec))

  pos_mat <-
    matrix(
      data = NA,
      ncol = length(condition_names),
      nrow = length(Y),
      dimnames = list(NULL, condition_names)
    )

  # Knock down to baseline
  for (i in 2:length(condition_names)) {
    pos_mat[, paste(condition_names[1])] <-
      Y - (assignment_vec == condition_names[i]) * sharp_hypothesis[i - 1]
  }

  # Build back up by conditions
  for (i in 2:length(condition_names)) {
    pos_mat[, paste(condition_names[i])] <-
      pos_mat[, paste(condition_names[1])] + sharp_hypothesis[i - 1]
  }

  return(pos_mat)
}
