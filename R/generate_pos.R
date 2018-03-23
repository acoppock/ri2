
generate_pos <- function(Y,
                         assignment_vec,
                         sharp_hypothesis) {
  condition_names <- sort(unique(assignment_vec))

  pos_mat <- matrix(NA, length(Y), length(condition_names))
  colnames(pos_mat) <- condition_names

  # Knock down to baseline
  pos_mat[,1] <- Y - c(0, sharp_hypothesis)[match(assignment_vec, condition_names)]

  # Build back up by conditions
  pos_mat[,-1] <- outer(pos_mat[,1], sharp_hypothesis, `+`)

  return(pos_mat)
}
