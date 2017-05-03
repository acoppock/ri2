
generate_pos <- function(Y, Z, sharp_hypothesis){

  condition_names <- sort(unique(Z))

  pos_mat <- matrix(data = NA, ncol = length(condition_names), nrow = length(Y),
                    dimnames = list(NULL, condition_names))
  pos_mat[,paste(condition_names[1])] <- Y - (Z == condition_names[2])*sharp_hypothesis
  pos_mat[,paste(condition_names[2])] <- pos_mat[,paste(condition_names[1])] + sharp_hypothesis

  return(pos_mat)
}
