switching_equation <- function(pos_mat, assignment_vec) {

  Y <- pos_mat[cbind(
         1:nrow(pos_mat),
         match(assignment_vec, colnames(pos_mat))
  )]

  return(Y)
}
