
#' @export
switching_equation <- function(pos_mat, Z){
  Y <- rep(NA, length(Z))

  condition_names <- sort(unique(Z))

  for(cond in condition_names){
    Y[Z == cond] <- pos_mat[Z == cond, paste(cond)]
  }

  return(Y)
}
