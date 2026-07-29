
#' @importFrom generics tidy
#' @export
generics::tidy

# randomizr moved from Depends to Imports in 0.5.0, which stopped library(ri2)
# from attaching it. These are re-exported so that declaring and conducting a
# random assignment still works from ri2 alone.

#' @importFrom randomizr declare_ra
#' @export
randomizr::declare_ra

#' @importFrom randomizr conduct_ra
#' @export
randomizr::conduct_ra

#' @importFrom randomizr simple_ra
#' @export
randomizr::simple_ra

#' @importFrom randomizr complete_ra
#' @export
randomizr::complete_ra

#' @importFrom randomizr block_ra
#' @export
randomizr::block_ra

#' @importFrom randomizr cluster_ra
#' @export
randomizr::cluster_ra

#' @importFrom randomizr block_and_cluster_ra
#' @export
randomizr::block_and_cluster_ra

#' @importFrom randomizr obtain_permutation_matrix
#' @export
randomizr::obtain_permutation_matrix

#' @importFrom randomizr obtain_num_permutations
#' @export
randomizr::obtain_num_permutations

#' @importFrom randomizr obtain_condition_probabilities
#' @export
randomizr::obtain_condition_probabilities
