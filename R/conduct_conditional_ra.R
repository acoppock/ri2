#' @importFrom randomizr simple_ra complete_ra block_ra cluster_ra block_and_cluster_ra
conduct_conditional_ra <-
  function(declaration, assignment_vec, conditions) {

    # Checks

    if (length(conditions) != 2) {
      stop("conditions must be of length 2.")
    }

    if (!all(conditions %in% assignment_vec)) {
      stop(
        "Both condition names supplied to the conditions argument must be present in the assigment vector supplied to assignment_vec."
      )
    }

    # Every branch below is an `if` with no `else`, and assignment_vec_new starts
    # as the observed assignment, so an unrecognised declaration class used to
    # fall through and return the observed assignment for every draw. That gives
    # a null distribution with one unique value and a p-value of 1 whatever the
    # data, with no error. Fail loudly instead.
    known <- c("ra_simple", "ra_complete", "ra_blocked", "ra_clustered",
               "ra_blocked_and_clustered")
    if (!inherits(declaration, known)) {
      stop(
        "Conditional random assignment is not implemented for a declaration of class ",
        paste(setdiff(class(declaration), "ra_declaration"), collapse = "/"),
        ".\nSupply `permutation_matrix` to conduct_ri() directly, or use a ",
        "declaration from declare_ra() of one of these kinds: ",
        paste(sub("^ra_", "", known), collapse = ", "), ".",
        call. = FALSE
      )
    }

    assignment_vec_new <- assignment_vec

    if (inherits(declaration, "ra_simple")) {
      prob_each_local <-
        declaration$probabilities_matrix[1, paste0("prob_", conditions)]
      prob_each_local <- prob_each_local / (sum(prob_each_local))


      assignment_vec_new[assignment_vec %in% conditions] <-
        simple_ra(
          N = sum(assignment_vec %in% conditions),
          prob_each = prob_each_local,
          conditions = conditions,
          check_inputs = FALSE
        )
    }

    if (inherits(declaration, "ra_complete")) {
      prob_each_local <-
        declaration$probabilities_matrix[1, paste0("prob_", conditions)]
      prob_each_local <- prob_each_local / (sum(prob_each_local))

      assignment_vec_new[assignment_vec %in% conditions] <-
        complete_ra(
          N = sum(assignment_vec %in% conditions),
          prob_each = prob_each_local,
          conditions = conditions,
          check_inputs = FALSE
        )
    }

    if (inherits(declaration, "ra_blocked")) {
      block_prob_each_local <- by(
        declaration$probabilities_matrix,
        INDICES = declaration$blocks,
        FUN = function(x) {
          x[1, paste0("prob_", conditions)]
        }
      )
      block_prob_each_local <- do.call("rbind", block_prob_each_local)
      block_prob_each_local <-
        block_prob_each_local / rowSums(block_prob_each_local)

      assignment_vec_new[assignment_vec %in% conditions] <-
        block_ra(
          blocks = declaration$blocks[assignment_vec %in% conditions],
          block_prob_each = block_prob_each_local,
          conditions = conditions,
          check_inputs = FALSE
        )
    }

    if (inherits(declaration, "ra_clustered")) {
      prob_each_local <-
        declaration$probabilities_matrix[1, paste0("prob_", conditions)]
      prob_each_local <- prob_each_local / (sum(prob_each_local))

      assignment_vec_new[assignment_vec %in% conditions] <-
        cluster_ra(
          clusters = declaration$clusters[assignment_vec %in% conditions],
          prob_each = prob_each_local,
          conditions = conditions,
          check_inputs = FALSE
        )
    }

    if (inherits(declaration, "ra_blocked_and_clustered")) {
      block_prob_each_local <- by(
        declaration$probabilities_matrix,
        INDICES = declaration$blocks,
        FUN = function(x) {
          x[1, paste0("prob_", conditions), drop = FALSE]
        },
        simplify = FALSE
      )
      block_prob_each_local <- do.call("rbind", block_prob_each_local)
      block_prob_each_local <- as.matrix(block_prob_each_local)
      block_prob_each_local <-
        block_prob_each_local / rowSums(block_prob_each_local)

      assignment_vec_new[assignment_vec %in% conditions] <-
        block_and_cluster_ra(
          blocks = declaration$blocks[assignment_vec %in% conditions],
          clusters = declaration$clusters[assignment_vec %in% conditions],
          block_prob_each = block_prob_each_local,
          conditions = conditions,
          check_inputs = FALSE
        )
    }

    return(assignment_vec_new)
  }
