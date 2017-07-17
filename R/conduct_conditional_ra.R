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

    assignment_vec_new <- assignment_vec

    if (declaration$ra_type == "simple")  {
      prob_each_local <-
        declaration$probabilities_matrix[1, paste0("prob_", conditions)]
      prob_each_local <- prob_each_local / (sum(prob_each_local))


      assignment_vec_new[assignment_vec %in% conditions] <-
        simple_ra(
          N = sum(assignment_vec %in% conditions),
          prob_each = prob_each_local,
          condition_names = conditions,
          check_inputs = FALSE
        )

    }

    if (declaration$ra_type == "complete") {
      prob_each_local <-
        declaration$probabilities_matrix[1, paste0("prob_", conditions)]
      prob_each_local <- prob_each_local / (sum(prob_each_local))

      assignment_vec_new[assignment_vec %in% conditions] <-
        complete_ra(
          N = sum(assignment_vec %in% conditions),
          prob_each = prob_each_local,
          condition_names = conditions,
          check_inputs = FALSE
        )
    }

    if (declaration$ra_type == "blocked") {
      block_prob_each_local <- by(
        declaration$probabilities_matrix,
        INDICES = declaration$block_var,
        FUN = function(x) {
          x[1, paste0("prob_", conditions)]
        }
      )
      block_prob_each_local <- do.call("rbind", block_prob_each_local)
      block_prob_each_local <-
        block_prob_each_local / rowSums(block_prob_each_local)

      assignment_vec_new[assignment_vec %in% conditions] <-
        block_ra(
          block_var = declaration$block_var[assignment_vec %in% conditions],
          block_prob_each = block_prob_each_local,
          condition_names = conditions,
          check_inputs = FALSE
        )
    }

    if (declaration$ra_type == "clustered") {
      prob_each_local <-
        declaration$probabilities_matrix[1, paste0("prob_", conditions)]
      prob_each_local <- prob_each_local / (sum(prob_each_local))

      assignment_vec_new[assignment_vec %in% conditions] <-
        cluster_ra(
          clust_var = declaration$clust_var[assignment_vec %in% conditions],
          prob_each = prob_each_local,
          condition_names = conditions,
          check_inputs = FALSE
        )
    }

    if (declaration$ra_type == "blocked_and_clustered") {
      block_prob_each_local <- by(
        declaration$probabilities_matrix,
        INDICES = declaration$block_var,
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
          block_var = declaration$block_var[assignment_vec %in% conditions],
          clust_var = declaration$clust_var[assignment_vec %in% conditions],
          block_prob_each = block_prob_each_local,
          condition_names = conditions,
          check_inputs = FALSE
        )
    }

    return(assignment_vec_new)

  }