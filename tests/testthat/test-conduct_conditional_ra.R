context("Conditional RA")


test_that("Conditional RA", {
  N <- 100


  declaration <-
    randomizr::declare_ra(N = N,
                          num_arms = 3,
                          simple = TRUE)

  # Simple ------------------------------------------------------------------


  Z <- randomizr::conduct_ra(declaration)

  Z2 <- conduct_conditional_ra(
    declaration = declaration,
    assignment_vec = Z,
    conditions = c("T1", "T2")
  )

  table(Z, Z2, useNA = "always")



  declaration <- randomizr::declare_ra(N = N, num_arms = 3)

  # Complete ----------------------------------------------------------------


  Z <- randomizr::conduct_ra(declaration)

  Z2 <- conduct_conditional_ra(
    declaration = declaration,
    assignment_vec = Z,
    conditions = c("T1", "T2")
  )

  table(Z)
  table(Z2)
  table(Z, Z2)


  # blocked -----------------------------------------------------------------
  block_var <- rep(c("A", "B", "C"), times = c(50, 100, 200))

  declaration <-
    declare_ra(block_var = block_var, prob_each = c(.1, .8, .1))
  Z <- conduct_ra(declaration)
  table(block_var, Z)


  Z2 <- conduct_conditional_ra(
    declaration = declaration,
    assignment_vec = Z,
    conditions = c("T1", "T2")
  )

  table(Z, block_var)
  table(Z2, block_var)
  table(Z, Z2)


  # Clustered ---------------------------------------------------------------
  clust_var <- rep(letters, times = 1:26)
  declaration <- declare_ra(clust_var = clust_var, num_arms = 3)
  Z <- conduct_ra(declaration)
  table(Z, clust_var)

  Z2 <- conduct_conditional_ra(
    declaration = declaration,
    assignment_vec = Z,
    conditions = c("T1", "T2")
  )

  table(Z, clust_var)
  table(Z2, clust_var)
  table(Z, Z2)

  # Blocked and Clustered ---------------------------------------------------

  clust_var <- rep(letters, times = 1:26)

  block_var <- rep(NA, length(clust_var))
  block_var[clust_var %in% letters[1:5]] <- "block_1"
  block_var[clust_var %in% letters[6:10]] <- "block_2"
  block_var[clust_var %in% letters[11:15]] <- "block_3"
  block_var[clust_var %in% letters[16:20]] <- "block_4"
  block_var[clust_var %in% letters[21:26]] <- "block_5"

  declaration <-
    declare_ra(clust_var = clust_var,
               block_var = block_var,
               num_arms = 3)
  Z <- conduct_ra(declaration)

  Z2 <- conduct_conditional_ra(
    declaration = declaration,
    assignment_vec = Z,
    conditions = c("T1", "T2")
  )

  table(Z, clust_var)
  table(Z2, clust_var)
  table(Z, Z2)

})



test_that("Conditional without conditions!", {
  N <- 100

  declaration <- declare_ra(N)

  Z <- randomizr::conduct_ra(declaration)

  Z2 <- conduct_conditional_ra(
    declaration = declaration,
    assignment_vec = Z,
    conditions = c(0, 1)
  )

  table(Z)
  table(Z2)
  table(Z, Z2)

})

