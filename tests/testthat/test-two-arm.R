context("Two Arm Trial")


test_that("Basic Two Arm Trial", {
  N <- 100
  declaration <- randomizr::declare_ra(N = N, m = 50)

  Z <- randomizr::conduct_ra(declaration)
  X <- rnorm(N)
  Y <- .9 * X + .2 * Z + rnorm(N)
  W <- runif(N)
  df <- data.frame(Y, X, Z, W)

  ri_out <-
    conduct_ri(
      formula = Y ~ Z,
      declaration = declaration,
      assignment = "Z",
      sharp_hypothesis = 0,
      data = df, sims = 100
    )

  plot(ri_out)
  summary(ri_out)

  ri_out <-
    conduct_ri(
      formula = Y ~ Z + X,
      declaration = declaration,
      assignment = "Z",
      sharp_hypothesis = 0,
      data = df, sims = 100
    )


  plot(ri_out)
  summary(ri_out)

  expect_true(TRUE)
})

test_that("A logical assignment column gives the same answer as 0/1", {
  # Z <- conduct_ra(declaration) == 1 gives conditions FALSE and TRUE, which a
  # default declaration (conditions 0 and 1) does not know; before the fix every
  # path stopped with "subscript out of bounds".
  N <- 60
  declaration <- randomizr::declare_ra(N = N)
  set.seed(343)
  Z <- randomizr::conduct_ra(declaration)
  X <- rnorm(N)
  Y <- 0.5 * Z + rnorm(N)
  df_num <- data.frame(Y, X, Z)
  df_lgl <- transform(df_num, Z = Z == 1)

  set.seed(1)
  out_lgl <- summary(conduct_ri(Y ~ Z, declaration = declaration, data = df_lgl, sims = 100))
  set.seed(1)
  out_num <- summary(conduct_ri(Y ~ Z, declaration = declaration, data = df_num, sims = 100))
  expect_equal(out_lgl, out_num)

  set.seed(1)
  out_lgl_int <- summary(conduct_ri(Y ~ Z * X, declaration = declaration, data = df_lgl, sims = 100))
  set.seed(1)
  out_num_int <- summary(conduct_ri(Y ~ Z * X, declaration = declaration, data = df_num, sims = 100))
  expect_equal(out_lgl_int, out_num_int)

  diff_means <- function(data) mean(data$Y[data$Z == 1]) - mean(data$Y[data$Z == 0])
  set.seed(1)
  tf_lgl <- summary(conduct_ri(test_function = diff_means, declaration = declaration,
                               assignment = "Z", outcome = "Y", data = df_lgl, sims = 100))
  set.seed(1)
  tf_num <- summary(conduct_ri(test_function = diff_means, declaration = declaration,
                               assignment = "Z", outcome = "Y", data = df_num, sims = 100))
  expect_equal(tf_lgl, tf_num)

  set.seed(1)
  f_lgl <- summary(conduct_ri(model_1 = Y ~ X, model_2 = Y ~ X + Z, declaration = declaration,
                              assignment = "Z", data = df_lgl, sims = 100))
  set.seed(1)
  f_num <- summary(conduct_ri(model_1 = Y ~ X, model_2 = Y ~ X + Z, declaration = declaration,
                              assignment = "Z", data = df_num, sims = 100))
  expect_equal(f_lgl, f_num)
})
