context("Multi Arm Trial")


test_that("Basic Multi Arm Trial", {
  N <- 100
  declaration <- randomizr::declare_ra(N = N, num_arms = 3)

  Z <- randomizr::conduct_ra(declaration)
  X <- rnorm(N)
  Y <- .9 * X + .5 * (Z == "T2") + rnorm(N)
  W <- runif(N)
  df <- data.frame(Y, X, Z, W)

  ri_out <-
    conduct_ri(
      formula = Y ~ Z,
      declaration = declaration,
      assignment = "Z",
      sharp_hypothesis = 0,
      data = df, sims = 50
    )

  ri_out
  summary(ri_out)
  plot(ri_out)

  summary(ri_out, p = "lower")
  plot(ri_out, p = "lower")

  summary(ri_out, p = "upper")
  plot(ri_out, p = "upper")

  expect_true(TRUE)
})

test_that("A character assignment column gives the same answer as a factor", {
  # A three-arm trial whose assignment column is read in as character, which is
  # what read.csv() gives. Before the fix, the conditional re-randomization for
  # each pairwise comparison returned a factor with only two levels, the design
  # matrix for that comparison was rebuilt with one column recycled into two,
  # and every simulated estimate came back NA, so every p-value was NA.
  N <- 90
  declaration <- randomizr::declare_ra(N = N, conditions = c("Control", "T1", "T2"))
  set.seed(42)
  Z <- as.character(randomizr::conduct_ra(declaration))
  X <- rnorm(N)
  Y <- 0.5 * (Z == "T2") + rnorm(N)
  df_chr <- data.frame(Y, X, Z, stringsAsFactors = FALSE)
  df_fac <- df_chr
  df_fac$Z <- factor(df_fac$Z, levels = c("Control", "T1", "T2"))

  set.seed(1)
  out_chr <- conduct_ri(Y ~ Z, declaration = declaration, data = df_chr, sims = 100)
  set.seed(1)
  out_fac <- conduct_ri(Y ~ Z, declaration = declaration, data = df_fac, sims = 100)

  s_chr <- summary(out_chr)
  s_fac <- summary(out_fac)
  expect_false(anyNA(s_chr$two_tailed_p_value))
  expect_equal(s_chr$estimate, s_fac$estimate)
  expect_equal(s_chr$two_tailed_p_value, s_fac$two_tailed_p_value)

  # The same for a formula that rebuilds the design matrix on every draw.
  set.seed(1)
  out_chr_int <- conduct_ri(Y ~ Z * X, declaration = declaration, data = df_chr, sims = 100)
  set.seed(1)
  out_fac_int <- conduct_ri(Y ~ Z * X, declaration = declaration, data = df_fac, sims = 100)
  expect_false(anyNA(summary(out_chr_int)$two_tailed_p_value))
  expect_equal(summary(out_chr_int)$two_tailed_p_value,
               summary(out_fac_int)$two_tailed_p_value)
})
