# GH #19: let users state the sharp hypothesis as a schedule of potential
# outcomes, which allows hypothesized effects that vary across units.

test_that("explicit PO columns reproduce a constant sharp_hypothesis exactly", {
  set.seed(70)
  N <- 120
  decl <- randomizr::declare_ra(N = N, m = 60)
  dat <- data.frame(Z = randomizr::conduct_ra(decl), X = rnorm(N))
  dat$Y <- 0.3 * dat$Z + rnorm(N)

  dat <- add_potential_outcomes(dat, "Y", "Z", sharp_hypothesis = 0.2)
  pm <- randomizr::obtain_permutation_matrix(decl, maximum_permutations = 30)

  via_scalar <- conduct_ri(Y ~ Z, declaration = decl, permutation_matrix = pm,
                           IPW = FALSE, sharp_hypothesis = 0.2, data = dat)
  via_columns <- conduct_ri(Y ~ Z, declaration = decl, permutation_matrix = pm,
                            IPW = FALSE,
                            potential_outcomes = c("Y_Z_0", "Y_Z_1"), data = dat)

  expect_equal(via_scalar$sims_df$est_sim, via_columns$sims_df$est_sim)
})

test_that("a heterogeneous sharp null matches a manual switching equation", {
  set.seed(71)
  N <- 120
  decl <- randomizr::declare_ra(N = N, m = 60)
  dat <- data.frame(Z = randomizr::conduct_ra(decl), X = rnorm(N))
  dat$Y <- 0.3 * dat$Z + rnorm(N)

  dat <- add_potential_outcomes(dat, "Y", "Z", sharp_hypothesis = 0.5 * dat$X)

  # the schedule must hold the observed outcome fixed
  expect_equal(ifelse(dat$Z == 1, dat$Y_Z_1, dat$Y_Z_0), dat$Y)

  pm <- randomizr::obtain_permutation_matrix(decl, maximum_permutations = 30)
  out <- conduct_ri(Y ~ Z, declaration = decl, permutation_matrix = pm,
                    IPW = FALSE, potential_outcomes = c("Y_Z_0", "Y_Z_1"),
                    data = dat)

  expected <- apply(pm, 2, function(z) {
    coef(lm(ifelse(z == 1, dat$Y_Z_1, dat$Y_Z_0) ~ z))[["z"]]
  })
  expect_equal(out$sims_df$est_sim, expected, ignore_attr = TRUE)
})

test_that("per-unit effects work in a multi-arm trial via a matrix", {
  set.seed(72)
  N <- 90
  decl <- randomizr::declare_ra(N = N, num_arms = 3)
  dat <- data.frame(Z = randomizr::conduct_ra(decl), X = rnorm(N))
  dat$Y <- rnorm(N)

  tau <- cbind(0.2 + 0.1 * dat$X, 0.6 - 0.2 * dat$X)
  dat <- add_potential_outcomes(dat, "Y", "Z", sharp_hypothesis = tau)

  expect_true(all(c("Y_Z_T1", "Y_Z_T2", "Y_Z_T3") %in% names(dat)))
  observed <- dat$Y_Z_T1
  observed[dat$Z == "T2"] <- dat$Y_Z_T2[dat$Z == "T2"]
  observed[dat$Z == "T3"] <- dat$Y_Z_T3[dat$Z == "T3"]
  expect_equal(observed, dat$Y)

  out <- conduct_ri(Y ~ Z, declaration = decl,
                    potential_outcomes = c("Y_Z_T1", "Y_Z_T2", "Y_Z_T3"),
                    data = dat, sims = 100)
  expect_equal(nrow(summary(out)), 2)
})

test_that("mis-ordered or inconsistent PO columns are caught", {
  set.seed(73)
  N <- 60
  decl <- randomizr::declare_ra(N = N, m = 30)
  dat <- data.frame(Z = randomizr::conduct_ra(decl))
  dat$Y <- 0.3 * dat$Z + rnorm(N)
  dat <- add_potential_outcomes(dat, "Y", "Z", sharp_hypothesis = 0.4)

  expect_error(
    conduct_ri(Y ~ Z, declaration = decl,
               potential_outcomes = c("Y_Z_1", "Y_Z_0"), data = dat, sims = 10),
    "do not reproduce the observed outcome"
  )
  expect_error(
    conduct_ri(Y ~ Z, declaration = decl, potential_outcomes = "Y_Z_0",
               data = dat, sims = 10),
    "one column per treatment condition"
  )
  expect_error(
    conduct_ri(Y ~ Z, declaration = decl,
               potential_outcomes = c("Y_Z_0", "nope"), data = dat, sims = 10),
    "not found in data"
  )
  expect_error(
    conduct_ri(Y ~ Z, declaration = decl,
               potential_outcomes = c("Y_Z_0", "Y_Z_1"),
               sharp_hypothesis = 0.3, data = dat, sims = 10),
    "not both"
  )
})

test_that("a sharp hypothesis without an outcome name is an error, not ignored", {
  set.seed(74)
  N <- 60
  decl <- randomizr::declare_ra(N = N, m = 30)
  dat <- data.frame(Z = randomizr::conduct_ra(decl))
  dat$Y <- 0.3 * dat$Z + rnorm(N)

  expect_error(
    conduct_ri(test_function = function(d) mean(d$Y), assignment = "Z",
               declaration = decl, sharp_hypothesis = 0.9, data = dat, sims = 10),
    "without naming the outcome variable"
  )
  # the sharp null of no effect needs no imputation, so this still runs
  expect_s3_class(
    conduct_ri(test_function = function(d) mean(d$Y[d$Z == 1]), assignment = "Z",
               declaration = decl, sharp_hypothesis = 0, data = dat, sims = 10),
    "ri2"
  )
})
