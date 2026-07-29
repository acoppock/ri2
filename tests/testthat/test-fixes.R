library(ri2)

# ---- setup ----
set.seed(42)
N <- 60
decl2 <- randomizr::declare_ra(N = N, m = 30)
Z2 <- randomizr::conduct_ra(decl2)
Y2 <- 0.4 * Z2 + rnorm(N)
dat2 <- data.frame(Y2, Z2)

# ---- #33: multi-arm permutation matrix — p-values must be sensible ----

test_that("multi-arm: null distributions for arms 2+ are centered near zero", {
  set.seed(1)
  N3 <- 90
  decl3 <- randomizr::declare_ra(N = N3, num_arms = 3)
  Z3 <- randomizr::conduct_ra(decl3)
  # No true effect on any arm
  Y3 <- rnorm(N3)
  dat3 <- data.frame(Y3, Z3)

  out <- conduct_ri(
    Y3 ~ Z3, declaration = decl3, assignment = "Z3",
    sharp_hypothesis = 0, data = dat3, sims = 200
  )

  sims_T2 <- out$sims_df[out$sims_df$term == "Z3T2", "est_sim"]
  sims_T3 <- out$sims_df[out$sims_df$term == "Z3T3", "est_sim"]

  # Each arm's null distribution should be centered near zero
  expect_lt(abs(mean(sims_T2)), 0.3)
  expect_lt(abs(mean(sims_T3)), 0.3)

  # SD of null distributions should be similar (both are DiM under the null)
  expect_lt(abs(sd(sims_T2) - sd(sims_T3)) / sd(sims_T2), 0.5)
})

test_that("multi-arm: arm with true effect has low p-value, null arm does not", {
  set.seed(7)
  N3 <- 120
  decl3 <- randomizr::declare_ra(N = N3, num_arms = 3)
  Z3 <- randomizr::conduct_ra(decl3)
  # Large effect on T2, no effect on T3
  Y3 <- 3.0 * (Z3 == "T2") + rnorm(N3)
  dat3 <- data.frame(Y3, Z3)

  out <- conduct_ri(
    Y3 ~ Z3, declaration = decl3, assignment = "Z3",
    sharp_hypothesis = 0, data = dat3, sims = 500
  )

  td <- tidy(out)
  p_T2 <- td$p.value[td$term == "Z3T2"]
  p_T3 <- td$p.value[td$term == "Z3T3"]

  expect_lt(p_T2, 0.05)   # large true effect should be detected
  expect_gt(p_T3, 0.05)   # no true effect should not be detected
})

# ---- #20: outcome transformation ----

test_that("log(Y) transformation is applied before computing the test statistic", {
  set.seed(3)
  Y_pos <- abs(rnorm(N)) + 0.5
  dat_log <- data.frame(Y_pos, Z2)

  out_raw <- conduct_ri(
    Y_pos ~ Z2, declaration = decl2, assignment = "Z2",
    sharp_hypothesis = 0, data = dat_log, sims = 50
  )
  out_log <- conduct_ri(
    log(Y_pos) ~ Z2, declaration = decl2, assignment = "Z2",
    sharp_hypothesis = 0, data = dat_log, sims = 50
  )

  # Observed estimates should differ (log vs raw outcome)
  expect_false(isTRUE(all.equal(
    unique(out_raw$sims_df$est_obs),
    unique(out_log$sims_df$est_obs)
  )))

  # log(Y) estimate should match manual calculation
  manual <- mean(log(dat_log$Y_pos[Z2 == 1])) - mean(log(dat_log$Y_pos[Z2 == 0]))
  expect_equal(unique(out_log$sims_df$est_obs), manual, tolerance = 1e-10)
})

# ---- #28: class name is "ri2", not "ri" ----

test_that("conduct_ri returns class ri2, not ri", {
  out <- conduct_ri(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
                    data = dat2, sims = 50)
  expect_s3_class(out, "ri2")
  expect_false(inherits(out, "ri"))
})

# ---- sampling weights ----

test_that("sampling_weights are applied to observed and permuted fits", {
  dat_sw <- dat2
  dat_sw$sw <- runif(N, 0.5, 1.5)

  out_plain <- conduct_ri(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
                          data = dat2, sims = 50, IPW = FALSE)
  out_sw <- conduct_ri(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
                       sampling_weights = "sw", data = dat_sw, sims = 50, IPW = FALSE)

  # Observed estimates differ when sampling weights are applied
  expect_false(isTRUE(all.equal(
    unique(out_plain$sims_df$est_obs),
    unique(out_sw$sims_df$est_obs)
  )))
})

# ---- clustered studentized SEs ----

test_that("studentize with clusters runs without error", {
  dat_cl <- dat2
  dat_cl$cl <- rep(1:10, each = 6)

  out <- conduct_ri(
    Y2 ~ Z2, declaration = decl2, assignment = "Z2",
    studentize = TRUE, clusters = "cl",
    data = dat_cl, sims = 50
  )
  expect_s3_class(out, "ri2")
  expect_true(grepl("studentized", unique(out$sims_df$term)))
})

test_that("studentize without clusters still works (HC2)", {
  out <- conduct_ri(
    Y2 ~ Z2, declaration = decl2, assignment = "Z2",
    studentize = TRUE, data = dat2, sims = 50
  )
  expect_s3_class(out, "ri2")
})

# ---- NA checks ----

test_that("NA in outcome raises an informative error", {
  dat_na <- dat2
  dat_na$Y2[5] <- NA
  expect_error(
    conduct_ri(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
               data = dat_na, sims = 10),
    "Missing values"
  )
})

test_that("NA in assignment raises an informative error", {
  dat_na <- dat2
  dat_na$Z2[5] <- NA
  expect_error(
    conduct_ri(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
               data = dat_na, sims = 10),
    "Missing values"
  )
})

# ---- rounding / tolerance: observed stat counted as extreme ----

test_that("observed statistic is always counted as extreme (two-tailed)", {
  # With exact enumeration, the observed Z appears in the permutation matrix.
  # The p-value must be >= 1/n_perms > 0.
  decl_small <- randomizr::declare_ra(N = 7, m = 2)
  table_2.2 <- data.frame(
    d = c(1, 0, 0, 0, 0, 0, 1),
    y = c(15, 15, 20, 20, 10, 15, 30)
  )
  out <- conduct_ri(y ~ d, declaration = decl_small, assignment = "d",
                    sharp_hypothesis = 0, data = table_2.2)
  p <- tidy(out)$p.value
  expect_gt(p, 0)
})

# ---- ri_ci ----

test_that("ri_ci returns a data frame with ci_lower < ci_upper", {
  out_ci <- ri_ci(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
                  data = dat2, sims = 200, n_grid = 20)
  expect_s3_class(out_ci, "data.frame")
  expect_true(all(c("term", "ci_lower", "ci_upper", "alpha") %in% names(out_ci)))
  expect_lt(out_ci$ci_lower, out_ci$ci_upper)
})

test_that("ri_ci interval contains the observed estimate", {
  set.seed(99)
  out <- conduct_ri(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
                    data = dat2, sims = 100)
  beta_hat <- unique(out$sims_df$est_obs)

  out_ci <- ri_ci(Y2 ~ Z2, declaration = decl2, assignment = "Z2",
                  data = dat2, sims = 200, n_grid = 30)
  expect_gte(beta_hat, out_ci$ci_lower)
  expect_lte(beta_hat, out_ci$ci_upper)
})

# GH PR #32 (mreece13): bins = nrow / 20 is fractional whenever the row count
# is not a multiple of 20, and geom_histogram() requires a whole number.
test_that("plot.ri2 builds when the simulation count is not a multiple of 20", {
  set.seed(101)
  N <- 60
  decl <- randomizr::declare_ra(N = N, m = 30)
  Z <- randomizr::conduct_ra(decl)
  dat <- data.frame(Y = 0.3 * Z + rnorm(N), Z = Z)

  out <- conduct_ri(Y ~ Z, declaration = decl, data = dat, sims = 1010)
  expect_equal(nrow(out$sims_df) %% 20, 10)

  built <- expect_no_warning(ggplot2::ggplot_build(plot(out)))
  expect_gt(nrow(built$data[[1]]), 0)
})

# The design matrix used to be patched by overwriting only the assignment
# columns, which left interaction columns at their observed values. The columns
# then disagreed about what the assignment was, narrowing the null distribution
# and making p-values anti-conservative.
test_that("terms involving the assignment are rebuilt under permutation", {
  set.seed(4)
  N <- 60
  decl <- randomizr::declare_ra(N = N, prob = 0.5)
  Z1 <- randomizr::conduct_ra(decl)
  Z2 <- rbinom(N, 1, 0.5)
  Y <- 0.5 * Z1 + 0.8 * Z2 + 1.2 * Z1 * Z2 + rnorm(N)
  dat <- data.frame(Y, Z1, Z2)

  pm <- randomizr::obtain_permutation_matrix(decl, maximum_permutations = 25)

  out <- conduct_ri(Y ~ Z1 * Z2, assignment = "Z1", declaration = decl,
                    permutation_matrix = pm, IPW = FALSE, data = dat)

  expected <- apply(pm, 2, function(z) coef(lm(Y ~ z * Z2))[["z"]])

  expect_equal(out$sims_df$est_sim, expected, tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("additive covariate adjustment is unaffected by the rebuild path", {
  set.seed(5)
  N <- 60
  decl <- randomizr::declare_ra(N = N, prob = 0.5)
  Z <- randomizr::conduct_ra(decl)
  X <- rnorm(N)
  Y <- 0.4 * Z + 0.7 * X + rnorm(N)
  dat <- data.frame(Y, Z, X)

  pm <- randomizr::obtain_permutation_matrix(decl, maximum_permutations = 25)

  out <- conduct_ri(Y ~ Z + X, assignment = "Z", declaration = decl,
                    permutation_matrix = pm, IPW = FALSE, data = dat)

  expected <- apply(pm, 2, function(z) coef(lm(Y ~ z + X))[["z"]])

  expect_equal(out$sims_df$est_sim, expected, tolerance = 1e-8,
               ignore_attr = TRUE)
})

# Only the assignment column is permuted, so a formula or model that does not
# mention it yields a degenerate null (every permutation reproduces the observed
# statistic) and a p-value of 1 regardless of the data. Guard rather than run.
test_that("conduct_ri errors when the assignment is absent from the model", {
  set.seed(6)
  N <- 60
  decl <- randomizr::declare_ra(N = N, conditions = c("00", "01", "10", "11"))
  cell <- randomizr::conduct_ra(decl)
  z1 <- as.numeric(substr(cell, 1, 1))
  z2 <- as.numeric(substr(cell, 2, 2))
  dat <- data.frame(Y = 0.5 * z1 + 0.8 * z2 + rnorm(N), cell, z1, z2)

  expect_error(
    conduct_ri(model_1 = Y ~ z1 + z2, model_2 = Y ~ z1 * z2,
               assignment = "cell", declaration = decl, data = dat, sims = 20),
    "does not appear in either model"
  )
  expect_error(
    conduct_ri(Y ~ z1 * z2, assignment = "cell", declaration = decl,
               data = dat, sims = 20),
    "does not appear in the formula"
  )

  # legitimate uses of the same design must still run
  expect_s3_class(
    conduct_ri(Y ~ cell, assignment = "cell", declaration = decl,
               data = dat, sims = 20),
    "ri2"
  )
})

test_that("sampling_weights with test_function errors instead of being ignored", {
  set.seed(7)
  N <- 40
  decl <- randomizr::declare_ra(N = N, m = 20)
  dat <- data.frame(Z = randomizr::conduct_ra(decl), sw = runif(N, 0.5, 2))
  dat$Y <- 0.3 * dat$Z + rnorm(N)

  expect_error(
    conduct_ri(test_function = function(d) mean(d$Y), assignment = "Z",
               outcome = "Y", declaration = decl, sampling_weights = "sw",
               data = dat, sims = 20),
    "not used with the test_function interface"
  )
})

# A factorial is a multi-arm trial over its cells; deriving the factors inside
# the test function gives exactly the right null distribution (#27).
test_that("factorial via cells and a test_function matches manual permutation", {
  set.seed(8)
  N <- 80
  decl <- randomizr::declare_ra(N = N, conditions = c("00", "01", "10", "11"))
  cell <- randomizr::conduct_ra(decl)
  z1 <- as.numeric(substr(cell, 1, 1))
  z2 <- as.numeric(substr(cell, 2, 2))
  Y <- 0.5 * z1 + 0.8 * z2 + 1.2 * z1 * z2 + rnorm(N)
  dat <- data.frame(Y, cell)

  interaction_coef <- function(data) {
    a <- as.numeric(substr(data$cell, 1, 1))
    b <- as.numeric(substr(data$cell, 2, 2))
    coef(lm(data$Y ~ a * b))[["a:b"]]
  }

  pm <- randomizr::obtain_permutation_matrix(decl, maximum_permutations = 25)
  out <- conduct_ri(test_function = interaction_coef, assignment = "cell",
                    outcome = "Y", declaration = decl,
                    permutation_matrix = pm, sharp_hypothesis = 0, data = dat)

  expected <- apply(pm, 2, function(cs) {
    a <- as.numeric(substr(cs, 1, 1))
    b <- as.numeric(substr(cs, 2, 2))
    coef(lm(Y ~ a * b))[["a:b"]]
  })

  expect_equal(out$sims_df$est_sim, expected, tolerance = 1e-8,
               ignore_attr = TRUE)
})

# The bin count came from the total row count, but plot.ri2() facets by term, so
# a multi-arm result got more bins per panel than a two-arm one over the same
# number of simulations per panel.
test_that("histogram bins do not scale with the number of arms", {
  set.seed(102)
  bins_drawn <- function(arms) {
    decl <- randomizr::declare_ra(N = 150, num_arms = arms)
    Z <- randomizr::conduct_ra(decl)
    dat <- data.frame(Y = rnorm(150), Z = Z)
    out <- conduct_ri(Y ~ Z, declaration = decl, data = dat, sims = 1000)
    built <- ggplot2::ggplot_build(plot(out))
    panel_1 <- built$data[[1]][built$data[[1]]$PANEL == 1, ]
    length(unique(panel_1$x))
  }

  two_arm <- bins_drawn(2)
  expect_equal(bins_drawn(3), two_arm)
  expect_equal(bins_drawn(4), two_arm)
})
