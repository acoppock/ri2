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
