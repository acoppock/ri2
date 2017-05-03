# context("Compare to ri")
#
#
# test_that("Compare to ri",{
#   N <- 100
#   declaration <- declare_ra(N = N, m = 50)
#
#   Z <- conduct_ra(declaration)
#   X <- rnorm(N)
#   Y <- .9*X + .2 * Z + rnorm(N)
#   W <- runif(N)
#   df <- data.frame(Y, X, Z, W)
#
#   ri_out <-
#   conduct_ri(formula = Y~Z,
#              declaration = declaration,
#              assignment = "Z",
#              sharp_hypothesis = 0,
#              data = df,sims = 100000)
#
#   summary(ri_out)
#
#   # compare to ri
#   library(ri)
#
#   perms <- genperms(Z, maxiter = 100000)
#   probs <- genprobexact(Z)
#   ate <- estate(Y, Z, prob = probs)
#   Ys <- genouts(Y, Z, ate = 0)
#   distout <- gendist(Ys, perms, prob = probs)
#   dispdist(distout, ate)
#
#
#   Ys <- genouts(y,Z,ate=ate) ## generate potential outcomes under tau = ATE
#   distout <- gendist(Ys,perms, prob=probs) # generate sampling dist. under tau = ATE
#   dispdist(distout, ate)  ## display characteristics of sampling dist. for inference
#
#
#
#
#
#   ri_out <-
#     conduct_ri(formula = Y ~ Z + X,
#                declaration = declaration,
#                assignment = "Z",
#                sharp_hypothesis = 0,
#                data = df)
#
#
#   plot(ri_out)
#   summary(ri_out)
#
#
#
# })
#
