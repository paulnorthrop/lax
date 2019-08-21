context("return levels")

# Check that the MLE of the m-observation return level returned by
# return_level() agrees with the value returned by evd::fgev() when
# prob = 1 - 1 / m.

got_evd <- requireNamespace("evd", quietly = TRUE)

if (got_evd) {
  library(evd)

  # An example from the evd::fgev documentation
  set.seed(4082019)
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  # Fir for GEV parameters
  M1 <- evd::fgev(uvdata)
  adj_M1 <- alogLik(M1)
  # Inference for 100-observation return level
  rl <- return_level(adj_M1, m = 100)
  # Fit for 99.9% quantile (100-observation return level)
  # Set the starting value based on the M1 fit to avoid differences owing
  # to optim's stopping criteria and the differing shape of the loglikelihoods
  my_start <- list(quantile = as.numeric(rl$rl_prof["mle"]),
                   scale = as.numeric(coef(M1)["scale"]),
                   shape = as.numeric(coef(M1)["shape"]))
  M1_prob <- evd::fgev(uvdata, start = my_start, prob = 0.01)

  test_that("return_level() vs evd::fgev", {
    testthat::expect_equal(as.numeric(rl$rl_prof["mle"]),
                           as.numeric(coef(M1_prob)["quantile"]))
  })
}
