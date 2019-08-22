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

  # Check that output from plot.retlev() and print.retlev() agree

  plot_output <- plot(rl)
  print_output <- print(rl)
  test_that("plot.retlev vs print.retlev", {
    testthat::expect_equal(plot_output, print_output$rl_prof)
  })

  # Check that when we reduce the CI level (from 95% to 75%) then
  # interval narrows

  plot_output_new <- plot(rl, level = 0.75)
  test_that("Lower CI level gives a narrower interval", {
    testthat::expect_lt(diff(range(plot_output_new)),
                        diff(range(plot_output)))
  })

  # Check that print.retlev and summary.retlev agree
  print_rl <- print(rl)
  summary_rl <- summary(rl)
  test_that("print.rl vs summary.rl, MLEs", {
    # Names differ
    testthat::expect_equivalent(summary_rl$matrix[, 1], print_rl$rl_sym["mle"])
  })
  test_that("print.rl vs summary.rl, EEs", {
    # Names differ
    testthat::expect_equivalent(summary_rl$matrix[, 2], print_rl$rl_se)
  })

}
