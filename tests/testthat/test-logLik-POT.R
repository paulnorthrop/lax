# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("POT", quietly = TRUE)) {
  library(POT)

  # POT::fitgpd

  # An example from the POT::fitgpd documentation.
  set.seed(4082019)
  x <- POT::rgpd(200, 1, 2, 0.25)
  fit <- POT::fitgpd(x, 1, "mle")
  temp <- fit
  adj_fit <- alogLik(fit)
  class(temp) <- "POT_fitgpd"

  test_that("POT::fitgpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("POT::fitgpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(adj_fit))
  })
  # Check logLik.POT_fitgpd: trivially correct
  test_that("POT::fitgpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

}
